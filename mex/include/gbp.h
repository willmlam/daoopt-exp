#ifndef __MEX_GBP_H
#define __MEX_GBP_H

#define GBP_USE_LOG

#include <assert.h>
#include <stdexcept>
#include <stdlib.h>
#include <stdint.h>
#include <cstdarg>
#include <cstring>

#include "graphmodel.h"
#include "regiongraph.h"
#include "alg.h"
#include "indexedHeap.h"

/*
*/

namespace mex {

// Factor graph algorithm specialization for loopy belief propagation
// 

class gbp : public gmAlg, public graphModel, virtual public mxObject {
public:
  typedef graphModel::findex       findex;        // factor index
  typedef graphModel::vindex       vindex;        // variable index
  typedef graphModel::flist        flist;         // collection of factor indices
  typedef vector<findex>           fvec;          // ordered collection of factor indices

  // No actual scheduling options as coded... ???
  MEX_ENUM( Schedule , Fixed,Random,Flood,Priority );

protected:
  graphModel     _gmo;

  regionGraph    _rg;
  vector<fvec>   _descend;
  vector<flist>  _descSet;
  vector<vector<vector<findex> > >   _ances;
  vector<findex> _rorder;
  //vector<Factor> _factors;			// TODO: inherited from graphmodel???
  double         _lnZ;

	double				 _lb;
	vector<index>  _best,_xhat;
	vector<Factor> _beliefs;

  Schedule         _sched;                         // schedule type
  Factor::Distance _dist;                           // message distance measure for priority
  double           _stopIter, _stopObj, _stopMsg, _stopTime;   // and stopping criteria
  size_t           _iter;
  uint32_t         _verbose;
  double           _damp;
  bool             _minimal;
  double           _dObj,_dObjCurr;

  indexedHeap      priority;

public:
  /// Constructors : from nothing, copy, list of factors, or input iterators
  gbp(const vector<Factor>& fs, const vector<VarSet>& regions=vector<VarSet>() ) : graphModel(fs), _rg(fs), _gmo(fs) {
    addRegions(regions);
    setProperties();
  }
  virtual gbp* clone() const { throw std::runtime_error("NOT IMPLEMENTED"); return NULL; }

  size_t nFactors() { return _factors.size(); }
  size_t nRegions() { return _rg.nregions(); }
  size_t nvar() { return _rg.nvar(); }

	size_t rorderSize() { return _rorder.size(); }

  const graphModel& gmOrig() const { return _gmo; }

  //void clearRegions() { _rg=regionGraph(); _rorder.clear(); }  //!!!TODO:DEBUG:UNSAFE
  //!!!TODO:DEBUG 
	void clearRegions() { _rg=regionGraph( _gmo.factors() ); _rorder.clear(); }  //!!!DEBUG:rorder
		
	struct myPair { findex r; size_t o; bool operator<(const myPair& b) const { return o < b.o; } }; //!!! TODO:DEBUG

  void addRegions( const vector<VarSet>& R ) { 
		//std::cout<<"AddRegions: called with "<<R.size()<<" regions\n";
		for (size_t r=0;r<R.size();++r) _rg.addRegion(R[r]);

		flist max = _rg.maximalRegions();
		size_t nMax = _rg.maximalRegions().size(), rMax=R.size();
		_rorder.clear(); _rorder.resize( 2*nMax ); 
		size_t r2=max.size();
		for (size_t r=0;r<R.size();++r) {  // find last-added maximal containing clique
			findex last = *((_rg.contains(R[rMax-r-1])&_rg.maximalRegions()).rbegin());  
			if (max.contains(last)) { _rorder[--r2]=last; max -= last; }
		}
		size_t M=max.size();
		for (size_t r=0;r<M;++r) { _rorder[--r2]=max[M-r-1]; } // add any leftovers?
		// Write them back in reverse order
		for (size_t r=0;r<nMax;++r) _rorder[2*nMax-1-r] = _rorder[r];

		//std::cout<<"Order:";
		//for (size_t r=0;r<_rorder.size();++r) std::cout<<" "<<_rorder[r];
		//std::cout<<"\n";
		//throw std::runtime_error("STOP");
	}


  double memory() const { 
    double mem=0.0;
		vector<double> memR(_rg.nregions(),0.0);
		for (size_t r=0;r<_rg.nregions();++r) memR[r] = _rg.region(r).nrStates()*sizeof(double)/1024.0/1024.0;
		for (size_t r=0;r<_rg.nregions();++r) mem+=memR[r];
		double memDampMax = 0.0;
    flist tmp = _rg.maximalRegions();
		for (vector<findex>::const_iterator it=tmp.begin();it!=tmp.end();++it) {
			double memDamp = 0.0;
      flist des = _rg.containedBy( _rg.region(*it) );
			for (vector<findex>::const_iterator d=des.begin();d!=des.end();++d) memDamp+=memR[*d];
			memDampMax = std::max(memDampMax,memDamp);
		}
		mem+=memDampMax;   //!!! TODO: fail to include memory for damping if damping off / exits...
		//if (_damp > 0) mem+=memDampMax;   //!!! TODO: fail to include memory for damping if damping off / exits...
/*
    double mem=0.0;
    for (size_t r=0;r<_rg.nregions();++r) 
      if (!_minimal || _rg.count(r)!=0) mem += _rg.region(r).nrStates()*sizeof(double)/1024.0/1024.0;
*/
    return mem;
  }

  void setFactors( const vector<Factor>& F ) { _gmo = graphModel(F); }

  const Factor& factor( findex r ) const { return _factors[r]; }
  
	Factor& __factor( findex r ) { return _factors[r]; }  // !!! non-safe!
 
  const regionGraph& getRG() const { return _rg; }
 
  findex regionWith( const VarSet& vs ) const {
    const flist& contains = _rg.contains(vs);
    assert( contains.size() > 0 );
    return _rg.sort(contains)[0];
  }

  Factor computeRegionBelief( findex r ) const {
    size_t nD = _descend[r].size();
#ifdef GBP_USE_LOG
    Factor bel = _factors[r];
#else
    Factor bel = _factors[r];
    double Zp = bel.sum(); bel/=Zp; 
#endif
    for (size_t i=0; i<nD-1; ++i) {
#ifdef GBP_USE_LOG
      bel += _factors[ _descend[r][i] ];
#else
      bel *= _factors[ _descend[r][i] ];
      double Zp = bel.sum(); bel/=Zp; 
#endif
    }
#ifdef GBP_USE_LOG
    double Zp = bel.logsumexp(); bel-=Zp; bel.exp();
#endif
    return bel;
  }

	void fixMaximal() { _rg.fixMaximal(); }
	void dumpRG() { _rg.dump(); }


/*
  lbp()                                 : factorGraph() { setProperties(); }
  lbp(const factorGraph& fg)            : factorGraph(fg) { setProperties(); }
  lbp(vector<Factor> fs)                : factorGraph(fs) { setProperties(); }
  template <class InputIterator>
  lbp(InputIterator f, InputIterator l) : factorGraph(f,l) { setProperties(); }

  virtual lbp* clone() const            { lbp* fg = new lbp(*this); return fg; }
*/

#ifdef MEX  
  // MEX Class Wrapper Functions //////////////////////////////////////////////////////////
  //void        mxInit();            // initialize any mex-related data structures, etc (private)
  //inline bool mxAvail() const;     // check for available matlab object
  bool        mxCheckValid(const mxArray*);   // check if matlab object is compatible with this object
  void        mxSet(mxArray*);     // associate with A by reference to data
  mxArray*    mxGet();             // get a pointer to the matlab object wrapper (creating if required)
  void        mxRelease();         // disassociate with a matlab object wrapper, if we have one
  void        mxDestroy();         // disassociate and delete matlab object
  void        mxSwap(lbp& gm);     // exchange with another object and matlab identity
  /////////////////////////////////////////////////////////////////////////////////////////
#endif


  Factor& belief(size_t f)              { throw std::runtime_error("Not implemented"); } //!!! const
  const Factor& belief(size_t f)  const { throw std::runtime_error("Not implemented"); }
  const Factor& belief(Var v)     const { //throw std::runtime_error("Not implemented"); }
		//if (iter() < 1) _belief[v] = _gbp.computeRegionBelief(_gbp.regionWith(v)).marginal(v);
		return _beliefs[v];
	}
  const Factor& belief(VarSet vs) const { throw std::runtime_error("Not implemented"); }
  //const Factor belief(Var v)     const { return computeRegionBelief(regionWith(v)).marginal(v);   }
  //const Factor belief(VarSet vs) const { return computeRegionBelief(regionWith(vs)).marginal(vs); }
  const vector<Factor>& beliefs() const { return _beliefs; } //throw std::runtime_error("Not implemented"); }

  // Only a bound-producing algorithm for max-elim 
   double lb() const { if (_task==Task::Max) return _lb; else throw std::runtime_error("Not available"); }
   double ub() const { if (_task==Task::Max) return _lnZ; else throw std::runtime_error("Not available"); }
  vector<index> best() const { if (_task==Task::Max) return _best; else throw std::runtime_error("Not available"); }
	void setBest(vector<index> best) { _best=best; _lb=_gmo.logP(_best); }

  // Gives an estimate of the partition function, but not a bound
  double logZ()   const { return _lnZ; }
  double logZub() const { throw std::runtime_error("Not available"); }
  double logZlb() const { throw std::runtime_error("Not available"); }
        
  MEX_ENUM( Task   , Max,Sum,Proximal,Mixed);
	Task _task;
	void setTask(Task m) { _task = m; } // !!!

  MEX_ENUM( ElimType   , Max,Sum);
  vector<ElimType> _elimType;							// !!! should be protected...
  void setElimType(ElimType e)            { for (int i=0;i<_elimType.size();++i) _elimType[i] = e; }
  void setElimType(Var v, ElimType e)     { _elimType[_vindex(v)] = e; }
  void setElimType(VarSet vs, ElimType e) { for (VarSet::const_iterator v=vs.begin();v!=vs.end();++v) _elimType[_vindex(*v)] = e; }


  MEX_ENUM( Property , Schedule,Distance,StopIter,StopObj,StopMsg,StopTime,Verbose,Damping,Minimal,Task );
  virtual void setProperties(std::string opt=std::string()) {
    if (opt.length()==0) {
      setProperties("Schedule=Fixed,Distance=HPM,StopIter=10,StopObj=-1,StopMsg=-1,StopTime=-1,Verbose=1,Damping=0.0,Task=Sum");
      setMinimal(false);
      return;
    }
    std::vector<std::string> strs = mex::split(opt,',');
    for (size_t i=0;i<strs.size();++i) {
      std::vector<std::string> asgn = mex::split(strs[i],'=');
      switch( Property(asgn[0].c_str()) ) {
        case Property::Schedule:  _sched  = Schedule(asgn[1].c_str()); break;
        case Property::Distance:  _dist   = Factor::Distance(asgn[1].c_str()); break;
        case Property::StopIter:  setStopIter( strtod(asgn[1].c_str(),NULL) ); break;
        case Property::StopObj:   setStopObj( strtod(asgn[1].c_str(),NULL) ); break;
        case Property::StopMsg:   setStopMsg( strtod(asgn[1].c_str(),NULL) ); break;
        case Property::StopTime:  setStopTime( strtod(asgn[1].c_str(),NULL) ); break;
        case Property::Verbose:   setVerbose( atol(asgn[1].c_str()) ); break;
        case Property::Damping:   setDamping( strtod(asgn[1].c_str(),NULL) ); break;
        case Property::Minimal:   setMinimal( atol(asgn[1].c_str()) ); break;
        case Property::Task:      setTask( Task(asgn[1].c_str()) ); break;
        default: break;
      }
    }
  }

  void setStopIter(double d) { _stopIter = d; }                   // stop when d*(# regions) updates done
  void setStopObj(double d)  { _stopObj = d;  }                   // stop when objective change is less than d
  void setStopMsg(double d)  { _stopMsg = d;  }                   // stop when message updates sare less than d
  void setStopTime(double d) { _stopTime = d;  }                  // stop when message updates sare less than d
  void setVerbose(bool d)    { _verbose = d; }
  void setDamping(double d)    { _damp = std::max(std::min(d,1.0),-1.0); } // negative values: never damp
	double getDamping()        { return _damp; }
  void setMinimal(bool d)    { _minimal = d; }

  /// Initialize the data structures
  virtual void init(const VarSet& vs) { init(); } // !!! inefficient


  void init() {
    _iter = 0;
    _dObjCurr = _dObj = infty();

    _factors.clear();    _descend.clear();    _ances.clear();    _descSet.clear();
		_beliefs.clear(); _beliefs.resize(nvar()); 	// !!! added saving belief support
		for (size_t v=0;v<nvar();++v) _beliefs[v]=Factor(var(v),1.0/var(v).states());

    if (_minimal) {
      _rg.cleanUp();      // !!!! move back into minimal
      //for (size_t f=0;f<_gmo.nFactors();++f) _rg.addRegion(_gmo.factor(f).vars()); // !!!! re-add minimal regions?
      /// !!! should change clean-up process to preserve these instead (if desired)
    }

    size_t nR = _rg.nregions();
    _factors.resize(nR); _descend.resize(nR); _ances.resize(nR); _descSet.resize(nR);
#ifdef SPARSE_FACTORS
    //for (size_t r=0;r<nR;++r) { Factor tmp(_rg.region(r),1.0); _factors[r].swap( tmp.log() ); }
    for (size_t r=0;r<nR;++r) { Factor tmp(1.0); _factors[r].swap( tmp.log() ); }
#else
#ifdef GBP_USE_LOG
    //_factors.resize(nR,Factor(0.0)); _descend.resize(nR); _ances.resize(nR); _descSet.resize(nR);
    for (size_t r=0;r<nR;++r) { Factor tmp(_rg.region(r),0.0); _factors[r].swap( tmp ); }
#else
    //_factors.resize(nR); _descend.resize(nR); _ances.resize(nR); _descSet.resize(nR);
    for (size_t r=0;r<nR;++r) { Factor tmp(_rg.region(r),1.0); _factors[r].swap( tmp ); }
#endif
#endif 

    // Associate factors in original model with some region
    flist added;
    _lnZ = 0.0;
    for (size_t f=0;f<_gmo.nFactors();++f) {
      if (_gmo.factor(f).nvar() == 0) _lnZ += std::log( _gmo.factor(f)[0] );
      else {
        vector<findex> tmp = _rg.sort( _rg.contains( _gmo.factor(f).vars() ) );
        assert( tmp.size() > 0 );
#ifdef GBP_USE_LOG
        _factors[ tmp[0] ] += log(_gmo.factor(f));
#else
        _factors[ tmp[0] ] *= _gmo.factor(f);
        double Z = _factors[tmp[0]].sum(); _factors[tmp[0]]/=Z; _lnZ+=std::log(Z);
#endif
        added += tmp[0];
      }
    }
		for (size_t v=0;v<nvar();++v) _withVariable(var(v)) = flist(); 
		for (size_t f=0;f<_factors.size();++f) {
			const VarSet& v = _factors[f].vars();
			for (VarSet::const_iterator i=v.begin();i!=v.end();++i) _withVariable(*i) |= f;
		}

		if (_task == Task::Max) {
			_best.resize( _gmo.nvar() , 0 ); 
			_xhat.resize( _gmo.nvar() , 0 ); 
			_lb = _gmo.logP(_best);
		}

    // Save descendant info from each maximal set, and ancestor subsets for use in updates
    for (size_t r=0;r<nR;++r) {  // !!! maximal only?
      flist des = _rg.containedBy( _rg.region(r) );
      _descSet[r] = des;
      _descend[r] = _rg.sort(des);
      _ances[r].resize( _descend[r].size() );
      for (size_t i=0;i<_descend[r].size();++i) {
        _ances[r][i] = _rg.sort( des & _rg.contains( _rg.region(_descend[r][i]) ) / _descend[r][i] );
      }
    }

    //entropy.resize(nRegions(),0.0);
    // Pre-execute one iteration on just the original factors (gets info from reparameterized loopy bp)
		double dampSave = _damp; setDamping(-1.0);
    for (size_t r=0;r<added.size();++r) update(added[r]); 
		setDamping(dampSave);

    if (_sched==Schedule::Priority) {                           // priority schedule =>
      flist tmp = _rg.maximalRegions();
      for (flist::const_iterator i=tmp.begin();i!=tmp.end();++i) 
        priority.insert(std::numeric_limits<double>::infinity(),*i);
    }

/* TODO:DEBUG: TRY NEW ORDERING
    // Order arbitrarily among maximal sets
    _rorder.clear(); 
    flist tmp = _rg.maximalRegions();
    // "Forward-backward" sweep (if cliques form a partial-ordered junction tree)
    _rorder.resize( 2*tmp.size() );
    vector<findex>::iterator next = std::copy( tmp.begin(), tmp.end(), _rorder.begin() );
    std::copy( tmp.rbegin(), tmp.rend(), next );   // duplicates some effort but aligned for convergence test
    // "Forward-backward" sweep minus endpoint
    //_rorder.resize( 2*tmp.size() -1 );
    //vector<findex>::iterator next = std::copy( tmp.begin(), tmp.end(), _rorder.begin() );
    //std::copy( tmp.rbegin(), tmp.rend(), --next );   // duplicates some effort but aligned for convergence test
    // "Random order"
    //_rorder.resize( tmp.size() );
    //std::copy( tmp.begin(), tmp.end(), _rorder.begin() );
    //for (size_t i=_rorder.size()-1;i>0;--i) std::swap( _rorder[i],_rorder[mex::randi(i+1)] );
*/

/*
std::cout<<"\n\n";
for (size_t r=0;r<nRegions();++r) std::cout<<r<<": "<<_rg.region(r)<<"\n";
for (size_t r=0;r<nRegions();++r) { std::cout<<r<<"=>"; for (size_t i=0;i<_descend[r].size();++i) std::cout<<_descend[r][i]<<", "; std::cout<<"\n"; }
for (size_t r=0;r<tmp.size();++r) std::cout<<tmp[r]<<" "; std::cout<<"\n";
throw std::runtime_error("QUIT");
*/
    //std::ofstream os("debug.txt");
    //for (size_t r=0;r<nR;++r) { os<<r<<" : "<<_rg.region(r)<<"\n"; } os<<"\n\n";
    //for (size_t r=0;r<nR;++r) { os<<r<<" "; for (double d=std::log10(0.9+r)-1;d<12;++d) os<<" "; }
    //os.close();

    _dObjStable = infty();
    _lnZStable = 0.0;
    _lnZStable = _lnZ;    // !!!! from "initial" updates
		if (_verbose) std::cout<<"GBP init(): "<<_lnZ<<"\n";
#ifdef DEBUG  
    std::ofstream ofs; ofs.open("updates.gbp.txt", std::ofstream::out | std::ofstream::app);
    ofs<<"END GBP INIT: "<<_lnZ<<"\n\n";
    ofs.close();
#endif

  }

  double iter() { return ((double) _iter)/_rg.maximalRegions().size(); }
  bool iter_boundary() { return (_iter % _rg.maximalRegions().size()) == 0; }
  double dObj() { return std::abs(_dObj); }

  double _dObjStable, _lnZStable;
  double logZStable() const { return _lnZStable; }

  /// Run the algorithm until one of the stopping criteria is reached
  virtual void run() {
    double startTime = timeSystem();
    double print=startTime+1;
    size_t stopIter = _stopIter * _rg.maximalRegions().size();
    size_t n=0; if (_rorder.size()) n=_iter % _rorder.size();
    double dMsg=infty();                            // initialize termination values

		/// While not done with GBP: ////////////////////////////////////////////////////////////////
    while ( (dMsg>=_stopMsg) && (_iter<stopIter) && (dObj()>=_stopObj) && (_dObjStable*(1-_damp))>_stopObj )  {
      if (_stopTime > 0 && _stopTime <= (timeSystem()-startTime)) break;       // time-out check
      if (iter_boundary()) _dObjCurr = 0.0;         // reset dObj each N iterations

      size_t r;
      if (_sched==Schedule::Priority) {                           // priority schedule =>
        //throw std::runtime_error("Priority schedule not implemented");
        r=priority.top().second; priority.pop();
      } else {                                                    // fixed schedule =>
        r=_rorder[n];                                             //   get next factor from list
        if (++n == _rorder.size()) n=0; 
      }

      double delta = _lnZ;            // keep track of differential in lnZ
      update(r);              				// update this region
      delta -= _lnZ; 
      _dObjCurr += std::abs(delta); 	// track total differential per iteration

      if (_sched==Schedule::Priority) {
        flist nbrs = _rg.maximalRegions() & _rg.intersects(_rg.region(r));
        double absdelta = std::abs(delta);
        for (flist::const_iterator i=nbrs.begin();i!=nbrs.end();++i) {
          double oldPriority = 0.0; if (priority.isIndexed(*i)) oldPriority = priority.byIndex(*i);
          priority.insert(oldPriority+absdelta,*i);
        }
      }

			//// Verbose outputs: ////
      if (_verbose && timeSystem()>print) { 
        print=timeSystem()+1; 
        std::cout<<"iter "<<iter()<<"; lnZ: "<<_lnZ<<" ("<<_dObj<<")\n"; 
      }

      _iter++;
			if (iter()<=1.0) _lnZStable = _lnZ;				// during first iteration, just take the current lnZ!

			//// Iteration boundary: ////
      if (iter_boundary()) { 
				if (_task == Task::Max) {
					vector<uint32_t> config = maxSequential(order(mex::graphModel::OrderMethod::Random));
					if (_gmo.logP(config) > _lb) { 
						_lb = _gmo.logP(config); 
						for (size_t v=0;v<config.size();++v) _best[v] = config[v]; 
					}
					std::cout<<"Best config found has value "<<_lb<<"\n";
				}
				if (_task == Task::Proximal) {
					for (size_t v=0;v<nvar();++v) if (_elimType[v]==ElimType::Max) { } //!!!
						//_factor[???v???] += log-belief[v];
				}
				double _dObjLast = _dObj;
        _dObj = _dObjCurr;                   // keep track of lagged objective
        if (dObj() <= _dObjStable || iter()<=2) {
					for (size_t v=0;v<nvar();++v) _beliefs[v] = computeRegionBelief(regionWith(var(v))).marginal(var(v));
				}
        if (dObj() <= _dObjStable || iter()<=1) {
        //if (dObj() <= _dObjStable * (1-_damp) || iter()<=1) { 
          _lnZStable = _lnZ;
          _dObjStable = dObj();
        }
				else if (_dObjLast < _dObjCurr && _damp>=0.0) {
					_iter = infty();		/// !!! TODO:DEBUG:FIX Exit without error if damping starts...
					setDamping( 1.0 - (_damp-1.0)/(_damp-2.0) );
					if (_verbose) std::cout << "  Shuffling order & damping: "<<_damp<<"... \n";
					_rorder.resize(_rg.maximalRegions().size());						// Reduce to just single set of max cliques
					std::copy(_rg.maximalRegions().begin(), _rg.maximalRegions().end(), _rorder.begin());
					std::random_shuffle( _rorder.begin(), _rorder.end() ); 	// Random shuffle update to avoid pathology
					n=0;  																									// and start back at the beginning of the order
				}
				//if (_damp>=0) _damp=1.0-1.0/(1.0+iter());
/* // Test: try shuffling update order proactively?
				if (iter()>=2.0 && _damp>=0.0) {
					if (_verbose) std::cout << "  Shuffling order...\n";
					_rorder.resize(_rg.maximalRegions().size());						// Reduce to just single set of max cliques
					std::copy(_rg.maximalRegions().begin(), _rg.maximalRegions().end(), _rorder.begin());
					std::random_shuffle( _rorder.begin(), _rorder.end() ); 	// Random shuffle update to avoid pathology
					n=0;  																									// and start back at the beginning of the order
				}
*/

				std::cout<<"end iter "<<iter()-1<<"; "<<_lnZ<<" ("<<dObj()<<") ; "<<_lnZStable<<" ("<<_dObjStable<<")\n";
#ifdef DEBUG  
		    std::ofstream ofs; ofs.open("updates.gbp.txt", std::ofstream::out | std::ofstream::app);
   			ofs<<"\n\n=== END ITER "<<iter()-1<<"; "<<_lnZ<<" ("<<dObj()<<") ; "<<_lnZStable<<" ("<<_dObjStable<<")\n\n";
    		ofs.close();
#endif
      }
      
      if (_lnZ == -infty()) { _dObj=0.0; break; }
    }
    if (_verbose) std::cout<<"GBP: "<<iter()<<"it, "<<timeSystem()-startTime<<"sec\n";

    //os.close();
  }



  regionGraph cce;
  flist cceRBase, cceRMax;
  size_t icce,jcce;   // j = index into rg-regions; i = index into cce-regions

  vector<double> cceLnZ;
  vector<double> countOld;
  
  void initCCE(double MemLimit) {
    cce = _rg;
    cceLnZ.clear(); cceLnZ.resize(cce.nregions(),0.0);
    icce=0; jcce=0;
    cceRMax  = cce.maximalRegions();   // if regions are removed, need to copy these off to avoid corruption
    cceRBase = _rg.maximalRegions() & _rg.intersects(cce.region(cceRMax[icce]));
    countOld.clear(); countOld.resize(cce.nregions(),0);
    for (size_t r=0;r<cce.nregions();++r) {
      double mem = cce.region(r).nrStates() * sizeof(double) / 1024 / 1024;  // required memory in MB
      if (mem <= MemLimit)
        try{
          cceLnZ[r] = cceTermMem( cce.region(r) );  // otherwise allocate memory for more speed
        } catch (std::exception& e) {
          MemLimit *= 0.9;
          cceLnZ[r] = cceTerm( cce.region(r) );     // if we can't afford it, call the slower version
        }
      else 
        cceLnZ[r] = cceTerm( cce.region(r) );     // if we can't afford it, call the slower version
    }
  }

  size_t getCCENR()   const { return cce.nregions(); }
  void cceStats() const {
    vector<size_t> cnt(40);
    for (size_t r=0;r<cce.nregions();++r) ++cnt[cce.region(r).nvar()];
    for (size_t r=0;r<27;++r) std::cout<<cnt[r]<<" "; std::cout<<"\n";
  }

  void cceDump() const {
    for (size_t r=0;r<cce.nregions();++r) {
      std::cout<<cce.region(r)<<" "<<cce.count(r)<<" "<<countOld[r]<<" "<<cceLnZ[r]/std::log(10)<<"\n";
    }
  }

/*
  void cceFillQueue() {
    cceQueue.clear();
    const flist& A = cce.maximalRegions();   // if regions are removed, need to copy these off to avoid corruption
    for (size_t i=0;i<A.size();++i) {
      flist B = _rg.maximalRegions() & _rg.intersects(cce.region(A[i]));
      for (size_t j=0;j<B.size();++j) {
        flist C = _rg.maximalRegions() & _rg.intersects( cce.region(A[i])-_rg.region(B[j]) )
                                        & _rg.intersects( _rg.region(B[j])-cce.region(A[i]) );
        for (size_t k=0;k<C.size();++k) {
          cceQueue.push_back(cce.region(A[i])+_rg.region(B[j])+_rg.region(C[k]));
        }
      }
    }
    std::cout<<"CCE: "<<cceQueue.size()<<" cliques queued\n";
  }
*/

  void cceUpdate() {
    countOld.resize(cce.nregions());
    for (size_t d=0;d<cce.nregions();++d) if (countOld[d]!=cce.count(d)) {
      _lnZ += (cce.count(d)-countOld[d])*cceLnZ[d];  // update lnZ using the change in counts
      countOld[d]=cce.count(d);                      // (and update counts)
    }
  }
  void cceIncrement() {
      // Increment regions to combine
      if (++jcce >= cceRBase.size()) { 
        jcce=0; 
        if (++icce >= cceRMax.size()) {
          cceRMax = cce.maximalRegions();   // get new set of maximal regions
          icce=0;
        }
        cceRBase = _rg.maximalRegions() & _rg.intersects(cce.region(cceRMax[icce]));
      }
  }

  size_t stepCCE(const VarSet& newRegion, double MemLimit) {
      flist added = cce.addRegion( newRegion );        // add new region to cce and track all added regions
      cceLnZ.resize(cce.nregions());                   // update lengths of lnZ, counts
      for (size_t a=0;a<added.size();++a) {
        //std::cout<<"    "<<cce.region(added[a])<<"\n";
        double mem = cce.region(added[a]).nrStates() * sizeof(double) / 1024 / 1024;  // required memory in MB
        if (mem > MemLimit)
          cceLnZ[added[a]] = cceTerm( cce.region(added[a]) );     // if we can't afford it, call the slower version
        else
          cceLnZ[added[a]] = cceTermMem( cce.region(added[a]) );  // otherwise allocate memory for more speed
      }
      cceUpdate();
      size_t redundant=0;
      for (size_t a=0;a<added.size();++a) 
        for (size_t r=0;r<_rg.maximalRegions().size();++r)
          if (_rg.region(_rg.maximalRegions()[r]) >> cce.region(added[a])) {
            ++redundant;
            break;
          }
      return added.size()-redundant;
  }

  VarSet nextCCE(void) {
    VarSet newRegion = _rg.region(cceRBase[jcce]) + cce.region(cceRMax[icce]);
    while (newRegion.size() == cce.region(cceRMax[icce]).size()) {   // simple check for uselessness
//std::cout<<icce<<","<<jcce<<" "<<newRegion<<"\n";
      cceIncrement();
      newRegion = _rg.region(cceRBase[jcce]) + cce.region(cceRMax[icce]);
    }
//std::cout<<icce<<","<<jcce<<" "<<newRegion<<"\n";
    cceIncrement();
    return newRegion;
  }

  void runCCE(double runTime, double MemLimit, bool verbose) { 
    double startTime = timeSystem();
    double print = timeSystem() + 1; 
    while ( (startTime+runTime) > timeSystem()) {
      // Get prospective new region and compute its contribution
      VarSet newRegion = _rg.region(cceRBase[jcce]) + cce.region(cceRMax[icce]);
if (newRegion.size() != cce.region(cceRMax[icce]).size()) {   // simple check for uselessness
      flist added = cce.addRegion( newRegion );        // add new region to cce and track all added regions
      cceLnZ.resize(cce.nregions());                   // update lengths of lnZ, counts
      for (size_t a=0;a<added.size();++a) {
        double mem = cce.region(added[a]).nrStates() * sizeof(double) / 1024 / 1024;  // required memory in MB
        if (mem > MemLimit)
          cceLnZ[added[a]] = cceTerm( cce.region(added[a]) );     // if we can't afford it, call the slower version
        else 
          cceLnZ[added[a]] = cceTermMem( cce.region(added[a]) );  // otherwise allocate memory for more speed
      }

      if (verbose && timeSystem()>print) {               // Update values and print progression
        cceUpdate();
        std::cout<<"CCE ("<<timeSystem()-startTime<<") "<<_lnZ<<"\n";
        cceStats();
        print = timeSystem()+1;
      }

}
      cceIncrement();
    } 
    cceUpdate();
    if (verbose) std::cout<<"Done CCE ("<<timeSystem()-startTime<<") "<<_lnZ<<"\n";
  }


/*
  void runCCE(double runTime, double MemLimit, bool verbose) { 
    double startTime = timeSystem();
    double print = timeSystem() + 1; 
    while ( (startTime+runTime) > timeSystem()) {
      if (iQueue==cceQueue.size()) { cceFillQueue(); iQueue=0; }
      if (cceQueue.size()==0) return;
      const VarSet& newRegion = cceQueue[iQueue];

      // Get prospective new region and compute its contribution
      flist added = cce.addRegion( newRegion );        // add new region to cce and track all added regions
      cceLnZ.resize(cce.nregions());                   // update lengths of lnZ, counts
      for (size_t a=0;a<added.size();++a) {
        double mem = cce.region(added[a]).nrStates() * sizeof(double) / 1024 / 1024;  // required memory in MB
        if (mem > MemLimit)
          cceLnZ[added[a]] = cceTerm( cce.region(added[a]) );     // if we can't afford it, call the slower version
        else 
          cceLnZ[added[a]] = cceTermMem( cce.region(added[a]) );  // otherwise allocate memory for more speed
      }

      if (verbose && timeSystem()>print) {               // Update values and print progression
        cceUpdate();
        std::cout<<"CCE ("<<timeSystem()-startTime<<") "<<_lnZ<<"\n";
        cceStats();
        print = timeSystem()+1;
      }
      ++iQueue;
    } 
    cceUpdate();
    if (verbose) std::cout<<"Done CCE ("<<timeSystem()-startTime<<") "<<_lnZ<<"\n";
  }
*/

  double cceTerm( const VarSet& vs ) {
    flist regions = _rg.containedBy(vs);
    std::vector<mex::subindex> idx; idx.reserve(regions.size());
    for (size_t i=0;i<regions.size();++i) idx.push_back( subindex(vs,_rg.region(regions[i])) );
    size_t nrStates = vs.nrStates();
    double lnZterm=-infty();
//std::cout<<"."; std::cout.flush();
    for (size_t n=0;n<nrStates;++n) {
//if (n%100000 == 0) { std::cout<<"."; std::cout.flush(); }
      double v=0.0;
      for (size_t i=0;i<regions.size();++i) { v+=_factors[regions[i]][idx[i]]; ++idx[i]; }
      if (lnZterm==-infty()) lnZterm  = v;
      else lnZterm += log1p(std::exp(v-lnZterm));
    }
    return lnZterm;
  }

  double cceTermMem( const VarSet& vs ) {
    //double mem = ((double)vs.nrStates())*sizeof(double)/1024.0/1024.0;
    //if (mem > MemLimit) return 0.0;

    flist regions = _rg.containedBy(vs);
    Factor F(vs,0.0);
    for (size_t i=0;i<regions.size();++i) {
      subindex idx(vs,_rg.region(regions[i]));
      for (size_t n=0;n<F.nrStates();++n,++idx) F[n]+=_factors[regions[i]][idx];
    }
    return F.logsumexp();
  }

  void cceTestMem(double MemLimit) {
    double time = timeSystem();
    for (size_t r=0;r<_rg.nregions();++r) {
      std::cout<<r<<" ("<<_rg.region(r).nrStates()<<") : "; std::cout.flush();
      if (MemLimit < ((double)_rg.region(r).nrStates())*sizeof(double)/1024/1024) std::cout<<"doesn't fit\n";
      else
        std::cout<<cceTermMem(_rg.region(r))<<"\n";
    }
    std::cout<<timeSystem()-time<<"\n";
  }

  void cceTest() {
    double time = timeSystem();
    for (size_t r=0;r<_rg.nregions();++r) {
      std::cout<<r<<" ("<<_rg.region(r).nrStates()<<") : "; std::cout.flush();
      std::cout<<cceTerm(_rg.region(r))<<"\n";
    }
    std::cout<<timeSystem()-time<<"\n";
  }


protected:	// Contained objects
	//vector<Factor>   _beliefs;                       // store calculated messages and beliefs
	//vector<Factor>   _msg;
	//vector<Factor>   _msgNew;

  void  logsumexpInto(Factor& src, Factor& target) {
		VarSet elim = src.vars() - target.vars();
		for (size_t t=0;t<target.nrStates();++t) {
			double mx = -infty(), val=0.0;
			for (superindex s(src.vars(),elim,t);s!=s.end();++s) mx = std::max(mx,src[s]);
			if (mx != -infty())
				for (superindex s(src.vars(),elim,t);s!=s.end();++s) val += std::exp(src[s]-mx);
			target[t] = mx + std::log(val);
		}
  }
  void  marginalInto(Factor& src, Factor& target) {
    target.fill(0.0);
    subindex s(src.vars(),target.vars());
    for (size_t i=0;i<src.nrStates(); ++i,++s) target[s]+=src[i];
  }
  void  maxmarginalInto(Factor& src, Factor& target) {
    target.fill(-infty());
    subindex s(src.vars(),target.vars());
    for (size_t i=0;i<src.nrStates(); ++i,++s) target[s]=(target[s]<src[i]) ? src[i] : target[s];
  }
  void  argmaxmarginalInto(Factor& src, Factor& target, const vector<uint32_t>& xhat) {
		// notation: (xA,xB) max vars; (xC,xD) sum vars; target is only over xB,xD (so eliminate xA,xC)
		// xhat is a configuration of (at least) xA chosen to maximize this factor's local marginal MAP
		VarSet xBD = target.vars(), xAC = src.vars()-xBD;
		VarSet xA,xB,xC,xD;
		for (VarSet::const_iterator v=xBD.begin();v!=xBD.end();++v) if (_elimType[*v]==ElimType::Max) xB|=*v;
		xD = xBD - xB;
		for (VarSet::const_iterator v=xAC.begin();v!=xAC.end();++v) if (_elimType[*v]==ElimType::Max) xA|=*v;
		xC = xAC - xA;
		if (xC.size()+xD.size()==0) maxmarginalInto(src,target);		// no sum variables at all => max-marginal
		else if (xA.size()==0)      marginalInto(src,target);				// no elim max variables => marginal
		else {
			superindex s(src.vars(), xB+xC+xD, sub2ind(xA, xhat));		// walk over src F(xA*,xBCD)
			subindex   t(xB+xC+xD,xB+xC);															// and simultaneously over target f(xBC)
			target.fill(0.0);
			for (;s!=s.end();++s,++t) target[t] += src[s];						// Marginalize this slice of src into target
		}
  }

	// find marginal MAP value of this factor
	double argmaxMMAP(Factor& F, vector<uint32_t>& xhat) {
		// Memory-inefficient way:
		VarSet xA;
		for (VarSet::const_iterator v=F.vars().begin();v!=F.vars().end();++v) if (_elimType[*v]==ElimType::Max) xA|=*v;
		Factor f = F.marginal(xA);
		size_t amax = f.argmax();
		ind2sub(xA, amax, xhat);
		return f[amax];

/*
		// Memory efficient way (!!! TODO)
		size_t bestx = 0;
		double bestv = -infty();
		for (size_t i=0;i<xA.nrStates();++i) {
			superindex j(F.vars(),xA,i);
			double val = 0.0;
			for (;j!=j.end();++j) val+=F[*j];
			if (val > bestv) { bestv=val; bestx = i; }
		}
		xA.ind2sub( bestx, xhat );
		return bestv;
*/
	}


/*
vector<double> entropy;

double recalcLogZ() {
  double lnZ=0.0;
  for (size_t f=0;f<_gmo.nFactors();++f) {
    if (_gmo.factor(f).nvar() == 0) lnZ += std::log( _gmo.factor(f)[0] );
    else {
      vector<findex> tmp = _rg.sort( _rg.contains( _gmo.factor(f).vars() ) );
      assert( tmp.size() > 0 );
#ifdef GBP_USE_LOG
      lnZ += (computeRegionBelief(tmp[0])*log0(_gmo.factor(f))).sum();
#else
      lnZ += (computeRegionBelief(tmp[0])*log0(_gmo.factor(f))).sum();
#endif
    }
  }
  for (size_t r=0;r<nRegions();++r) lnZ += _rg.count(r)*entropy[r];
  return lnZ;
}
*/

  /// Perform a reparameterization update on some maximal clique r, renormalizing after
  void update(findex r) { 
    assert( *_descend[r].rbegin() == r );
    size_t nD = _descend[r].size();

    vector<Factor> fOld; // !!! damp
    if (_stopMsg > 0) { /* save all descendant's factors */ }
    if (_damp > 0) { fOld.resize(nD); for (size_t i=0;i<nD;++i) fOld[i]=_factors[_descend[r][i]]; } // !!! damp

    // Collect: push all beliefs up to maximal clique
#ifdef GBP_USE_LOG
#else
    for (size_t i=0; i<nD; ++i) _factors[_descend[r][i]].log();
#endif

    for (size_t i=0; i<nD-1; ++i) {
      findex n = _descend[r][i], p=_ances[r][i][0];
      _factors[ _ances[r][i][0] ] += _factors[ _descend[r][i] ];
    }

    // Update objective value 
		double logZ;
		if (_task==Task::Max) {
			//size_t mx = _factors[r].argmax();
			//ind2sub(_factors[r].vars(),mx,_xtmp);
			//_xhat=_xtmp; _lb=_gmo.logP(_xhat); 
			//if (_gmo.logP(_xtmp) > _lb) { _xhat=_xtmp; _lb=_gmo.logP(_xhat); }
			//logZ = _factors[r][mx];  // =_factors[r].max()
			logZ = _factors[r].max();  // TODO: CHANGE?
		} else {
			logZ = _factors[r].logsumexp();
		}
    //double logZ = _factors[r].logsumexp(); 
	
#ifdef DEBUG	
		std::ofstream ofs; ofs.open("updates.gbp.txt", std::ofstream::out | std::ofstream::app);
		ofs<<"  dZ="<<logZ<<" Region "<<r<<"="<<_factors[r].vars();
		if (std::abs(logZ)>.01) ofs<<" (!!!!) ";  // !!! DEBUG
		ofs<<"\n";
		ofs.close();
#endif

    _factors[r]-=logZ; 

    if (_damp > 0) _lnZ+=(1.0-_damp)*logZ; else _lnZ+=logZ; // !!! damp
    //_lnZ+=logZ;  // !!! damp (if no damping)
    //_factors[r].exp(); // !!! logsumexpInto instead
		// TODO: don't do most of this if maximizing?  don't care?
		// TODO: if MMAP, v=argmaxMMAP(_factors[r],xhat); _factors[r]-=v; lnZ+=v;
		// !!! TODO: if mixed, compute maximal region's argmax into vector/map structure
  
    //entropy[r] = _factors[r].entropy(); 
    // Eliminate: calculate beliefs at each (sub) region by elimination
    for (size_t i=2; i<nD+1; ++i) {   // !!! ugly for backward indexed loop
      findex n = _descend[r][nD-i], p=_ances[r][nD-i][0];
#ifdef SPARSE_FACTORS
			if (_task==Task::Max) _factors[n] = _factors[p].max(_factors[p].vars() - _factors[n].vars());
			else                  _factors[n] = _factors[p].logsumexp(_factors[p].vars() - _factors[n].vars());
#else
			if (_task==Task::Max) maxmarginalInto( _factors[p] , _factors[n]); // memory-efficient elim
			else                  logsumexpInto( _factors[p] , _factors[n]);  	
#endif
      //entropy[n] = _factors[n].entropy();										// if entropy calc used?
    }
#ifdef GBP_USE_LOG
    //for (size_t i=0; i<nD; ++i) _factors[_descend[r][i]].log();	// !!! logsumexpInto instead
#else
#endif

		if (_task==Task::Max) for (size_t i=0;i<nD;++i) {
			size_t n=_descend[r][i]; _factors[n]*=((double)_factors[n].nvar())/_factors[r].nvar();
		}


    // Reparameterize: pull out subsets' beliefs from their ancestors
    vector<flist> done(_descend[r].size());         // make a list for each descendant
    for (size_t ii=1; ii<nD+1; ++ii) {              // Down pass: try to divide by large beliefs
      size_t i=nD-ii; findex n = _descend[r][i];    //   but make sure we don't divide by any vars twice
      for (fvec::reverse_iterator d=++_descend[n].rbegin();d!=_descend[n].rend();++d) {
        if ( !done[i].intersects(_descSet[*d]) ) {  // if nothing below d has been divided away yet,
#ifdef GBP_USE_LOG
          _factors[n] -= _factors[*d];              //   divide by d's belief (below n, factors are beliefs)
#else
          _factors[n] /= _factors[*d];              //   divide by d's belief (below n, factors are beliefs)
#endif
          done[i] += _descSet[*d];                  //   and this will remove everything below d
        }
      }
    }
    for (size_t i=0; i<nD; ++i) {                   // Up pass: divide by residual factor of any cliques
      findex n=_descend[r][i];                      //   that we missed in the downward pass
      flist todo = _descSet[n] - done[i] - n;
#ifdef GBP_USE_LOG
      for (flist::const_iterator d=todo.begin();d!=todo.end();++d) _factors[n] -= _factors[*d];
#else
      for (flist::const_iterator d=todo.begin();d!=todo.end();++d) _factors[n] /= _factors[*d];
#endif
    }

    if (_damp > 0) {  // !!! damp
      for (size_t i=0;i<nD;++i) { size_t n=_descend[r][i];
        _factors[n] *= (1.0 - _damp); fOld[i] *= _damp; _factors[n] += fOld[i];
      }  
			fOld.clear();
    }

  }
    
protected:  // Contained objects

};


#ifdef MEX
//////////////////////////////////////////////////////////////////////////////////////////////
// MEX specific functions, and non-mex stubs for compatibility
//////////////////////////////////////////////////////////////////////////////////////////////
bool gbp::mxCheckValid(const mxArray* GM) { throw std::runtime_error("NOT IMPLEMENTED"); return false; }
void gbp::mxSet(mxArray* GM) { throw std::runtime_error("NOT IMPLEMENTED"); }
mxArray* gbp::mxGet() { throw std::runtime_error("NOT IMPLEMENTED"); return NULL; }
void gbp::mxRelease() {  throw std::runtime_error("NOT IMPLEMENTED"); }
void gbp::mxDestroy() { throw std::runtime_error("NOT IMPLEMENTED"); }
void gbp::mxSwap(lbp& gm) { throw std::runtime_error("NOT IMPLEMENTED"); }
#endif


//////////////////////////////////////////////////////////////////////////////////////////////
}       // namespace mex
#endif  // re-include
