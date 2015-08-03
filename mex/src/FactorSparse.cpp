/////////////////////////////////////////////////////////////////////////////////////
// Factor.h  --  class definition for matlab-compatible factor class
//
// A few functions are defined only for MEX calls (construction & load from matlab)
// Most others can be used more generally.
// 
//////////////////////////////////////////////////////////////////////////////////////

//#include <assert>
#include <cmath>
#include <cstring>
#include <iostream>
#include <limits>
#include <algorithm>
#include <numeric>
#include <float.h>
#include <vector>

#include "mxUtil.h"
#include "VarSet.h"
#include "subindex.h"
#include "vector.h"

#include "FactorSparse.h"

#if defined _WINDOWS || WINDOWS
#define strcasecmp(s1,s2) stricmp(s1,s2)
#endif 

namespace mex {


#ifdef __FACTOR_H_MEMORY
size_t Factor::memused = 0;
size_t Factor::mmax = 0;
#endif

/************************************************************************************
 ************************************************************************************
 IMPLEMENTATION
 ************************************************************************************
************************************************************************************/
  // Functors defined for binary operations on the factor table : Op(a,b) and Op.IP(a,b) (in-place version)
  struct binOpPlus   {
    Factor::value  operator()(Factor::value  a, const Factor::value b) { return a+b; };
    Factor::value&         IP(Factor::value& a, const Factor::value b) { return a+=b;};
  };
  struct binOpMinus  {
    //value  operator()(value  a, const value b) { return a-b; };
    //value&         IP(value& a, const value b) { return a-=b;};
    Factor::value  operator()(Factor::value  a, const Factor::value b) { return (b!=-infty()) ? a-b : b; };
    Factor::value&         IP(Factor::value& a, const Factor::value b) { return (b!=-infty()) ? a-=b: a=b; };
    //value  operator()(value  a, const value b) { return (isnan(a-=b)) ? a=-infty() : a; };
    //value&         IP(value& a, const value b) { return (isnan(a-=b)) ? a=-infty() : a; };
    //value  operator()(value  a, const value b) { value c=a-b; if (c!=c) c=a; return c; };
    //value&         IP(value& a, const value b) { value c=a-b; if (c==c) a=c; return a; };
  };
  struct binOpTimes  {
    Factor::value  operator()(Factor::value  a, const Factor::value b) { return a*b; };
    Factor::value&         IP(Factor::value& a, const Factor::value b) { return a*=b;};
  };
  struct binOpDivide {
    Factor::value  operator()(Factor::value  a, const Factor::value b) { return (b) ? a/b  : 0;  };
    Factor::value&         IP(Factor::value& a, const Factor::value b) { return (b) ? a/=b : a=0;};
  };
  struct binOpPower {
    Factor::value  operator()(Factor::value  a, const Factor::value b) { return std::pow(a,b);  };
    Factor::value&         IP(Factor::value& a, const Factor::value b) { return a=std::pow(a,b);};
  };

  Factor  Factor::operator+ (const Factor& B) const  { return binaryOp(  B, binOpPlus()  ); };
  Factor& Factor::operator+=(const Factor& B)        { return binaryOpIP(B, binOpPlus()  ); };
  Factor  Factor::operator- (const Factor& B) const  { return binaryOp(  B, binOpMinus() ); };
  Factor& Factor::operator-=(const Factor& B)        { return binaryOpIP(B, binOpMinus() ); };
  Factor  Factor::operator* (const Factor& B) const  { return binaryOp(  B, binOpTimes() ); };
  Factor& Factor::operator*=(const Factor& B)        { return binaryOpIP(B, binOpTimes() ); };
  Factor  Factor::operator/ (const Factor& B) const  { return binaryOp(  B, binOpDivide()); };
  Factor& Factor::operator/=(const Factor& B)        { return binaryOpIP(B, binOpDivide()); };

  Factor  Factor::operator+ (const value B) const    { return binaryOp(  B, binOpPlus()  ); };
  Factor& Factor::operator+=(const value B)          { return binaryOpIP(B, binOpPlus()  ); };
  Factor  Factor::operator- (const value B) const    { return binaryOp(  B, binOpMinus() ); };
  Factor& Factor::operator-=(const value B)          { return binaryOpIP(B, binOpMinus() ); };
  Factor  Factor::operator* (const value B) const    { return binaryOp(  B, binOpTimes() ); };
  Factor& Factor::operator*=(const value B)          { return binaryOpIP(B, binOpTimes() ); };
  Factor  Factor::operator/ (const value B) const    { return binaryOp(  B, binOpDivide()); };
  Factor& Factor::operator/=(const value B)          { return binaryOpIP(B, binOpDivide()); };
  Factor  Factor::operator^ (const value B) const    { return binaryOp(  B, binOpPower() ); };
  Factor& Factor::operator^=(const value B)          { return binaryOpIP(B, binOpPower() ); };


  // Helper:  r1v2 = unionVals( v1, r1, v2, r2, &v1v2)  ; r1^2=isectVals(v1,v1^2,r1)
//TODO:BIG:  
  // if v1 < v2 then  // (use pointers so can swap if reversed)
  //    t12 = t2; v1^2=v1&v2; // copy & use v2's map:
  //    for (iterate over tuples i2 in t12), 
  //      i1=isectVals(v2,v1^2,i2); 
  //      if (~in t1) remove from t12; else t12[i12]op=t1(i1);  // (hint form?)
  // else:
  //   v1^2 = v1&v2; 
  //   map find1(i12 -> vector(i1)), find2(i12 -> i2);
  //   for (i1 in t1), i12=isectVals(v1,v1^2,i1); find1(i12)+=i1;
  //   for (i2 in t2), i12=isectVals(v2,v1^2,i2); find2(i12)+=i2;
  //   for (i12 in find1),
  //     for i1 in find1(i12), for i2 in find2(i12), 
  //       i1v2 = unionVals(v1,v2,i1,i2,v1^2); 
  //       t12(i1v2) = t1(i1) OP t2(i2);   //(use hint form?)
  //

  size_t Factor::unionVals(const VarSet& v1, size_t r1, const VarSet& v2, size_t r2, const VarSet& vUnion) {
    assert(v1 << vUnion); assert(v2 << vUnion);
    if (vUnion.size()==0) return 0;
    std::vector<vsize> tuple(vUnion.last().label()+1);  // allocate storage for joint tuple
    ind2sub(v1,r1,tuple);
    ind2sub(v2,r2,tuple);
    return sub2ind(vUnion,tuple);
  }
  size_t Factor::isectVals(const VarSet& v1, const VarSet& v12, size_t r1) {
    assert(v12 << v1);
    if (v1.size()==0) return 0;
    std::vector<vsize> tuple(v1.last().label()+1);  // allocate storage for tuple
    ind2sub(v1,r1,tuple);
    return sub2ind(v12,tuple);
  }

  // Binary operations (eg A + B); returns new object
  template<typename Function> Factor Factor::binaryOp( const Factor& B, Function Op) const {
    VarSet vUnion = v_ + B.v_;              // expand scope to union
    // TODO: something else if v_ << B.v_ ...
    VarSet vIsect = v_ & B.v_;
    Factor F; F.v_ = vUnion; F.t_.clear(); // valid? use constructor? !!!
    std::map<vsize,std::vector<vsize> > findA, findB;
    // find non-empty configurations of the intersection variables in both A & B:
    for (_tableCIter i=t_.begin(); i!=t_.end(); ++i) {
      //std::cout<< "Pushing " << i->first <<"=> "<<isectVals(v_,vIsect,i->first)<<"\n";
      findA[isectVals(v_,vIsect,i->first)].push_back(i->first);
    }
    //std::cout<<"-------------------------------------------------\n";
    for (_tableCIter i=B.t_.begin(); i!=B.t_.end(); ++i) {
      //std::cout<< "Pushing " << i->first <<"=> "<<isectVals(v_,vIsect,i->first)<<"\n";
      findB[isectVals(B.v_,vIsect,i->first)].push_back(i->first);
    }
    //std::cout<<"=================================================\n";
    // for all such intersect configurations, and all completions in A & B, insert an element:
    for (std::map<vsize,std::vector<vsize> >::iterator i=findA.begin(); i!=findA.end(); ++i) {
      for (std::vector<vsize>::iterator jA=findA[i->first].begin(); jA!=findA[i->first].end(); ++jA) {
        for (std::vector<vsize>::iterator jB=findB[i->first].begin(); jB!=findB[i->first].end(); ++jB) {
          vsize idx = unionVals(v_,*jA,B.v_,*jB,vUnion);
          //std::cout<<"Setting "<<idx<<" <= "<<*jA<<","<<*jB<<"\n";
          F[idx] = Op(get(*jA), B[*jB]);
        }
      }
    }
    return F;
  }

  template<typename Function> Factor Factor::binaryOp( const value B, Function Op) const {
    Factor F=*this; F.binaryOpIP(B,Op); return F; // for scalar args, define with an in-place operator
  }
  // Binary in-place operations (eg A += B); returns reference to modified A
  template<typename Function> Factor& Factor::binaryOpIP( const Factor& B, Function Op) {
    if (!(v_ >> B.v_)) { Factor F=binaryOp(B,Op); swap(F); }  // if A's scope is too small, call non-in-place version
    else {
      for (_tableType::iterator i=t_.begin(); i!=t_.end(); ++i) {   // otherwise iterate over the entries of A
				_tableType::const_iterator bi = B.t_.find(isectVals(v_,B.v_,i->first)); // project to B's scope and find
				if (bi == B.t_.end()) t_.erase(i);						// if "zero" in B (missing), remove from this too
				else i->second=Op.IP(i->second, bi->second);  //   otherwise do the inplace operation
        //i->second=Op.IP(i->second, B[isectVals(v_,B.v_,i->first)]); //  project to B's scope & operate
			}
    }
    return *this;
  }
  template<typename Function> Factor& Factor::binaryOpIP( const value B, Function Op) {
    for (_tableType::iterator i=t_.begin(); i!=t_.end(); ++i) Op.IP(i->second, B);
    return *this;
  }

  Factor::value Factor::entropy(void) const {
    value H=0, Z=0;
    for (_tableCIter i=t_.begin();i!=t_.end();++i) {
      Z += i->second;
      if (i->second != 0) H -= (i->second)*std::log(i->second);
    }
    H/=Z; H+=std::log(Z);
    return H;
  }


  Factor Factor::sumPower(VarSet const& sumOut,Factor::value pow)  const {
    if (pow==1.0)          return sum(sumOut);
    else if(pow==-infty()) return min(sumOut);
    else if(pow== infty()) return max(sumOut);
    else {
      Factor F=*this; F.log(); F*=pow; F=F.logsumexp(sumOut); F/=pow; F.exp();
      return F;
    }
  }

// TODO: fix
  Factor Factor::logsumexp(const VarSet& sumOut) const {
    VarSet target = v_ - sumOut;
    Factor mx = maxmarginal(target);
    Factor B(target,0.0);      // TODO: FIX
    for (_tableCIter i=t_.begin(); i!=t_.end(); ++i) { size_t s2=isectVals(v_,target,i->first); if (mx[s2]!=-infty()) B[s2]+=std::exp(i->second - mx[s2]); }
    for (_tableCIter i=B.t_.begin(); i!=B.t_.end(); ++i) mx[i->first] += std::log(i->second); //=B[s]
    return mx;
  }

  double Factor::logsumexp() const {
    double r=0, mx=max();
    if (mx == -infty()) return mx;
    for (_tableCIter i=t_.begin();i!=t_.end();++i) r+=std::exp(i->second - mx);
    return std::log(r)+mx;
  }

// TODO: fix
  size_t Factor::argmax(const VarSet& vCond, Factor::vsize vState) const {
    if (vCond.size()==0) return argmax();
    size_t mxi=0; double mx=-infty();
    for (_tableCIter i=t_.begin();i!=t_.end();++i) { size_t s=isectVals(v_,vCond,i->first); if (s==vState) if (i->second > mx) { mx=i->second; mxi=i->first; } };
    return mxi;
  }
  // TODO: Implement more efficient version:
  //size_t Factor::argmax(const VarSet& vCond, Factor::vsize vState) const {
  //  VarSet vKeep = vars()-vCond;
  //  superindex sup(vars(),vKeep,vState); size_t N=vKeep.nrStates();
  //  size_t mxi=0; double mx=-infty();
  //  for (size_t i=0;i<N;++i,++sup) if (t_[sup] > mx) { mx=t_[sup]; mxi=sup; }
  //  return mxi;
  //}


  Factor Factor::condition(const VarSet& vRem, Factor::vsize vState) const {
    VarSet vKeep = vars() - vRem;
    Factor F(vKeep,0.0);
    for (_tableCIter i=t_.begin();i!=t_.end();++i) {size_t src=isectVals(v_,vRem,i->first), dst=isectVals(v_,vKeep,i->first); if (src==vState) F[dst]=i->second; }
    return F;
  }

  // TODO: implement more efficient version:
  //Factor condition(const VarSet& vRem, vsize vState) const {
  //  VarSet vKeep = vars() - vRem;
  //  Factor F(vKeep,0.0);
  //  superindex sup(vars(),vKeep,vState); size_t N=vKeep.nrStates();
  //  for (size_t i=0;i<N;++i,++sup) F[i]=t_[sup];
  //  for (size_t i=0;sup!=sup.end();++i,++sup) F[i]=t_[sup];
  //  return F;
  //}



  size_t Factor::sample() const {
    double x=0.0, z = sum();
    if (z==0.0) return mex::randi(nrStates());
    double y = mex::randu() * z;
    for (_tableCIter i=t_.begin();i!=t_.end();++i)
      if ((x+=i->second) > y) return i->first;
    return t_.size()-1;
  }

  Factor Factor::marginal(VarSet const& target) const {
    Factor F(target&vars(),0.0); // !!! sparse TODO
    for (_tableCIter i=t_.begin(); i!=t_.end(); ++i)
      F[isectVals(vars(),F.vars(),i->first)] += i->second;  // TODO: assumes default constructor to 0.0
    return F;
  }
  void Factor::marginalInto(VarSet const& target, Factor& F) const {
    assert(F.vars()==(target&vars()) && "marginalInto: target factor has incorrect variables");
    F.fill(0.0);  // TODO check
    for (_tableCIter i=t_.begin(); i!=t_.end(); ++i)
      F[isectVals(vars(),F.vars(),i->first)] += i->second;
  }

  Factor Factor::maxmarginal(VarSet const& target) const {
    Factor F(target&vars(),-infty()); // !!! sparse TODO
    for (_tableCIter i=t_.begin(); i!=t_.end(); ++i) {
      vsize s = isectVals(vars(),F.vars(),i->first);
      _tableIter Fs = F.t_.find(s);
      if (Fs == F.t_.end()) F[s]=i->second;       // can't find it => take current value
      else Fs->second=(Fs->second < i->second) ? i->second : Fs->second; // otherwise take max
    }
    return F;
  }

  Factor Factor::minmarginal(VarSet const& target) const {
    Factor F(target&vars(),infty()); // !!! sparse TODO
    for (_tableCIter i=t_.begin(); i!=t_.end(); ++i) {
      vsize s = isectVals(vars(),F.vars(),i->first);
      _tableIter Fs = F.t_.find(s);
      if (Fs == F.t_.end()) F[s]=i->second;       // can't find it => take current value
      else Fs->second=(Fs->second > i->second) ? i->second : Fs->second; // otherwise take min
    }
    return F;
  }

////////////////////////////////////////////////////////////////////////////////
// Misc other functions
////////////////////////////////////////////////////////////////////////////////

  // MEX_ENUM( Distance , L1,L2,LInf,KL,HPM,MAS,OptGap );

  double Factor::distance(Factor const& F2, Distance type) const {
    assert( vars() == F2.vars() && "distance: factors scopes do not match");
    Factor F(*this),Ftmp;               // make a copy for manipulation
    double dist=-1.0;                   // local variables
    value Z;
    switch (type) {
      case Distance::L2:                // L2, sum of squared errors
        F-=F2; F*=F; dist=F.sum();
        break;
      case Distance::L1:                // L1, sum of absolute errors  (=TV !!)
        F-=F2; dist=F.abs().sum();
        break;
      case Distance::LInf:              // L-infinity, max absolute error
        F-=F2; dist=F.abs().max();
        break;
      case Distance::KL:                // KL-divergence (relative entropy)
        Z=sum(); F/=F2; F*=F2.sum()/Z; F.log(); F*=*this; dist=F.sum()/Z;
        break;
      case Distance::HPM:               // Hilbert's projective metric
        F/=F2; F.log(); dist=F.max()-F.min(); //   aka "dynamic range"
        break;
      case Distance::MAS:               // "MAS" error value (not a metric)
        F.log(); Ftmp=F2; F/=Ftmp.log();
        dist = std::max( F.max(), 1.0/F.min() )-1.0;
        break;
      case Distance::OptGap:            // "Primal/Dual Gap"-like
        double mx1,mx2,gap1,gap2;
        mx1=mx2=gap1=gap2=0.0;
        for (size_t i=0;i<nrStates();++i) {  // TODO: fix
          if (mx1<F[i]) { gap2=F2[i]; mx1=F[i];} else if (mx1==F[i]) gap2 = std::min(gap2,F2[i]);
          if (mx2<F2[i]){ gap1=F[i]; mx2=F2[i];} else if (mx2==F2[i]) gap1= std::min(gap1,F[i]);
        }
        return (mx1-gap1)+(mx2-gap2);
        break;
      //case Distance::Hellinger:   (!!)
      //  F^=0.5; F-=F2^0.5; F*=F; dist=(0.5*F.sum())^0.5;  // straightforward computation
      //  F*=F2; F^=0.5; dist=(1-F.sum())^0.5;              // alternate computation if F,F2 normalized
      //  break;
      default: throw std::runtime_error("Invalid distance type");
    }
    return dist;
  }

  double Factor::norm(Distance type) const {
  Factor F(*this);                // make a copy for manipulation
  double dist=-1.0;               // 
  switch (type) {
    case Distance::L2:                // L2, sum of squared errors
      F*=F; dist=F.sum();
      break;
    case Distance::L1:                // L1, sum of absolute errors
      dist=F.abs().sum();
      break;
    case Distance::LInf:                // L-infinity, max absolute error
      dist=F.abs().max();
      break;
    case Distance::KL:                // KL-divergence (relative entropy => entropy?)
      return entropy();
      break;
    case Distance::HPM:               // Hilbert's projective metric
      F.log(); dist=F.max()-F.min();      //   aka "dynamic range"
      break;
    default:
      throw std::runtime_error("Invalid norm type");
  }
  return dist;
  }

  std::vector<Factor> Factor::decompSum(std::vector<VarSet> vlist, Factor::Decomp method) const {
    int nF=vlist.size();
    double mx,mn;
    std::vector<Factor> Flist(nF);

    Factor tmp,F=*this;
    switch (method) {
      case Decomp::L2: //L2
        double Cn,Cd;
        Cd=F.numel(); Cn=F.sum(); // /Cd*(1-1.0/nF);
        for (int j=0;j<nF;j++) {
          Flist[j] = F.marginal( vlist[j] );
          double D = Cd/Flist[j].numel();
          Flist[j]/= D;
          Flist[j]-= Cn/Cd*(1.0-1.0/(nF-j));
          F  -= Flist[j];
          Cn -= Flist[j].sum()*D;
        }
        break;
      case Decomp::L2_HPM: //L2+HPM
        Flist = decompSum(vlist,Decomp::L2);
        for (int j=0;j<nF;j++) F-=Flist[j];
        mx=F.max(); mn=F.min();
        for (int j=0;j<nF;j++) Flist[j]+=(mx+mn)/2/nF;
        break;
      case Decomp::L2_MAS: //L2+MAS
        Flist = decompSum(vlist,Decomp::L2);
        F=Flist[0]; for (int j=1;j<nF;j++) F+=Flist[j];
        F /= *this; F.log();
        mx=F.max(); mn=F.min();
        for (int j=0;j<nF;j++) Flist[j]*=std::exp(-(mx+mn)/2/nF);
        break;
      }
    return Flist;
  }

  std::vector<Factor> Factor::decompProd(std::vector<VarSet> vlist, Factor::Decomp method) const {
    Factor F=*this; F.log();
    std::vector<Factor> Flist = F.decompSum(vlist,method);
    for (size_t j=0;j<vlist.size();j++) Flist[j].exp();
    return Flist;
  }

  std::ostream& operator<<( std::ostream& out, const mex::Factor& F) {
    out << "Factor over " << F.variables() << ":";
    size_t last=-1;
    for (Factor::_tableCIter i=F.t_.begin();i!=F.t_.end();++i) {
      if (last+1!=i->first) out<<" ("<<i->first - last - 1<<")";
      out<<" "<<i->second;  // TODO:fix
      last = i->first;
    }
    if (last+1 != F.nrStates()) out<<" ("<<F.nrStates()-last-1<<")";
    return out;
  };


  Factor Factor::delta(const VarSet& vs, size_t idx) { Factor F(vs,0.0); F[idx]=1.0; return F; };



vector<Factor> Factor::readUai10(std::istream& is) {
  size_t nvar, ncliques, csize, v, nval;
	bool _sparse = false, _bayes = false;
  char* st; st = new char[20];
  is >> st;
  if ( strcasecmp(st,"MARKOV")==0 ) 
     ;               // 2010 Markov format (standard)
  else if ( strcasecmp(st,"BAYES")==0 )
     _bayes = true;  // 2008 Bayes format OK (TODO:check?) without modification
  else if ( strcasecmp(st,"SPARSE")==0 ) 
     _sparse = true; // New "sparse", run-oriented format
  else
    throw std::runtime_error("Unsupported UAI format type");

  // Read in # of variables and their dimensions:
  is >> nvar;
  vector<size_t> dims(nvar);
  for (size_t i=0;i<nvar;i++) is>>dims[i]; 
  for (size_t i=0;i<nvar;i++) if (dims[i] == 0) { dims[i]=1; std::cout<<"Warning: dimension of "<<i<<"=0; set to 1\n"; }
  
  // Read in # of cliques, then each clique: (size) v1 v2 ... vC 
  is >> ncliques;
  std::vector<mex::vector<Var> > cliques(ncliques);
  std::vector<VarSet > sets(ncliques);
  std::vector<bool> used(nvar, false);
  for (size_t i=0;i<ncliques;i++) {
    is >> csize;
    cliques[i].reserve(csize);
    for (size_t j=0;j<csize;j++) { 
      is>>v; Var V(v,dims[v]);
      cliques[i].push_back(V); 
      sets[i] |= V;
      used[v] = true;
    }   
  }
  // Check that model is well-defined for all variables:
  for (size_t i=0;i<nvar;i++) if (!used[i]) { std::cout<<"Warning: variable "<<i<<" does not appear in any factor\n"; }
  
  // Now read the factors themselves: 
  vector<Factor> tables(ncliques);
  for (size_t i=0;i<ncliques;i++) {
    is >> nval;                            // number of entries in the table
    assert(nval == sets[i].nrStates());    //   (should equal # of states of the clique)
    tables[i] = Factor(sets[i],0.0);        // preallocate memory and convert from given order, bigEndian
		if (nval==1) { is>>tables[i][0]; continue; }  // special case for constant factor
    permuteIndex pi( cliques[i], true); pi=pi.inverse();   // to our order
    for (size_t j=0;j<nval;j++) {
      if (_sparse && (is>>std::ws).peek() == '(') {
        char open,mid,close;
        size_t run;
        double val;
//        if ( is >> open >> run >> mid >> val >> close  && open=='(' && mid==':' && close==')' ) { }
        if ( is >> open >> val >> mid >> run >> close  && open=='(' && mid==':' && close==')' ) { }
        else {
          throw std::runtime_error("Error in parsing sparse entry, UAI file"); // TODO: make helpful
          // e.g., std::string("Error in parse: "+... )
        }
        if (val) for (size_t r=0;r<run;++r,++j) tables[i][pi.convert(j)] = val;  // TODO: Sparse
        else j+=run;
        --j; // subtract before loop re-adds
      }
      else {
        //if ( is>>tables[i][pi.convert(j)] )  { }
        double val;
        if ( is>>val ) { if (val || !_sparse) tables[i][pi.convert(j)]=val; }
        else {
          throw std::runtime_error("Error in parsing UAI file");  // TODO: make helpful (line/table#/etc)
        }
      }
    }
  }
  //TODO: for debugging, output fraction of tables used 
  //for (size_t i=0;i<ncliques;++i) std::cout<<"("<<tables[i].nrStates()<<","<<tables[i].numel()<<") "; std::cout<<"\n";
  //for (size_t i=0;i<ncliques;++i) std::cout<<tables[i]<<"\n"; 

  delete[] st;
  return tables;
}


void Factor::writeUai10(std::ostream& os, const vector<Factor>& fs) {
  size_t nvar=0;

  for (size_t f=0;f<fs.size();++f) if (fs[f].nvar()>0)          // find maximum variable label
    nvar=std::max(nvar,(size_t)(fs[f].vars().rbegin()->label()+1));
  
  vector<uint32_t> dims(nvar,1);                                // collect all variable dimensions
  for (size_t f=0,nvar=0;f<fs.size();++f) 
    for (size_t v=0;v<fs[f].nvar();++v) { Var V=fs[f].vars()[v]; dims[V.label()] = V.states(); }

  os << "MARKOV\n";                                             // Markov Random Field:
  os << nvar << "\n";                                           // # variables
  for (size_t v=0;v<nvar;++v) os << dims[v] << " ";  os<<"\n";  // dimensions
  os << fs.size() << "\n";                                      // # cliques
  for (size_t f=0;f<fs.size();++f) {
    os << fs[f].nvar() << "   ";                                // each clique: # var, var ids
    for (size_t v=0;v<fs[f].nvar();++v) os << fs[f].vars()[v].label() << " ";
    os << "\n";
  }
  for (size_t f=0;f<fs.size();++f) {                            // each function
    os << fs[f].nrStates() << "   ";                            // # states, then values
		if (fs[f].nrStates()==1) { os <<std::scientific << fs[f][0] <<"\n"; continue; } // special case: no vars
    vector<Var> clique( fs[f].vars().begin(), fs[f].vars().end() );
    permuteIndex pi( clique, true); pi=pi.inverse();
    for (size_t j=0;j<fs[f].nrStates();++j) os << std::scientific << fs[f][pi.convert(j)] << " ";
    os << "\n";
  }
}



//////////////////////////////////////////////////////////////////////////////////////////////
// MEX specific functions
//////////////////////////////////////////////////////////////////////////////////////////////
#ifdef MEX
bool Factor::mxCheckValid(const mxArray* M) {
  return mxIsStruct(M) && t_.mxCheckValid(mxGetField(M,0,"t")) && v_.mxCheckValid(mxGetField(M,0,"v"));
}

void Factor::mxSet(mxArray* M) {
  if (mxAvail()) {
    //mxDestroyArray(M_);            //   destroy old matlab object if it existed
  } // !!! what to do about v_ and t_?
  M_=M;
  if (M) {                            // if we got a matlab object, set to it
    t_.mxSetND(mxGetField(M_,0,"t"));
    v_.mxSet(mxGetField(M_,0,"v"), mxGetDimensions(t_.mxGet()));
  } else {                            // null pointer => clear 
    *this = Factor();
  }
}

mxArray* Factor::mxGet() {
  if (!mxAvail()) {
    mxArray* m; 
    int retval = mexCallMATLAB(1,&m,0,NULL,"factor");           // create new empty factor object
    if (retval) mexErrMsgTxt("Error creating new factor");      //   associated with a matlab mxArray
    Factor f; f.mxSet(m);                                        // use the copy constructor to put our data
    f = *this;                                                  //   there (in matlab-allocated memory)
    mxSwap(f);                                                  // then swap everything for the new object
    setDims();
  }
  return M_;
};

void Factor::mxRelease() {                                 // Disassociate with a given matlab object
  if (mxAvail()) {                                         //   without deallocating / destroying it
    t_.mxRelease(); v_.mxRelease(); 
  }
}
void Factor::mxDestroy() {                                  //  Disassociate with a given matlab object
  mxArray* m=M_;                                            //    and also destroy it
  mxRelease();
  if (m) mxDestroyArray(m);
}
void Factor::mxSwap(Factor& f) {
  v_.mxSwap(f.v_);
  t_.mxSwap(f.t_);
  std::swap(M_,f.M_);
}
#endif
//////////////////////////////////////////////////////////////////////////////////////////////



} // end namespace mex

