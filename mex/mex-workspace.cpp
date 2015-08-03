// Wrapper for UAI competition code
// 

#include <cstdio>
#include <iostream>
#include <fstream>

#include "factorgraph.h"

#include "mplp.h"
#include "mbe.h"

#include "lbp.h"
#include "gbp.h"

#include "mex/mex-workspace.h"

using namespace std;
using mex::mxObject;
using mex::Var;
using mex::VarSet;
using mex::Factor;
using mex::vector;
using mex::graphModel;
using mex::factorGraph;
using mex::mplp;

//MEX_ENUM(Task, MPE, PR, MAR) ;

using mex::timeSystem;

#define c_log10 std::log(10)

//double MemLimit;
double lbpTime, lbpIter, lbpObj, lbpErr;
double gbpTime, gbpIter;
double dt;
int nOrders,nExtra;
double timeOrder;
double memUseRandom = mex::infty(); // if memory *way* out of reach, just use random ordering...

/*
int main__(int argc, char* argv[])
{
  double timeStart = timeSystem();

  //if (argc < 4) { cout<<"Usage: "<<argv[0]<<" <file.uai> <seed> <PR|MPE|MAR>\n"; return 0; }
  const char* probName; // = argv[1];
  const char* taskName;//  = argv[3];
  mex::Task task;
  mex::vector<Factor> bel;

  po::options_description desc("Available options");
  desc.add_options()
    ("help", "print help message")
    ("file,f", po::value<std::string>(), "input problem filename")
    ("evidence,e", po::value<std::string>(), "input evidence filename")
    ("seed,S", po::value<int>(),         "random number initial seed")
    ("task,T", po::value<std::string>(), "inference task string")
    ("ibound,i", po::value<int>(),       "initial i-bound")
    ("orders,o",    po::value<int>(&nOrders)->default_value(1),      "number of variable orderings to try")
    ("order-time,t",po::value<double>(&timeOrder)->default_value(1), "max time spend on variable orderings")
    ("order-rand",  po::value<int>(&nExtra)->default_value(0),   "var order randomness; n=-1 (none), or among best+n")
		("order-file", po::value<std::string>(), "problem elimination ordering filename")
    ("memory,m", po::value<double>(&MemLimit)->default_value(2*1024.0),    "memory bound (MB)")
    ("simulate", "simulate only for memory info")
  ;

  po::variables_map vm;
  po::store(po::parse_command_line(argc,argv,desc),vm);
  po::notify(vm);

  if (vm.count("help")) { std::cout<<desc<<"\n"; return 1; }
  if (vm.count("file")) { probName=vm["file"].as<std::string>().c_str(); }
  else { std::cout<<"Missing input problem file!\n"; return 1; }
  if (vm.count("seed")) { mex::randSeed( vm["seed"].as<int>() ); }
  if (vm.count("task")) { taskName = vm["task"].as<std::string>().c_str(); task=Task(taskName); }
  else { std::cout<<"Missing task!\n"; return 1; }


  //*** READ IN PROBLEM FILE **********************************************************
	std::cout<<"Reading model file: "<<probName<<"\n";
  ifstream is; is.open(probName);
  if (!is.is_open()) throw std::runtime_error("Failed to open problem file");
  mex::vector<Factor> forig = Factor::readUai10(is);
	mex::graphModel gmo( forig );
  size_t nvar = gmo.nvar();
  bel.resize(nvar); 
	mex::vector<size_t> evid(nvar, 0);

	//*** READ IN EVIDENCE FILE *********************************************************
  VarSet evVar;
  ifstream is2;
  if (vm.count("evidence")) { is2.open(vm["evidence"].as<std::string>().c_str()); }
  if (is2.is_open()) {
    std::cout<<"Got evidence file\n";
    //int nEvid; is2 >> nEvid;
		int nEvid=1;      // 2014 format: single evidence, no # of evidences entry
		std::cout<<"Got "<<nEvid<<" evidences?\n";
    if (nEvid > 0) {
      int nEvidVar; is2 >> nEvidVar;
      for (size_t i=0;i<nEvidVar;i++) {
        uint32_t vid; size_t vval; is2>>vid>>vval; 
        evid[vid]=vval; evVar |= gmo.var(vid);
				bel[vid] = Factor::delta(gmo.var(vid),vval);
      }
			std::cout<<"Evidence on variables "<<evVar<<"\n";
			gmo.condition(evVar,evid);
    }
  } else std::cout<<"Evidence file not specified or not found\n";


  //*** PREPARE REQUESTED TASK ********************************************************
  std::cout<<"Task is "<<task<<"\n";
  //std::cout<<"Task is "<<(const char*)task<<"\n";
  if (task==Task::MPE) { std::cout<<"MPE task not supported\n"; return 1; }


  //*** PREPARE OUTPUT FILE ***********************************************************
  std::string outfiles(probName); outfiles += '.'; outfiles += taskName;
  std::string::size_type start = outfiles.find_last_of('/');
  if (start==std::string::npos) start=0; else ++start;
  outfiles = outfiles.substr(start,std::string::npos);
  const char* outfile = outfiles.c_str();
  std::cout<<"Writing to "<<outfile<<"\n";

  double ln10 = std::log(10);


  //*** LOOPY BELIEF PROPAGATION ******************************************************
  mex::lbp fg(gmo.factors()); 

    double mbCutoff = MemLimit/sizeof(double)*1024*1024;     // translate memory into MBE cutoff
    mex::mbe mb(gmo.factors());  mb.setProperties("ElimOp=SumUpper,sBound=inf,DoMplp=0,DoMatch=0,DoFill=0");
    if (task==Task::PR)       mb.setProperties("DoJG=0");
    else if (task==Task::MAR) mb.setProperties("DoJG=1");

    // Calculate elimination order(s) ////////////////////////////
    double startOrder = timeSystem();
    mex::VarOrder order;
    size_t InducedWidth, iOrder = 1;
    double score = fg.order( mex::graphModel::OrderMethod::Random, order, 0, memUseRandom );

    const char *orderFile = NULL;       // Check for pre-specified elimination order file
    if (vm.count("order-file")) { orderFile = vm["order-file"].as<std::string>().c_str(); }
    ifstream orderIStream; if (orderFile!=NULL) orderIStream.open(orderFile);
    if (orderIStream.is_open()) {
      // If we were given an input file with an elimination ordering, just use that
      std::cout << "Reading elimination order from "<<orderFile<<"\n";
      size_t ordersize;  orderIStream>>ordersize; assert(ordersize == order.size());
      for (size_t i=0;i<order.size();++i) { size_t tmp;  orderIStream>>tmp; order[i]=tmp;  };
      orderIStream.close();
			
    } else {
      // Otherwise, calculate elimination order(s) ////////////////////////////////////
      double startOrder = timeSystem();
      size_t iOrder = 0;
      // Try to build new orders until time or count limit reached ////////////////////
      while (iOrder < nOrders && (timeSystem()-startOrder < timeOrder)) {
        score = fg.order(mex::graphModel::OrderMethod::WtMinFill, order, nExtra, score);
        ++iOrder;
      }

      // If we were given an ordering file name but no file, write our order out 
      ofstream orderOStream; if (orderFile!=NULL) orderOStream.open(orderFile);
      if (orderOStream.is_open()) {
        std::cout << "Writing elimination order to "<<orderFile<<"\n";
        orderOStream<<order.size();
        for (size_t i=0;i<order.size();++i) orderOStream<<" "<<order[i];
        orderOStream<<"\n";
        orderOStream.close();
      }
		}

    InducedWidth = fg.inducedWidth(order);
    std::cout<<"Best order of "<<iOrder<<" has induced width "<<InducedWidth<<", score "<<score<<"\n";
    mb.setOrder(order);

    mb.setIBound(InducedWidth);

    VarSet cond;
    mex::vector<VarSet> cliques;
    mex::vector<Factor> blank;
    mex::gbp _gbp2(blank);
    double mem = -1;
    if (task==Task::MAR) mbCutoff *= 2;

    double  mbMem = mb.simulateMemory(&cliques,NULL,mbCutoff);
    if (task==Task::MAR) { _gbp2.addRegions(cliques); mem = _gbp2.memory(); }
    std::cout<<"Starting: "<<mbMem*sizeof(double)/1024/1024<<"MB\n";
    while (mbMem > mbCutoff || mem > MemLimit) {
      cond += fg.bestConditioner(order,cond);
      cliques.clear();
      mbMem = mb.simulateMemory(&cliques,&cond,mbCutoff);
      if (task==Task::MAR) { _gbp2.clearRegions(); _gbp2.addRegions(cliques); mem=_gbp2.memory(); }
std::cout<<_gbp2.nRegions()<<" => "<<cliques.size()<<" cliques\n";
      std::cout<<"Conditioning "<<cond<<" => "<<mbMem*sizeof(double)/1024/1024<<"MB, "<<mem<<"MB\n";
    }
    std::cout<<"Conditioner has "<<cond.nrStates()<<" values\n";


    if (vm.count("simulate")) { std::cout<<"Simulation Done\n"; return 0; }



    Factor lnZ(cond);
    mex::vector<mex::vector<Factor> > condMarginals;
    //mex::vector<mex::gbp::findex> regions(fg.nvar());
    //if (task==Task::MAR) for (size_t v=0;v<fg.nvar();++v) {
    //  if (!evVar.contains(Var(v,0))&& !cond.contains(Var(v,0))) regions[v]=_gbp.regionWith(Var(v,0));
    //}

    for (size_t i=0;i<lnZ.nrStates();++i) {
      std::map<Var,size_t> val;  ind2sub(cond,i,val);
			mex::graphModel gmcond( gmo.factors() ); 
			ind2sub(cond,i,evid); gmcond.condition(cond, evid);	// !!! use evid vector as temp storage for conditioning
      if (task==Task::PR) {
        mex::mbe mbc(gmcond.factors()); mbc.setOrder(order); mbc.setIBound(InducedWidth); mbc.setProperties("ElimOp=SumUpper");
        mbc.init();
        lnZ[i] = mbc.logZ();
      } else if (task == Task::MAR) {
    		mex::gbp _gbp(gmcond.factors()); _gbp.setProperties("Schedule=Fixed");
				mex::vector<mex::gbp::findex> regions(nvar);
        //_gbp.setFactors(gmcond.factors()); _gbp.setProperties("Schedule=Fixed");
				_gbp.clearRegions(); _gbp.addRegions(cliques);
        _gbp.init();
        _gbp.setStopIter(-1); _gbp.setStopObj(1e-7); _gbp.setStopMsg(-1.0); _gbp.setStopTime( 10000 ); _gbp.setVerbose(0);
        _gbp.run();
        lnZ[i] = _gbp.logZ();
        if (task==Task::MAR) {
					for (size_t v=0;v<fg.nvar();++v) regions[v]=_gbp.regionWith(Var(v,0));
					//for (size_t v=0;v<fg.nvar();++v) if (!evVar.contains(Var(v,0))&& !cond.contains(Var(v,0))) regions[v]=_gbp.regionWith(Var(v,0));
          condMarginals.push_back(bel);
          for (size_t v=0;v<fg.nvar();++v)
            if (!evVar.contains(Var(v,0)) && !cond.contains(Var(v,0)))
              condMarginals[i][v]=_gbp.computeRegionBelief(regions[v]).marginal(Var(v,0));
        }
      }

      for (size_t v=0;v<cond.size();++v) std::cout<<((mex::VarSet const&)cond[v])<<"="<<val[cond[v]]<<" "; std::cout<<lnZ[i]<<"\n";
    }
    std::cout<<"lnZ scores "<<lnZ<<"\n";
    double lnZtot = lnZ.logsumexp();
    std::cout<<"Final lnZ "<<lnZtot<<"\n";
    if (task==Task::PR) writePR(outfile,lnZtot);
    if (task==Task::MAR) {
      Factor probs = (lnZ - lnZtot).exp();
      for (size_t v=0;v<fg.nvar();++v) {
        if (evVar.contains(Var(v,0))) { } // evidence variables not updated
        else if (cond.contains(Var(v,0)))  { bel[v] = probs.marginal(Var(v,0)); }
        else {
          bel[v] = condMarginals[0][v] * probs[0];
          for (size_t i=1;i<lnZ.nrStates();++i) bel[v] += condMarginals[i][v] * probs[i];
        }
      }
      writeMAR(outfile, bel);
    }

  return 0;

}
*/

/*
static std::string getFileContents(const char* filename)
{
	std::ifstream in(filename, std::ios::in | std::ios::binary);
	std::string contents;
	in.seekg(0, std::ios::end);
	contents.resize(in.tellg());
	in.seekg(0, std::ios::beg);
	in.read(&contents[0], contents.size());
	in.close();
	return(contents);
}
*/

int mex::Workspace::Solve(void)
{
	size_t nvar = _gmo->nvar() ;

	/*** LOOPY BELIEF PROPAGATION ******************************************************/
	mex::lbp fg(_gmo->factors()); 

	if (_MemLimit < 1024) 
		_MemLimit = 1024 ;
	else if (_MemLimit > 4294967296) 
		_MemLimit = 4294967296 ;

	double mbCutoff = _MemLimit/sizeof(double)*1024*1024;     // translate memory into MBE cutoff
	mex::mbe mb(_gmo->factors());  mb.setProperties("ElimOp=SumUpper,sBound=inf,DoMplp=0,DoMatch=0,DoFill=0");
	if (_task==Task::PR)       mb.setProperties("DoJG=0");
	else if (_task==Task::MAR) mb.setProperties("DoJG=1");

	// Calculate elimination order(s) ////////////////////////////
	double startOrder = timeSystem();
	mex::VarOrder order;
	size_t InducedWidth, iOrder = 1;
	double score = fg.order( mex::graphModel::OrderMethod::Random, order, 0, memUseRandom );

	// var elim order must be given
	if (_VarElimOrder.size() != nvar) 
		return 1 ;
	for (int i = 0 ; i < nvar ; i++) order[i] = _VarElimOrder[i] ;

	InducedWidth = fg.inducedWidth(order);
//	std::cout<<"Best order of "<<iOrder<<" has induced width "<<InducedWidth<<", score "<<score<<"\n";
	mb.setOrder(order);
	mb.setIBound(InducedWidth);

	VarSet cond;
	mex::vector<VarSet> cliques;
	mex::vector<Factor> blank;
	mex::gbp _gbp2(blank);
	double mem = -1;
	if (_task==Task::MAR) mbCutoff *= 2;

	double  mbMem = mb.simulateMemory(&cliques,NULL,mbCutoff);
	if (_task==Task::MAR) { _gbp2.addRegions(cliques); mem = _gbp2.memory(); }
//	std::cout<<"Starting: "<<mbMem*sizeof(double)/1024/1024<<"MB\n";
	while (mbMem > mbCutoff || mem > _MemLimit) {
		cond += fg.bestConditioner(order,cond);
		cliques.clear();
		mbMem = mb.simulateMemory(&cliques,&cond,mbCutoff);
		if (_task==Task::MAR) { _gbp2.clearRegions(); _gbp2.addRegions(cliques); mem=_gbp2.memory(); }
//		std::cout<<_gbp2.nRegions()<<" => "<<cliques.size()<<" cliques\n";
//		std::cout<<"Conditioning "<<cond<<" => "<<mbMem*sizeof(double)/1024/1024<<"MB, "<<mem<<"MB\n";
		}
//	std::cout<<"Conditioner has "<<cond.nrStates()<<" values\n";

//	if (vm.count("simulate")) { std::cout<<"Simulation Done\n"; return 0; }

	Factor lnZ(cond);
	mex::vector<mex::vector<Factor> > condMarginals;
	//mex::vector<mex::gbp::findex> regions(fg.nvar());
	//if (task==Task::MAR) for (size_t v=0;v<fg.nvar();++v) {
	//  if (!evVar.contains(Var(v,0))&& !cond.contains(Var(v,0))) regions[v]=_gbp.regionWith(Var(v,0));
	//}

	for (size_t i=0;i<lnZ.nrStates();++i) {
		std::map<Var,size_t> val;  ind2sub(cond,i,val);
		mex::graphModel gmcond( _gmo->factors() ); 
		ind2sub(cond,i,_evid); gmcond.condition(cond, _evid);	// !!! use evid vector as temp storage for conditioning
		if (_task==Task::PR) {
			mex::mbe mbc(gmcond.factors()); mbc.setOrder(order); mbc.setIBound(InducedWidth); mbc.setProperties("ElimOp=SumUpper");
			mbc.init();
			lnZ[i] = mbc.logZ();
			} 
		else if (_task == Task::MAR) {
			mex::gbp _gbp(gmcond.factors()); _gbp.setProperties("Schedule=Fixed");
			mex::vector<mex::gbp::findex> regions(nvar);
			//_gbp.setFactors(gmcond.factors()); _gbp.setProperties("Schedule=Fixed");
			_gbp.clearRegions(); _gbp.addRegions(cliques);
			_gbp.init();
			_gbp.setStopIter(-1); _gbp.setStopObj(1e-7); _gbp.setStopMsg(-1.0); _gbp.setStopTime( 10000 ); _gbp.setVerbose(0);
			_gbp.run();
			lnZ[i] = _gbp.logZ();
			if (_task==Task::MAR) {
				for (size_t v=0;v<fg.nvar();++v) regions[v]=_gbp.regionWith(Var(v,0));
				//for (size_t v=0;v<fg.nvar();++v) if (!evVar.contains(Var(v,0))&& !cond.contains(Var(v,0))) regions[v]=_gbp.regionWith(Var(v,0));
				condMarginals.push_back(_bel);
				for (size_t v=0;v<fg.nvar();++v)
					if (!_evVar.contains(Var(v,0)) && !cond.contains(Var(v,0)))
						condMarginals[i][v]=_gbp.computeRegionBelief(regions[v]).marginal(Var(v,0));
				}
			}

		for (size_t v=0;v<cond.size();++v) std::cout<<((mex::VarSet const&)cond[v])<<"="<<val[cond[v]]<<" "; std::cout<<lnZ[i]<<"\n";
		}

//	std::cout<<"lnZ scores "<<lnZ<<"\n";
	_lnZtot = lnZ.logsumexp();
//	std::cout<<"Final lnZ "<<_lnZtot<<"\n";
//	if (_task==Task::PR) writePR(outfile,_lnZtot);
	if (_task==Task::MAR) {
		Factor probs = (lnZ - _lnZtot).exp();
		for (size_t v=0;v<fg.nvar();++v) {
			if (_evVar.contains(Var(v,0))) { } // evidence variables not updated
			else if (cond.contains(Var(v,0)))  { _bel[v] = probs.marginal(Var(v,0)); }
			else {
				_bel[v] = condMarginals[0][v] * probs[0];
				for (size_t i=1;i<lnZ.nrStates();++i) _bel[v] += condMarginals[i][v] * probs[i];
				}
			}
//		writeMAR(outfile, _bel);
		}

/*	if (false) {
		std::string uaifile("C:\\UCI\\problems\\handbuilt_burglary.uai") ;
		std::string problem = getFileContents(uaifile.c_str()) ;
		std::string evidence ;

		mex::Workspace mexws ;
		mexws.SetTask("MAR") ;
		mexws.LoadUAI(problem.c_str(), problem.length()) ;

		std::vector<int> order ;
		int n = 4 ;
		for (int i = n-1 ; i >= 0 ; i--) order.push_back(i) ;
		mexws.SetVarElimOrdering(order, -1) ;

		mexws._lnZtot = DBL_MAX ;
		int res = mexws.Solve() ;

		return 0 ;
		}*/

	return 0 ;
}

