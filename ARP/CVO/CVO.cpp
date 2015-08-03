// CVO.cpp : Defines the entry point for the console application.
//

#include <stdlib.h>
#include <stdio.h>
#include <string>
#include <time.h>

#if defined WINDOWS || _WINDOWS
#include <process.h>    /* _beginthread, _endthread */
#endif // WINDOWS

#ifdef LINUX
#define stricmp strcasecmp
#define strnicmp strncasecmp
#endif

#include "AVLtreeSimple.hxx"
#include "ARPall.hxx"

static std::string logfile ;
static std::string uaifile ;
static std::string vofile ;
static std::string alg ;
static ARE::VarElimOrderComp::ObjectiveToMinimize objCode = ARE::VarElimOrderComp::Width; // 0=width, 1=state space size
static ARE::VarElimOrderComp::NextVarPickCriteria algCode = ARE::VarElimOrderComp::MinFill ; // 0=MinFill, 1=MinInducedWidth(MinDegree), 2=MinComplexity
static int nThreads = 0 ;
static int nRunsToDo = 0 ;
static INT64 timelimit = 0 ;
static int nRandomPick = 4 ;
static double eRandomPick = 0.5 ;
static bool PerformSingletonConsistencyChecking = false ;
static bool EarlyTerminationOfBasic_W = false, EarlyTerminationOfBasic_C = false ;

int main(int argc, char *argv[])
{
	ARE::Initialize() ;

	int nParams = argc ;
	if (14 != nParams) {
		printf("\nUSAGE:\n $ %s logfile uaifile vofile objective={w,s} algorithm={MinFill,MinDegree,MinComplexity} nRandomPick eRandomPick #tries timelimit[msec] nThreads doSingletonConsCheck={0,1} earlyTerminationW={0,1} earlyTerminationC={0,1}\n", argv[0]) ;
		printf("\nExample:\n $ %s log.txt problem.uai problem.vo w MinFill 4 0.5 500 60 2 1 1 0\n\n", argv[0]) ;
		return 0 ;
		}

	logfile = (NULL != argv[1] ? argv[1] : "") ;
	uaifile = (NULL != argv[2] ? argv[2] : "") ;
	vofile = (NULL != argv[3] ? argv[3] : "") ;
	if (NULL != argv[4]) objCode = 's' == argv[4][0] || 'S' == argv[4][0] ? ARE::VarElimOrderComp::StateSpaceSize : ARE::VarElimOrderComp::Width;
	alg = (NULL != argv[5]) ? argv[5] : 0 ;
	nRandomPick = (NULL != argv[6]) ? atoi(argv[6]) : 0 ;
	eRandomPick = (NULL != argv[7]) ? atof(argv[7]) : 0 ;
	nRunsToDo = (NULL != argv[8]) ? atoi(argv[8]) : 0 ;
	timelimit = (NULL != argv[9]) ? _atoi64(argv[9]) : 0 ;
	nThreads = (NULL != argv[10]) ? atoi(argv[10]) : 0 ;
	PerformSingletonConsistencyChecking = '1' == *argv[11] ? true : false ;
	EarlyTerminationOfBasic_W = '1' == *argv[12] ? true : false ;
	EarlyTerminationOfBasic_C = '1' == *argv[13] ? true : false ;

	ARE::fpLOG = fopen(logfile.c_str(), "w") ;

	if (uaifile.length() < 1) {
		printf("\nNo UAI file; will exit ...") ;
		return 0 ;
		}

	if (alg.length() < 1) 
		alg = "MinFill" ;
	if (0 == stricmp("MinFill", alg.c_str())) 
		algCode = ARE::VarElimOrderComp::MinFill ;
	else if (0 == stricmp("MinDegree", alg.c_str())) 
		algCode = ARE::VarElimOrderComp::MinDegree ;
	else if (0 == stricmp("MinComplexity", alg.c_str())) 
		algCode = ARE::VarElimOrderComp::MinStateSpaceSize ;
	else 
		algCode = ARE::VarElimOrderComp::MinFill ;

#if defined WINDOWS || _WINDOWS
	SYSTEM_INFO sysinfo ;
	GetSystemInfo(&sysinfo) ;
	int nThreadsMax = sysinfo.dwNumberOfProcessors ;
#elif defined (LINUX)
	int nThreadsMax = sysconf(_SC_NPROCESSORS_ONLN) ;
#endif
	if (nThreadsMax > 128) 
		nThreadsMax = 128;
	else if (nThreadsMax < 1) {
		printf("\nNO THREADS ... will exit ...") ;
		return 0 ;
		}
	if (nThreads > nThreadsMax) 
		nThreads = nThreadsMax ;
	else if (nThreads < 1) 
		nThreads = 1 ;
	if (nThreads > nRunsToDo-1) 
		nThreads = nRunsToDo-1 ;
// DEBUGGG
//nThreads = 1 ;
//	if (nThreads > 1 && nThreads == nThreadsMax) 
//		// if there is lotsa threads, leave 1 for others
//		nThreads -= 1 ;

	// at least one of nRunsToDo/timelimit has to be given
	if (nRunsToDo < 1 && timelimit < 1) 
		nRunsToDo = 1 ;

	ARE::VarElimOrderComp::Order BestOrder ;
	ARE::VarElimOrderComp::CVOcontext *cvocontext = new ARE::VarElimOrderComp::CVOcontext ;

	printf("\nobjCode=%d algCode=%d nRunsToDo=%d timelimit=%I64d nThreads=%d/%d practicalorder=%c ...", (int) objCode, (int) algCode, (int) nRunsToDo, (INT64) timelimit, (int) nThreads, (int) nThreadsMax, cvocontext->_FindPracticalVariableOrder ? 'Y' : 'N') ;
	printf("\nnRandomPick=%d eRandomPick=%g ...", (int) nRandomPick, (double) eRandomPick) ;
	printf("\nEarlyTerminationOfBasic_W=%c EarlyTerminationOfBasic_C=%c ...", EarlyTerminationOfBasic_W ? 'Y' : 'N', EarlyTerminationOfBasic_C ? 'Y' : 'N') ;
	if (NULL != ARE::fpLOG) {
		fprintf(ARE::fpLOG, "\nobjCode=%d algCode=%d nRunsToDo=%d timelimit=%I64d nThreads=%d/%d practicalorder=%c ...", (int) objCode, (int) algCode, (int) nRunsToDo, (INT64) timelimit, (int) nThreads, (int) nThreadsMax, cvocontext->_FindPracticalVariableOrder ? 'Y' : 'N') ;
		fprintf(ARE::fpLOG, "\nnRandomPick=%d eRandomPick=%g ...", (int) nRandomPick, (double) eRandomPick) ;
		fprintf(ARE::fpLOG, "\nEarlyTerminationOfBasic_W=%c EarlyTerminationOfBasic_C=%c ...", EarlyTerminationOfBasic_W ? 'Y' : 'N', EarlyTerminationOfBasic_C ? 'Y' : 'N') ;
		fflush(ARE::fpLOG) ;
		}

	int ret = ARE::VarElimOrderComp::Compute(uaifile, objCode, algCode, nThreads, nRunsToDo, timelimit, nRandomPick, eRandomPick, PerformSingletonConsistencyChecking, EarlyTerminationOfBasic_W, EarlyTerminationOfBasic_C, BestOrder, cvocontext) ;

	if (0 == ret) {
		INT64 runtime = cvocontext->_tEnd - cvocontext->_tStart ;
		printf("\n\nBEST ORDER w width=%d; TotalVarElimComplexity=%g, TotalNewFunctionStorageAsNumOfElements=%g, nFillEdges=%d, nImprvmnts=%d, nRuns=%d/%d, runtime=%I64dmsec, file=%s", BestOrder._Width, BestOrder._Complexity, BestOrder._TotalNewFunctionStorageAsNumOfElements, BestOrder._nFillEdges, cvocontext->_nImprovements, cvocontext->_nRunsCompleted, cvocontext->_nRunsDone, runtime, vofile.c_str()) ;
		if (NULL != ARE::fpLOG) {
			fprintf(ARE::fpLOG, "\n\nBEST ORDER w width=%d; TotalVarElimComplexity=%g, TotalNewFunctionStorageAsNumOfElements=%g, nFillEdges=%d, nImprvmnts=%d, nRuns=%d/%d, runtime=%I64dmsec, file=%s", BestOrder._Width, BestOrder._Complexity, BestOrder._TotalNewFunctionStorageAsNumOfElements, BestOrder._nFillEdges, cvocontext->_nImprovements, cvocontext->_nRunsCompleted, cvocontext->_nRunsDone, runtime, vofile.c_str()) ;
			fprintf(ARE::fpLOG, "\nImprovements list :") ;
			for (int i = 0 ; i < cvocontext->_nImprovements ; i++) 
				fprintf(ARE::fpLOG, " %d", cvocontext->_Improvements[i]._dt) ;
			fprintf(ARE::fpLOG, "\nWidth2StatMap :") ;
			for (int i = 0 ; i < 255 ; i++) {
				if (0 == cvocontext->_Width2CountMap[i]) 
					continue ;
				fprintf(ARE::fpLOG, "\n w=%d count=%I64d min=%g max=%g", i, cvocontext->_Width2CountMap[i], cvocontext->_Width2MinComplexityMap[i], cvocontext->_Width2MaxComplexityMap[i]) ;
				}
			fprintf(ARE::fpLOG, "\nWidth2StatMap : done") ;
			}
		FILE *fpVO_w = fopen(vofile.c_str(), "w") ;
		if (NULL != fpVO_w) {
			fprintf(fpVO_w, "#filename=%s, width=%d, TotalVarElimComplexity=%g, TotalNewFunctionStorageAsNumOfElements=%g", uaifile.c_str(), BestOrder._Width, BestOrder._Complexity, BestOrder._TotalNewFunctionStorageAsNumOfElements) ;
			fprintf(fpVO_w, "\n%d", (int) BestOrder._nVars) ;
			for (int i = 0 ; i < BestOrder._nVars ; i++) 
				fprintf(fpVO_w, "\n%d", (int) BestOrder._VarListInElimOrder[i]) ;
			fclose(fpVO_w) ;
			}

		std::string vofile_w_rss(vofile + "_rss.txt") ;
		FILE *fpVO_w_rss = fopen(vofile_w_rss.c_str(), "w") ;
		if (NULL != fpVO_w_rss) {
			fprintf(fpVO_w_rss, "#filename=%s, vofile=%s, nRandomPick=%g, eRandomPick=%g, nRunsToDo=%d, timelimit=%d, nThreads=%d, doSCC=%c", uaifile.c_str(), vofile.c_str(), (double) nRandomPick, (double) eRandomPick, (int) nRunsToDo, (int) timelimit, (int) nThreads, PerformSingletonConsistencyChecking ? 'Y' : 'N') ;
			fprintf(fpVO_w_rss, "\n%d", (int) cvocontext->_nImprovements) ;
			for (int i = 0 ; i < cvocontext->_nImprovements ; i++) 
				fprintf(fpVO_w_rss, "\n%d %d %g", cvocontext->_Improvements[i]._dt, cvocontext->_Improvements[i]._width, cvocontext->_Improvements[i]._complexity) ;
			fclose(fpVO_w_rss) ;
			}

/*
		// transform from elimination order to BE-order.
		for (i = 0 ; i < BestOrderN >> 1 ; i++) {
			int u = BestOrder[i] ;
			int v = BestOrder[BestOrderN-i-1] ;
			BestOrder[i] = v ;
			BestOrder[BestOrderN-i-1] = u ;
			}
		for (i = 0 ; i < BestOrderN >> 1 ; i++) {
			int u = MinTotalComplexity_Order[i] ;
			int v = MinTotalComplexity_Order[BestOrderN-i-1] ;
			MinTotalComplexity_Order[i] = v ;
			MinTotalComplexity_Order[BestOrderN-i-1] = u ;
			}

		// TEST THE ORDER WIDTH
		p.ComputeGraph() ;
		{
		BucketElimination::BEworkspace *bews = new BucketElimination::BEworkspace(NULL) ;
		if (0 != bews->Initialize(p, BestOrder)) {
			if (NULL != ARE::fpLOG) 
				fprintf(ARE::fpLOG, "\nbews.Initialize() failed ...") ;
			}
		int maxNumChildrenBT = bews->ComputeMaxNumChildren() ;
		int maxBFWidth = bews->GetMaxBucketFunctionWidth() ;
		int maxNumVarsInBucket = bews->ComputeMaxNumVarsInBucket() ;
		__int64 newFnSpace = bews->ComputeNewFunctionSpace() ;
		__int64 newFnComplexity = bews->ComputeNewFunctionComputationComplexity() ;
		printf("\n\n Verified w : BE workspace maxNumVarsInBucket=%d newFnComplexity=%I64d newFnSpace(bytes)=%I64d", maxNumVarsInBucket, newFnComplexity, newFnSpace) ;
		}
		{
		BucketElimination::BEworkspace *bews = new BucketElimination::BEworkspace(NULL) ;
		if (0 != bews->Initialize(p, MinTotalComplexity_Order)) {
			if (NULL != ARE::fpLOG) 
				fprintf(ARE::fpLOG, "\nbews.Initialize() failed ...") ;
			}
		int maxNumChildrenBT = bews->ComputeMaxNumChildren() ;
		int maxBFWidth = bews->GetMaxBucketFunctionWidth() ;
		int maxNumVarsInBucket = bews->ComputeMaxNumVarsInBucket() ;
		__int64 newFnSpace = bews->ComputeNewFunctionSpace() ;
		__int64 newFnComplexity = bews->ComputeNewFunctionComputationComplexity() ;
		printf("\n\n Verified c : BE workspace maxNumVarsInBucket=%d newFnComplexity=%I64d newFnSpace(bytes)=%I64d", maxNumVarsInBucket, newFnComplexity, newFnSpace) ;
		}
*/

		}

	return ret ;
}

