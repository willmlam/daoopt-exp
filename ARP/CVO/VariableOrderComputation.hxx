#ifndef ARE_VariableOrderComputation_HXX_INCLUDED
#define ARE_VariableOrderComputation_HXX_INCLUDED

#include <stdlib.h>
#include <stdio.h>
#include <string>

#include "Graph.hxx"

namespace BucketElimination { class MBEworkspace ; }

namespace ARE
{
	
namespace VarElimOrderComp
{

enum ObjectiveToMinimize
{
	None,
	Width, 
	StateSpaceSize, 
	BTheight
} ;

enum NextVarPickCriteria
	{
	MinFill,
	MinDegree, 
	MinStateSpaceSize
	};

#define TempAdjVarSpaceSize 1000000
#define TempAdjVarSpaceSizeExtraArraySize 1000

class ResultSnapShot
{
public :
	int64_t _dt ; // in milliseconds
	int _width ;
	double _complexity ;
public :
	void operator=(const ResultSnapShot & O)
	{
		_dt = O._dt ;
		_width = O._width ;
		_complexity = O._complexity ;
	}
public :
	ResultSnapShot(void)
		:
		_dt(0), 
		_width(0), 
		_complexity(0.0)
	{
	}
} ;

class Order
{
public :
	int _nVars ;
	int *_VarListInElimOrder ;
	int _Width ;
	int _WidthLowerBound ;
	double _Complexity_Log10 ; // log of
	double _TotalNewFunctionStorageAsNumOfElements_Log10 ; // log of
	double _MaxSingleVarElimComplexity ;
	int _nFillEdges ;
public :
	int SerializeAsElimOrder(const char *FN)
	{
		FILE *fp_eo = fopen(FN, "w") ;
		if (NULL == fp_eo) 
			return 1 ;
		fprintf(fp_eo, "%d", _nVars) ;
		for (int i = 0 ; i < _nVars ; i++) 
			fprintf(fp_eo, "\n%d", _VarListInElimOrder[i]) ;
		fflush(fp_eo) ;
		fclose(fp_eo) ;
		return 0 ;
	}
	inline int Initialize(int n)
	{
		if (NULL != _VarListInElimOrder) { delete [] _VarListInElimOrder ; _VarListInElimOrder = NULL ; }
		_nVars = 0 ;
		_Width = INT_MAX ;
		_WidthLowerBound = -1 ;
		_Complexity_Log10 = 0.0 ;
		_TotalNewFunctionStorageAsNumOfElements_Log10 = 0.0 ;
		_MaxSingleVarElimComplexity = DBL_MAX ;
		_nFillEdges = 0 ;
		if (n > 0) {
			_VarListInElimOrder = new int[n] ;
			if (NULL == _VarListInElimOrder) 
				return 1 ;
			_nVars = n ;
			}
		return 0 ;
	}
	int32_t SerializeTreeDecomposition(ARE::ARP & P, BucketElimination::MBEworkspace & bews, bool one_based_indexing, bool ConnectedComponents, std::string & sOutput) ;
	void Destroy(void)
	{
		_nVars = 0 ;
		_Width  = -1 ;
		_WidthLowerBound = -1 ;
		_Complexity_Log10 = -1.0 ;
		_TotalNewFunctionStorageAsNumOfElements_Log10 = -1.0 ;
		_MaxSingleVarElimComplexity = DBL_MAX ;
		_nFillEdges = 0 ;
		if (NULL != _VarListInElimOrder) 
			{ delete [] _VarListInElimOrder ; _VarListInElimOrder = NULL ; }
	}
public :
	Order(void)
		: 
		_nVars(0), 
		_VarListInElimOrder(NULL), 
		_Width(-1), 
		_WidthLowerBound(-1), 
		_Complexity_Log10(-1.0), 
		_TotalNewFunctionStorageAsNumOfElements_Log10(-1.0), 
		_MaxSingleVarElimComplexity(DBL_MAX), 
		_nFillEdges(0) 
	{
	}
	~Order(void)
	{
		if (NULL != _VarListInElimOrder) 
			delete [] _VarListInElimOrder ;
	}
} ;

class CVOcontext
{
public :
	// IN
	ARE::ARP *_Problem ;
	// OPTIONS
	ARE::VarElimOrderComp::ObjectiveToMinimize _ObjCode ;
	ARE::VarElimOrderComp::NextVarPickCriteria _AlgCode ;
	ARE::VarElimOrderComp::ObjectiveToMinimize _SecondaryObjCode ;
	int _nThreads ;
	int _nRunsToDoMin, _nRunsToDoMax ;
	int64_t _TimeLimitInMilliSeconds ;
	int _nRandomPick ;
	double _eRandomPick ;
	bool _EarlyTerminationOfBasic_W ;
	bool _EarlyTerminationOfBasic_C ;
	int _LogIncrement ;
	bool _FindPracticalVariableOrder ; // this will cut off computation when we know that the result will not be practical
	int _PracticalOrderLimit_W ; // default is 42 (2^42 = 4TB)
	double _PracticalOrderLimit_C ; // default is 13.0 (10^13 = 10TB)
	// OUT
	ARE::utils::RecursiveMutex _BestOrderMutex ;
	int _ret ;
	ARE::VarElimOrderComp::Order *_BestOrder ;
	// CONTROL
	FILE *_fpLOG ;
	unsigned long _RandomGeneratorSeed ;
#if defined WINDOWS || _WINDOWS
	LONG volatile _StopAndExit ;
	uintptr_t _ThreadHandle ;
#elif defined (LINUX)
	int64_t volatile _StopAndExit ;
	pthread_t _ThreadHandle ;
#endif 
	// INTERNALS
	// _OriginalGraph = original copy of the problem
	// _MasterGraph = global copy of the problem; used by all threads as a starting point
	ARE::Graph _OriginalGraph, _MasterGraph ;
	int64_t _tStart, _tEnd, _tToStop ;
	// AdjVar space is allocated it blocks (each size is TempAdjVarSpaceSize) and here we store ptrs to each block.
	int _TempAdjVarSpaceSizeExtraArrayN ;
	ARE::AdjVar *_TempAdjVarSpaceSizeExtraArray[TempAdjVarSpaceSizeExtraArraySize] ;
	// STATISTICS
	volatile long _nRunsStarted ;
	int _nRunsCompleted ;
	int64_t _Width2CountMap[1024] ; // for widths [0,1023], how many times it was obtained
	double _Width2MinComplexityMap[1024] ; // for widths [0,1023], log of smallest complexity
	double _Width2MaxComplexityMap[1024] ; // for widths [0,1023], log of largest complexity
	int _nImprovements ;
	ARE::VarElimOrderComp::ResultSnapShot _Improvements[1024] ;
public :
	int NoteVarOrderComputationCompletion(int w_IDX, Graph & G) ;
	int CreateCVOthread(void) ;
	int RequestStopCVOthread(void) ; // ret=0 means stopped; 1=stop requested, but still running; -1=stop requested before, but still running.
	int StopCVOthread(int64_t TimeoutInMilliseconds = 10000) ;
	int Reset(void)
	{
		_nRunsStarted = 0 ;
		_nRunsCompleted = 0 ;
		_nImprovements = 0 ;
		for (int i = 0 ; i < 1024 ; i++) {
			_Width2CountMap[i] = 0 ;
			_Width2MinComplexityMap[i] = DBL_MAX ;
			_Width2MaxComplexityMap[i] = 0 ;
			}
		return 0 ;
	}
	int Destroy(void)
	{
		StopCVOthread() ;
		for (int i = 0 ; i < _TempAdjVarSpaceSizeExtraArrayN ; i++) {
			delete [] _TempAdjVarSpaceSizeExtraArray[i] ;
			}
		_TempAdjVarSpaceSizeExtraArrayN = 0 ;
		if (NULL != _BestOrder) 
			_BestOrder->Destroy() ;
		Reset() ;
		_Problem = NULL ;
		return 0 ;
	}
public :
	CVOcontext(void) 
		:
		_Problem(NULL), 
		_ObjCode(ARE::VarElimOrderComp::Width), 
		_AlgCode(ARE::VarElimOrderComp::MinFill), 
		_SecondaryObjCode(ARE::VarElimOrderComp::None), 
		_nThreads(1), 
		_nRunsToDoMin(1), 
		_nRunsToDoMax(1), 
		_TimeLimitInMilliSeconds(3600000), 
		_nRandomPick(8), 
		_eRandomPick(0.5), 
		_EarlyTerminationOfBasic_W(true), 
		_EarlyTerminationOfBasic_C(false), 
		_LogIncrement(1), 
		_FindPracticalVariableOrder(true), 
		_PracticalOrderLimit_W(42), 
		_PracticalOrderLimit_C(13.0), 
		_ret(-1), 
		_BestOrder(NULL), 
		_fpLOG(NULL), 
		_RandomGeneratorSeed(0), 
		_StopAndExit(0), 
		_ThreadHandle(0), 
		_tStart(0), _tEnd(0), _tToStop(0), 
		_TempAdjVarSpaceSizeExtraArrayN(0), 
		_nRunsStarted(0), 
		_nRunsCompleted(0), 
		_nImprovements(0)
	{
		for (int i = 0 ; i < 1024 ; i++) {
			_Width2CountMap[i] = 0 ;
			_Width2MinComplexityMap[i] = DBL_MAX ;
			_Width2MaxComplexityMap[i] = 0 ;
			}
	}
	~CVOcontext(void)
	{
		Destroy() ;
	}
} ;

class Worker
{
public :
	CVOcontext *_CVOcontext ;
	int _IDX ;
	Graph *_G ;
public :
//	ARE::VarElimOrderComp::NextVarPickCriteria _Algorithm ; // 0=MinFill, 1=MinDegree(MinInducedWidth), 2=MinComplexity
//	int _nRandomPick ;
#if defined WINDOWS || _WINDOWS
	uintptr_t _ThreadHandle ;
#elif defined (LINUX)
	pthread_t _ThreadHandle ;
#endif 
	bool _ThreadIsRunning ;
	bool _ThreadStop ; // signal the thread to stop
	int _nRunsDone ;
	// AdjVar space is allocated it blocks (each size is TempAdjVarSpaceSize) and here we store ptrs to each block.
	int _TempAdjVarSpaceSizeExtraArrayN ;
	ARE::AdjVar *_TempAdjVarSpaceSizeExtraArray[TempAdjVarSpaceSizeExtraArraySize] ;
public :
	Worker(CVOcontext *CVOcontext = NULL, int idx = -1, Graph *G = NULL)
		:
		_CVOcontext(CVOcontext), 
		_IDX(idx), 
		_G(G), 
//		_Algorithm(ARE::VarElimOrderComp::MinFill), 
//		_nRandomPick(1), 
		_ThreadHandle(0), 
		_ThreadIsRunning(false), 
		_ThreadStop(true), 
		_nRunsDone(0), 
		_TempAdjVarSpaceSizeExtraArrayN(0)
	{
	}
	~Worker(void)
	{
		for (int i = 0 ; i < _TempAdjVarSpaceSizeExtraArrayN ; i++) {
			delete [] _TempAdjVarSpaceSizeExtraArray[i] ;
			}
		if (NULL != _G) 
			delete _G ;
	}
} ;

int Compute(
	// IN
	const std::string & ProblemInputFile, 
	ARE::VarElimOrderComp::ObjectiveToMinimize objCode,
	ARE::VarElimOrderComp::NextVarPickCriteria algCode,
	ARE::VarElimOrderComp::ObjectiveToMinimize objCodeSecondary,
	int nThreads, 
	int nRunsToDo, 
	int64_t TimeLimitInMilliSeconds, 
	int nRandomPick, 
	double eRandomPick, 
	bool PerformSingletonConsistencyChecking, 
	bool EliminateSingletonDomainVariables, 
	bool EarlyTerminationOfBasic_W, 
	bool EarlyTerminationOfBasic_C, 
	bool FindPracticalVariableOrder,
	unsigned long random_seed,
	// OUT
	Order & BestOrder, 
	CVOcontext * & CVOcontext
	) ;

inline void DeleteNewAdjVarList(int & n, ARE::AdjVar *TempAdjVarSpaceSizeExtraArray[]) 
{
	for (int i = n-1 ; i >=0 ; i--) 
		delete [] TempAdjVarSpaceSizeExtraArray[i] ;
	n = 0 ;
}

} // namespace VarElimOrderComp
} // namespace ARE

#if defined DEFINE_PACE16_MAIN_FN
int main(int argc, char* argv[]) ;
#endif

#endif // ARE_VariableOrderComputation_HXX_INCLUDED
