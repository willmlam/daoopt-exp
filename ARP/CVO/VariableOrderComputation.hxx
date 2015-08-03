#ifndef ARE_VariableOrderComputation_HXX_INCLUDED
#define ARE_VariableOrderComputation_HXX_INCLUDED

#include <stdlib.h>
#include <stdio.h>
#include <string>

#include "Graph.hxx"

namespace ARE
{

namespace VarElimOrderComp
{

enum ObjectiveToMinimize
{
	Width, 
	StateSpaceSize
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
	INT64 _dt ; // in milliseconds
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
	double _Complexity ; // log of
	double _TotalNewFunctionStorageAsNumOfElements ; // log of
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
		_Complexity = 0.0 ;
		_TotalNewFunctionStorageAsNumOfElements = 0.0 ;
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
public :
	Order(void)
		: 
		_nVars(0), 
		_VarListInElimOrder(NULL), 
		_Width(-1), 
		_Complexity(-1.0), 
		_TotalNewFunctionStorageAsNumOfElements(-1.0), 
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
	int _nThreads ;
	int _nRunsToDoMin, _nRunsToDoMax ;
	INT64 _TimeLimitInMilliSeconds ;
	int _nRandomPick ;
	double _eRandomPick ;
	bool _EarlyTerminationOfBasic_W ;
	bool _EarlyTerminationOfBasic_C ;
	int _LogIncrement ;
	bool _FindPracticalVariableOrder ; // this will cut off computation when we know that the result will not be practical
	int _PracticalOrderLimit_W ; // default is 32 (2^32 = 4GB)
	double _PracticalOrderLimit_C ; // default is 10.0 (10^10 = 10GB)
	// OUT
	int _ret ;
	ARE::VarElimOrderComp::Order *_BestOrder ;
	// CONTROL
	FILE *_fpLOG ;
#if defined WINDOWS || _WINDOWS
	LONG volatile _StopAndExit ;
	uintptr_t _ThreadHandle ;
#elif defined (LINUX)
  INT64 volatile _StopAndExit ;
	pthread_t _ThreadHandle ;
#endif 
	// INTERNALS
	ARE::utils::RecursiveMutex _MasterGraphMutex ;
	// _OriginalGraph = original copy of the problem
	// _MasterGraph = global copy of the problem; used by all threads as a starting point
	ARE::Graph _OriginalGraph, _MasterGraph ;
	INT64 _tStart, _tEnd, _tToStop ;
	// AdjVar space is allocated it blocks (each size is TempAdjVarSpaceSize) and here we store ptrs to each block.
	int _TempAdjVarSpaceSizeExtraArrayN ;
	ARE::AdjVar *_TempAdjVarSpaceSizeExtraArray[TempAdjVarSpaceSizeExtraArraySize] ;
	// STATISTICS
	volatile long _nRunsDone ;
	int _nRunsCompleted ;
	INT64 _Width2CountMap[1024] ; // for widths [0,1023], how many times it was obtained
	double _Width2MinComplexityMap[1024] ; // for widths [0,1023], log of smallest complexity
	double _Width2MaxComplexityMap[1024] ; // for widths [0,1023], log of largest complexity
	int _nImprovements ;
	ARE::VarElimOrderComp::ResultSnapShot _Improvements[1024] ;
public :
	int NoteVarOrderComputationCompletion(int w_IDX, Graph & G) ;
	int CreateCVOthread(void) ;
	int StopCVOthread(INT64 TimeoutInMilliseconds = 10000) ;
	int Reset(void)
	{
		_nRunsDone = 0 ;
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
		_PracticalOrderLimit_W(32), 
		_PracticalOrderLimit_C(10.0), 
		_ret(-1), 
		_BestOrder(NULL), 
		_fpLOG(NULL), 
		_StopAndExit(0), 
		_ThreadHandle(0), 
		_tStart(0), _tEnd(0), _tToStop(0), 
		_TempAdjVarSpaceSizeExtraArrayN(0), 
		_nRunsDone(0), 
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
	const std::string & uaifile, 
	ARE::VarElimOrderComp::ObjectiveToMinimize objCode, // 0=width, 1=space size (complexity)
	ARE::VarElimOrderComp::NextVarPickCriteria algCode, // 0=MinFill, 1=MinDegree, 2=MinComplexity
	int nThreads, 
	int nRunsToDo, 
	INT64 TimeLimitInMilliSeconds, 
	int nRandomPick, 
	double eRandomPick, 
	bool PerformSingletonConsistencyChecking, 
	bool EarlyTerminationOfBasic_W, 
	bool EarlyTerminationOfBasic_C, 
	// OUT
	Order & BestOrder, 
	CVOcontext *CVOcontext
	) ;

} // namespace VarElimOrderComp
} // namespace ARE

#endif // ARE_VariableOrderComputation_HXX_INCLUDED
