#ifndef BEworkspace_HXX_INCLUDED
#define BEworkspace_HXX_INCLUDED

#include "Function.hxx"
#include "Problem.hxx"
#include "Bucket.hxx"
#include "Workspace.hxx"

namespace ARE { class Workspace ; }

namespace BucketElimination
{

class BEworkspace : public ARE::Workspace
{
public:
#if defined (LINUX)
  static pthread_mutex_t stopSignalMutex;
#endif

protected :

	int _nVars ; // number of variables
	int *_VarOrder ; // order of variables, as an array of variable indeces; _VarOrder[0] is first bucket, _VarOrder[N-1] is last bucket.
	int *_VarPos ; // positions of variables in the order; the first variable to be eliminated is at the end; the last variable to be eliminated is at [0].
	Bucket **_Var2BucketMapping ; // bucket for each variable

	bool _DeleteUsedTables ; // delete tables used (child tables when parent bucket is computed)
	FILE *_fpLOG ;

public :

	inline int N(void) const { return _nVars ; }
	inline const int *VarOrder(void) const { return _VarOrder ; }
	inline const int *VarPos(void) const { return _VarPos ; }
	inline Bucket *MapVar2Bucket(int v) const { return _Var2BucketMapping[v] ; }

	inline bool DeleteUsedTables(void) const { return _DeleteUsedTables ; }
	inline void SetDeleteUsedTables(bool v) { _DeleteUsedTables = v ; }

	inline FILE * & logFile(void) { return _fpLOG ; }

public :

	// stop/exit computation thread
	long volatile _StopAndExit ;

#if defined WINDOWS || _WINDOWS
	uintptr_t _ThreadHandle ;
#elif defined (LINUX)
	pthread_t _ThreadHandle ;
#endif 

	int CreateThread(void) ;
	int StopThread(void) ;

	// ****************************************************************************************
	// Query.
	// ****************************************************************************************

protected :

	ARE_Function_TableType _AnswerFactor ;
	int _FnCombinationType ; // 0 = undef, 1 = product, 2 = sum
	int _VarEliminationType ; // 0 = undef, 1 = sum, 2 = max, 3 = min

	// this is the answer when all variables are eliminated
	ARE_Function_TableType _CompleteEliminationResult ;

	// this is the answer when distribution on 1 variable is computed
	ARE_Function_TableType _MarginalSingleVariableDistribution[MAX_NUM_VALUES_PER_VAR_DOMAIN] ;
	int _MarginalSingleVariableDistributionK ;
	int _MarginalSingleVariableDistributionVar ;

	INT64 _tStart, _tEnd, _tToStop, _RunTimeInMilliseconds ;

public :

	inline ARE_Function_TableType AnswerFactor(void) const { return _AnswerFactor ; }
	inline ARE_Function_TableType CompleteEliminationResult(void) const { return _CompleteEliminationResult ; }
	inline ARE_Function_TableType MarginalSingleVariableDistribution(int idx) const { return _MarginalSingleVariableDistribution[idx] ; }
	inline ARE_Function_TableType *MarginalSingleVariableDistribution(void) { return _MarginalSingleVariableDistribution ; }
	inline int & MarginalSingleVariableDistributionK(void) { return _MarginalSingleVariableDistributionK ; }
	inline int & MarginalSingleVariableDistributionVar(void) { return _MarginalSingleVariableDistributionVar ; }
	inline int FnCombinationType(void) const { return _FnCombinationType ; }
	inline int VarEliminationType(void) const { return _VarEliminationType ; }
	inline void AddAnswerFactor(ARE_Function_TableType & F)
	{
		if (FN_COBINATION_TYPE_PROD == _FnCombinationType) 
			_AnswerFactor *= F ;
		else if (FN_COBINATION_TYPE_SUM == _FnCombinationType) 
			_AnswerFactor += F ;
	}
	inline void ApplyFnCombinationOperator(ARE_Function_TableType & V, ARE_Function_TableType v)
	{
		if (FN_COBINATION_TYPE_PROD == _FnCombinationType) 
			V *= v ;
		else if (FN_COBINATION_TYPE_SUM == _FnCombinationType) 
			V += v ;
	}
	inline void ApplyVarEliminationOperator(ARE_Function_TableType & V, ARE_Function_TableType v)
	{
		if (VAR_ELIMINATION_TYPE_SUM == _VarEliminationType) 
			V += v ;
		else if (VAR_ELIMINATION_TYPE_MAX == _VarEliminationType) 
			{ if (v > V) V = v ; }
		else if (VAR_ELIMINATION_TYPE_MIN == _VarEliminationType) 
			{ if (v < V) V = v ; }
	}
	inline ARE_Function_TableType FnCombinationNeutralValue(void) const
	{
		if (FN_COBINATION_TYPE_PROD == _FnCombinationType) 
			return 1.0 ;
		else if (FN_COBINATION_TYPE_SUM == _FnCombinationType) 
			return 0.0 ;
		return DBL_MAX ;
	}
	inline ARE_Function_TableType VarEliminationDefaultValue(void) const
	{
		if (VAR_ELIMINATION_TYPE_SUM == _VarEliminationType) 
			return 0.0 ;
		else if (VAR_ELIMINATION_TYPE_MAX == _VarEliminationType) 
			return -DBL_MAX ;
		else if (VAR_ELIMINATION_TYPE_MIN == _VarEliminationType) 
			return DBL_MAX ;
		return DBL_MAX ;
	}
	inline void SetFnCombinationType(int v) { _FnCombinationType = v ; }
	inline void SetVarEliminationType(int v) { _VarEliminationType = v ; }
	inline void SetCompleteEliminationResult(ARE_Function_TableType v) { _CompleteEliminationResult = v ; }

	inline INT64 & tStart(void) { return _tStart ; }
	inline INT64 & tEnd(void) { return _tEnd ; }
	inline INT64 & tToStop(void) { return _tToStop ; }
	inline INT64 & RunTimeInMilliseconds(void) { return _RunTimeInMilliseconds ; }

	// ****************************************************************************************
	// Bucket-tree structure.
	// ****************************************************************************************

	int _nBuckets ;
	Bucket *_Buckets[MAX_NUM_BUCKETS] ;

	// an ordered list of buckets to compute; typically this list is ordered in the order of decreasing height (of the bucket), 
	// and the computation is carried out from last-to-first.
	// the allocated size of this array is _nVars.
	long *_BucketOrderToCompute ;

	int _nBucketsWithSingleChild_initial ; // number of buckets with a single child; these buckets can be eliminated by using superbuckets.
	int _nBucketsWithSingleChild_final ; // number of buckets with a single child; these buckets can be eliminated by using superbuckets.
	int _nBuckets_initial ;
	int _nBuckets_final ;
	int _nVarsWithoutBucket ;
	int _nConstValueFunctions ;

	int _MaxNumChildren ; // maximum number of children of any bucket
	int _MaxNumVarsInBucket ; // maximum number of variables in any bucket
	int _MaxBucketFunctionWidth ; // maximum size of any bucket output function
	int _nBucketsWithSingleChild ; // number of buckets with single child
	__int64 _TotalOriginalFunctionSize ; // size of original functions in existance during any point during BE execution
	__int64 _TotalOriginalFunctionSpace ; // space of original functions in existance during any point during BE execution
	__int64 _TotalNewFunctionSize ; // return number of entries over all tables generated during BE execution
	__int64 _TotalNewFunctionSpace ; // return space (in bytes) over all tables generated during BE execution
	__int64 _TotalNewFunctionComputationComplexity ; // return bucket_width*bucket_nFNs over all buckets; this is the number of operations it takes to compute BE execution.
	__int64 _MaxSimultaneousNewFunctionSize ; // max size of new functions in existance during any point during BE execution; this is determined by DeleteUsedTables flag.
	__int64 _MaxSimultaneousNewFunctionSpace ; // max space of new functions in existance during any point during BE execution; this is determined by DeleteUsedTables flag.
	__int64 _MaxTotalFunctionSize ; // max size of new functions in existance during any point during BE execution; this is determined by DeleteUsedTables flag.
	__int64 _MaxTotalFunctionSpace ; // max space of new functions in existance during any point during BE execution; this is determined by DeleteUsedTables flag.

	__int64 _TotalNewFunctionSizeComputed ; // statistics computed during the execution

public :

	inline int nBuckets(void) const { return _nBuckets ; }
	inline BucketElimination::Bucket *getBucket(int IDX) const { return _Buckets[IDX] ; }
	inline int BucketInComputationOrder(int IDX) const { return _BucketOrderToCompute[IDX] ; }

	int CheckBucketTreeIntegrity(void) ;

	inline int MaxNumChildren(void) const { return _MaxNumChildren ; }
	inline int MaxNumVarsInBucket(void) const { return _MaxNumVarsInBucket ; }
	inline int nBucketsWithSingleChild_initial(void) const { return _nBucketsWithSingleChild_initial ; }
	inline int nBucketsWithSingleChild_final(void) const { return _nBucketsWithSingleChild_final ; }
	inline int nBuckets_initial(void) const { return _nBuckets_initial ; }
	inline int nBuckets_final(void) const { return _nBuckets_final ; }
	inline int nVarsWithoutBucket(void) const { return _nVarsWithoutBucket ; }
	inline int nConstValueFunctions(void) const { return _nConstValueFunctions ; }
	inline int MaxBucketFunctionWidth(void) const { return _MaxBucketFunctionWidth ; }
	inline int nBucketsWithSingleChild(void) const { return _nBucketsWithSingleChild ; }
	inline __int64 TotalOriginalFunctionSize(void) const { return _TotalOriginalFunctionSize ; }
	inline __int64 TotalOriginalFunctionSpace(void) const { return _TotalOriginalFunctionSpace ; }
	inline __int64 TotalNewFunctionSize(void) const { return _TotalNewFunctionSize ; }
	inline __int64 TotalNewFunctionSpace(void) const { return _TotalNewFunctionSpace ; }
	inline __int64 TotalNewFunctionComputationComplexity(void) const { return _TotalNewFunctionComputationComplexity ; }
	inline __int64 MaxSimultaneousNewFunctionSize(void) const { return _MaxSimultaneousNewFunctionSize ; }
	inline __int64 MaxSimultaneousNewFunctionSpace(void) const { return _MaxSimultaneousNewFunctionSpace ; }
	inline __int64 MaxTotalFunctionSize(void) const { return _MaxTotalFunctionSize ; }
	inline __int64 MaxTotalFunctionSpace(void) const { return _MaxTotalFunctionSpace ; }

	inline __int64 & TotalNewFunctionSizeComputed(void) { return _TotalNewFunctionSizeComputed ; }
	inline double GetSolutionCompletionPercentage(void) { return _TotalNewFunctionSize > 0 ? (100.0*((double) _TotalNewFunctionSizeComputed))/_TotalNewFunctionSize : -1.0 ; }

	int ComputeMaxNumChildren(void) ;
	int ComputeMaxNumVarsInBucket(void) ;
	int ComputeNBucketsWithSingleChild(void) ;
	int ComputeMaxBucketFunctionWidth(void) ;
	__int64 ComputeTotalOriginalFunctionSizeAndSpace(void) ;
	__int64 ComputeTotalNewFunctionSizeAndSpace(void) ;
	__int64 ComputeTotalNewFunctionComputationComplexity(void) ;
	// simulate execution and compute the minimum amount of size/space required (for new functions).
	// this depends on the DeleteUsedTables flag.
	__int64 SimulateComputationAndComputeMinSpace(bool IgnoreInMemoryTables) ;

	virtual int CreateBuckets(bool CreateSuperBuckets) ;

	// create an order in which buckets should be processed
	// algorithm=0 means from leaves to root in terms of uniform height, i.e. one level must be finished before starting next level. this is default.
	// algorithm=1 means in the order that minimizes space, i.e. finish computing one bucket before starting next sibling.
	virtual int CreateComputationOrder(int algorithm) ;

public :

	// run regular BE on the bucket-tree
	virtual int RunSimple(void) ;

	// **************************************************************************************************
	// Problem generation
	// **************************************************************************************************

public :

	// this function generates a random uniform Bayesian network for the given parameters.
	// it will not fill in any function tables.
	int GenerateRandomBayesianNetworkStructure(
		int N, // # of variables
		int K, // same domain size for all variables
		int P, // # of parents per CPT
		int C,  // # of CPTs; variables that are not a child in a CPT will get a prior.
		int ProblemCharacteristic // 0 = totally random, 1 = 1 leaf node
		) ;

public :

	virtual int Initialize(ARE::ARP & Problem, const int *VarOrderingAsVarList, int DeleteUsedTables) ;
	virtual int Destroy(void) ;
	BEworkspace(const char *BEEMDiskSpaceDirectory) ;
	virtual ~BEworkspace(void)
	{
		Destroy() ;
	}
} ;

// this function generates given number of random uniform Bayesian networks for the given parameters, 
// within specified complexity bounds.
// problems will be saved in uai format.
// return value is number of problem generated/saved.
int GenerateRandomBayesianNetworksWithGivenComplexity(
	int nProblems, 
	int N, // # of variables
	int K, // same domain size for all variables
	int P, // # of parents per CPT
	int C,  // # of CPTs; variables that are not a child in a CPT will get a prior.
	int ProblemCharacteristic, // 0 = totally random, 1 = 1 leaf node
	__int64 MinSpace, 
	__int64 MaxSpace
	) ;

} // namespace BucketElimination

#endif // BEworkspace_HXX_INCLUDED
