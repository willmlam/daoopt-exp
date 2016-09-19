#ifndef MBEworkspace_HXX_INCLUDED
#define MBEworkspace_HXX_INCLUDED

#include <inttypes.h>

#include "Function.hxx"
#include "Problem.hxx"
#include "Bucket.hxx"
#include "MiniBucket.hxx"
#include "Workspace.hxx"

namespace ARE { class Workspace ; }

namespace BucketElimination
{

class MBEworkspace : public ARE::Workspace
{
public:
#if defined (LINUX)
  static pthread_mutex_t stopSignalMutex;
#endif

protected :

	int32_t _nVars ; // number of variables
	int32_t *_VarOrder ; // order of variables, as an array of variable indeces; _VarOrder[0] is first bucket, _VarOrder[N-1] is last bucket.
	int32_t *_VarPos ; // positions of variables in the order; the first variable to be eliminated is at the end; the last variable to be eliminated is at [0].
	Bucket **_Var2BucketMapping ; // bucket for each variable

	// when keeping track of the children of each buckets, we need space; all the space for keeping track of all children of all buckets, 
	// takes no more than _n space. allocate it at once.
	int32_t *_BTchildlistStorage = NULL ;

	bool _DeleteUsedTables ; // delete tables used (child tables when parent bucket is computed)
	FILE *_fpLOG ;

	int32_t _InducedWidth ; // induced width of the given ordering

	int32_t _iBound ; // iBound is the max num of variables in a mini-bucket, including ones being eliminated and ones remaining; default is 1000000 = meaning infinite
	int32_t _nBucketsWithPartitioning ; // number of buckets with >1 minibuckets

	// memory bounds to run MBE
	double _MaxSpaceAllowed_Log10 ;

public :

	inline int32_t N(void) const { return _nVars ; }
	inline const int32_t *VarOrder(void) const { return _VarOrder ; }
	inline const int32_t *VarPos(void) const { return _VarPos ; }
	inline Bucket *MapVar2Bucket(int32_t v) const { return v >= 0 && v < _nVars ? _Var2BucketMapping[v] : NULL ; }

	inline int32_t InducedWidth(void) const { return _InducedWidth ; }
	inline int32_t & iBound(void) { return _iBound ; }
	inline double & MaxSpaceAllowed_Log10(void) { return _MaxSpaceAllowed_Log10 ; }
	inline int32_t nBucketsWithPartitioning(void) const { return _nBucketsWithPartitioning ; }

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

	int32_t CreateThread(void) ;
	int32_t StopThread(void) ;

	// ****************************************************************************************
	// Query.
	// ****************************************************************************************

protected :

	ARE_Function_TableType _AnswerFactor ; // it is in the same normal/log scale as the problem.
	int32_t _FnCombinationType ; // 0 = undef, 1 = product, 2 = sum
	int32_t _VarEliminationType ; // 0 = undef, 1 = sum, 2 = max, 3 = min

	// this is the answer when all variables are eliminated; it is in the same normal/log scale as the problem.
	ARE_Function_TableType _CompleteEliminationResult ;

	// this is the answer when distribution on 1 variable is computed; it is in the same normal/log scale as the problem.
	ARE_Function_TableType _MarginalSingleVariableDistribution[MAX_NUM_VALUES_PER_VAR_DOMAIN] ;
	int32_t _MarginalSingleVariableDistributionK ;
	int32_t _MarginalSingleVariableDistributionVar ;

	int64_t _tStart, _tEnd, _tToStop, _RunTimeInMilliseconds ;

public :

	inline bool ProblemIsLogScale(void) const { return _Problem->FunctionsAreConvertedToLogScale() ; }
	inline ARE_Function_TableType AnswerFactorEx(void) const { return _Problem->FunctionsAreConvertedToLogScale() ? pow(10.0, _AnswerFactor) : _AnswerFactor ; }
	inline ARE_Function_TableType AnswerFactor(void) const { return _AnswerFactor ; }
	inline ARE_Function_TableType CompleteEliminationResultEx(void) const { return _Problem->FunctionsAreConvertedToLogScale() ? pow(10.0, _CompleteEliminationResult) : _CompleteEliminationResult ; }
	inline ARE_Function_TableType CompleteEliminationResult(void) const { return _CompleteEliminationResult ; }
	inline ARE_Function_TableType MarginalSingleVariableDistributionEx(int32_t idx) const { return _Problem->FunctionsAreConvertedToLogScale() ? pow(10.0, _MarginalSingleVariableDistribution[idx]) : _MarginalSingleVariableDistribution[idx] ; }
	inline ARE_Function_TableType MarginalSingleVariableDistribution(int32_t idx) const { return _MarginalSingleVariableDistribution[idx] ; }
	inline ARE_Function_TableType *MarginalSingleVariableDistribution(void) { return _MarginalSingleVariableDistribution ; }
	inline int32_t & MarginalSingleVariableDistributionK(void) { return _MarginalSingleVariableDistributionK ; }
	inline int32_t & MarginalSingleVariableDistributionVar(void) { return _MarginalSingleVariableDistributionVar ; }
	inline int32_t FnCombinationType(void) const { return _FnCombinationType ; }
	inline int32_t VarEliminationType(void) const { return _VarEliminationType ; }
	inline void AddAnswerFactor(ARE_Function_TableType & F)
	{
		ApplyFnCombinationOperator(_AnswerFactor, F) ;
	}
	inline void ApplyFnCombinationOperator(ARE_Function_TableType & V, ARE_Function_TableType v)
	{
		if (FN_COBINATION_TYPE_PROD == _FnCombinationType) {
			if (_Problem->FunctionsAreConvertedToLogScale()) 
				V += v ;
			else 
				V *= v ;
			}
		else if (FN_COBINATION_TYPE_SUM == _FnCombinationType) {
			if (_Problem->FunctionsAreConvertedToLogScale()) 
				LOG_OF_SUM_OF_TWO_NUMBERS_GIVEN_AS_LOGS(V, V, v)
			else 
				V += v ;
			}
	}
	inline void ApplyFnDivisionOperator(ARE_Function_TableType & V, ARE_Function_TableType v)
	{
		if (FN_COBINATION_TYPE_PROD == _FnCombinationType) {
			if (_Problem->FunctionsAreConvertedToLogScale()) 
				V -= v ;
			else 
				V /= v ;
			}
		else if (FN_COBINATION_TYPE_SUM == _FnCombinationType) {
			if (_Problem->FunctionsAreConvertedToLogScale()) 
				LOG_OF_SUB_OF_TWO_NUMBERS_GIVEN_AS_LOGS(V, V, v)
			else 
				V -= v ;
			}
	}
	inline void ApplyVarEliminationOperator(ARE_Function_TableType & V, ARE_Function_TableType v)
	{
		if (VAR_ELIMINATION_TYPE_SUM == _VarEliminationType) {
			if (_Problem->FunctionsAreConvertedToLogScale()) 
				LOG_OF_SUM_OF_TWO_NUMBERS_GIVEN_AS_LOGS(V, V, v)
			else 
				V += v ;
			}
		else if (VAR_ELIMINATION_TYPE_MAX == _VarEliminationType) 
			{ if (v > V) V = v ; }
		else if (VAR_ELIMINATION_TYPE_MIN == _VarEliminationType) 
			{ if (v < V) V = v ; }
	}
	inline ARE_Function_TableType FnCombinationNeutralValue(void) const
	{
		if (FN_COBINATION_TYPE_PROD == _FnCombinationType) {
			if (_Problem->FunctionsAreConvertedToLogScale()) 
				return 0.0 ;
			return 1.0 ;
			}
		else if (FN_COBINATION_TYPE_SUM == _FnCombinationType) {
			if (_Problem->FunctionsAreConvertedToLogScale()) 
				return ARP_nInfinity ;
			return 0.0 ;
			}
		return DBL_MAX ;
	}
	inline ARE_Function_TableType VarEliminationDefaultValue(void) const
	{
		if (VAR_ELIMINATION_TYPE_SUM == _VarEliminationType) {
			if (_Problem->FunctionsAreConvertedToLogScale()) 
				return ARP_nInfinity ;
			return 0.0 ;
			}
		else if (VAR_ELIMINATION_TYPE_MAX == _VarEliminationType) {
			return ARP_nInfinity ;
			}
		else if (VAR_ELIMINATION_TYPE_MIN == _VarEliminationType) {
			return ARP_pInfinity ;
			}
		return DBL_MAX ;
	}
	inline void SetFnCombinationType(int32_t v) { _FnCombinationType = v ; }
	inline void SetVarEliminationType(int32_t v) { _VarEliminationType = v ; }
	inline void SetCompleteEliminationResult(ARE_Function_TableType v) { _CompleteEliminationResult = v ; }

	inline int64_t & tStart(void) { return _tStart ; }
	inline int64_t & tEnd(void) { return _tEnd ; }
	inline int64_t & tToStop(void) { return _tToStop ; }
	inline int64_t & RunTimeInMilliseconds(void) { return _RunTimeInMilliseconds ; }

	// ****************************************************************************************
	// Bucket-tree structure.
	// ****************************************************************************************

protected :

	int32_t _nBuckets ;
	Bucket **_Buckets ; // note that _Buckets[] array size may be smaller than _Var2BucketMapping[], but the order must be the same.

	// an ordered list of buckets to compute; typically this list is ordered in the order of decreasing height (of the bucket), 
	// and the computation is carried out from last-to-first.
	// the allocated size of this array is _nVars.
	int32_t *_BucketOrderToCompute ;

	int32_t _nBucketsWithSingleChild_initial ; // number of buckets with a single child; these buckets can be eliminated by using superbuckets.
	int32_t _nBucketsWithSingleChild_final ; // number of buckets with a single child; these buckets can be eliminated by using superbuckets.
	int32_t _nBuckets_initial ;
	int32_t _nBuckets_final ;
	int32_t _nVarsWithoutBucket ;
	int32_t _nConstValueFunctions ;

	int32_t _MaxNumChildren ; // maximum number of children of any bucket
	int32_t _MaxNumVarsInBucket ; // maximum number of variables in any bucket
	int32_t _MaxTreeHeight ; // maximum height of the bucket tree
	int32_t _MaxNumMiniBucketsPerBucket ; // maximum number of minibuckets in a bucket
	int32_t _MaxBucketFunctionWidth ; // maximum size of any bucket output function
	int32_t _nBucketsWithSingleChild ; // number of buckets with single child
	int32_t _nBucketsWithNoChildren ; // number of leaf buckets
	int32_t _nRoots ; // number of roots in the bucket tree
	int64_t _TotalOriginalFunctionSize ; // size of original functions in existance during any point during BE execution
	int64_t _TotalOriginalFunctionSpace ; // space of original functions in existance during any point during BE execution
	double _TotalNewFunctionSize_Log10 ; // return number of entries over all tables generated during BE execution
	double _TotalNewFunctionSpace_Log10 ; // return space (in bytes) over all tables generated during BE execution
	double _TotalNewFunctionComputationComplexity_Log10 ; // return bucket_width*bucket_nFNs over all buckets; this is the number of operations it takes to compute BE execution.
	double _MaxSimultaneousNewFunctionSize_Log10 ; // max size of new functions in existance during any point during BE execution; this is determined by DeleteUsedTables flag.
	double _MaxSimultaneousNewFunctionSpace_Log10 ; // max space of new functions in existance during any point during BE execution; this is determined by DeleteUsedTables flag.
	double _MaxSimultaneousTotalFunctionSize_Log10 ; // max size of all (old/new) functions in existance during any point during BE execution; this is determined by DeleteUsedTables flag.
	double _MaxSimultaneousTotalFunctionSpace_Log10 ; // max space of all (old/new) functions in existance during any point during BE execution; this is determined by DeleteUsedTables flag.

	double _TotalNewFunctionSizeComputed_Log10 ; // statistics computed during the execution

public :

	inline int32_t nBuckets(void) const { return _nBuckets ; }
	inline BucketElimination::Bucket *getBucket(int32_t IDX) const { return NULL != _Buckets ? _Buckets[IDX] : NULL ; }
	inline int32_t BucketOrderToCompute(int32_t IDX) const { return _BucketOrderToCompute[IDX] ; }
	inline int32_t Height(void)
	{
		if (_nBuckets <= 0 || NULL == _Buckets) 
			return -1 ;
		return _Buckets[0]->Height() ;
	}

	inline int32_t MaxNumChildren(void) const { return _MaxNumChildren ; }
	inline int32_t MaxNumVarsInBucket(void) const { return _MaxNumVarsInBucket ; }
	inline int32_t MaxTreeHeight(void) const { return _MaxTreeHeight ; }
	inline int32_t MaxNumMiniBucketsPerBucket(void) const { return _MaxNumMiniBucketsPerBucket ; }
	inline int32_t nBucketsWithSingleChild_initial(void) const { return _nBucketsWithSingleChild_initial ; }
	inline int32_t nBucketsWithSingleChild_final(void) const { return _nBucketsWithSingleChild_final ; }
	inline int32_t nBuckets_initial(void) const { return _nBuckets_initial ; }
	inline int32_t nBuckets_final(void) const { return _nBuckets_final ; }
	inline int32_t nVarsWithoutBucket(void) const { return _nVarsWithoutBucket ; }
	inline int32_t nConstValueFunctions(void) const { return _nConstValueFunctions ; }
	inline int32_t MaxBucketFunctionWidth(void) const { return _MaxBucketFunctionWidth ; }
	inline int32_t nBucketsWithSingleChild(void) const { return _nBucketsWithSingleChild ; }
	inline int32_t nBucketsWithNoChildren(void) const { return _nBucketsWithNoChildren ; }
	inline int32_t nRoots(void) const { return _nRoots ; }
	inline int64_t TotalOriginalFunctionSize(void) const { return _TotalOriginalFunctionSize ; }
	inline int64_t TotalOriginalFunctionSpace(void) const { return _TotalOriginalFunctionSpace ; }
	inline double TotalNewFunctionSize_Log10(void) const { return _TotalNewFunctionSize_Log10 ; }
	inline double TotalNewFunctionSpace_Log10(void) const { return _TotalNewFunctionSpace_Log10 ; }
	inline double TotalNewFunctionComputationComplexity_Log10(void) const { return _TotalNewFunctionComputationComplexity_Log10 ; }
	inline double MaxSimultaneousNewFunctionSize_Log10(void) const { return _MaxSimultaneousNewFunctionSize_Log10 ; }
	inline double MaxSimultaneousNewFunctionSpace_Log10(void) const { return _MaxSimultaneousNewFunctionSpace_Log10 ; }
	inline double MaxSimultaneousTotalFunctionSize_Log10(void) const { return _MaxSimultaneousTotalFunctionSize_Log10 ; }
	inline double MaxSimultaneousTotalFunctionSpace_Log10(void) const { return _MaxSimultaneousTotalFunctionSpace_Log10 ; }

	inline double & TotalNewFunctionSizeComputed_Log10(void) { return _TotalNewFunctionSizeComputed_Log10 ; }
	inline double GetSolutionCompletionPercentage(void) { return _TotalNewFunctionSize_Log10 >= 0 ? 100.0*pow(10.0, _TotalNewFunctionSizeComputed_Log10 - _TotalNewFunctionSize_Log10) : -1.0 ; }

	// these functions can be executed immediately, as soon as ws is initialled.
	int64_t ComputeTotalOriginalFunctionSizeAndSpace(void) ;

	// this fn can be run any time, whether MB partitioning is done or not
	int32_t ComputeMaxNumVarsInBucket(void) ;

	// Compute num of roots
	int32_t ComputeNumRoots(void)
	{
		_nRoots = 0 ;
		for (int i = 0 ; i < _nBuckets ; i++) {
			BucketElimination::Bucket *b = _Buckets[i] ;
			if (NULL == b) continue ;
			if (NULL == b->ParentBucket()) 
				_nRoots++ ;
			}
		return 0 ;
	}

	// these functions can be executed once bucket tree structure is set up; no MB partitioning is required.
	// these values are wrt original (full) bucket tree.
	int32_t ComputeMaxNumChildren(void) ;
	int32_t ComputeNBucketsWithSingleChild(void) ;

	// these functions can be executed once MB partitioning is done.
	int32_t ComputeMaxBucketFunctionWidth(void) ;
	void ComputeTotalNewFunctionSizeAndSpace(void) ;
	void ComputeTotalNewFunctionComputationComplexity(void) ;

	// **************************************************************************************************
	// Setting up workspace/bucket-tree
	// **************************************************************************************************

public :

	virtual int32_t Initialize(ARE::ARP & Problem, bool UseLogScale, const int32_t *VarOrderingAsVarList, int32_t DeleteUsedTables) ;

	// allocate buckets; assign given/input/original functions to their buckets.
	// this function does no MB partitioning!
	virtual int32_t CreateBuckets(bool KeepBTsignature, bool SimplifyBTstructure, bool CreateSuperBuckets) ;

	// apply BE for all buckets, bottom up, recursively, for which the width is 3. i.e. resulting output function is over 2 variables.
	// at the end, all leaves of the resulting BT have width > 3.
	// save this reduced (but equivalent to the original problem) in a uai format file; also new_var->old_var mapping.
	// note : this fn will destroy current MB partitioning.
	int32_t SaveReducedProblem(int32_t & nNewVariables, std::vector<int32_t> & Old2NewVarMap, std::vector<ARE::Function*> & ReducedProblemFunctions) ;

	// find largest i-bound so that new function space is within given space limit
	int32_t FindIBoundForSpaceAllowed(int32_t MinIBound, int32_t & bestIboundFound, double & NewFnSpaceUsed_Log10, int32_t & nBucketsPartitioned) ;

	// compute output functions of all buckets
	int32_t ComputeOutputFunctions(bool DoMomentMatching) ;

	// delete all MB generated tables
	int32_t DeleteMBEgeneratedTables(void) ;

	// create an order in which buckets should be processed
	// algorithm=0 means from leaves to root in terms of uniform height, i.e. one level must be finished before starting next level. this is default.
	// algorithm=1 means in the order that minimizes space, e.g. finish computing one bucket before starting next sibling.
	// MB partitioning should be done when this fn is called.
	virtual int32_t CreateComputationOrder(int32_t algorithm) ;

	// create/destroy MB partitionin.
	int32_t CreateMBPartitioning(bool CreateTables, bool doMomentMatching, int32_t ComputeComputationOrder) ;
	int32_t DestroyMBPartitioning(void) ;

	int32_t CheckBucketTreeIntegrity(void) ;
	
	// assuming entire MBE elimination is completed, compute some data
	int32_t PostComputationProcessing(void) ;

	// **************************************************************************************************
	// Executing (M)BE.
	// **************************************************************************************************

public :

	// simulate execution and compute the minimum amount of size/space required (for new functions).
	// this depends on the DeleteUsedTables flag.
	void SimulateComputationAndComputeMinSpace(bool IgnoreInMemoryTables) ;

	// run regular (M)BE on the bucket-tree
	virtual int32_t RunSimple(void) ;

	// extract solution assignment from current MBE execution; operator (min/max) will be obtained from the problem.
	// solution will be stored in the problem.
	int32_t BuildSolution(void) ;

	// **************************************************************************************************
	// Problem generation
	// **************************************************************************************************

public :

	// this function generates a random uniform Bayesian network for the given parameters.
	// it will not fill in any function tables.
	int32_t GenerateRandomBayesianNetworkStructure(
		int32_t N, // # of variables
		int32_t K, // same domain size for all variables
		int32_t P, // # of parents per CPT
		int32_t C,  // # of CPTs; variables that are not a child in a CPT will get a prior.
		int32_t ProblemCharacteristic // 0 = totally random, 1 = 1 leaf node
		) ;

public :

	virtual int32_t Destroy(void) ;
	virtual int32_t DestroyBucketPartitioning(void) ;

	MBEworkspace(const char *BEEMDiskSpaceDirectory = NULL) ;
	virtual ~MBEworkspace(void)
	{
		Destroy() ;
	}
} ;

// this function generates given number of random uniform Bayesian networks for the given parameters, 
// within specified complexity bounds.
// problems will be saved in uai format.
// return value is number of problem generated/saved.
int32_t GenerateRandomBayesianNetworksWithGivenComplexity(
	int32_t nProblems, 
	int32_t N, // # of variables
	int32_t K, // same domain size for all variables
	int32_t P, // # of parents per CPT
	int32_t C,  // # of CPTs; variables that are not a child in a CPT will get a prior.
	int32_t ProblemCharacteristic, // 0 = totally random, 1 = 1 leaf node
	int64_t MinSpace, 
	int64_t MaxSpace
	) ;

} // namespace BucketElimination

#endif // MBEworkspace_HXX_INCLUDED
