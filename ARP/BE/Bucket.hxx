#ifndef BUCKET_HXX_INCLUDED
#define BUCKET_HXX_INCLUDED

#include <inttypes.h>

#include "Utils/MiscUtils.hxx"
#include "Function.hxx"
#include "MiniBucket.hxx"

namespace BucketElimination
{

class MBEworkspace ;
class MiniBucket ;

class Bucket
{
protected :
	MBEworkspace *_Workspace ;
	int32_t _IDX ; // index of this bucket in the workspace
	int32_t _V ; // in case the bucket is created for a specic variable in bucket elimination, this is the var. -1 otherwise (e.g. in case of superbuckets).
public :
	inline MBEworkspace *Workspace(void) const { return _Workspace ; }
	inline int32_t IDX(void) const { return _IDX ; }
	inline int32_t V(void) const { return _V ; }
	inline void SetIDX(int32_t IDX) { _IDX = IDX ; }

	// width/signature of this bucket; this includes Original/Augmented functions, but not Intermediate functions.
	// when processing a bucket (bottom-up over the bucket tree) we combine original/augmented functions, ignoring intermediate functions.
	// note that when computing the heuristic, we combine augmented/intermediate functions, ignoring original functions.
protected :
	int32_t _Width ; // cardinality of signature of this bucket; this includes variables eliminated in this bucket.
	int32_t *_Signature ; // a union of the scopes of all functions (original/augmented) in this bucket, including variables eliminated in this bucket.
public :
	inline int32_t Width(void) const { return _Width ; }
	inline const int32_t *Signature(void) const { return _Signature ; }
	inline void InvalidateSignature(void) { _Width = -1 ; if (NULL != _Signature) { delete [] _Signature ; _Signature = NULL ; }}
	int32_t ComputeSignature(void) ;

	// variable(s) of the bucket that are eliminated when a bucket output fn is computed from functions of this bucket.
	// normally this is 1 variable. in case of superbuckets, this is more.
protected :
	int32_t _nVars ;
	int32_t *_Vars ;
	int32_t _VarsSpace ;
public :
	inline int32_t nVars(void) const { return _nVars ; }
	inline int32_t Var(int32_t IDX) const { return _Vars[IDX] ; }
	inline int32_t *VarsArray(void) { return _Vars ; }
	inline int32_t AddVar(int32_t Var)
	{
		int32_t i ;
		for (i = 0 ; i < _nVars ; i++) {
			if (_Vars[i] == Var) 
				return 0 ;
			}
		if (_nVars >= _VarsSpace) {
			int32_t *v = new int32_t[_VarsSpace + 2] ;
			if (NULL == v) 
				return 1 ;
			for (i = 0 ; i < _nVars ; i++) 
				v[i] = _Vars[i] ;
			_VarsSpace += 2 ;
			delete [] _Vars ;
			_Vars = v ;
			}
		_Vars[_nVars++] = Var ;
		return 0 ;
	}

	// ****************************************************************************************
	// Bucket-tree structure.
	// ****************************************************************************************

protected :

	// parent bucket is determined by _OutputFunction.
	Bucket *_ParentBucket ;
	// root of the bucket tree this bucket belongs to; note that a BE workspace may have a bucket-forest.
	Bucket *_RootBucket ;
	// distance to the root; defined as number of edges to travel to get to the root.
	int32_t _DistanceToRoot ;
	// distance to the farthest leaf node from this bucket.
	// this is useful, e.g., when computing BEEM computation order.
	int32_t _Height ;
	// maximum number of variables in any descendant bucket; not including this bucket.
	int32_t _MaxDescendantNumVars ;

	// size of the computation of all mb output functions
	int64_t _ComputationNewFunctionSize ;

	// max ComputationNewFunctionSize() for all descendants of this bucket; it does not include this bucket.
	int64_t _MaxDescendantComputationNewFunctionSize ;

	// child buckets of this bucket
	int32_t _nChildren ;
	int32_t *_ChildVars ; // this space belongs to MBE workspace

public :

	inline Bucket *ParentBucket(void) const { return _ParentBucket ; }
	inline Bucket *RootBucket(void) const { return _RootBucket ; }
	inline void SetParentBucket(Bucket *B) { _ParentBucket = B ; }
	inline void SetRootBucket(Bucket *B) { _RootBucket = B ; }
	inline void SetDistanceToRoot(int32_t D) { _DistanceToRoot = D ; }
	inline void SetHeight(int32_t H) { _Height = H ; }
	inline void SetMaxDescendantNumVars(int32_t v) { _MaxDescendantNumVars = v ; }
	inline void SetComputationNewFunctionSize(int64_t v) { _ComputationNewFunctionSize = v ; }
	inline void SetMaxDescendantComputationNewFunctionSize(int64_t v) { _MaxDescendantComputationNewFunctionSize = v ; }
	inline int32_t DistanceToRoot(void) const { return _DistanceToRoot ; }
	inline int32_t Height(void) const { return _Height ; }
	inline int32_t MaxDescendantNumVars(void) const { return _MaxDescendantNumVars ; }
	inline int32_t MaxDescendantNumVarsEx(void) const { return _MaxDescendantNumVars > _Width ? _MaxDescendantNumVars : _Width ; }
	inline int64_t ComputationNewFunctionSize(void) const { return _ComputationNewFunctionSize ; }
	inline int64_t MaxDescendantComputationNewFunctionSize(void) const { return _MaxDescendantComputationNewFunctionSize ; }
	inline int64_t MaxDescendantComputationNewFunctionSizeEx(void) const { return _ComputationNewFunctionSize > _MaxDescendantComputationNewFunctionSize ? _ComputationNewFunctionSize : _MaxDescendantComputationNewFunctionSize ; }

	inline int32_t nChildren(void) const { return _nChildren ; }
	inline int32_t ChildVar(int32_t idx) const { return _ChildVars[idx] ; }
	inline void SetChildren(int32_t n, int32_t *C) { _nChildren = n ; _ChildVars = C ; }

	// functions part of the original problem, assigned to this bucket
	// they functions don't belong to the bucket; normally they belong to the problem.
protected :
	int32_t _nOriginalFunctions ;
	ARE::Function **_OriginalFunctions ;
	int32_t _OriginalWidth ; // cardinality of the original signature of this bucket; this includes _V; if <0, then unknown, should be computed
	// a union of the scopes of all original functions, including _V.
	// 2016-02-23 KK : these variables are in no particular order.
	int32_t *_OriginalSignature ;
public :
	inline int32_t nOriginalFunctions(void) const { return _nOriginalFunctions ; }
	inline ARE::Function *OriginalFunction(int32_t IDX) const { return _OriginalFunctions[IDX] ; }
	inline ARE::Function **OriginalFunctionsArray(void) { return _OriginalFunctions ; }
	inline int32_t OriginalWidth(void) const { return _OriginalWidth ; }
	int32_t SetOriginalFunctions(int32_t N, ARE::Function *FNs[]) ;
	int32_t AddOriginalFunctions(int32_t N, ARE::Function *FNs[]) ;

	// functions generated by other buckets (higher in the ordering; i.e. below this bucket in the bucket-tree), assigned to this bucket.
	// i.e. functions generated during (M)BE that contain _V.
	// these functions don't belong to this bucket; they belong to the (mini) bucket that generated them (i.e. e.g. this bucket should not delete them).
protected :
	int32_t _nAugmentedFunctions ;
	ARE::Function **_AugmentedFunctions ;
	int32_t _AugmentedFunctionsArraySize ;
public :
	inline int32_t nAugmentedFunctions(void) const { return _nAugmentedFunctions ; }
	inline void ResetnAugmentedFunctions(void) { _nAugmentedFunctions = 0 ; }
	inline ARE::Function *AugmentedFunction(int32_t IDX) const { return _AugmentedFunctions[IDX] ; }
	int32_t AddAugmentedFunction(ARE::Function & F) ;
	int32_t RemoveAugmentedFunction(ARE::Function & F, bool InvalidateSignature) ;

	// functions generated by other buckets (higher in the ordering; i.e. below this bucket in the bucket-tree), assigned to an ancestor of this bucket.
	// i.e. functions generated during (M)BE that do not contain _V, but that come from below and contain a variable that is ancestor of _V (in the bucket tree).
	// these functions don't belong to this bucket; they belong to the bucket that generated them (i.e. e.g. this bucket should not delete them).
protected :
	int32_t _nIntermediateFunctions ;
	ARE::Function **_IntermediateFunctions ;
	int32_t _IntermediateFunctionsArraySize ;
public :
	inline int32_t nIntermediateFunctions(void) const { return _nIntermediateFunctions ; }
	inline void ResetnIntermediateFunctions(void) { _nIntermediateFunctions = 0 ; }
	inline ARE::Function *IntermediateFunction(int32_t IDX) const { return _IntermediateFunctions[IDX] ; }
	int32_t AddIntermediateFunction(ARE::Function & F) ;
	int32_t RemoveIntermediateFunction(ARE::Function & F) ;

protected :
	std::vector<MiniBucket *> _MiniBuckets ;
public :
	inline int32_t nMiniBuckets(void) const { return _MiniBuckets.size() ; }
	inline std::vector<MiniBucket *> & MiniBuckets(void) { return _MiniBuckets ; }

public :

	// note : 
	// 1) scope(bucketfuncion) = _Signature - _V.
	// 2) when this function is called, we assume that scope of _OutputFunction is already ordered wrt the parent-bucket, 
	// since _OutputFunction belongs to the parent-bucket, and it is supposed to be already sorted.
	// this function will reorder the scopes of all functions in this bucket so that 
	// 1) _V is the last variable, 
	// 2) order of other variables agrees with the order of scope(bucketfuncion).
// 2016-02-20 KK : commented out; not using this fn right now.
//	int32_t ReorderFunctionScopesWrtParentBucket(bool IncludeOriginalFunctions, bool IncludeNewFunctions) ;

	// this creates the output function of this bucket and its scope; it does not fill in the table.
	// arguments of the buckets are sorted in the increasing order of distance to root.
	int32_t ComputeOutputFunctionWithScopeWithoutTable(
		// IN 
		int32_t * & TempSpaceForArglist, 
		int32_t TempSpaceForArglistSize,
		// OUT
		ARE::Function * & FN, 
		int32_t & fnMaxVar) ;

	// output function has been computed.
	// to stuff, e.g. cleanup (release FTBs of all child buckets).
	int32_t NoteOutputFunctionComputationCompletion(void) ;

	// processing complexity is the size of output table
	int64_t ComputeProcessingComplexity(void) ;

	// do mini-bucket partitioning for this bucket.
	int32_t CreateMBPartitioning(
		// IN
		int32_t iBound, bool CreateTables, bool doMomentMatching, 
		bool & AbandonIfActualPartitioning, 
		// IN : helper arrays; passed in so that they can be created by a calling fn and reused for all buckets
		std::vector<int32_t> & key, std::vector<int64_t> & data, std::vector<int32_t> & helperArray) ;

	// compute output functions of all minibuckets
	int32_t ComputeOutputFunctions(bool DoMomentMatching) ;

	// compute distribution on the first variable; this is used when a marginal distribution is required.
	// this eliminates all variables of this bucket, except for the given variable, combining original/augmented/intermediate functions.
	int32_t ComputeFirstVariableDistribution(
		// OUT
		ARE_Function_TableType *dist) ;

	// assuming this bucket has 1 variable, distribution on this variable.
	// values of all variables, except for the variable(s) of this bucket, are given in ContextValues.
	// this fn is used, e.g., when computing a cost-to-go for each configuration of this bucket's variable(s).
	int32_t ComputeFirstVariableDistributionEx(
		// IN : values to all context variables (= ancestors of this variable in the bucket tree)
		int32_t *ContextValues, 
		// OUT
		ARE_Function_TableType *dist) ;

	// delete all minibuckets
	int32_t DestroyPartitioning(void) ;

protected :
	// this var is used temporarily when BEworkspace is computing computation order
	Bucket *_NextInOrderComputationGenList ;
public :
	inline Bucket * & NextInOrderComputationGenList(void) { return _NextInOrderComputationGenList ; }

public :

	void Destroy(void) ;
	Bucket(void) ;
	Bucket(MBEworkspace & WS, int32_t IDX, int32_t V) ;
	virtual ~Bucket(void) ;
} ;

} // namespace BucketElimination

#endif // BUCKET_HXX_INCLUDED
