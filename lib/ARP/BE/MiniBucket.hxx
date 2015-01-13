#ifndef MINIBUCKET_HXX_INCLUDED
#define MINIBUCKET_HXX_INCLUDED

#include <Function.hxx>

namespace BucketElimination
{

class BEworkspace ;
class Bucket ;

class MiniBucket
{
protected :
	BucketElimination::Bucket *_B ; // the bucket this minibucket belongs to
public :
	inline BucketElimination::Bucket *Bucket(void) { return _B ; }

	// width/signature of this minibucket
protected :
	int _Width ; // cardinality of signature of this minibucket; this includes variables eliminated in this minibucket.
	int *_Signature ; // a union of the scopes of all functions in this minibucket, including variables eliminated in this minibucket.
public :
	inline int Width(void) const { return _Width ; }
	inline const int *Signature(void) const { return _Signature ; }
	int ComputeSignature(void) ;

	// variable(s) of the minibucket that are eliminated when a minibucket output fn is computed from functions of this minibucket
protected :
	int _nVars ;
	int *_Vars ;
	int _VarsSpace ;
public :
	inline int nVars(void) const { return _nVars ; }
	inline int Var(int IDX) const { return _Vars[IDX] ; }
	inline int *VarsArray(void) { return _Vars ; }
	inline int AddVar(int Var)
	{
		int i ;
		for (i = 0 ; i < _nVars ; i++) {
			if (_Vars[i] == Var) 
				return 0 ;
			}
		if (_nVars >= _VarsSpace) {
			int *v = new int[_VarsSpace + 2] ;
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

	// functions part of the original problem, assigned to this bucket
	// they functions don't belong to the bucket; normally they belong to the problem.
protected :
	int _nOriginalFunctions ;
	ARE::Function **_OriginalFunctions ;
	int _OriginalWidth ; // cardinality of the original signature of this bucket; this includes _V; if <0, then unknown, should be computed
	int *_OriginalSignature ; // a union of the scopes of all original functions, including _V
public :
	inline int nOriginalFunctions(void) const { return _nOriginalFunctions ; }
	inline ARE::Function *OriginalFunction(int IDX) const { return _OriginalFunctions[IDX] ; }
	inline ARE::Function **OriginalFunctionsArray(void) { return _OriginalFunctions ; }
	inline int OriginalWidth(void) const { return _OriginalWidth ; }
	int SetOriginalFunctions(int N, ARE::Function *FNs[]) ;
	int AddOriginalFunctions(int N, ARE::Function *FNs[]) ;

	// functions generated by other buckets (higher in the ordering), assigned to this bucket
	// these functions don't belong to this bucket; they belong to the bucket that generated them (i.e. e.g. this bucket should not delete them).
protected :
	int _nChildBucketFunctions ;
	ARE::Function **_ChildBucketFunctions ;
	int _ChildBucketFunctionArraySize ;
public :
	inline int nChildBucketFunctions(void) const { return _nChildBucketFunctions ; }
	inline ARE::Function *ChildBucketFunction(int IDX) const { return _ChildBucketFunctions[IDX] ; }
	int AddChildBucketFunction(ARE::Function & F) ;
	int RemoveChildBucketFunction(ARE::Function & F) ;

	// ****************************************************************************************
	// Function generated by this bucket; its scope is _Signature-_V.
	// This function is generated by this bucket and belongs to this bucket.
	// However, after it is generated, it is placed in the parent-bucket.
	// Note that bucket function (as any other funcion), may have its table broken into pieces (blocks).
	// If the table is small, it may be kept in memory in its entirety (as 1 block); if it is large it may be broken into pieces (blocks) 
	// and each block stored as a table on the disk.
	// ****************************************************************************************

protected :

	ARE::Function _OutputFunction ;

public :

	inline ARE::Function & OutputFunction(void) { return _OutputFunction ; }
	int ComputeOutputFunctionWithScopeWithoutTable(void) ; // this creates bucket-function of this bucket and its scope; it does not fill in the table.

	// processing complexity is the size of output table
	__int64 ComputeProcessingComplexity(void) ;

	// here we keep track of which OutputFunction blocks have been computed.
	// this allows us to 1) tell if a block needed has already been computed, 2) stop programs execution in the middle and store current state, 
	// and later load the state and continue computation where we left it off.
	// the size of the array is |scope(OutputFunction)|/8.
	// each byte is a bit-vector; compution status of a block with an idx 'idx' can be determined by checking 
	// _OutputFunctionBlockComputationResult[idx >> 3] & (1 << (idx & 7))
protected :
	int _OutputFunctionBlockComputationResultSize ;
	__int64 _nOutputFunctionBlocks ;
	unsigned char *_OutputFunctionBlockComputationResult ;
	__int64 _nOutputFunctionBlocksComputed ;
public :
	int AllocateOutputFunctionBlockComputationResult(__int64 MaxBlockSizeInNumberOfCells, int nComputingThreads) ;
	bool IsOutputFunctionBlockComputed(__int64 IDX) ;
	__int64 nOutputFunctionBlocksComputed(void) const ;
	void MarkOutputFunctionBlockComputed(__int64 IDX) ;

public :

	// note : 
	// 1) scope(bucketfuncion) = _Signature - _V.
	// 2) when this function is called, we assume that scope of _OutputFunction is already ordered wrt the parent-bucket, 
	// since _OutputFunction belongs to the parent-bucket, and it is supposed to be already sorted.
	// this function will reorder the scopes of all functions in this bucket so that 
	// 1) _V is the last variable, 
	// 2) order of other variables agrees with the order of scope(bucketfuncion).
	int ReorderFunctionScopesForExternalMemory(bool IncludeOriginalFunctions, bool IncludeNewFunctions) ;

	// output function has been computed.
	// to stuff, e.g. cleanup (release FTBs of all child buckets).
	int NoteOutputFunctionComputationCompletion(void) ;

	// serialize this bucket as XML; append it to the given string.
	virtual int SaveXMLString(const char *prefixspaces, const std::string & Dir, std::string & S) ;

	// compute output function from scratch completely.
	// this fn is usually used when regular bucket elimination is used.
	// when a more sophistifacted version of BE is used (e.g. BEEM), this fn should not be used.
	virtual int ComputeOutputFunction_1Block(void) ;
	virtual int ComputeOutputFunction_EliminateAllVars_1Block(void) ;

	// compute distribution on the first variable; this is used when a marginal distribution is required
	virtual int ComputeFirstVariableDistribution_1Block(ARE_Function_TableType *dist) ;

	// this var is used temporarily when BEworkspace is computing computation order
protected :
	Bucket *_NextInOrderComputationGenList ;
public :
	inline Bucket * & NextInOrderComputationGenList(void) { return _NextInOrderComputationGenList ; }

public :

	void Destroy(void) ;
	Bucket(void) ;
	Bucket(BEworkspace & WS, int IDX, int V) ;
	virtual ~Bucket(void) ;
} ;

} // namespace BucketElimination

#endif // MINIBUCKET_HXX_INCLUDED