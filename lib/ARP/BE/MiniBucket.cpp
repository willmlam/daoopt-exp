#include <stdlib.h>
#include <memory.h>

#include <Function.hxx>
#include <Bucket.hxx>
#include <BEworkspace.hxx>


BucketElimination::Bucket::Bucket(void)
	:
	_Workspace(NULL),
	_IDX(-1), 
	_nVars(0), 
	_Vars(NULL), 
	_VarsSpace(0), 
	_ParentBucket(NULL), 
	_RootBucket(NULL), 
	_DistanceToRoot(-1), 
	_Height(-1), 
	_MaxDescendantNumVars(-1), 
	_ComputationNewFunctionSize(-1), 
	_MaxDescendantComputationNewFunctionSize(-1), 
	_nOriginalFunctions(0), 
	_OriginalFunctions(NULL), 
	_OriginalWidth(-1), 
	_OriginalSignature(NULL), 
	_nChildBucketFunctions(0), 
	_ChildBucketFunctions(NULL), 
	_ChildBucketFunctionArraySize(0), 
	_Width(-1), 
	_Signature(NULL), 
	_OutputFunctionBlockComputationResultSize(0), 
	_nOutputFunctionBlocks(-1), 
	_OutputFunctionBlockComputationResult(NULL), 
	_nOutputFunctionBlocksComputed(0)
{
}


BucketElimination::Bucket::~Bucket(void)
{
	Destroy() ;
}


BucketElimination::Bucket::Bucket(BEworkspace & WS, int IDX, int V)
	:
	_Workspace(&WS),
	_IDX(IDX), 
	_nVars(0), 
	_Vars(NULL), 
	_VarsSpace(0), 
	_ParentBucket(NULL), 
	_RootBucket(NULL), 
	_DistanceToRoot(-1), 
	_Height(-1), 
	_MaxDescendantNumVars(-1), 
	_ComputationNewFunctionSize(-1), 
	_MaxDescendantComputationNewFunctionSize(-1), 
	_nOriginalFunctions(0), 
	_OriginalFunctions(NULL), 
	_OriginalWidth(-1), 
	_OriginalSignature(NULL), 
	_nChildBucketFunctions(0), 
	_ChildBucketFunctions(NULL), 
	_ChildBucketFunctionArraySize(0), 
	_Width(-1), 
	_Signature(NULL), 
	_OutputFunctionBlockComputationResultSize(0), 
	_nOutputFunctionBlocks(-1), 
	_OutputFunctionBlockComputationResult(NULL), 
	_nOutputFunctionBlocksComputed(0)
{
	// if Var is given, add it
	if (V >= 0) 
		AddVar(V) ;
	// fix up output function of the bucket
	_OutputFunction.Initialize(_Workspace, _Workspace->Problem(), -1) ;
	_OutputFunction.SetOriginatingBucket(this) ;
}


void BucketElimination::Bucket::Destroy(void)
{
	if (NULL != _OriginalFunctions) {
		delete [] _OriginalFunctions ;
		_OriginalFunctions = NULL ;
		}
	if (NULL != _OriginalSignature) {
		delete [] _OriginalSignature ;
		_OriginalSignature = NULL ;
		}
	if (NULL != _ChildBucketFunctions) {
		delete [] _ChildBucketFunctions ;
		_ChildBucketFunctions = NULL ;
		_ChildBucketFunctionArraySize = 0 ;
		}
	if (NULL != _Signature) {
		delete [] _Signature ;
		_Signature = NULL ;
		}
	_OutputFunction.Destroy() ;
	if (NULL != _OutputFunctionBlockComputationResult) {
		delete [] _OutputFunctionBlockComputationResult ;
		_OutputFunctionBlockComputationResult = NULL ;
		}
	_OutputFunctionBlockComputationResultSize = 0 ;
	_Width = -1 ;
	_nChildBucketFunctions = 0 ;
	_OriginalWidth = -1 ;
	_nOriginalFunctions = 0 ;
	_nVars = 0 ;
	if (NULL != _Vars) {
		delete [] _Vars ;
		_Vars = NULL ;
		}
	_VarsSpace = 0 ;
	_MaxDescendantNumVars = -1 ;
	_ComputationNewFunctionSize = _MaxDescendantComputationNewFunctionSize = -1 ;
}


int BucketElimination::Bucket::SaveXMLString(const char *prefixspaces, const std::string & Dir, std::string & S)
{
	char s[1024] ;
	std::string temp ;
	int i, j ;
	sprintf(s, "%s<bucket IDX=\"%d\" nVars=\"%d\" Vars=\"", prefixspaces, _IDX, _nVars) ;
	S += s ;
	for (i = 0 ; i < _nVars ; i++) {
		sprintf(s, "%d", _Vars[i]) ;
		if (i > 0) 
			S += ';' ;
		S += s ;
		}
	S += '"' ;
	if (NULL !=_ParentBucket) {
		sprintf(s, " parentbucket=\"%d\"", _ParentBucket->IDX()) ;
		S += s ;
		}
	S += ">" ;

	// save original functions
	if (_nOriginalFunctions > 0 && NULL != _OriginalFunctions) {
		sprintf(s, "\n%s <originalfunctions n=\"%d\" list=\"", prefixspaces, _nOriginalFunctions) ;
		S += s ;
		for (i = 0 ; i < _nOriginalFunctions ; i++) {
			ARE::Function *f = _OriginalFunctions[i] ;
			if (NULL == f) continue ;
			sprintf(s, "%d", f->IDX()) ;
			if (i > 0) 
				S += ';' ;
			S += s ;
			}
		S += "\"/>" ;
		}

	// save incoming child bucket functions
	if (_nChildBucketFunctions > 0 && NULL != _ChildBucketFunctions) {
		sprintf(s, "%s ", prefixspaces) ;
		for (i = 0 ; i < _nChildBucketFunctions ; i++) {
			ARE::Function *f = _ChildBucketFunctions[i] ;
			if (NULL == f) continue ;
			S += "\n" ;
			// bucketfunctions don't have IDX; set the IDX of the originating bucket as the IDX of this function
			int idx = f->IDX() ;
			Bucket *b = f->OriginatingBucket() ;
			if (NULL != b) 
				f->SetIDX(b->IDX()) ;
			f->SaveXMLString(s, "childbucketfn", Dir, S) ;
			f->SetIDX(idx) ;
			}
		}

	// serialize _OutputFunction
	if (_OutputFunction.N() > 0) {
		sprintf(s, "%s ", prefixspaces) ;
		S += "\n" ;
		_OutputFunction.SaveXMLString(s, "ownbucketfn", Dir, S) ;
		}

	// serialize _OutputFunctionBlockComputationResult
	if (_OutputFunction.N() > 0 && NULL != _OutputFunctionBlockComputationResult) {
		__int64 nTB = _OutputFunction.nTableBlocks() ;
		int size = (7 + nTB) >> 3 ;
		sprintf(s, "\n%s <ownbucketfncomputationresult nComputed=\"%I64d/%I64d\" bits=\"", prefixspaces, _nOutputFunctionBlocksComputed, nTB) ;
		S += s ;
		int n = 0 ;
		for (i = 0 ; i < size && n < nTB ; i++) {
			for (j = 0 ; j < 8 && n < nTB ; j++, n++) {
				int bit = _OutputFunctionBlockComputationResult[i] & (1 << j) ;
				if (0 != bit) 
					S += '1' ;
				else 
					S += '0' ;
				}
			}
		S += "\"/>" ;
		}

	sprintf(s, "\n%s</bucket>", prefixspaces) ;
	S += s ;

	return 0 ;
}


int BucketElimination::Bucket::SetOriginalFunctions(int N, ARE::Function *FNs[]) 
{
	if (NULL != _OriginalFunctions) {
		delete [] _OriginalFunctions ;
		_OriginalFunctions = NULL ;
		}
	if (NULL != _OriginalSignature) {
		delete [] _OriginalSignature ;
		_OriginalSignature = NULL ;
		}
	_OriginalWidth = 0 ;

	if (N < 1) 
		return 0 ;

	int i, j, k ;

	_OriginalFunctions = new ARE::Function*[N] ;
	if (NULL == _OriginalFunctions) 
		{ Destroy() ; return 1 ; }
	_nOriginalFunctions = N ;
	for (i = 0 ; i < _nOriginalFunctions ; i++) 
		_OriginalFunctions[i] = FNs[i] ;

	// compute approx width
	if (0 == _nOriginalFunctions) 
		return 0 ;
	int n = 0 ;
	for (i = 0 ; i < _nOriginalFunctions ; i++) 
		n += _OriginalFunctions[i]->N() ;
	// n is an upper bound on the width

	// compute width/signature
	_OriginalSignature = new int[n] ;
	if (NULL == _OriginalSignature) 
		goto failed ;
	_OriginalWidth = 0 ;
	for (i = 0 ; i < _nOriginalFunctions ; i++) {
		ARE::Function *f = _OriginalFunctions[i] ;
		f->SetBEBucket(this) ;
		for (j = 0 ; j < f->N() ; j++) {
			int v = f->Argument(j) ;
			for (k = 0 ; k < _OriginalWidth ; k++) 
				{ if (_OriginalSignature[k] == v) break ; }
			if (k < _OriginalWidth) 
				continue ;
			_OriginalSignature[_OriginalWidth++] = v ;
			}
		}

	return 0 ;
failed :
	Destroy() ;
	return 1 ;
}


int BucketElimination::Bucket::AddOriginalFunctions(int N, ARE::Function *FNs[]) 
{
	if (N < 1) 
		return 0 ;

	int i, j, k ;

	// check if functions in FNs are already in this bucket
	for (i = N-1 ; i >= 0 ; i--) {
		ARE::Function *f = FNs[i] ;
		for (j = 0 ; j < _nOriginalFunctions ; j++) {
			if (_OriginalFunctions[j] == f) {
				FNs[i] = FNs[--N] ;
				break ;
				}
			}
		}
	if (N < 1) 
		return 0 ;

	// reallocate original functions array
	int space = _nOriginalFunctions + N ;
	ARE::Function **fnlist = new ARE::Function*[space] ;
	if (NULL == fnlist) 
		return 1 ;
	for (i = 0 ; i < _nOriginalFunctions ; i++) 
		fnlist[i] = _OriginalFunctions[i] ;
	for (; i < space ; i++) 
		fnlist[i] = FNs[i-_nOriginalFunctions] ;
	delete [] _OriginalFunctions ;
	_OriginalFunctions = fnlist ;
	_nOriginalFunctions = space ;

	// compute approx width
	int n = 0 ;
	for (i = 0 ; i < _nOriginalFunctions ; i++) 
		n += _OriginalFunctions[i]->N() ;
	// n is an upper bound on the width

	// compute width/signature
	if (NULL != _OriginalSignature) {
		delete [] _OriginalSignature ;
		_OriginalSignature = NULL ;
		}
	_OriginalSignature = new int[n] ;
	if (NULL == _OriginalSignature) 
		goto failed ;
	_OriginalWidth = 0 ;
	for (i = 0 ; i < _nOriginalFunctions ; i++) {
		ARE::Function *f = _OriginalFunctions[i] ;
		f->SetBEBucket(this) ;
		for (j = 0 ; j < f->N() ; j++) {
			int v = f->Argument(j) ;
			for (k = 0 ; k < _OriginalWidth ; k++) 
				{ if (_OriginalSignature[k] == v) break ; }
			if (k < _OriginalWidth) 
				continue ;
			_OriginalSignature[_OriginalWidth++] = v ;
			}
		}

	return 0 ;
failed :
	Destroy() ;
	return 1 ;
}


int BucketElimination::Bucket::AddChildBucketFunction(ARE::Function & F)
{
	// check if we have enough space
	if (_nChildBucketFunctions+1 > _ChildBucketFunctionArraySize) {
		int newsize = _ChildBucketFunctionArraySize + 8 ;
		ARE::Function **newspace = new ARE::Function*[newsize] ;
		if (NULL == newspace) 
			return 1 ;
		if (_nChildBucketFunctions > 0) 
			memcpy(newspace, _ChildBucketFunctions, sizeof(ARE::Function *)*_nChildBucketFunctions) ;
		delete []_ChildBucketFunctions ;
		_ChildBucketFunctions = newspace ;
		_ChildBucketFunctionArraySize = newsize ;
		}

	_ChildBucketFunctions[_nChildBucketFunctions++] = &F ;

	return 0 ;
}


int BucketElimination::Bucket::RemoveChildBucketFunction(ARE::Function & F)
{
	int i ;
	for (i = _nChildBucketFunctions - 1 ; i >= 0 ; i--) {
		if (&F == _ChildBucketFunctions[i]) {
			_ChildBucketFunctions[i] = _ChildBucketFunctions[--_nChildBucketFunctions] ;
			F.SetBEBucket(NULL) ;
			}
		}
	return 0 ;
}


int BucketElimination::Bucket::ComputeSignature(void)
{
	if (NULL != _Signature) {
		delete [] _Signature ;
		_Signature = NULL ;
		}
	_Width = 0 ;

	// compute approx width
	if (_OriginalWidth < 0) {
		_Width = -1 ;
		return 1 ;
		}
	if (0 == _OriginalWidth && 0 == _nChildBucketFunctions) 
		return 0 ;
	int i, n = _OriginalWidth ;
	for (i = 0 ; i < _nChildBucketFunctions ; i++) 
		n += _ChildBucketFunctions[i]->N() ;
	// n is an upper bound on the width

	// OriginalSignature is part of Signature
	_Signature = new int[n] ;
	if (NULL == _Signature) 
		goto failed ;
	if (_OriginalWidth > 0) {
		memcpy(_Signature, _OriginalSignature, _OriginalWidth*sizeof(int)) ;
		_Width = _OriginalWidth ;
		}

	// add scopes of bucketfunctions to the signature
	int j, k ;
	for (i = 0 ; i < _nChildBucketFunctions ; i++) {
		ARE::Function *f = _ChildBucketFunctions[i] ;
		for (j = 0 ; j < f->N() ; j++) {
			int v = f->Argument(j) ;
			for (k = 0 ; k < _Width ; k++) 
				{ if (_Signature[k] == v) break ; }
			if (k < _Width) 
				continue ;
			_Signature[_Width++] = v ;
			}
		}

	return 0 ;
failed :
	return 1 ;
}


__int64 BucketElimination::Bucket::ComputeProcessingComplexity(void)
{
	__int64 n = 1 ;
	for (int i = 0 ; i < _Width ; i++) {
		n *= _Workspace->Problem()->K(_Signature[i]) ;
		}
	return n ;
}


int BucketElimination::Bucket::ComputeOutputFunctionWithScopeWithoutTable(void)
{
	INT64 newFNsize;
	// dump current bucket fn; cleanup.
	if (NULL != _ParentBucket) {
		_ParentBucket->RemoveChildBucketFunction(_OutputFunction) ;
		_ParentBucket = NULL ;
		}
	_OutputFunction.Destroy() ;

	if (_Width < 0) {
		if (0 != ComputeSignature()) 
			return 1 ;
		}
	if (_Width <= _nVars) 
		// _Width=0 can only be when this bucket has no functions
		// _Width=_nVars means all functions in this bucket have only _Vars in their scope; this bucket has a function with const-value; in this case
		// we will still compute the const-value, but the output function does not go anywhere.
		return 0 ;
	if (_Width > MAX_NUM_VARIABLES_PER_BUCKET) 
		return 1 ;

	// make a local copy, while ignoring bucket variables
	int vlist[MAX_NUM_VARIABLES_PER_BUCKET] ;
	int i, j, n = 0 ;
	for (i = 0 ; i < _Width ; i++) {
		int v = _Signature[i] ;
		for (j = 0 ; j < _nVars ; j++) {
			if (_Vars[j] == v) 
				break ;
			}
		if (j < _nVars) 
			continue ;
		vlist[n++] = v ;
		}

	// create and initialize bucket function
	if (0 != _OutputFunction.SetArguments(n, vlist)) 
		goto failed ;
	_OutputFunction.ComputeTableSize() ;
	if (_nOriginalFunctions > 0) 
		_OutputFunction.SetType(_OriginalFunctions[0]->Type()) ;
	else if (_nChildBucketFunctions > 0) 
		_OutputFunction.SetType(_ChildBucketFunctions[0]->Type()) ;

	// find appropriate parent bucket and assign the bucket function to it
	{
	const int *varpos = _Workspace->VarPos() ;
	int v = _OutputFunction.GetHighestOrderedVariable(varpos) ;
	Bucket *parentbucket = _Workspace->MapVar2Bucket(v) ;
	if (NULL == parentbucket) 
		// this is not supposed to happen
		goto failed ;
	if (parentbucket->IDX() >= _IDX) 
		// this is not supposed to happen
		goto failed ;
	if (0 != parentbucket->AddChildBucketFunction(_OutputFunction)) 
		goto failed ;
	_OutputFunction.SetBEBucket(parentbucket) ;

	// make sure this bucket knows its parent bucket
	_ParentBucket = parentbucket ;
	}

	// compute new fn size of computing this bucket; we assume child buckets are set up.
	newFNsize = 0 ;
	for (i = 0 ; i < nChildBucketFunctions() ; i++) {
		ARE::Function *f = ChildBucketFunction(i) ;
		if (f->TableSize() < 0) {
			f->ComputeTableSize() ;
			if (f->TableSize() < 0) 
				{ newFNsize = -1 ; break ; }
			}
		newFNsize += ChildBucketFunction(i)->TableSize() ;
		}
	if (newFNsize >= 0) 
		newFNsize += _OutputFunction.TableSize() ;
	if (newFNsize >= 0) 
		_ComputationNewFunctionSize = newFNsize ;

	return 0 ;
failed :
	if (_ParentBucket) {
		_ParentBucket->RemoveChildBucketFunction(_OutputFunction) ;
		_ParentBucket = NULL ;
		}
	_OutputFunction.Destroy() ;
	return 1 ;
}


int BucketElimination::Bucket::ComputeFirstVariableDistribution_1Block(ARE_Function_TableType *dist)
{
	int i, j, k, ret = 0 ;

	BEworkspace *bews = dynamic_cast<BEworkspace*>(_Workspace) ;
	if (NULL == bews) 
		return ERRORCODE_generic ;
	ARE::ARP *problem = bews->Problem() ;
	if (NULL == problem) 
		return ERRORCODE_generic ;
	ARE_Function_TableType nv = bews->FnCombinationNeutralValue() ;
	const int w = Width() ;
	const int *signature = Signature() ;
	const int nOF = nOriginalFunctions() ;
	const int nCF = nChildBucketFunctions() ;
	int nTotalFunctions = nOF + nCF ;
	if (w < 0 || nTotalFunctions < 1) {
		return 0 ;
		}

	if (w > MAX_NUM_VARIABLES_PER_BUCKET) 
		return ERRORCODE_too_many_variables ;
	if (nTotalFunctions > MAX_NUM_FUNCTIONS_PER_BUCKET) 
		return ERRORCODE_too_many_functions ;

	int values[MAX_NUM_VARIABLES_PER_BUCKET] ; // this is the current value combination of the arguments of this function (table block).
	ARE::Function *flist[MAX_NUM_FUNCTIONS_PER_BUCKET] ; // this is a list of input functions
	ARE::FunctionTableBlock *ftblist[MAX_NUM_FUNCTIONS_PER_BUCKET] ; // this is a list of input function table blocks, regardless of whether they are in-memory or external(disk)-memory.

	ARE_Function_TableType const_factor = bews->FnCombinationNeutralValue() ;
	int nFNs = 0 ;
	for (j = 0 ; j < nOF ; j++) {
		ARE::Function *f = OriginalFunction(j) ;
		if (NULL == f) return ERRORCODE_fn_ptr_is_NULL ;
		if (0 == f->N()) bews->ApplyFnCombinationOperator(const_factor, f->ConstValue()) ;
		else  { flist[nFNs] = f ; ftblist[nFNs++] = f->Table() ; f->ComputeArgumentsPermutationList(w, signature) ; }
		}
	for (; j < nTotalFunctions ; j++) {
		ARE::Function *f = ChildBucketFunction(j - nOF) ;
		if (NULL == f) return ERRORCODE_fn_ptr_is_NULL ;
		if (0 == f->N()) bews->ApplyFnCombinationOperator(const_factor, f->ConstValue()) ;
		else  { flist[nFNs] = f ; ftblist[nFNs++] = f->Table() ; f->ComputeArgumentsPermutationList(w, signature) ; }
		}

	INT64 ElimSize = 1 ;
	for (j = 1 ; j < w ; j++) 
		ElimSize *= problem->K(signature[j]) ;

	ARE::Function *MissingFunction = NULL ;
	INT64 MissingBlockIDX = -1 ;

	for (j = 0 ; j < w ; j++) 
		values[j] = 0 ;
	int v0 = _Vars[0] ;
	for (i = 0 ; i < problem->K(v0) ; i++) {
		ARE_Function_TableType V = bews->VarEliminationDefaultValue() ;
		for (INT64 ElimIDX = 0 ; ElimIDX < ElimSize ; ElimIDX++) {
			ARE_Function_TableType value = nv ;
			for (j = 0 ; j < nFNs ; j++) { // note : it should be that 0 != flist[j]->N(). note : it is assumed that flist[j] has 1 block. note : it is assumed that order of flist[j] arguments is the same as signature.
				INT64 adr = flist[j]->ComputeFnTableAdr_wrtLocalPermutation(w, values, problem->K()) ;
				bews->ApplyFnCombinationOperator(value, ftblist[j]->Entry(adr)) ;
				}
			bews->ApplyVarEliminationOperator(V, value) ;
			// go to next argument value combination
			ARE::EnumerateNextArgumentsValueCombination(w-1, signature+1, values+1, problem->K()) ;
			}
		bews->ApplyFnCombinationOperator(V, const_factor) ;
		dist[i] = V ;
		values[0]++ ;
		}
done :
	return ret ;
}


int BucketElimination::Bucket::ComputeOutputFunction_EliminateAllVars_1Block(void)
{
	int j, k, ret = 0 ;

	BEworkspace *bews = dynamic_cast<BEworkspace*>(_Workspace) ;
	if (NULL == bews) 
		return ERRORCODE_generic ;
	ARE::ARP *problem = bews->Problem() ;
	if (NULL == problem) 
		return ERRORCODE_generic ;
	ARE_Function_TableType nv = bews->FnCombinationNeutralValue() ;
	ARE::Function & f = OutputFunction() ;
	ARE_Function_TableType & V = f.ConstValue() ;
	const int w = Width() ;
	const int *signature = Signature() ;
	const int nOF = nOriginalFunctions() ;
	const int nCF = nChildBucketFunctions() ;
	int nTotalFunctions = nOF + nCF ;
	if (w < 0 || nTotalFunctions < 1) {
		V = nv ;
		return 0 ;
		}

	if (w > MAX_NUM_VARIABLES_PER_BUCKET) 
		return ERRORCODE_too_many_variables ;
	if (nTotalFunctions > MAX_NUM_FUNCTIONS_PER_BUCKET) 
		return ERRORCODE_too_many_functions ;

	int values[MAX_NUM_VARIABLES_PER_BUCKET] ; // this is the current value combination of the arguments of this function (table block).
	ARE::Function *flist[MAX_NUM_FUNCTIONS_PER_BUCKET] ; // this is a list of input functions
	ARE::FunctionTableBlock *ftblist[MAX_NUM_FUNCTIONS_PER_BUCKET] ; // this is a list of input function table blocks, regardless of whether they are in-memory or external(disk)-memory.

	ARE_Function_TableType const_factor = bews->FnCombinationNeutralValue() ;
	int nFNs = 0 ;
	for (j = 0 ; j < nOF ; j++) {
		ARE::Function *f = OriginalFunction(j) ;
		if (NULL == f) return ERRORCODE_fn_ptr_is_NULL ;
		if (0 == f->N()) bews->ApplyFnCombinationOperator(const_factor, f->ConstValue()) ;
		else  { flist[nFNs] = f ; ftblist[nFNs++] = f->Table() ; f->ComputeArgumentsPermutationList(w, signature) ; }
		}
	for (; j < nTotalFunctions ; j++) {
		ARE::Function *f = ChildBucketFunction(j - nOF) ;
		if (NULL == f) return ERRORCODE_fn_ptr_is_NULL ;
		if (0 == f->N()) bews->ApplyFnCombinationOperator(const_factor, f->ConstValue()) ;
		else  { flist[nFNs] = f ; ftblist[nFNs++] = f->Table() ; f->ComputeArgumentsPermutationList(w, signature) ; }
		}

	INT64 ElimSize = 1 ;
	for (j = 0 ; j < w ; j++) 
		ElimSize *= problem->K(signature[j]) ;

	ARE::Function *MissingFunction = NULL ;
	INT64 MissingBlockIDX = -1 ;

	V = bews->VarEliminationDefaultValue() ;
	for (j = 0 ; j < w ; j++) 
		values[j] = 0 ;
	for (INT64 ElimIDX = 0 ; ElimIDX < ElimSize ; ElimIDX++) {
		ARE_Function_TableType value = nv ;
		for (j = 0 ; j < nFNs ; j++) { // note : it should be that 0 != flist[j]->N(). note : it is assumed that flist[j] has 1 block. note : it is assumed that order of flist[j] arguments is the same as signature.
			INT64 adr = flist[j]->ComputeFnTableAdr_wrtLocalPermutation(w, values, problem->K()) ;
			bews->ApplyFnCombinationOperator(value, ftblist[j]->Entry(adr)) ;
			}
		bews->ApplyVarEliminationOperator(V, value) ;
		// go to next argument value combination
		ARE::EnumerateNextArgumentsValueCombination(w, signature, values, problem->K()) ;
		}
	bews->ApplyFnCombinationOperator(V, const_factor) ;
done :
	return ret ;
}


int BucketElimination::Bucket::ComputeOutputFunction_1Block(void)
{
	int i, j, k ;

	BEworkspace *bews = dynamic_cast<BEworkspace*>(_Workspace) ;
	if (NULL == bews) 
		return ERRORCODE_generic ;
	ARE::ARP *problem = bews->Problem() ;
	if (NULL == problem) 
		return ERRORCODE_generic ;
	if (nVars() < 1) 
		return 1 ;

	ARE::Function & f = OutputFunction() ;
	AllocateOutputFunctionBlockComputationResult(1000000, 0) ;
	f.ComputeTableSize() ;
	f.SetnTableBlocks(0) ;
	f.AllocateInMemoryAsSingleTableBlock() ;
	int nA = f.N() ;
	if (0 == nA) 
		return ComputeOutputFunction_EliminateAllVars_1Block() ;
	ARE::FunctionTableBlock *ftb = f.Table() ;
	if (NULL == ftb) 
		{ return 1 ; }

	const int *refFNarguments = f.Arguments() ;
	const int w = Width() ;
	const int *signature = Signature() ;
	const int nOF = nOriginalFunctions() ;
	const int nCF = nChildBucketFunctions() ;
	int nTotalFunctions = nOF + nCF ;
	if (w < 0 || nTotalFunctions < 1) {
		return 1 ;
		}

	if (0 != ftb->AllocateData())
		return ERRORCODE_memory_allocation_failure ;

	int values[MAX_NUM_VARIABLES_PER_BUCKET] ; // this is the current value combination of the arguments of this function (table block).
	ARE::Function *flist[MAX_NUM_FUNCTIONS_PER_BUCKET] ; // this is a list of input functions
	ARE::FunctionTableBlock *ftblist[MAX_NUM_FUNCTIONS_PER_BUCKET] ; // this is a list of input function table blocks, regardless of whether they are in-memory or external(disk)-memory.

	int vars[MAX_NUM_VARIABLES_PER_BUCKET] ; // this is the list of variables : vars2Keep + vars2Eliminate
	int nAtotal = nA + nVars() ;
	for (j = 0 ; j < nA ; j++) {
		vars[j] = refFNarguments[j] ;
		values[j] = 0 ;
		}
	for (; j < nAtotal ; j++) 
		vars[j] = Var(j - nA) ;

	ARE_Function_TableType const_factor = bews->FnCombinationNeutralValue() ;
	int nFNs = 0 ;
	for (j = 0 ; j < nOF ; j++) {
		ARE::Function *f = OriginalFunction(j) ;
		if (NULL == f) return ERRORCODE_fn_ptr_is_NULL ;
		if (0 == f->N()) bews->ApplyFnCombinationOperator(const_factor, f->ConstValue()) ;
		else  { flist[nFNs] = f ; ftblist[nFNs++] = f->Table() ; f->ComputeArgumentsPermutationList(nAtotal, vars) ; }
		}
	for (; j < nTotalFunctions ; j++) {
		ARE::Function *f = ChildBucketFunction(j - nOF) ;
		if (NULL == f) return ERRORCODE_fn_ptr_is_NULL ;
		if (0 == f->N()) bews->ApplyFnCombinationOperator(const_factor, f->ConstValue()) ;
		else  { flist[nFNs] = f ; ftblist[nFNs++] = f->Table() ; f->ComputeArgumentsPermutationList(nAtotal, vars) ; }
		}

	INT64 ElimIDX, ElimSize = 1 ;
	for (j = 0 ; j < nVars() ; j++) 
		ElimSize *= problem->K(Var(j)) ;

	__int64 KeepIDX ;
	ARE_Function_TableType *data = ftb->Data() ;
	for (KeepIDX = 0 ; KeepIDX < ftb->Size() ; KeepIDX++) {
		for (j = nA ; j < nAtotal ; j++) 
			values[j] = 0 ;
		data[KeepIDX] = bews->VarEliminationDefaultValue() ;
		for (ElimIDX = 0 ; ElimIDX < ElimSize ; ElimIDX++) {
			ARE_Function_TableType value = bews->FnCombinationNeutralValue() ;
			for (j = 0 ; j < nFNs ; j++) {
				__int64 adr = flist[j]->ComputeFnTableAdr_wrtLocalPermutation(nAtotal, values, problem->K()) ;
				bews->ApplyFnCombinationOperator(value, ftblist[j]->Entry(adr)) ;
				}
			bews->ApplyVarEliminationOperator(data[KeepIDX], value) ;
			// go to next argument value combination
			ARE::EnumerateNextArgumentsValueCombination(nAtotal, vars, values, problem->K()) ;
			}
		bews->ApplyFnCombinationOperator(data[KeepIDX], const_factor) ;
		}

	return 0 ;
}


int BucketElimination::Bucket::NoteOutputFunctionComputationCompletion(void)
{
	int i ;

	// child bucket output functions are no longer needed; release them.
	for (i = 0 ; i < _nChildBucketFunctions ; i++) {
		ARE::Function *f = _ChildBucketFunctions[i] ;
		if (NULL == f) continue ;
		// if needed, delete input function (tables)
		f->ReleaseAllFTBs(_Workspace->DeleteUsedTables()) ;
		}

	return 0 ;
}


int BucketElimination::Bucket::ReorderFunctionScopesForExternalMemory(bool IncludeOriginalFunctions, bool IncludeNewFunctions)
{
	// destroy table of all childbucketfunctions of this bucket;
	// reordering of scopes should be done in the beginning when these tables don't exist yet, 
	// so this should not be a problem.
	int i ;
	if (IncludeNewFunctions) {
		for (i = 0 ; i < _nChildBucketFunctions ; i++) {
			ARE::Function *f = _ChildBucketFunctions[i] ;
			f->DestroyTable() ;
			f->DestroyFTBList() ;
			}
		}

	// reorder scopes
	if (_OutputFunction.N() < 1) {
		if (IncludeNewFunctions) {
			for (i = 0 ; i < _nChildBucketFunctions ; i++) {
				if (0 != _ChildBucketFunctions[i]->ReOrderArguments(_nVars, _Vars, 0, NULL)) 
					return 1 ;
				}
			}
		if (IncludeOriginalFunctions) {
			for (i = 0 ; i < _nOriginalFunctions ; i++) {
				if (0 != _OriginalFunctions[i]->ReOrderArguments(_nVars, _Vars, 0, NULL)) 
					return 1 ;
				}
			}
		}
	else {
		if (IncludeNewFunctions) {
			for (i = 0 ; i < _nChildBucketFunctions ; i++) {
				if (0 != _ChildBucketFunctions[i]->ReOrderArguments(_OutputFunction.N(), _OutputFunction.Arguments(), _nVars, _Vars)) 
					return 1 ;
				}
			}
		if (IncludeOriginalFunctions) {
			for (i = 0 ; i < _nOriginalFunctions ; i++) {
				if (0 != _OriginalFunctions[i]->ReOrderArguments(_OutputFunction.N(), _OutputFunction.Arguments(), _nVars, _Vars)) 
					return 1 ;
				}
			}
		}

	return 0 ;
}


int BucketElimination::Bucket::AllocateOutputFunctionBlockComputationResult(__int64 MaxBlockSize, int nComputingThreads)
{
	ARE::utils::AutoLock lock(_Workspace->FTBMutex()) ;
	if (_OutputFunction.nTableBlocks() < 0) {
		if (0 != _OutputFunction.ComputeTableBlockSize(MaxBlockSize, nComputingThreads)) 
			return 1 ;
		}
	_nOutputFunctionBlocks = _OutputFunction.nTableBlocks() ;
	if (_nOutputFunctionBlocks < 0) 
		// this means fn has no table 
		return 0 ;
	int size = _nOutputFunctionBlocks > 0 ? (7 + _nOutputFunctionBlocks) >> 3 : 1 ;
	if (_OutputFunctionBlockComputationResultSize == size) {
		memset(_OutputFunctionBlockComputationResult, 0, _OutputFunctionBlockComputationResultSize) ;
		_nOutputFunctionBlocksComputed = 0 ;
		return 0 ;
		}
	if (NULL != _OutputFunctionBlockComputationResult) 
		delete [] _OutputFunctionBlockComputationResult ;
	_OutputFunctionBlockComputationResult = new unsigned char[size] ;
	if (NULL == _OutputFunctionBlockComputationResult) 
		return 0 ;
	_OutputFunctionBlockComputationResultSize = size ;
	memset(_OutputFunctionBlockComputationResult, 0, _OutputFunctionBlockComputationResultSize) ;
	_nOutputFunctionBlocksComputed = 0 ;
	return 0 ;
}


bool BucketElimination::Bucket::IsOutputFunctionBlockComputed(__int64 IDX)
{
	ARE::utils::AutoLock lock(_Workspace->FTBMutex()) ;
	if (0 == IDX) {
		// this means function table has 1 in-memory block
		return _nOutputFunctionBlocksComputed > 0 ;
		}
	--IDX ; // disk memeory block indeces run [1,nBlocks], but we need [0,nBlocks)
	return 0 != (_OutputFunctionBlockComputationResult[IDX >> 3] & (1 << (IDX & 7))) ;
}


__int64 BucketElimination::Bucket::nOutputFunctionBlocksComputed(void) const
{
	ARE::utils::AutoLock lock(_Workspace->FTBMutex()) ;
	return _nOutputFunctionBlocksComputed ;
}


void BucketElimination::Bucket::MarkOutputFunctionBlockComputed(__int64 IDX)
{
	ARE::utils::AutoLock lock(_Workspace->FTBMutex()) ;
	if (0 == IDX) {
		// this means function table has 1 in-memory block
		if (0 == _nOutputFunctionBlocksComputed) 
			_nOutputFunctionBlocksComputed = 1 ;
		NoteOutputFunctionComputationCompletion() ;
		return ;
		}
	if (! IsOutputFunctionBlockComputed(IDX)) {
		--IDX ; // disk memeory block indeces run [1,nBlocks], but we need [0,nBlocks)
		_OutputFunctionBlockComputationResult[IDX >> 3] |= (1 << (IDX & 7)) ;
		++_nOutputFunctionBlocksComputed ;
		}
	if (_nOutputFunctionBlocksComputed >= _nOutputFunctionBlocks) 
		NoteOutputFunctionComputationCompletion() ;
}
