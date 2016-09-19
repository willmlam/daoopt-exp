#include <stdlib.h>
#include <memory.h>

#include <Function.hxx>
#include <Bucket.hxx>
#include <MBEworkspace.hxx>
#include <Bucket.hxx>
#include <MiniBucket.hxx>
#include <Sort.hxx>


BucketElimination::MiniBucket::MiniBucket(void)
	:
	_Workspace(NULL),
	_IDX(-1), 
	_Width(-1), 
	_Signature(NULL), 
	_SortedSignature(NULL), 
	_SignatureArraySize(0), 
	_nVars(0), 
	_Vars(NULL), 
	_VarsSpace(0), 
	_nFunctions(0), 
	_Functions(NULL), 
	_FunctionsArraySize(0), 
	_ComputationNewFunctionSize(-1) 
{
}


BucketElimination::MiniBucket::~MiniBucket(void)
{
	Destroy() ;
}


BucketElimination::MiniBucket::MiniBucket(MBEworkspace & WS, int32_t IDX, int32_t V)
	:
	_Workspace(&WS),
	_IDX(IDX), 
	_V(V), 
	_Width(-1), 
	_Signature(NULL), 
	_SortedSignature(NULL), 
	_SignatureArraySize(0), 
	_nVars(0), 
	_Vars(NULL), 
	_VarsSpace(0), 
	_nFunctions(0), 
	_Functions(NULL), 
	_FunctionsArraySize(0), 
	_ComputationNewFunctionSize(-1) 
{
	BucketElimination::Bucket *b = WS.getBucket(V) ;
	// if Var is given, add it
	if (V >= 0) 
		AddVar(V) ;
	// fix up output function of the bucket
	_OutputFunction.Initialize(_Workspace, _Workspace->Problem(), -1) ;
	_OutputFunction.SetOriginatingBucket(b) ;
	_OutputFunction.SetOriginatingMiniBucket(this) ;
}


void BucketElimination::MiniBucket::Destroy(void)
{
	if (NULL != _Functions) {
		delete [] _Functions ;
		_Functions = NULL ;
		_FunctionsArraySize = 0 ;
		}
	if (NULL != _Signature) {
		delete [] _Signature ;
		_Signature = NULL ;
		}
	_SortedSignature = NULL ; // don't delete this; it is allocated (and deleted) as second half of _Signature.
	_SignatureArraySize = 0 ;
	_OutputFunction.Destroy() ;
	_Width = -1 ;
	_nFunctions = 0 ;
	_nVars = 0 ;
	if (NULL != _Vars) {
		delete [] _Vars ;
		_Vars = NULL ;
		}
	_VarsSpace = 0 ;
	_ComputationNewFunctionSize = -1 ;
}


void BucketElimination::MiniBucket::Initalize(BucketElimination::Bucket & B, int32_t IDX)
{
	_Bucket = &B ;
	_Workspace = B.Workspace() ;
	_IDX = IDX ;
	_V = B.Var(0) ;
	_OutputFunction.SetProblem(_Workspace->Problem()) ;
	_OutputFunction.SetWS(B.Workspace()) ;
	_OutputFunction.SetIDX(-(_V+1)) ;
}


int32_t BucketElimination::MiniBucket::ComputeMissingVariables(ARE::Function & F, std::vector<int32_t> & MissingVars)
{
	if (F.N() <= 0) 
		return 0 ;
	int32_t i, j ;
	int32_t *fSortedArgs = F.SortedArgumentsList(true) ;
	if (NULL == fSortedArgs) {
		for (i = F.N()-1 ; i >= 0 ; i--) {
			int32_t v = F.Argument(i) ;
			for (j = _Width-1 ; j >= 0 ; j--) {
				if (_Signature[j] == v) 
					break ;
				}
			if (j < 0) 
				MissingVars.push_back(v) ;
			}
		}
	else {
		i = j = 0 ;
		while (i < _Width && j < F.N()) {
			if (_SortedSignature[i] == fSortedArgs[j]) 
				{ i++ ; j++ ; }
			else if (_SortedSignature[i] < fSortedArgs[j]) 
				{ i++ ; }
			else 
				{ MissingVars.push_back(fSortedArgs[j++]) ; }
			}
		for (; j < F.N() ; j++) 
			MissingVars.push_back(fSortedArgs[j]) ;
		}

	return 0 ;
}


int32_t BucketElimination::MiniBucket::AllowsFunction(ARE::Function & F, int32_t MBmaxsize, std::vector<int32_t> & HelperArray)
{
	if (0 == F.N()) 
		return 1 ; // const fn can placed in any minibucket

	// iBound is max number of variables in the minibucket.

	// must have current width of the minibucket.
	if (_Width < 0) {
		if (0 != ComputeSignature()) 
			return -1 ;
		if (_Width < 0) 
			return -1 ;
		}
	if (_Width <= 0) 
		return 1 ; // empty minibucket always allows a fn

	if (F.N() > MBmaxsize) 
		return 0 ; // this minibucket is not empty (_Width>0), and since this F signature is larger than MBmaxsize, so we know we will go over limit.
	if ((_Width + F.N()) <= MBmaxsize) 
		return 1 ; // _Width + F.N() is upper bound on the new width if we added this fn.

	if (HelperArray.capacity() < F.N()) {
		HelperArray.reserve(F.N()) ;
		if (HelperArray.capacity() < F.N()) 
			return -1 ;
		}
	HelperArray.clear() ;
	ComputeMissingVariables(F, HelperArray) ;
	int32_t nMissing = HelperArray.size() ;
	if (0 == nMissing) 
		return 1 ; // all arguments of the function are already in the mini-bucket
	return (_Width + nMissing) <= MBmaxsize ? 1 : 0 ;
}


int32_t BucketElimination::MiniBucket::AddFunction(ARE::Function & F, std::vector<int32_t> & HelperArray)
{
	if (_Width < 0) {
		ComputeSignature() ;
		if (_Width < 0) 
			return 1 ;
		}

	// check if we have enough space
	if (_nFunctions+1 > _FunctionsArraySize) {
		int32_t newsize = _FunctionsArraySize + 8 ;
		ARE::Function **newspace = new ARE::Function*[newsize] ;
		if (NULL == newspace) 
			return 1 ;
		if (_nFunctions > 0) 
			memcpy(newspace, _Functions, sizeof(ARE::Function *)*_nFunctions) ;
		if (NULL == _Functions) 
			delete [] _Functions ;
		_Functions = newspace ;
		_FunctionsArraySize = newsize ;
		}

	_Functions[_nFunctions++] = &F ;

	// fix up Signature/SortedSignature
	if (HelperArray.capacity() < F.N()) {
		HelperArray.reserve(F.N()) ;
		if (HelperArray.capacity() < F.N()) 
			return -1 ;
		}
	HelperArray.clear() ;
	ComputeMissingVariables(F, HelperArray) ;
	int32_t nMissing = HelperArray.size() ;
	if (nMissing > 0) {
		if (_Width + nMissing > _SignatureArraySize) {
			int32_t newsize = _Width + nMissing + 8 ;
			int32_t *newspace = new int32_t[2*newsize] ;
			if (NULL == newspace) 
				return 1 ;
			if (_Width > 0) {
				memcpy(newspace, _Signature, sizeof(int32_t)*_Width) ;
				memcpy(newspace+newsize, _SortedSignature, sizeof(int32_t)*_Width) ;
				}
			if (NULL != _Signature) 
				delete [] _Signature ;
			_Signature = newspace ;
			_SortedSignature = newspace+newsize ;
			_SignatureArraySize = newsize ;
			}
		for (int32_t i = 0 ; i < nMissing ; i++) {
			_Signature[_Width] = HelperArray[i] ;
			_SortedSignature[_Width] = HelperArray[i] ;
			_Width++ ;
			}
		int32_t left[32], right[32] ;
		QuickSortLong2((int32_t*) _SortedSignature, _Width, left, right) ;
		}

	return 0 ;
}


int32_t BucketElimination::MiniBucket::RemoveFunction(ARE::Function & F)
{
	int32_t i, n = 0 ;
	for (i = _nFunctions - 1 ; i >= 0 ; i--) {
		if (&F == _Functions[i]) {
			_Functions[i] = _Functions[--_nFunctions] ;
			++n ;
			}
		}
	// if any removals, update signature
	if (n > 0) 
		ComputeSignature() ;
	return 0 ;
}


int32_t BucketElimination::MiniBucket::ComputeSignature(void)
{
// DEBUGGGGG
/*if (NULL != ARE::fpLOG) {
	fprintf(ARE::fpLOG, "\n         ComputeSignature var=%d", (int32_t) V()) ;
	fflush(ARE::fpLOG) ;
	}*/
	if (NULL != _Signature) {
		delete [] _Signature ;
		_Signature = NULL ;
		}
	_SignatureArraySize = 0 ;
	_SortedSignature = NULL ;
	_Width = 0 ;

	// compute approx width
	if (0 == _nFunctions) 
		return 0 ;
	int32_t i, n = 0 ;
	for (i = 0 ; i < _nFunctions ; i++) 
		n += _Functions[i]->N() ;
	// n is an upper bound on the width
	if (n <= 0) 
		return 0 ;

	// OriginalSignature is part of Signature
	_Signature = new int32_t[2*n] ;
	if (NULL == _Signature) 
		goto failed ;
	_SignatureArraySize = n ;
	_SortedSignature = _Signature + n ;

	// add scopes of bucketfunctions to the signature
	int32_t j, k ;
	for (i = 0 ; i < _nFunctions ; i++) {
		ARE::Function *f = _Functions[i] ;
		for (j = 0 ; j < f->N() ; j++) {
			int32_t v = f->Argument(j) ;
			for (k = 0 ; k < _Width ; k++) 
				{ if (_Signature[k] == v) break ; }
			if (k < _Width) 
				continue ;
			_Signature[_Width] = v ;
			_SortedSignature[_Width++] = v ;
			}
		}

	if (_Width > 1) {
		int32_t left[32], right[32] ;
		QuickSortLong2((int32_t*) _SortedSignature, _Width, left, right) ;
		}

	return 0 ;
failed :
	return 1 ;
}


int64_t BucketElimination::MiniBucket::ComputeProcessingComplexity(void)
{
	int64_t n = 1 ;
	for (int32_t i = 0 ; i < _Width ; i++) {
		n *= _Workspace->Problem()->K(_Signature[i]) ;
		}
	return n ;
}


int32_t BucketElimination::MiniBucket::ComputeOutputFunctionWithScopeWithoutTable(void)
{
	if (NULL == _Bucket) 
		return 1 ;

// DEBUGGGGG
/*if (NULL != ARE::fpLOG) {
	fprintf(ARE::fpLOG, "\n         ComputeOutputFunctionWithScopeWithoutTable var=%d", (int32_t) V()) ;
	fflush(ARE::fpLOG) ;
	}*/

	// dump current bucket fn; cleanup.
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

#ifdef _DEBUG
	for (int32_t jjj = 0 ; jjj < _Width ; jjj++) {
		if (_Signature[jjj] < 0 || _Signature[jjj] >= _Bucket->Workspace()->Problem()->N()) {
			int32_t error = 1 ;
			}
		}
#endif // _DEBUG

	// make a local copy, while ignoring bucket variables
	std::vector<int32_t> vlist ;
	vlist.reserve(_Width) ;
	int32_t i, j ;
	for (i = 0 ; i < _Width ; i++) {
		int32_t v = _Signature[i] ;
if (v < 0 || v >= _Bucket->Workspace()->Problem()->N()) {
	int32_t error = 1 ;
	}
		for (j = 0 ; j < _nVars ; j++) {
			if (_Vars[j] == v) 
				break ;
			}
		if (j < _nVars) 
			continue ;
		vlist.push_back(v) ;
		}
	int32_t n = vlist.size() ;

#ifdef _DEBUG
	for (j = 0 ; j < n ; j++) {
		if (vlist[j] < 0 || vlist[j] >= _Bucket->Workspace()->Problem()->N()) {
			int32_t error = 1 ;
			}
		}
#endif // _DEBUG
	double newFNsize;
  double temp_sum;
  double of_log;
  double log_max;

	// create and initialize bucket function
	if (0 != _OutputFunction.SetArguments(n, vlist.data())) 
		goto failed ;
	of_log = _OutputFunction.GetTableSize_Log10() ;
	log_max = log10((double) _I64_MAX) ;
	if (of_log < log_max) 
		_OutputFunction.ComputeTableSize() ;
	else {
		// output fn size too large; continue since we can still compute the output fn scope, which is what this fn is supposed to do.
		}
	if (_nFunctions > 0) 
		_OutputFunction.SetType(_Functions[0]->Type()) ;

	// find appropriate parent bucket and assign the bucket function to it
	{
	const int32_t *varpos = _Workspace->VarPos() ;
	int32_t v = _OutputFunction.GetHighestOrderedVariable(varpos) ;
	BucketElimination::Bucket *parentbucket = _Workspace->MapVar2Bucket(v) ;
	if (NULL == parentbucket) 
		// this is not supposed to happen
		goto failed ;
	if (parentbucket->IDX() >= _Bucket->IDX()) 
		// this is not supposed to happen
		goto failed ;
	_OutputFunction.SetBucket(parentbucket) ;
	_OutputFunction.SetOriginatingBucket(_Bucket) ;
	_OutputFunction.SetOriginatingMiniBucket(this) ;
	}

	// compute new fn size of computing this bucket; we assume child buckets are set up.
	// compute in log space since individual fn sizes may be huge.
	/*
		from http://en.wikipedia.org/wiki/List_of_logarithmic_identities
		log10(SUM a_i) = log10(a_0) + log10(1 + SUM (a_i/a_0)) = log10(a_0) + log10(1 + SUM 10^(log10(a_i) - log10(a_0)))
	*/
	newFNsize = _OutputFunction.GetTableSize_Log10();
  temp_sum = 1.0 ;
	if (newFNsize >= 0.0) {
		for (i = 0 ; i < _nFunctions ; i++) {
			ARE::Function *f = _Functions[i] ;
			if (NULL == f->OriginatingMiniBucket()) 
				continue ; // this is original fn, not MBE generated function.
			double fnsize = f->GetTableSize_Log10() ;
			if (fnsize < 0.0) 
				{ newFNsize = -1.0 ; break ; }
			temp_sum += pow(10.0, fnsize - newFNsize) ;
			}
		}
	_ComputationNewFunctionSize = 0 ;
	if (newFNsize >= 0.0) {
		newFNsize += log10(temp_sum) ;
		// if size overflows int64_t, set it to max.
		double log_max = log10((double) _I64_MAX) ;
		if (newFNsize < log_max) 
			_ComputationNewFunctionSize = pow(10.0, newFNsize) ;
		else 
			_ComputationNewFunctionSize = _I64_MAX ;
		}
// DEBUGGGGG
/*if (NULL != ARE::fpLOG) {
	int64_t tNOW = ARE::GetTimeInMilliseconds() ;
	fprintf(ARE::fpLOG, "\n         ComputeOutputFunctionWithScopeWithoutTable var=%d newfncompsize=%lld", (int32_t) V(), _ComputationNewFunctionSize) ;
	fflush(ARE::fpLOG) ;
	}*/

	return 0 ;
failed :
	_OutputFunction.Destroy() ;
	return 1 ;
}


int32_t BucketElimination::MiniBucket::ComputeOutputFunction_EliminateAllVars(void)
{
	int32_t j, k, ret = 0 ;

	MBEworkspace *bews = dynamic_cast<MBEworkspace*>(_Workspace) ;
	if (NULL == bews) 
		return ERRORCODE_generic ;
	ARE::ARP *problem = bews->Problem() ;
	if (NULL == problem) 
		return ERRORCODE_generic ;
	ARE_Function_TableType nv = bews->FnCombinationNeutralValue() ;
	ARE::Function & f = OutputFunction() ;
	ARE_Function_TableType & V = f.ConstValue() ;
	const int32_t w = Width() ;
	if (w < 0 || _nFunctions < 1) {
		V = nv ;
		return 0 ;
		}
	const int32_t *signature = Signature() ;

	if (w > MAX_NUM_VARIABLES_PER_BUCKET) 
		return ERRORCODE_too_many_variables ;
	if (_nFunctions > MAX_NUM_FUNCTIONS_PER_BUCKET) 
		return ERRORCODE_too_many_functions ;

	int32_t values[MAX_NUM_VARIABLES_PER_BUCKET] ; // this is the current value combination of the arguments of this function (table block).
	ARE::Function *flist[MAX_NUM_FUNCTIONS_PER_BUCKET] ; // this is a list of input functions

	ARE_Function_TableType const_factor = bews->FnCombinationNeutralValue() ;
	int32_t nFNs = 0 ;
	for (j = 0 ; j < _nFunctions ; j++) {
		ARE::Function *f = _Functions[j] ;
		if (NULL == f) return ERRORCODE_fn_ptr_is_NULL ;
		if (0 == f->N()) bews->ApplyFnCombinationOperator(const_factor, f->ConstValue()) ;
		else { flist[nFNs++] = f ; f->ComputeArgumentsPermutationList(w, signature) ; }
		}

	int64_t ElimSize = 1 ;
	for (j = 0 ; j < w ; j++) 
		ElimSize *= problem->K(signature[j]) ;

	ARE::Function *MissingFunction = NULL ;
	int64_t MissingBlockIDX = -1 ;

	V = bews->VarEliminationDefaultValue() ;
	for (j = 0 ; j < w ; j++) 
		values[j] = 0 ;
	for (int64_t ElimIDX = 0 ; ElimIDX < ElimSize ; ElimIDX++) {
		ARE_Function_TableType value = nv ;
		for (j = 0 ; j < nFNs ; j++) { // note : it should be that 0 != flist[j]->N(). note : it is assumed that flist[j] has 1 block. note : it is assumed that order of flist[j] arguments is the same as signature.
			int64_t adr = flist[j]->ComputeFnTableAdr_wrtLocalPermutation(w, values, problem->K()) ;
			bews->ApplyFnCombinationOperator(value, flist[j]->TableEntry(adr)) ;
			}
		bews->ApplyVarEliminationOperator(V, value) ;
		// go to next argument value combination
		ARE::EnumerateNextArgumentsValueCombination(w, signature, values, problem->K()) ;
		}
	bews->ApplyFnCombinationOperator(V, const_factor) ;
done :
	return ret ;
}


int32_t BucketElimination::MiniBucket::ComputeOutputFunction(ARE::Function *fMaxMarginal, ARE::Function *fAvgMaxMarginal)
{
	int32_t i, j, k ;

	MBEworkspace *bews = dynamic_cast<MBEworkspace*>(_Workspace) ;
	if (NULL == bews) 
		return ERRORCODE_generic ;
	ARE::ARP *problem = bews->Problem() ;
	if (NULL == problem) 
		return ERRORCODE_generic ;
	if (_Width < 0) {
		if (0 != ComputeSignature()) 
			return ERRORCODE_generic ;
		if (_Width < 0) 
			return ERRORCODE_generic ;
		}
	if (_nFunctions < 1) {
		return 0 ;
		}
	if (nVars() < 1) 
		// nothing to do; should not happen; however, allow this to pass as ok; calling fn should be able to handle this.
		return 0 ;

	ARE::Function & f = OutputFunction() ;
	f.ComputeTableSize() ;
	int32_t nA = f.N() ;
	if (0 == nA) 
		return ComputeOutputFunction_EliminateAllVars() ;

	if (0 != f.AllocateTableData())
		return ERRORCODE_memory_allocation_failure ;

	int32_t values[MAX_NUM_VARIABLES_PER_BUCKET] ; // this is the current value combination of the arguments of this function (table block).
	int32_t vars[MAX_NUM_VARIABLES_PER_BUCKET] ; // this is the list of variables : vars2Keep + vars2Eliminate
	const int32_t *refFNarguments = f.Arguments() ;
	for (j = 0 ; j < nA ; j++) {
		vars[j] = refFNarguments[j] ;
		values[j] = 0 ;
		}
	for (; j < _Width ; j++) 
		vars[j] = Var(j - nA) ;

	ARE_Function_TableType const_factor = bews->FnCombinationNeutralValue() ;
	int32_t nFNs = 0 ;
	ARE::Function *flist[MAX_NUM_FUNCTIONS_PER_BUCKET] ; // this is a list of input functions
	for (j = 0 ; j < _nFunctions ; j++) {
		ARE::Function *f = _Functions[j] ;
		if (NULL == f) return ERRORCODE_fn_ptr_is_NULL ;
		if (0 == f->N()) bews->ApplyFnCombinationOperator(const_factor, f->ConstValue()) ;
		else  { flist[nFNs++] = f ; f->ComputeArgumentsPermutationList(_Width, vars) ; }
		}
	if (NULL != fMaxMarginal && NULL != fAvgMaxMarginal) {
		fMaxMarginal->ComputeArgumentsPermutationList(_Width, vars) ;
		fAvgMaxMarginal->ComputeArgumentsPermutationList(_Width, vars) ;
		}

	int64_t ElimIDX, ElimSize = 1 ;
	for (j = 0 ; j < nVars() ; j++) 
		ElimSize *= problem->K(Var(j)) ;

	ARE_Function_TableType *data = f.TableData() ;
	int64_t tablesize = f.TableSize(), adr ;
	for (int64_t KeepIDX = 0 ; KeepIDX < tablesize ; KeepIDX++) {
		for (j = nA ; j < _Width ; j++) 
			values[j] = 0 ;
		data[KeepIDX] = bews->VarEliminationDefaultValue() ;
		for (ElimIDX = 0 ; ElimIDX < ElimSize ; ElimIDX++) {
			ARE_Function_TableType value = bews->FnCombinationNeutralValue() ;
			for (j = 0 ; j < nFNs ; j++) {
				adr = flist[j]->ComputeFnTableAdr_wrtLocalPermutation(_Width, values, problem->K()) ;
				bews->ApplyFnCombinationOperator(value, flist[j]->TableEntry(adr)) ;
				}
			if (NULL != fMaxMarginal && NULL != fAvgMaxMarginal) {
				adr = fAvgMaxMarginal->ComputeFnTableAdr_wrtLocalPermutation(_Width, values, problem->K()) ;
				bews->ApplyFnCombinationOperator(value, fAvgMaxMarginal->TableEntry(adr)) ;
				adr = fMaxMarginal->ComputeFnTableAdr_wrtLocalPermutation(_Width, values, problem->K()) ;
				bews->ApplyFnDivisionOperator(value, fMaxMarginal->TableEntry(adr)) ;
				}
			bews->ApplyVarEliminationOperator(data[KeepIDX], value) ;
			// go to next argument value combination
			ARE::EnumerateNextArgumentsValueCombination(_Width, vars, values, problem->K()) ;
			}
		bews->ApplyFnCombinationOperator(data[KeepIDX], const_factor) ;
		}

	return 0 ;
}


int32_t BucketElimination::MiniBucket::ComputeOutputFunction(ARE::Function & f, const int32_t *ElimVars, int32_t nElimVars, int32_t *TempSpaceForVars)
{
	int32_t i, j, k ;

	MBEworkspace *bews = dynamic_cast<MBEworkspace*>(_Workspace) ;
	if (NULL == bews) 
		return ERRORCODE_generic ;
	ARE::ARP *problem = bews->Problem() ;
	if (NULL == problem) 
		return ERRORCODE_generic ;
	if (nElimVars < 1) 
		// nothing to do; should not happen; however, allow this to pass as ok; calling fn should be able to handle this.
		return 0 ;
	if (_Width < 0) {
		if (0 != ComputeSignature()) 
			return 1 ;
		if (_Width < 0) 
			return 1 ;
		}
	if (_nFunctions < 1) 
		return 0 ;

	// compute output function signature
	if (nElimVars < _Width) {
		const int32_t *sorted_scope = SortedSignature() ;
		for (i = 0 ; i < _Width ; i++) TempSpaceForVars[i] = sorted_scope[i] ; j = _Width ;
		ARE::SetMinus(TempSpaceForVars, j, ElimVars, nElimVars) ;
		if (0 != f.SetArguments(j, TempSpaceForVars)) 
			return 1 ;
		}

	// set upp output fn table
	f.ComputeTableSize() ;
	int32_t nA = f.N() ;
	if (nA > 0) {
		if (0 != f.AllocateTableData()) 
			return 1 ;
		}
	const int32_t *refFNarguments = f.Arguments() ;
	const int32_t *signature = Signature() ;

	int32_t values[MAX_NUM_VARIABLES_PER_BUCKET] ; // this is the current value combination of the arguments of this function (table block).
	int32_t vars[MAX_NUM_VARIABLES_PER_BUCKET] ; // this is the list of variables : vars2Keep + vars2Eliminate
	for (j = 0 ; j < nA ; j++) {
		vars[j] = refFNarguments[j] ;
		values[j] = 0 ;
		}
	for (; j < _Width ; j++) 
		vars[j] = ElimVars[j - nA] ;

	ARE_Function_TableType const_factor = bews->FnCombinationNeutralValue() ;
	int32_t nFNs = 0 ;
	ARE::Function *flist[MAX_NUM_FUNCTIONS_PER_BUCKET] ; // this is a list of input functions
	for (j = 0 ; j < _nFunctions ; j++) {
		ARE::Function *f = _Functions[j] ;
		if (NULL == f) return ERRORCODE_fn_ptr_is_NULL ;
		if (0 == f->N()) bews->ApplyFnCombinationOperator(const_factor, f->ConstValue()) ;
		else  { flist[nFNs++] = f ; f->ComputeArgumentsPermutationList(_Width, vars) ; }
		}

	int64_t ElimIDX, ElimSize = 1 ;
	for (j = 0 ; j < nElimVars ; j++) 
		ElimSize *= problem->K(ElimVars[j]) ;

	ARE_Function_TableType *data = f.TableData() ;
	int64_t tablesize = f.TableSize() ;
	for (int64_t KeepIDX = 0 ; KeepIDX < tablesize ; KeepIDX++) {
		for (j = nA ; j < _Width ; j++) 
			values[j] = 0 ;
		data[KeepIDX] = bews->VarEliminationDefaultValue() ;
		for (ElimIDX = 0 ; ElimIDX < ElimSize ; ElimIDX++) {
			ARE_Function_TableType value = bews->FnCombinationNeutralValue() ;
			for (j = 0 ; j < nFNs ; j++) {
				int64_t adr = flist[j]->ComputeFnTableAdr_wrtLocalPermutation(_Width, values, problem->K()) ;
				bews->ApplyFnCombinationOperator(value, flist[j]->TableEntry(adr)) ;
				}
			bews->ApplyVarEliminationOperator(data[KeepIDX], value) ;
			// go to next argument value combination
			ARE::EnumerateNextArgumentsValueCombination(_Width, vars, values, problem->K()) ;
			}
		bews->ApplyFnCombinationOperator(data[KeepIDX], const_factor) ;
		}

	return 0 ;
}


int32_t BucketElimination::MiniBucket::NoteOutputFunctionComputationCompletion(void)
{
	int32_t i ;

	MBEworkspace *bews = dynamic_cast<MBEworkspace*>(_Workspace) ;
	bool deleteMBfunctions = NULL != bews ? bews->DeleteUsedTables() : false ;

	if (deleteMBfunctions) {
		for (i = 0 ; i < _nFunctions ; i++) {
			ARE::Function *f = _Functions[i] ;
			if (NULL == f->OriginatingMiniBucket()) 
				continue ; // this is original fn, not MBE generated function.
			f->DestroyTableData() ;
			}
		}

	return 0 ;
}

