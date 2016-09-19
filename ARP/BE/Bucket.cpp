#include <stdlib.h>
#include <memory.h>

#include <Function.hxx>
#include <Bucket.hxx>
#include <MBEworkspace.hxx>


BucketElimination::Bucket::Bucket(void)
	:
	_Workspace(NULL),
	_IDX(-1), 
	_Width(-1), 
	_Signature(NULL), 
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
	_nChildren(-1), 
	_ChildVars(NULL), 
	_nOriginalFunctions(0), 
	_OriginalFunctions(NULL), 
	_OriginalWidth(-1), 
	_OriginalSignature(NULL), 
	_nAugmentedFunctions(0), 
	_AugmentedFunctions(NULL), 
	_AugmentedFunctionsArraySize(0), 
	_nIntermediateFunctions(0), 
	_IntermediateFunctions(NULL), 
	_IntermediateFunctionsArraySize(0), 
	_NextInOrderComputationGenList(NULL)
{
}


BucketElimination::Bucket::~Bucket(void)
{
	Destroy() ;
}


BucketElimination::Bucket::Bucket(MBEworkspace & WS, int32_t IDX, int32_t V)
	:
	_Workspace(&WS),
	_IDX(IDX), 
	_V(V),
	_Width(-1), 
	_Signature(NULL), 
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
	_nChildren(-1), 
	_ChildVars(NULL), 
	_nOriginalFunctions(0), 
	_OriginalFunctions(NULL), 
	_OriginalWidth(-1), 
	_OriginalSignature(NULL), 
	_nAugmentedFunctions(0), 
	_AugmentedFunctions(NULL), 
	_AugmentedFunctionsArraySize(0), 
	_nIntermediateFunctions(0), 
	_IntermediateFunctions(NULL), 
	_IntermediateFunctionsArraySize(0), 
	_NextInOrderComputationGenList(NULL)
{
	// if Var is given, add it
	if (V >= 0) 
		AddVar(V) ;
}


void BucketElimination::Bucket::Destroy(void)
{
	for (MiniBucket *mb : _MiniBuckets) {
		mb->Destroy() ;
		delete mb ;
		}
	_MiniBuckets.clear() ;

	if (NULL != _OriginalFunctions) {
		delete [] _OriginalFunctions ;
		_OriginalFunctions = NULL ;
		}
	if (NULL != _OriginalSignature) {
		delete [] _OriginalSignature ;
		_OriginalSignature = NULL ;
		}
	if (NULL != _AugmentedFunctions) {
		delete [] _AugmentedFunctions ;
		_AugmentedFunctions = NULL ;
		_AugmentedFunctionsArraySize = 0 ;
		}
	if (NULL != _IntermediateFunctions) {
		delete [] _IntermediateFunctions ;
		_IntermediateFunctions = NULL ;
		_IntermediateFunctionsArraySize = 0 ;
		}
	if (NULL != _Signature) {
		delete [] _Signature ;
		_Signature = NULL ;
		}
	_Width = -1 ;
	_nAugmentedFunctions = 0 ;
	_nIntermediateFunctions = 0 ;
	_OriginalWidth = -1 ;
	_nOriginalFunctions = 0 ;
	_nVars = 0 ;
	if (NULL != _Vars) {
		delete [] _Vars ;
		_Vars = NULL ;
		}
	_VarsSpace = 0 ;
	_MaxDescendantNumVars = -1 ;
	_ComputationNewFunctionSize = -1 ;
	_MaxDescendantComputationNewFunctionSize = -1 ;
	_nChildren = -1 ;
	_ChildVars = NULL ;
}

/*
int32_t BucketElimination::Bucket::SaveXMLString(const char *prefixspaces, const std::string & Dir, std::string & S)
{
	char s[1024] ;
	std::string temp ;
	int32_t i, j ;
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
			int32_t idx = f->IDX() ;
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
		int32_t size = (7 + nTB) >> 3 ;
		sprintf(s, "\n%s <ownbucketfncomputationresult nComputed=\"%I64d/%I64d\" bits=\"", prefixspaces, _nOutputFunctionBlocksComputed, nTB) ;
		S += s ;
		int32_t n = 0 ;
		for (i = 0 ; i < size && n < nTB ; i++) {
			for (j = 0 ; j < 8 && n < nTB ; j++, n++) {
				int32_t bit = _OutputFunctionBlockComputationResult[i] & (1 << j) ;
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
*/

int32_t BucketElimination::Bucket::SetOriginalFunctions(int32_t N, ARE::Function *FNs[]) 
{
	if (NULL != _OriginalFunctions) {
		delete [] _OriginalFunctions ;
		_OriginalFunctions = NULL ;
		}
	if (NULL != _OriginalSignature) {
		delete [] _OriginalSignature ;
		_OriginalSignature = NULL ;
		}
	_OriginalWidth = -1 ;
	_nOriginalFunctions = 0 ;
	return AddOriginalFunctions(N, FNs) ;
}


int32_t BucketElimination::Bucket::AddOriginalFunctions(int32_t N, ARE::Function *FNs[]) 
{
	if (N < 1) 
		return 0 ;

	int32_t i, j, k ;

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
	int32_t space = _nOriginalFunctions + N ;
	ARE::Function **fnlist = new ARE::Function*[space] ;
	if (NULL == fnlist) 
		return 1 ;
	for (i = 0 ; i < _nOriginalFunctions ; i++) 
		fnlist[i] = _OriginalFunctions[i] ;
	for (; i < space ; i++) {
		fnlist[i] = FNs[i-_nOriginalFunctions] ;
		fnlist[i]->SetBucket(this) ;
		}
	delete [] _OriginalFunctions ;
	_OriginalFunctions = fnlist ;
	_nOriginalFunctions = space ;

	if (NULL != _OriginalSignature) { delete [] _OriginalSignature ; _OriginalSignature = NULL ; }
	_OriginalWidth = -1 ;

	if (0 != ARE::ComputeSignature(N, FNs, _OriginalWidth, _OriginalSignature)) 
		goto failed ;

	return 0 ;

failed :
	Destroy() ;
	return 1 ;
}


int32_t BucketElimination::Bucket::AddAugmentedFunction(ARE::Function & F)
{
	// check if we have enough space
	if (_nAugmentedFunctions+1 > _AugmentedFunctionsArraySize) {
		// 2016-05-06 KK : make it 16 for reallocations (used to be 8) so that there are fewer reallocations
		int32_t newsize = 0 == _AugmentedFunctionsArraySize ? 4 : _AugmentedFunctionsArraySize + 16 ;
		ARE::Function **newspace = new ARE::Function*[newsize] ;
		if (NULL == newspace) 
			return 1 ;
		if (_nAugmentedFunctions > 0) 
			memcpy(newspace, _AugmentedFunctions, sizeof(ARE::Function *)*_nAugmentedFunctions) ;
		delete []_AugmentedFunctions ;
		_AugmentedFunctions = newspace ;
		_AugmentedFunctionsArraySize = newsize ;
		}

	_AugmentedFunctions[_nAugmentedFunctions++] = &F ;

	// if fn has arguments, current signature (may be) is invalid
	if (F.N() > 0) {
		_Width = -1 ;
		if (NULL != _Signature) { delete [] _Signature ; _Signature = NULL ; }
		}

	return 0 ;
}


int32_t BucketElimination::Bucket::RemoveAugmentedFunction(ARE::Function & F, bool InvalidateSignature)
{
	int32_t i ;
	int32_t n = 0 ;
	for (i = _nAugmentedFunctions - 1 ; i >= 0 ; i--) {
		if (&F == _AugmentedFunctions[i]) {
			_AugmentedFunctions[i] = _AugmentedFunctions[--_nAugmentedFunctions] ;
			F.SetBucket(NULL) ;
			n++ ;
			}
		}

	// if fn has arguments, current signature (may be) is invalid
	if (InvalidateSignature && n > 0 && F.N() > 0) {
		_Width = -1 ;
		if (NULL != _Signature) { delete [] _Signature ; _Signature = NULL ; }
		}

	return 0 ;
}


int32_t BucketElimination::Bucket::AddIntermediateFunction(ARE::Function & F)
{
	// check if we have enough space
	if (_nIntermediateFunctions+1 > _IntermediateFunctionsArraySize) {
		// 2016-05-06 KK : make reallcation size larger (used to be 8) so that there are fewer reallocations
		int32_t newsize = 5 ; 
		if (_IntermediateFunctionsArraySize >= 65) 
			newsize = _IntermediateFunctionsArraySize + 60 ;
		else if (_IntermediateFunctionsArraySize >= 25) 
			newsize = _IntermediateFunctionsArraySize + 40 ;
		else if (_IntermediateFunctionsArraySize > 0) 
			newsize = _IntermediateFunctionsArraySize + 20 ;
//		int32_t newsize = (0 == _IntermediateFunctionsArraySize) ? 8 : _IntermediateFunctionsArraySize + 32 ;
		ARE::Function **newspace = new ARE::Function*[newsize] ;
		if (NULL == newspace) 
			return 1 ;
		if (_nIntermediateFunctions > 0) 
			memcpy(newspace, _IntermediateFunctions, sizeof(ARE::Function *)*_nIntermediateFunctions) ;
		delete []_IntermediateFunctions ;
		_IntermediateFunctions = newspace ;
		_IntermediateFunctionsArraySize = newsize ;
		}

	_IntermediateFunctions[_nIntermediateFunctions++] = &F ;

	return 0 ;
}


int32_t BucketElimination::Bucket::RemoveIntermediateFunction(ARE::Function & F)
{
	int32_t i ;
	for (i = _nIntermediateFunctions - 1 ; i >= 0 ; i--) {
		if (&F == _IntermediateFunctions[i]) {
			_IntermediateFunctions[i] = _IntermediateFunctions[--_nIntermediateFunctions] ;
			F.SetBucket(NULL) ;
			}
		}
	return 0 ;
}


int32_t BucketElimination::Bucket::ComputeSignature(void)
{
	if (NULL != _Signature) {
		delete [] _Signature ;
		_Signature = NULL ;
		}
	_Width = -1 ;

	if (_OriginalWidth < 0) {
		if (NULL != _OriginalSignature) {
			delete [] _OriginalSignature ;
			_OriginalSignature = NULL ;
			}
		if (_nOriginalFunctions > 0) {
			if (0 != ARE::ComputeSignature(_nOriginalFunctions, _OriginalFunctions, _OriginalWidth, _OriginalSignature)) 
				return 1 ;
			}
		else if (_nVars > 0) {
			_OriginalSignature = new int32_t[_nVars] ;
			if (NULL == _OriginalSignature) 
				return 1 ;
			for (int32_t i = 0 ; i < _nVars ; i++) 
				_OriginalSignature[i] = _Vars[i] ;
			_OriginalWidth = _nVars ;
			}
		if (_OriginalWidth < 0) 
			return 1 ;
		}

	if (_OriginalWidth <= 0 && _nAugmentedFunctions <= 0) 
		{ _Width = 0 ; return 0 ; }

	int nF = _nOriginalFunctions + _nAugmentedFunctions ;
	if (nF > 0) {
		std::vector<ARE::Function *> FNs ; 
		FNs.reserve(nF) ;
		if (FNs.capacity() != nF) 
			return 1 ;
		for (int32_t i = 0 ; i < _nOriginalFunctions ; i++) 
			FNs.push_back(_OriginalFunctions[i]) ;
		for (int32_t i = 0 ; i < _nAugmentedFunctions ; i++) 
			FNs.push_back(_AugmentedFunctions[i]) ;

		if (NULL != _Signature) { delete [] _Signature ; _Signature = NULL ; }
		_Width = -1 ;

		if (0 != ARE::ComputeSignature(FNs.size(), FNs.data(), _Width, _Signature)) 
			return 1 ;
		}
	else if (_nVars > 0) {
		_Signature = new int32_t[_nVars] ;
		if (NULL == _Signature) 
			return 1 ;
		for (int32_t i = 0 ; i < _nVars ; i++) 
			_Signature[i] = _Vars[i] ;
		_Width = _nVars ;
		}

	return 0 ;
}


int32_t BucketElimination::Bucket::ComputeOutputFunctionWithScopeWithoutTable(int32_t * & TempSpaceForArglist, int32_t TempSpaceForArglistSize, ARE::Function * & FN, int32_t & max_var)
{
	max_var = -1 ;
	if (NULL != FN) {
		FN->Destroy() ;
		FN->SetIDX(-(_V+1)) ;
		}
	if (NULL == _Workspace) 
		return 1 ;

	if (_Width < 0) {
		int32_t res = ComputeSignature() ;
		if (0 != res) {
//			if (createdFN && NULL != FN) 
//				{ delete FN ; FN = NULL ; }
			return 1 ;
			}
		}
	if (_Width <= 0) 
		return 0 ; // bucket has no variables; that should not happen
	// if _Width - _nVars <= 0, then all variables are eliminated; output is a const fn.

	bool createdFN = false ;
	if (NULL == FN) {
		FN = new ARE::Function(_Workspace, NULL != _Workspace ? _Workspace->Problem() : NULL, -(_V+1)) ;
		if (NULL == FN) 
			return 1 ;
		FN->SetOriginatingBucket(this) ;
		createdFN = true ;
		}

	// check memory is ok for sorting
	if (TempSpaceForArglistSize < _Width) {
		if (NULL != TempSpaceForArglist) 
			{ delete [] TempSpaceForArglist ; TempSpaceForArglist = NULL ; TempSpaceForArglistSize = 0 ; }
		TempSpaceForArglist = new int32_t[_Width] ;
		if (NULL == TempSpaceForArglist) 
			return 1 ;
		TempSpaceForArglistSize = _Width ;
		}

	// construct args list; find max argument.
	const int32_t *varpos = _Workspace->VarPos() ;
	int32_t nArgs = 0 ;
	for (int32_t i = 0 ; i < _Width ; i++) {
		int32_t v = _Signature[i], j ;
		for (j = 0 ; j < _nVars ; j++) 
			{ if (v == _Vars[j]) break ; }
		if (j < _nVars) 
			continue ;
		TempSpaceForArglist[nArgs++] = v ;
		if (max_var < 0) 
			max_var = v ;
		else if (varpos[max_var] < varpos[v]) 
			max_var = v ;
		}

	// set up fn; nArgs should be _Width - _nVars
	FN->SetArguments(nArgs, TempSpaceForArglist) ;

	return 0 ;
}


/*
int32_t BucketElimination::Bucket::ReorderFunctionScopesWrtParentBucket(bool IncludeOriginalFunctions, bool IncludeNewFunctions)
{
	// destroy table of all childbucketfunctions of this bucket;
	// reordering of scopes should be done in the beginning when these tables don't exist yet, 
	// so this should not be a problem.
	int32_t i ;
	if (IncludeNewFunctions) {
		for (i = 0 ; i < _nAugmentedFunctions ; i++) {
			ARE::Function *f = _AugmentedFunctions[i] ;
			f->DestroyTableData() ;
			}
		}

	// reorder scopes
	int32_t nOutputVars = _Width - _nVars ;
	if (nOutputVars < 1) {
		if (IncludeNewFunctions) {
			for (i = 0 ; i < _nAugmentedFunctions ; i++) {
				if (0 != _AugmentedFunctions[i]->ReOrderArguments(_nVars, _Vars, 0, NULL)) 
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
		for (MiniBucket & mb : _MiniBuckets) {
			}
		if (IncludeNewFunctions) {
			for (i = 0 ; i < _nAugmentedFunctions ; i++) {
				if (0 != _AugmentedFunctions[i]->ReOrderArguments(nOutputVars, _OutputFunction.Arguments(), _nVars, _Vars)) 
					return 1 ;
				}
			}
		if (IncludeOriginalFunctions) {
			for (i = 0 ; i < _nOriginalFunctions ; i++) {
				if (0 != _OriginalFunctions[i]->ReOrderArguments(nOutputVars, _OutputFunction.Arguments(), _nVars, _Vars)) 
					return 1 ;
				}
			}
		}

	return 0 ;
}
*/

int32_t BucketElimination::Bucket::NoteOutputFunctionComputationCompletion(void)
{
	for (MiniBucket *mb : _MiniBuckets) {
		mb->NoteOutputFunctionComputationCompletion() ;
		}
	return 0 ;
}


__int64 BucketElimination::Bucket::ComputeProcessingComplexity(void)
{
	__int64 N = 0 ;
	for (MiniBucket *mb : _MiniBuckets) {
		__int64 n = mb->ComputeProcessingComplexity() ;
		if (n > 0) 
			N += n ;
		}
	return N ;
}


int32_t BucketElimination::Bucket::ComputeOutputFunctions(bool DoMomentMatching)
{
	// this is the avg max-marginals fn
	ARE::Function fAvgMM(_Workspace, _Workspace->Problem(), _MiniBuckets.size()) ;
	double *average_mm_table = NULL ;
	// this is max-marginals for each mb
	ARE::Function **max_marginals = NULL ;
	// this is joint-scope of all mb output functions
	int32_t *joint_scope = NULL, joint_scope_size = 0 ;

	int32_t res = 1 ;
	int32_t idx ;
	if (_MiniBuckets.size() > 1 && DoMomentMatching) {
		// make sure width is computed
		int32_t min_width = INT_MAX, max_width = -INT_MAX ;
		for (MiniBucket *mb : _MiniBuckets) {
			if (mb->Width() < 0) {
				if (0 != mb->ComputeSignature()) 
					return 1 ;
				if (mb->Width() < 0) 
					return 1 ;
				}
			if (min_width > mb->Width()) 
				min_width = mb->Width() ;
			if (max_width < mb->Width()) 
				max_width = mb->Width() ;
			}
		if (min_width <= 0 || max_width <= 0) 
			return 1 ;

		// compute joint scope of all MBs; next line, allocate some more space for other arrays
		joint_scope = new int32_t[min_width + max_width + max_width] ; if (NULL == joint_scope) goto done_MM ;
		int32_t *temp_scope = joint_scope + min_width, temp_scope_size ; 
		int32_t *temp2_scope = temp_scope + max_width ;
		bool first_mb = true ;
		for (MiniBucket *mb : _MiniBuckets) {
			const int32_t *sig = mb->SortedSignature() ;
			if (first_mb) 
				{ for (int32_t i = 0 ; i < mb->Width() ; i++) joint_scope[i] = sig[i] ; joint_scope_size = mb->Width() ; first_mb = false ; }
			else 
				{ ARE::SetIntersection(joint_scope, joint_scope_size, sig, mb->Width()) ; }
			if (0 == joint_scope_size) 
				goto done_MM ; // error; this is not supposed to happen
			}

		// compute max-marginals for each MB.
		max_marginals = new ARE::Function*[_MiniBuckets.size()] ; if (NULL == max_marginals) goto done_MM ;
		memset(max_marginals, 0, sizeof(ARE::Function*)*_MiniBuckets.size()) ;
		int32_t idx = 0 ;
		for (MiniBucket *mb : _MiniBuckets) {
			const int32_t *sig = mb->SortedSignature() ;
			for (int32_t i = 0 ; i < mb->Width() ; i++) temp_scope[i] = sig[i] ; temp_scope_size = mb->Width() ;
			ARE::SetMinus(temp_scope, temp_scope_size, joint_scope, joint_scope_size) ;
			max_marginals[idx] = new ARE::Function(_Workspace, _Workspace->Problem(), idx) ;
			if (NULL == max_marginals[idx]) 
				goto done_MM ;
			if (0 != mb->ComputeOutputFunction(*(max_marginals[idx]), temp_scope, temp_scope_size, temp2_scope)) 
				goto done_MM ;
			++idx ;
			}
		// NOTE : all max_marginals[*] have same scope and the order of arguments is the same (i.e. argument lists are identical).

		// compute avg max-marginals
		INT64 table_size = max_marginals[0]->TableSize() ;
		if (table_size <= 0) 
			goto done_MM ;
		average_mm_table = new double[table_size] ;
		if (NULL == average_mm_table) 
			goto done_MM ;
		average_mm_table[0] = _Workspace->FnCombinationNeutralValue() ;
		for (int32_t i = 1 ; i < table_size ; i++) average_mm_table[i] = average_mm_table[0] ;
		for (int32_t j = _MiniBuckets.size() - 1 ; j >= 0 ; j--) {
			for (int32_t i = 0; i < table_size ; i++) 
				_Workspace->ApplyFnCombinationOperator(average_mm_table[i], max_marginals[j]->TableEntry(i)) ;
			}
		double N_ = (double) _MiniBuckets.size() ;
		for (int32_t i = 0; i < table_size ; i++) 
			average_mm_table[i] /= N_ ;

		if (0 != fAvgMM.SetArguments(max_marginals[0]->N(), max_marginals[0]->Arguments())) 
			goto done_MM ;
		if (0 != fAvgMM.SetTableData(table_size, average_mm_table)) 
			goto done_MM ;
		average_mm_table = NULL ;
		}

	idx = 0 ;
	for (MiniBucket *mb : _MiniBuckets) {
		if (NULL != max_marginals) 
			mb->ComputeOutputFunction(max_marginals[idx], &fAvgMM) ;
		else 
			mb->ComputeOutputFunction(NULL, NULL) ;
		++idx ;
		}

	// done with MM; delete stuff.
	res = 0 ;
done_MM :
	if (NULL != average_mm_table) delete [] average_mm_table ;
	if (NULL != max_marginals) {
		for (int32_t j = 0 ; j < _MiniBuckets.size() ; j++) { if (NULL != max_marginals[j]) delete max_marginals[j] ; }
		delete [] max_marginals ;
		}
	if (NULL != joint_scope) 
		delete [] joint_scope ;
	return res ;
}


int32_t BucketElimination::Bucket::CreateMBPartitioning(int32_t iBound, bool CreateTables, bool doMomentMatching, bool & AbandonIfActualPartitioning, std::vector<int32_t> & key, std::vector<int64_t> & data, std::vector<int32_t> & helperArray)
{
	bool abandonIfActualPartitioning = AbandonIfActualPartitioning ;
	AbandonIfActualPartitioning = false ;

//INT64 tB = ARE::GetTimeInMilliseconds() ;

	// sort functions in the order of decreasing scope size
	int32_t nF = nOriginalFunctions() + nAugmentedFunctions() ;
	if (key.capacity() < nF) { key.reserve(nF) ; if (key.capacity() < nF) goto failed ; }
	if (data.capacity() < nF) { data.reserve(nF) ; if (data.capacity() < nF) goto failed ; }
	key.clear() ; data.clear() ;
	for (int32_t j = _nOriginalFunctions - 1 ; j >= 0 ; j--) {
		ARE::Function *f = _OriginalFunctions[j] ;
		key.push_back(-f->N()) ;
		data.push_back((INT64) f) ;
		if (data.size() != key.size()) {
			fprintf(ARE::fpLOG, "\n   BucketElimination::Bucket::CreateMBPartitioning ERROR : data_size!=key_size(1); i=%d, v=%d, nMBs=%d", (int32_t) _IDX, (int32_t) _V, (int32_t) _MiniBuckets.size()) ;
			::fflush(ARE::fpLOG) ;
			goto failed ;
			}
		}
	for (int32_t j = _nAugmentedFunctions - 1 ; j >= 0 ; j--) {
		ARE::Function *f = _AugmentedFunctions[j] ;
		key.push_back(-f->N()) ;
		data.push_back((INT64) f) ;
		if (data.size() != key.size()) {
			fprintf(ARE::fpLOG, "\n   BucketElimination::Bucket::CreateMBPartitioning ERROR : data_size!=key_size(2); i=%d, v=%d, nMBs=%d", (int32_t) _IDX, (int32_t) _V, (int32_t) _MiniBuckets.size()) ;
			::fflush(ARE::fpLOG) ;
			goto failed ;
			}
		}
	int32_t left[32], right[32] ;
	QuickSortLong_i64(key.data(), key.size(), data.data(), left, right) ;

/*
INT64 tS = ARE::GetTimeInMilliseconds() ;
if (_IDX < 2000) {
INT64 dt = tS - tB ;
fprintf(ARE::fpLOG, "\n   BucketElimination::Bucket::CreateMBPartitioning; i=%d, v=%d, iBound=%d, nFNs=%d, sorting time=%lld[msec]", (int32_t) _IDX, (int32_t) _V, (int32_t) iBound, (int32_t) data.size(), dt) ;
::fflush(ARE::fpLOG) ; }*/

#ifdef _DEBUG
		if (_MiniBuckets.size() > 0) {
			int32_t error = 1 ;
			}
		_MiniBuckets.clear() ;
#endif

//INT64 dtSUM1 = 0, dtSUM2 = 0, dtSUM3 = 0 ; ;

		// process functions, one at a time; for each fn, try to add to an existing MB.
		for (int32_t j = 0 ; j < key.size() ; j++) {
			ARE::Function *f = (ARE::Function*) data[j] ;
			bool placed = false ;
			for (MiniBucket *mb : _MiniBuckets) {
//INT64 t1 = ARE::GetTimeInMilliseconds() ;
				int32_t res = mb->AllowsFunction(*f, iBound, helperArray) ;
//INT64 t2 = ARE::GetTimeInMilliseconds() ;
//dtSUM1 += t2 - t1 ;
				if (res < 0) 
					goto failed ;
				if (res > 0) {
					placed = true ;
//INT64 t1_ = ARE::GetTimeInMilliseconds() ;
					if (0 != mb->AddFunction(*f, helperArray)) 
						goto failed ;
//INT64 t2_ = ARE::GetTimeInMilliseconds() ;
//dtSUM2 += t2_ - t1_ ;
					break ;
					}
				}
			if (! placed) {
				if (_MiniBuckets.size() > 0 && abandonIfActualPartitioning) {
					AbandonIfActualPartitioning = true ;
					goto failed ;
					}
//INT64 t3 = ARE::GetTimeInMilliseconds() ;
				MiniBucket *mb = new MiniBucket ;
				if (NULL == mb) 
					goto failed ;
				int32_t existing_size = _MiniBuckets.size() ;
				_MiniBuckets.push_back(mb) ;
				if (_MiniBuckets.size() != (1+existing_size)) 
					{ delete mb ; goto failed ; }
				mb->Initalize(*this, existing_size) ;
				for (int32_t idxV = 0 ; idxV < _nVars ; idxV++) 
					mb->AddVar(_Vars[idxV]) ;
				if (0 != mb->AddFunction(*f, helperArray)) 
					goto failed ;
//INT64 t4 = ARE::GetTimeInMilliseconds() ;
//dtSUM3 += t4 - t3 ;
				}
			}

/*INT64 tMB = ARE::GetTimeInMilliseconds() ;
int32_t div = _IDX%1000 ;
if (_IDX < 2000) {
	INT64 dt = tMB - tS ;
	fprintf(ARE::fpLOG, "\n   BucketElimination::Bucket::CreateMBPartitioning; i=%d, v=%d, nMBs=%d, nFNs=%d, MB partitioning time=%lld(check=%lld, addfn=%lld, new=%lld)", (int32_t) _IDX, (int32_t) _V, (int32_t) _MiniBuckets.size(), (int32_t) data.size(), dt, dtSUM1, dtSUM2, dtSUM3) ;
	::fflush(ARE::fpLOG) ;
	}*/

	// minibuckets for current bucket are now ready, process each and place resulting function
//dtSUM1 = dtSUM2 = dtSUM3 = 0 ; ;
	for (MiniBucket *mb : _MiniBuckets) {
		ARE::Function & f = mb->OutputFunction() ;
//INT64 t1 = ARE::GetTimeInMilliseconds() ;
		if (0 != mb->ComputeOutputFunctionWithScopeWithoutTable()) {
			fprintf(ARE::fpLOG, "\n   BucketElimination::Bucket::CreateMBPartitioning ERROR : mb->ComputeOutputFunctionWithScopeWithoutTable() failed ...") ;
			::fflush(ARE::fpLOG) ;
			goto failed ;
			}
//INT64 t2 = ARE::GetTimeInMilliseconds() ;
//dtSUM1 += t2 - t1 ;
		Bucket *target_bucket = f.Bucket() ; // note target_bucket may be NULL (i.e. when 0==f.N())
		for (Bucket *mb_b = ParentBucket() ; mb_b != target_bucket && NULL != mb_b ; mb_b = mb_b->ParentBucket()) {
//t1 = ARE::GetTimeInMilliseconds() ;
//dtSUM1 += t2 - t1 ;
			if (0 != mb_b->AddIntermediateFunction(f)) {
				fprintf(ARE::fpLOG, "\n   BucketElimination::Bucket::CreateMBPartitioning ERROR : mb_b->AddIntermediateFunction() failed ...") ;
				::fflush(ARE::fpLOG) ;
				goto failed ;
				}
//t2 = ARE::GetTimeInMilliseconds() ;
//dtSUM2 += t2 - t1 ;
			}
		if (NULL != target_bucket) {
//t1 = ARE::GetTimeInMilliseconds() ;
			if (0 != target_bucket->AddAugmentedFunction(f)) {
				fprintf(ARE::fpLOG, "\n   BucketElimination::Bucket::CreateMBPartitioning ERROR : target_bucket->AddAugmentedFunction() failed ...") ;
				::fflush(ARE::fpLOG) ;
				goto failed ;
				}
//t2 = ARE::GetTimeInMilliseconds() ;
//dtSUM3 += t2 - t1 ;
			}
		}

/*INT64 tE = ARE::GetTimeInMilliseconds() ;
if (_IDX < 2000) {
INT64 dt = tE - tMB ;
fprintf(ARE::fpLOG, "\n   BucketElimination::Bucket::CreateMBPartitioning; i=%d, v=%d, MB placing  time=%lld[msec] compOFN=%lld +IFN=%lld +AFN=%lld", (int32_t) _IDX, (int32_t) _V, dt, dtSUM1, dtSUM2, dtSUM3) ;
::fflush(ARE::fpLOG) ; }*/

	// if requested, do moment-matching
	if (CreateTables) {
//INT64 t1 = ARE::GetTimeInMilliseconds() ;
		if (0 != ComputeOutputFunctions(doMomentMatching)) {
			fprintf(ARE::fpLOG, "\n   BucketElimination::Bucket::CreateMBPartitioning ERROR : mb->ComputeOutputFunction() failed ...") ;
			::fflush(ARE::fpLOG) ;
			goto failed ;
			}
/*INT64 t2 = ARE::GetTimeInMilliseconds() ;
INT64 dt = t2 - t1 ;
if (_IDX < 2000) {
INT64 dt = t2 - t1 ;
fprintf(ARE::fpLOG, "\n   BucketElimination::Bucket::CreateMBPartitioning; i=%d, v=%d, ComputeOutputFunctions doMomentMatching=%c time=%lld[msec]", (int32_t) _IDX, (int32_t) _V, doMomentMatching ? 'Y' : 'N', dt) ;
::fflush(ARE::fpLOG) ; }*/
		}

	return 0 ;
failed :
	fprintf(ARE::fpLOG, "\n   BucketElimination::Bucket::CreateMBPartitioning ERROR : failed label reached (AbandonIfActualPartitioning=%c) ...", AbandonIfActualPartitioning ? 'Y' : 'N') ;
	::fflush(ARE::fpLOG) ;
	DestroyPartitioning() ;
	return 1 ;
}


int32_t BucketElimination::Bucket::ComputeFirstVariableDistribution(ARE_Function_TableType *dist)
{
	int32_t i, j, k, ret = 0 ;

	MBEworkspace *bews = dynamic_cast<MBEworkspace*>(_Workspace) ;
	if (NULL == bews) 
		return ERRORCODE_generic ;
	ARE::ARP *problem = bews->Problem() ;
	if (NULL == problem) 
		return ERRORCODE_generic ;
	ARE_Function_TableType nv = bews->FnCombinationNeutralValue() ;
	int32_t w = Width() ;
	if (w < 0) {
		ComputeSignature() ;
		w = Width() ;
		if (w < 0) 
			return 0 ; // should not happen; means bucket has no functions/variables.
		}
	const int32_t *signature = Signature() ;
	const int32_t nOF = nOriginalFunctions() ;
	const int32_t nAF = nAugmentedFunctions() ;
	const int32_t nIF = nIntermediateFunctions() ;
	int32_t nFunctions_OA = nOF + nAF ;
	int32_t nTotalFunctions = nOF + nAF + nIF ;
	if (w < 0 || nTotalFunctions < 1) {
		return 0 ;
		}

	if (w > MAX_NUM_VARIABLES_PER_BUCKET) 
		return ERRORCODE_too_many_variables ;
	if (nTotalFunctions > MAX_NUM_FUNCTIONS_PER_BUCKET) 
		return ERRORCODE_too_many_functions ;

	int32_t values[MAX_NUM_VARIABLES_PER_BUCKET] ; // this is the current value combination of the arguments of this function (table block).
	ARE::Function *flist[MAX_NUM_FUNCTIONS_PER_BUCKET] ; // this is a list of input functions

	ARE_Function_TableType const_factor = bews->FnCombinationNeutralValue() ;
	int32_t nFNs = 0 ;
	for (j = 0 ; j < nOF ; j++) {
		ARE::Function *f = OriginalFunction(j) ;
		if (NULL == f) continue ;
		if (0 == f->N()) bews->ApplyFnCombinationOperator(const_factor, f->ConstValue()) ;
		else { flist[nFNs++] = f ; f->ComputeArgumentsPermutationList(w, signature) ; }
		}
	for (; j < nTotalFunctions ; j++) {
		ARE::Function *f = AugmentedFunction(j - nOF) ;
		if (NULL == f) continue ;
		if (0 == f->N()) bews->ApplyFnCombinationOperator(const_factor, f->ConstValue()) ;
		else { flist[nFNs++] = f ; f->ComputeArgumentsPermutationList(w, signature) ; }
		}
	for (; j < nTotalFunctions ; j++) {
		ARE::Function *f = IntermediateFunction(j - nFunctions_OA) ;
		if (NULL == f) continue ;
		if (0 == f->N()) bews->ApplyFnCombinationOperator(const_factor, f->ConstValue()) ;
		else { flist[nFNs++] = f ; f->ComputeArgumentsPermutationList(w, signature) ; }
		}

	INT64 ElimSize = 1 ;
	for (j = 1 ; j < w ; j++) 
		ElimSize *= problem->K(signature[j]) ;

	for (j = 0 ; j < w ; j++) 
		values[j] = 0 ;
	int32_t v0 = _Vars[0] ;
	for (i = 0 ; i < problem->K(v0) ; i++) {
		ARE_Function_TableType V = bews->VarEliminationDefaultValue() ;
		for (INT64 ElimIDX = 0 ; ElimIDX < ElimSize ; ElimIDX++) {
			ARE_Function_TableType value = nv ;
			for (j = 0 ; j < nFNs ; j++) { // note : it should be that 0 != flist[j]->N(). note : it is assumed that flist[j] has 1 block. note : it is assumed that order of flist[j] arguments is the same as signature.
				INT64 adr = flist[j]->ComputeFnTableAdr_wrtLocalPermutation(w, values, problem->K()) ;
				bews->ApplyFnCombinationOperator(value, flist[j]->TableEntry(adr)) ;
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


int32_t BucketElimination::Bucket::ComputeFirstVariableDistributionEx(int32_t *ContextValues, ARE_Function_TableType *dist)
{
	int32_t j ;

	if (1 != _nVars) 
		return ERRORCODE_generic ;

	MBEworkspace *bews = dynamic_cast<MBEworkspace*>(_Workspace) ;
	if (NULL == bews) 
		return ERRORCODE_generic ;
	ARE::ARP *problem = bews->Problem() ;
	if (NULL == problem) 
		return ERRORCODE_generic ;
	ARE_Function_TableType nv = bews->FnCombinationNeutralValue() ;

	const int32_t nOF = nOriginalFunctions() ;
	const int32_t nAF = nAugmentedFunctions() ;
	const int32_t nIF = nIntermediateFunctions() ;
	int32_t nTotalFunctions = nOF + nAF + nIF ;
	if (nTotalFunctions < 1) {
		return 0 ;
		}

	std::vector<ARE::Function *> flist ; // this is a list of input functions
	flist.reserve(nTotalFunctions) ;
	if (flist.capacity() != nTotalFunctions) 
		return ERRORCODE_out_of_memory ;

	ARE_Function_TableType const_factor = bews->FnCombinationNeutralValue() ;
	for (j = 0 ; j < nOF ; j++) {
		ARE::Function *f = OriginalFunction(j) ;
		if (NULL == f) continue ;
		if (0 == f->N()) bews->ApplyFnCombinationOperator(const_factor, f->ConstValue()) ;
		else { flist.push_back(f) ; }
		}
	for (j = 0 ; j < nAF ; j++) {
		ARE::Function *f = AugmentedFunction(j) ;
		if (NULL == f) continue ;
		if (0 == f->N()) bews->ApplyFnCombinationOperator(const_factor, f->ConstValue()) ;
		else { flist.push_back(f) ; }
		}
	for (j = 0 ; j < nIF ; j++) {
		ARE::Function *f = IntermediateFunction(j) ;
		if (NULL == f) continue ;
		if (0 == f->N()) bews->ApplyFnCombinationOperator(const_factor, f->ConstValue()) ;
		else {
			// all variables in the scope should be assigned
			INT64 adr = f->ComputeFnTableAdr(ContextValues, problem->K()) ;
			bews->ApplyFnCombinationOperator(const_factor, f->TableEntry(adr)) ;
			}
		}

	int32_t val_bucket_var = ContextValues[_V] ;
	for (int32_t i = 0 ; i < problem->K(_V) ; i++) {
		ContextValues[_V] = i ;
		ARE_Function_TableType V = const_factor ;
		for (j = flist.size() - 1 ; j >= 0 ; j--) { // note : it should be that 0 != flist[j]->N(). note : it is assumed that flist[j] has 1 block. note : it is assumed that order of flist[j] arguments is the same as signature.
			INT64 adr = flist[j]->ComputeFnTableAdr(ContextValues, problem->K()) ;
			bews->ApplyFnCombinationOperator(V, flist[j]->TableEntry(adr)) ;
			}
		dist[i] = V ;
		}
	ContextValues[_V] = val_bucket_var ;

	return 0 ;
}


int32_t BucketElimination::Bucket::DestroyPartitioning(void)
{
	for (MiniBucket *mb : _MiniBuckets) {
		ARE::Function & f = mb->OutputFunction() ;
		BucketElimination::Bucket *B = f.Bucket() ;
		if (NULL == B) 
			continue ; // only way this is possible is that this fn is const fn
		// remove output function from their buckets.
		for (BucketElimination::Bucket *b = ParentBucket() ; NULL != b && B != b ; b = b->ParentBucket()) {
			b->RemoveIntermediateFunction(f) ;
			}
		B->RemoveAugmentedFunction(f, true) ;
		// destroy output functio; so that its table can be released.
		f.Destroy() ;
		// destroy mb
		mb->Destroy() ;
		delete mb ;
		}
	_MiniBuckets.clear() ;

	return 0 ;
}

