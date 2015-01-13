#include <stdlib.h>
#include <time.h>
#include <math.h>

#if defined WINDOWS || _WINDOWS
#include <process.h>    /* _beginthread, _endthread */
#endif // WINDOWS

#include "Utils/Sort.hxx"

#include <Problem/Globals.hxx>
#include <Utils/MiscUtils.hxx>
#include <Function.hxx>
#include <BEworkspace.hxx>

#if defined WINDOWS || _WINDOWS
typedef unsigned int (*pBEWSThreadFn)(void *X) ;
static unsigned int __stdcall BEWSThreadFn(void *X) 
#elif defined (LINUX)
typedef void *(*pBEWSThreadFn)(void *X) ;
static void *BEWSThreadFn(void *X)
#endif 
{
	BucketElimination::BEworkspace *bews = (BucketElimination::BEworkspace *)(X) ;

	bews->SetCompleteEliminationResult(DBL_MAX) ;
	bews->MarginalSingleVariableDistributionK() = 0 ;
	bews->MarginalSingleVariableDistributionVar() = -1 ;

	long stop_signalled = 0 ;
	int n = bews->nBuckets() - 1 ;
	bews->tStart() = ARE::GetTimeInMilliseconds() ;
	while (bews->IsValid()) {
		stop_signalled = InterlockedCompareExchange(&(bews->_StopAndExit), 1, 1) ;
		if (0 != stop_signalled) {
			break ;
			}
		BucketElimination::Bucket *b = bews->getBucket(n) ;
		if (NULL == b) {
			if (NULL != bews->logFile()) {
				fprintf(bews->logFile(), "\n   BE elimination n=%d, WARNING : bucket==NULL", (int) n) ;
				}
			continue ;
			}
		int res = b->ComputeOutputFunction_1Block() ;
/*		if (NULL != bews->logFile()) {
			ARE::Function & f = b->OutputFunction() ;
			if (0 == f.N()) 
				fprintf(bews->logFile(), "\n   BE elimination var=%d, const-output=%g", b->Var(0), f.ConstValue()) ;
			else if (1 == f.N()) 
				fprintf(bews->logFile(), "\n   BE elimination var=%d, output=%g,%g", b->Var(0), f.Table()->Data()[0], f.Table()->Data()[1]) ;
			else if (2 == f.N()) 
				fprintf(bews->logFile(), "\n   BE elimination var=%d, output=%g,%g,%g,%g", b->Var(0), f.Table()->Data()[0], f.Table()->Data()[1], f.Table()->Data()[2], f.Table()->Data()[3]) ;
			}*/
		if (0 != res) {
			// failed; abandon.
			if (NULL != bews->logFile()) {
				fprintf(bews->logFile(), "\n   BE elimination n=%d v=%d, ERROR : ComputeOutputFunction_1Block() returned error=%d; will quit ...", (int) n, (int) b->Var(0), (int) res) ;
				}
			break ;
			}
		INT64 table_size = b->OutputFunction().TableSize() ;
		bews->TotalNewFunctionSizeComputed() += table_size >= 0 ? table_size : 0 ;
		if (--n < 0) {// all computed
			// compose complete_elimination_answer
			ARE_Function_TableType v = bews->AnswerFactor() ;
			ARE_Function_TableType v_all_except_first = v ;
			for (int i = bews->nBuckets() - 1 ; i >= 0 ; i--) {
				BucketElimination::Bucket *b = bews->getBucket(i) ;
				if (0 == b->DistanceToRoot()) {
					ARE::Function & f = b->OutputFunction() ;
					bews->ApplyFnCombinationOperator(v, f.ConstValue()) ;
					if (0 != i) 
						v_all_except_first = v ;
					}
				}
			bews->SetCompleteEliminationResult(v) ;
			// compose var[0] distribution
			ARE_Function_TableType *dist = bews->MarginalSingleVariableDistribution() ;
			BucketElimination::Bucket *b_0 = bews->getBucket(0) ;
			b_0->ComputeFirstVariableDistribution_1Block(dist) ;
			int v_0 = *(b_0->Signature()) ;
			int k_0 = bews->Problem()->K(v_0) ;
			for (int i = 0 ; i < k_0 ; i++) {
				bews->ApplyFnCombinationOperator(dist[i], v_all_except_first) ;
				}
			bews->MarginalSingleVariableDistributionK() = k_0 ;
			bews->MarginalSingleVariableDistributionVar() = v_0 ;
			break ;
			}
		}

done :
	bews->tEnd() = ARE::GetTimeInMilliseconds() ;
	bews->RunTimeInMilliseconds() = bews->tEnd() - bews->tStart() ;
	bews->_ThreadHandle = 0 ;
#if defined WINDOWS || _WINDOWS
	_endthreadex(0) ;
	return 0  ;
#else
	return NULL ;
#endif
}


BucketElimination::BEworkspace::BEworkspace(const char *BEEMDiskSpaceDirectory)
	:
	ARE::Workspace(BEEMDiskSpaceDirectory), 
	_VarOrder(NULL), 
	_VarPos(NULL), 
	_Var2BucketMapping(NULL), 
	_DeleteUsedTables(false), 
	_fpLOG(NULL), 
	_AnswerFactor(1.0), 
	_CompleteEliminationResult(DBL_MAX), 
	_tStart(0), 
	_tEnd(0), 
	_tToStop(0), 
	_RunTimeInMilliseconds(-1), 
	_StopAndExit(0), 
	_ThreadHandle(0), 
	_FnCombinationType(FN_COBINATION_TYPE_NONE), 
	_VarEliminationType(VAR_ELIMINATION_TYPE_NONE), 
	_nBuckets(0), 
	_MaxNumChildren(-1), 
	_MaxNumVarsInBucket(-1),
	_nBucketsWithSingleChild_initial(-1), 
	_nBucketsWithSingleChild_final(-1), 
	_nBuckets_initial(-1), 
	_nBuckets_final(-1), 
	_nVarsWithoutBucket(-1), 
	_nConstValueFunctions(-1), 
	_MaxBucketFunctionWidth(-1), 
	_nBucketsWithSingleChild(-1), 
	_TotalOriginalFunctionSize(-1), 
	_TotalOriginalFunctionSpace(-1), 
	_TotalNewFunctionSize(-1), 
	_TotalNewFunctionSpace(-1), 
	_TotalNewFunctionComputationComplexity(-1), 
	_MaxSimultaneousNewFunctionSize(-1), 
	_MaxSimultaneousNewFunctionSpace(-1), 
	_MaxTotalFunctionSize(-1), 
	_MaxTotalFunctionSpace(-1), 
	_TotalNewFunctionSizeComputed(0), 
	_BucketOrderToCompute(NULL) 
{
	if (! _IsValid) 
		return ;

	for (int i = 0 ; i < MAX_NUM_BUCKETS ; i++) 
		_Buckets[i] = NULL ;

	for (int i = 0 ; i < MAX_NUM_VALUES_PER_VAR_DOMAIN ; i++) 
		_MarginalSingleVariableDistribution[i] = DBL_MAX ;
	_MarginalSingleVariableDistributionK = 0 ;
	_MarginalSingleVariableDistributionVar = -1 ;

	_IsValid = true ;
}


int BucketElimination::BEworkspace::Destroy(void)
{
	StopThread() ;

	if (_VarOrder) {
		delete [] _VarOrder ;
		_VarOrder = NULL ;
		}
	if (_VarPos) {
		delete [] _VarPos ;
		_VarPos = NULL ;
		}
	if (_Var2BucketMapping) {
		delete [] _Var2BucketMapping ;
		_Var2BucketMapping = NULL ;
		}
	if (NULL != _BucketOrderToCompute) {
		delete [] _BucketOrderToCompute ;
		_BucketOrderToCompute = NULL ;
		}
	_nVars = 0 ;
	for (int i = 0 ; i < MAX_NUM_BUCKETS ; i++) {
		if (NULL != _Buckets[i]) {
			delete _Buckets[i] ;
			_Buckets[i] = NULL ;
			}
		}
	_nBuckets = 0 ;
	_nVarsWithoutBucket = -1 ;

	Workspace::Destroy() ;

	_tStart = _tEnd = _tToStop = 0 ;
	_RunTimeInMilliseconds = -1 ;

	_CompleteEliminationResult = DBL_MAX ;
	_MarginalSingleVariableDistributionK = 0 ;
	_MarginalSingleVariableDistributionVar = -1 ;

	return 0 ;
}


int BucketElimination::BEworkspace::CreateThread(void)
{
#if defined WINDOWS || _WINDOWS
	_ThreadHandle = _beginthreadex(NULL, 0, BEWSThreadFn, this, 0, NULL) ;
#else
	pthread_create(&_ThreadHandle, NULL, BEWSThreadFn, this) ; // TODO third argument
#endif
	return 0 != _ThreadHandle ? 0 : 1 ;
}


int BucketElimination::BEworkspace::StopThread(void)
{
	if (0 == _ThreadHandle) {
		if (NULL != ARE::fpLOG) {
			INT64 tNowLog = ARE::GetTimeInMilliseconds() ;
			fprintf(ARE::fpLOG, "\n%I64d BEWS_th : stop variable order computation; already stopped ...", tNowLog) ;
			fflush(ARE::fpLOG) ;
			}
		return 0 ;
		}
	InterlockedCompareExchange(&_StopAndExit, 1, 0) ;
	_tToStop = ARE::GetTimeInMilliseconds() ;
	if (NULL != ARE::fpLOG) {
		fprintf(ARE::fpLOG, "\n%I64d BEWS_th : stop variable order computation; stop signalled, will wait ...", _tToStop) ;
		fflush(ARE::fpLOG) ;
		}
	while (true) {
		SLEEP(50) ;
		if (0 == _ThreadHandle) 
			break ;
		INT64 tNow = ARE::GetTimeInMilliseconds() ;
		INT64 dt = tNow - _tToStop ;
		if (dt > 10000) {
			// we asked the thread to stop and waited for it to stop, but it won't stop, so kill the thread.
#if defined WINDOWS || _WINDOWS
			TerminateThread((HANDLE) _ThreadHandle, 0) ;
			CloseHandle((HANDLE) _ThreadHandle) ;
			_ThreadHandle = 0 ;
#else
			// TODO : handle linux
#endif
			if (NULL != ARE::fpLOG) {
				INT64 tNowLog = ARE::GetTimeInMilliseconds() ;
				fprintf(ARE::fpLOG, "\n%I64d BEWS_th : stop variable order computation, hard kill ...", tNowLog) ;
				fflush(ARE::fpLOG) ;
				}
			break ;
			}
		}
	return 0 ;
}


int BucketElimination::BEworkspace::Initialize(ARE::ARP & Problem, const int *VarOrderingAsVarList, int DeleteUsedTables)
{
	int i ;

	Destroy() ;

	if (NULL == VarOrderingAsVarList) {
		VarOrderingAsVarList = Problem.VarOrdering_VarList() ;
		if (NULL == VarOrderingAsVarList) 
			return 1 ;
		}

	if (! _IsValid || 0 != ARE::Workspace::Initialize(Problem)) 
		{ _IsValid = false ; return 1 ; }
	if (NULL == _Problem) 
		{ _IsValid = false ; return 2 ; }

	_FnCombinationType = Problem.FnCombinationType() ;
	_VarEliminationType = Problem.VarEliminationType() ;

	if (FN_COBINATION_TYPE_PROD == _FnCombinationType) 
		_AnswerFactor = 1.0 ;
	else if (FN_COBINATION_TYPE_SUM == _FnCombinationType) 
		_AnswerFactor = 0.0 ;
	else 
		_AnswerFactor = 0.0 ;

	if (0 == DeleteUsedTables) 
		SetDeleteUsedTables(false) ;
	else if (DeleteUsedTables > 0) 
		SetDeleteUsedTables(true) ;

	_nVars = _Problem->N() ;
	if (_nVars <= 0) 
		return 0 ;

	// allocate space; initialize

	_VarOrder = new int[_nVars] ;
	_VarPos = new int[_nVars] ;
	_Var2BucketMapping = new BucketElimination::Bucket*[_nVars] ;
	_BucketOrderToCompute = new long[_nVars] ;
	if (NULL == _VarOrder || NULL == _VarPos || NULL == _Var2BucketMapping || NULL == _BucketOrderToCompute) {
		Destroy() ;
		_IsValid = false ;
		return 1 ;
		}
	for (i = 0 ; i < _nVars ; i++) 
		_VarOrder[i] = _VarPos[i] = -1 ;
	for (i = 0 ; i < _nVars ; i++) {
		_VarOrder[i] = VarOrderingAsVarList[i] ;
		_VarPos[_VarOrder[i]] = i ;
		_Var2BucketMapping[i] = NULL ;
		}
	for (i = 0 ; i < _nVars ; i++) {
		if (_VarOrder[i] < 0 || _VarOrder[i] >= _nVars) {
			Destroy() ;
			_IsValid = false ;
			return 1 ;
			}
		if (_VarPos[i] < 0 || _VarPos[i] >= _nVars) {
			Destroy() ;
			_IsValid = false ;
			return 1 ;
			}
		}

	// create buckets
	if (0 != CreateBuckets(false)) {
		_IsValid = false ;
		return 1 ;
		}

	// create computation order
	CreateComputationOrder(1) ;

#ifdef _DEBUG
	// verify integrity of functions
	for (i = 0 ; i < _Problem->nFunctions() ; i++) {
		ARE::Function *f = _Problem->getFunction(i) ;
		if (0 != f->CheckIntegrity()) 
			{ _IsValid = false ; return 1 ; }
		}
	for (i = 0 ; i < _nBuckets ; i++) {
		BucketElimination::Bucket *b = _Buckets[i] ;
		ARE::Function & f = b->OutputFunction() ;
		if (0 != f.CheckIntegrity()) 
			{ _IsValid = false ; return 1 ; }
		}
#endif // _DEBUG

	ComputeTotalOriginalFunctionSizeAndSpace() ;
	ComputeMaxNumChildren() ;
	ComputeMaxNumVarsInBucket() ;
	ComputeNBucketsWithSingleChild() ;
	ComputeMaxBucketFunctionWidth() ;
	ComputeTotalNewFunctionSizeAndSpace() ;
	ComputeTotalNewFunctionComputationComplexity() ;
	SimulateComputationAndComputeMinSpace(false) ;

	_TotalNewFunctionSizeComputed = 0 ;

	return 0 ;
}


int BucketElimination::BEworkspace::CreateBuckets(bool CreateSuperBuckets)
{
	if (_nVars <= 0) 
		return 0 ;

	int i, j, n, ret = 1 ;

	// ***************************************************************************************************
	// create initial bucket tree structure
	// ***************************************************************************************************

	for (i = 0 ; i < _nVars ; i++) {
		BucketElimination::Bucket *b = _Buckets[i] = new BucketElimination::Bucket(*this, i, _VarOrder[i]) ;
		if (NULL == b) {
			Destroy() ;
			return 1 ;
			}
		if (1 != b->nVars()) {
			Destroy() ;
			return 1 ;
			}
		_Var2BucketMapping[_VarOrder[i]] = b ;
		}
	_nBuckets = _nVars ;

	// ***************************************************************************************************
	// create initial function assignment
	// ***************************************************************************************************

	int NF = _Problem->nFunctions() ;
	if (0 == NF) 
		return 0 ;
	for (i = 0 ; i < NF ; i++) {
		ARE::Function *f = _Problem->getFunction(i) ;
		f->SetBEBucket(NULL) ;
		f->SetOriginatingBucket(NULL) ;
		}

	ARE::Function **fl = new ARE::Function*[2*NF] ;
	if (NULL == fl) 
		return 1 ;
	// for debugging purposes, we want to check which functions get assigned to a bucket
	ARE::Function **fl_assigned = fl + NF ;

	// process buckets, from last to first
	int nfPlaced = 0 ;
	int nfQirrelevant  = 0 ;
	for (i = 0 ; i < NF ; i++) {
		ARE::Function *f = _Problem->getFunction(i) ;
		if (f->IsQueryIrrelevant()) 
			++nfQirrelevant ;
		fl_assigned[i] = f ;
		}
	for (i = _nBuckets - 1 ; i >= 0 ; i--) {
		BucketElimination::Bucket *b = _Buckets[i] ;
		int v = _VarOrder[i] ;
		// find functions whose highest indexed variable is b->Var(0)
		n = 0 ;
		for (j = 0 ; j < _Problem->nAdjFunctions(v) ; j++) {
			ARE::Function *f = _Problem->AdjFunction(v, j) ;
			if (NULL == f) 
				continue ;
			if (f->IsQueryIrrelevant()) 
				continue ;
			if (NULL != f->BEBucket()) 
				// this means f was placed in a later bucket; i.e. v is not the highest-ordered variable in f.
				continue ;
			fl[n++] = f ;
			}
		nfPlaced += n ;
		b->SetOriginalFunctions(n, fl) ;
		// mark off the assigned functions
		for (j = 0 ; j < n ; j++) {
			ARE::Function *f = fl[j] ;
			int idxFN = f->IDX() ;
			if (idxFN >= 0 && idxFN < NF) {
				if (NULL == fl_assigned[idxFN]) {
					int error_FN_assigned_more_than_once = 1 ;
					}
				else 
					fl_assigned[idxFN] = NULL ;
				}
			}
		}

	// collect all const-functions; add as factor to the workspace
	{
	_nConstValueFunctions = 0 ;
	for (i = 0 ; i < NF ; i++) {
		ARE::Function *f = _Problem->getFunction(i) ;
		if (0 == f->N()) {
			AddAnswerFactor(f->ConstValue()) ;
			++_nConstValueFunctions ;
			}
		}

// DEBUGGGG
if (NULL != ARE::fpLOG) {
fprintf(ARE::fpLOG, "\nBEws::CreateBuckets _FnCombinationType=%d _VarElimType=%d AnswerFactor=%g", (int) _FnCombinationType, (int) _VarEliminationType, (double) _AnswerFactor) ;
fflush(ARE::fpLOG) ;
}
	}

	// test all functions are processed
	if ((nfPlaced + nfQirrelevant + _nConstValueFunctions) != NF) {
		for (i = 0 ; i < NF ; i++) {
			ARE::Function *f = _Problem->getFunction(i) ;
			int break_point_here = 1 ;
			}
		goto failed ;
		}

	// ***************************************************************************************************
	// create a bucket-function for each bucket, from last to first (bottom up on the bucket tree)
	// ***************************************************************************************************

	BucketElimination::Bucket *b ;
	for (i = _nBuckets - 1 ; i >= 0 ; i--) {
		b = _Buckets[i] ;
		if (b->ComputeSignature()) 
			goto failed ;
		// generate output fn signature; this will also compute width/signature of this bucket.
		if (0 != b->ComputeOutputFunctionWithScopeWithoutTable()) 
			goto failed ;
		}

	// eliminate buckets that have no functions in them
	_nVarsWithoutBucket = 0 ;
	for (i = _nBuckets - 1 ; i >= 0 ; i--) {
		BucketElimination::Bucket *b = _Buckets[i] ;
		if (0 == b->nOriginalFunctions() + b->nChildBucketFunctions()) {
			for (j = 0 ; j < b->nVars() ; j++) {
				_Var2BucketMapping[b->Var(j)] = NULL ;
				++_nVarsWithoutBucket ;
				}
			--_nBuckets ;
			for (j = b->IDX() ; j < _nBuckets ; j++) {
				_Buckets[j] = _Buckets[j+1] ;
				_Buckets[j]->SetIDX(j) ;
				}
			_Buckets[_nBuckets] = NULL ;
			delete b ;
			}
		}

#ifdef _DEBUG
	if (0 != CheckBucketTreeIntegrity()) 
		goto failed ;
#endif // _DEBUG

	_nBuckets_initial = _nBuckets ;
	_nBucketsWithSingleChild_initial = ComputeNBucketsWithSingleChild() ;

	// ***************************************************************************************************
	// traverse the bucket-tree bottom up; at each step, check if buckets can be merged.
	// in particular, we work with candidate buckets, where a candidate is something 
	// that is either a leaf (in the bucket tree) or has more than 1 child.
	// if a parent of a candidate has no other children (than the candidate), 
	// then the parent can be merged into the candidate.
	// ***************************************************************************************************

	n = 0 ;
	BucketElimination::Bucket *buckets[MAX_NUM_BUCKETS] ;
	if (CreateSuperBuckets) {
		for (i = 0 ; i < _nBuckets ; i++) {
			BucketElimination::Bucket *b = _Buckets[i] ;
			if (1 == b->nChildBucketFunctions()) 
				continue ;
			BucketElimination::Bucket *B = b->ParentBucket() ;
			if (NULL == B) 
				continue ;
			if (1 != B->nChildBucketFunctions()) 
				continue ;
			buckets[n++] = b ;
			}
		}

//int nMergers = 0 ;

	while (n > 0) {
		BucketElimination::Bucket *b = buckets[--n] ;
		BucketElimination::Bucket *B = b->ParentBucket() ;
		if (NULL == B) 
			continue ;
		// if B has no other children, other than b, B can be eliminated, by merging its original functions into b.
// TODO : check that b does not get too big (in terms of #Vars and PROD-of-domain-sizes-of-Vars).
		if (1 != B->nChildBucketFunctions()) 
			// B has other children
			continue ;
//fprintf(ARE::fpLOG, "\nSUPERBUCKETS adding B=%d to b=%d", B->IDX(), b->IDX()) ;
//fflush(ARE::fpLOG) ;
		// check that adding B to b does not increase complexity too much
		__int64 cmplxty_b = b->ComputeProcessingComplexity() ;
		__int64 cmplxty = cmplxty_b ;
		for (i = 0 ; i < B->Width() ; i++) {
			int vB = (B->Signature())[i] ;
			for (j = 0 ; j < b->Width() ; j++) {
				int vb = (b->Signature())[j] ;
				if (vB == vb) break ;
				}
			if (j < b->Width()) continue ;
			cmplxty *= _Problem->K(vB) ;
			}
		if (cmplxty >= (1 << 20)) {
			__int64 cmplxty_B = B->ComputeProcessingComplexity() ;
			if (cmplxty > (cmplxty_B + cmplxty_b)) {
				// merging would increase complexity too much, above a threshold;
				// don't merge; check B instead.
				if (NULL != ARE::fpLOG) 
					fprintf(ARE::fpLOG, "\nCANNOT merge bucket B=%d into b=%d, complexities would grow from %I64d+%I64d to %I64d", 
					B->IDX(), b->IDX(), cmplxty_B, cmplxty_b, cmplxty) ;
				buckets[n++] = B ;
				continue ;
				}
			}
		// add variables of B to b
		for (i = 0 ; i < B->nVars() ; i++) {
			j = B->Var(i) ;
			if (0 != b->AddVar(j)) 
				goto failed ;
			_Var2BucketMapping[j] = b ;
			}
		// remove b's bucketfunction from B
		B->RemoveChildBucketFunction(b->OutputFunction()) ;
		b->SetParentBucket(NULL) ;
		// add original functions of B to b
		if (b->AddOriginalFunctions(B->nOriginalFunctions(), B->OriginalFunctionsArray())) 
			goto failed ;
		if (b->ComputeSignature()) 
			goto failed ;
		// recompute bucketfunction of b
		if (0 != b->ComputeOutputFunctionWithScopeWithoutTable()) 
			goto failed ;
		// run a check
		if (b->Width() != b->OutputFunction().N() + b->nVars()) 
			{ int error = 1 ; }
		// detach B from its parent
		BucketElimination::Bucket *pB = B->ParentBucket() ;
		if (NULL != pB) 
			pB->RemoveChildBucketFunction(B->OutputFunction()) ;
		B->SetParentBucket(NULL) ;
		// shrink the _Buckets array
		j = _nBuckets - 1 ;
// DEBUGGG
//int BIDX = B->IDX() ;
		for (i = B->IDX() ; i < j ; i++) {
			_Buckets[i] = _Buckets[i+1] ;
			_Buckets[i]->SetIDX(i) ;
			}
		_nBuckets = j ;
		_Buckets[j] = NULL ;
		// delete B
		delete B ;
		// b needs to be rechecked
		buckets[n++] = b ;
// merger; DEBUGGGG
//fprintf(ARE::fpLOG, "\nELIM B=%d", BIDX) ;
//if (++nMergers >= 31) 
//	break ;
		}

	// ***************************************************************************************************
	// check the bucket tree
	// ***************************************************************************************************

#ifdef _DEBUG
	if (0 != CheckBucketTreeIntegrity()) 
		goto failed ;
#endif // _DEBUG

	// ***************************************************************************************************
	// set root ptr; set distance2root
	// ***************************************************************************************************

	for (i = 0 ; i < _nBuckets ; i++) {
		BucketElimination::Bucket *b = _Buckets[i] ;
		BucketElimination::Bucket *B = b->ParentBucket() ;
		if (NULL == B) {
			b->SetDistanceToRoot(0) ;
			b->SetRootBucket(b) ;
			}
		else {
			b->SetDistanceToRoot(B->DistanceToRoot() + 1) ;
			b->SetRootBucket(B->RootBucket()) ;
			}
		b->SetHeight(-1) ;
		}

	// ***************************************************************************************************
	// set height
	// ***************************************************************************************************

	for (i = 0 ; i < _nBuckets ; i++) {
		BucketElimination::Bucket *b = _Buckets[i] ;
		b->SetMaxDescendantNumVars(-1) ;
		b->SetMaxDescendantComputationNewFunctionSize(-1) ;
		}
	for (i = _nBuckets - 1 ; i >= 0 ; i--) {
		BucketElimination::Bucket *b = _Buckets[i] ;
		BucketElimination::Bucket *B = b->ParentBucket() ;
		if (b->Height() < 0) 
			b->SetHeight(0) ;
		if (NULL != B) {
			int h = b->Height() + 1 ;
			if (B->Height() < h) 
				B->SetHeight(h) ;
			int mdnv = b->MaxDescendantNumVarsEx() ;
			if (B->MaxDescendantNumVars() < mdnv) 
				B->SetMaxDescendantNumVars(mdnv) ;
			INT64 mdnfns = b->MaxDescendantComputationNewFunctionSizeEx() ;
			if (B->MaxDescendantComputationNewFunctionSize() < mdnfns) 
				B->SetMaxDescendantComputationNewFunctionSize(mdnfns) ;
			}
		}

	_nBuckets_final = _nBuckets ;
	_nBucketsWithSingleChild_final = ComputeNBucketsWithSingleChild() ;

	ret = 0 ;

failed :
	if (NULL != fl) 
		delete [] fl ;
	return ret ;
}


int BucketElimination::BEworkspace::CheckBucketTreeIntegrity(void)
{
	int i, j, k ;

	for (i = 0 ; i < _nBuckets ; i++) {
		BucketElimination::Bucket *b = _Buckets[i] ;
		BucketElimination::Bucket *B = b->ParentBucket() ;
		// check indeces are within range
		if (b->IDX() < 0 || b->IDX() >= _nBuckets) 
			goto failed ;
		if (NULL != B ? B->IDX() < 0 || B->IDX() >= _nBuckets : false) 
			goto failed ;
		// check width is computed
		if (b->Width() < 0) 
			goto failed ;
		// check output function
		ARE::Function & f = b->OutputFunction() ;
		if (f.BEBucket() != B) 
			goto failed ;
		if (f.N() > 0 && NULL == B) 
			goto failed ;
		if (0 == f.N() && NULL != B) 
			goto failed ;
		// check bVar, Signature, f->args are unique
		for (j = 0 ; j < b->nVars() ; j++) {
			int u = b->Var(j) ;
			for (k = j+1 ; k < b->nVars() ; k++) {
				if (u == b->Var(k)) 
					goto failed ;
				}
			}
		const int *bSig = b->Signature() ;
		for (j = 0 ; j < b->Width() ; j++) {
			int u = bSig[j] ;
			for (k = j+1 ; k < b->Width() ; k++) {
				if (u == bSig[k]) 
					goto failed ;
				}
			}
		const int *fArgs = f.Arguments() ;
		for (j = 0 ; j < f.N() ; j++) {
			int u = fArgs[j] ;
			for (k = j+1 ; k < f.N() ; k++) {
				if (u == fArgs[k]) 
					goto failed ;
				}
			}
		// check f->args and bVars are disjoint
		for (j = 0 ; j < f.N() ; j++) {
			int u = fArgs[j] ;
			for (k = 0 ; k < b->nVars() ; k++) {
				if (u == b->Var(k)) 
					goto failed ;
				}
			}
		for (j = 0 ; j < b->nVars() ; j++) {
			int u = b->Var(j) ;
			for (k = 0 ; k < f.N() ; k++) {
				if (u == fArgs[k]) 
					goto failed ;
				}
			}
		// check all f->args are in signature
		for (j = 0 ; j < f.N() ; j++) {
			int u = fArgs[j] ;
			for (k = 0 ; k < b->Width() ; k++) {
				if (u == bSig[k]) 
					break ;
				}
			if (k >= b->Width()) 
				goto failed ;
			}
		// bVar + nFN + width much match
		if (0 == b->nOriginalFunctions() + b->nChildBucketFunctions()) {
			if (0 != f.N()) 
				goto failed ;
			if (0 != b->Width()) 
				goto failed ;
			}
		else if (b->Width() != f.N() + b->nVars()) 
			goto failed ;
		// check children
		for (j = 0 ; j < b->nChildBucketFunctions() ; j++) {
			ARE::Function *cf = b->ChildBucketFunction(j) ;
			if (NULL == cf) 
				goto failed ;
			if (b != cf->BEBucket()) 
				goto failed ;
			if (NULL == cf->OriginatingBucket()) 
				goto failed ;
			if (cf->OriginatingBucket()->ParentBucket() != b) 
				goto failed ;
			}
		}

	return 0 ;
failed :
	fprintf(ARE::fpLOG, "\nERROR : bucket tree integrity check") ;
	fflush(ARE::fpLOG) ;
	return 1 ;
}


int BucketElimination::BEworkspace::CreateComputationOrder(int algorithm)
{
	int i ;

	if (NULL == _BucketOrderToCompute || _nBuckets < 1) 
		return 1 ;

	if (1 == algorithm) {
		BucketElimination::Bucket *BLhead = NULL, *BLtail = NULL ;
		int nAdded = 0 ;
		// collect all top-level buckets (i.e. buckets that have no parent).
		for (i = 0 ; i < _nBuckets; i++) {
			BucketElimination::Bucket *b = _Buckets[i] ;
			if (0 == b->DistanceToRoot()) {
				if (NULL == BLhead) 
					{ BLhead = BLtail = b ; }
				else 
					{ BLtail->NextInOrderComputationGenList() = b ; BLtail = b ; }
				b->NextInOrderComputationGenList() = NULL ;
				++nAdded ;
				}
			}
		// go through the list; for each bucket, compute all its children before buckets lates
		BucketElimination::Bucket *B = BLhead ;
		long bucket_fn_num_vars[256] ;
		BucketElimination::Bucket *child_buckets[256] ;
		long left[32], right[32] ;
		int n = 0 ;
		while (NULL != B) {
			_BucketOrderToCompute[n++] = B->IDX() ;
			// sort children of b
			long nChildren = 0 ;
			for (int j = 0 ; j < B->nChildBucketFunctions() ; j++) {
				ARE::Function *f = B->ChildBucketFunction(j) ;
				child_buckets[nChildren] = f->OriginatingBucket() ;
//				bucket_fn_num_vars[nChildren] = -f->N() ; // order from largest to smallest
//				bucket_fn_num_vars[nChildren] = child_buckets[nChildren]->MaxDescendantNumVarsEx() ; // order from smallest to largest
				bucket_fn_num_vars[nChildren] = child_buckets[nChildren]->MaxDescendantComputationNewFunctionSizeEx() ; // order from smallest to largest
				++nChildren ;
				}
			if (nChildren > 0) {
				QuickSortLong_i64(bucket_fn_num_vars, nChildren, (INT64 *) child_buckets, left, right) ;
				// add all children of B right after B
				BucketElimination::Bucket *Bnext = B->NextInOrderComputationGenList(), *b = B ;
				for (int j = 0 ; j < nChildren ; j++) {
					Bucket *bChild = child_buckets[j] ;
					b->NextInOrderComputationGenList() = bChild ;
					b = bChild ;
					++nAdded ;
					}
				b->NextInOrderComputationGenList() = Bnext ;
				}
			B = B->NextInOrderComputationGenList() ;
			}
		if (n != _nBuckets) {
			int bug_here = 1 ;
			}
		}
	else {
		// sort by height descending
		long *keys = new long[_nBuckets] ;
		if (NULL == keys)
			return 1 ;
		for (i = 0; i < _nBuckets; i++) {
			BucketElimination::Bucket *b = _Buckets[i] ;
			keys[i] = b->Height() ;
			_BucketOrderToCompute[i] = b->IDX() ;
			}
		long left[32], right[32] ;
		QuickSortLong_Descending(keys, _nBuckets, _BucketOrderToCompute, left, right) ;
		delete [] keys ;
		}

	return 0 ;
}


__int64 BucketElimination::BEworkspace::ComputeTotalOriginalFunctionSizeAndSpace(void)
{
	_TotalOriginalFunctionSize = _TotalOriginalFunctionSpace = 0 ;
	if (NULL == _Problem) 
		return 0 ;
	int NF = _Problem->nFunctions() ;
	for (int i = 0 ; i < NF ; i++) {
		ARE::Function *f = _Problem->getFunction(i) ;
		_TotalOriginalFunctionSize += f->ComputeTableSize() ;
		_TotalOriginalFunctionSpace += f->ComputeTableSpace() ;
		}
	return _TotalOriginalFunctionSize ;
}


__int64 BucketElimination::BEworkspace::ComputeTotalNewFunctionSizeAndSpace(void)
{
	_TotalNewFunctionSize = _TotalNewFunctionSpace = 0 ;
	for (int i = 0 ; i < _nBuckets ; i++) {
		BucketElimination::Bucket *b = _Buckets[i] ;
		ARE::Function & f = b->OutputFunction() ;
		_TotalNewFunctionSize += f.ComputeTableSize() ;
		_TotalNewFunctionSpace += f.ComputeTableSpace() ;
		}
	return _TotalNewFunctionSize ;
}


__int64 BucketElimination::BEworkspace::ComputeTotalNewFunctionComputationComplexity(void)
{
	_TotalNewFunctionComputationComplexity = 0 ;
	for (int i = 0 ; i < _nBuckets ; i++) {
		BucketElimination::Bucket *b = _Buckets[i] ;
		int w = b->Width() ;
		__int64 n = 1 ;
		const int *vars = b->Signature() ;
		for (int j = 0 ; j < w ; j++) 
			n *= _Problem->K(vars[j]) ;
		__int64 N = n * (b->nChildBucketFunctions() + b->nOriginalFunctions()) ;
		_TotalNewFunctionComputationComplexity += N ;
		}
	return _TotalNewFunctionComputationComplexity ;
}


__int64 BucketElimination::BEworkspace::SimulateComputationAndComputeMinSpace(bool IgnoreInMemoryTables)
{
	int i ;

	// compute table sizes
	for (i = _nBuckets - 1 ; i >= 0 ; i--) {
		BucketElimination::Bucket *b = _Buckets[i] ;
		// add the space of the output function of this bucket
		ARE::Function & f = b->OutputFunction() ;
		f.ComputeTableSize() ;
		for (int j = 0 ; j < b->nChildBucketFunctions() ; j++) {
			ARE::Function *cf = b->ChildBucketFunction(j) ;
			if (NULL != cf) 
				cf->ComputeTableSize() ;
			}
		}

	// compute maximum space for any bucket
	__int64 size = 0 ;
	_MaxSimultaneousNewFunctionSize = 0 ;
	for (i = _nBuckets - 1 ; i >= 0 ; i--) {
		int idx = _BucketOrderToCompute[i] ;
		BucketElimination::Bucket *b = _Buckets[idx] ;
		// add the space of the output function of this bucket
		ARE::Function & f = b->OutputFunction() ;
		__int64 s = f.TableSize() ; if (s < 0) s = 0 ;
		size += s ;
		if (size > _MaxSimultaneousNewFunctionSize) 
			_MaxSimultaneousNewFunctionSize = size ;
		if (_DeleteUsedTables) {
			for (int j = 0 ; j < b->nChildBucketFunctions() ; j++) {
				ARE::Function *cf = b->ChildBucketFunction(j) ;
				if (NULL == cf) continue ;
				if (cf->TableSize() <= 0) continue ;
				size -= cf->TableSize() ;
				}
			}
		}

	_MaxSimultaneousNewFunctionSpace = sizeof(ARE_Function_TableType)*_MaxSimultaneousNewFunctionSize ;
	_MaxTotalFunctionSize = _MaxSimultaneousNewFunctionSize + _TotalOriginalFunctionSize ;
	_MaxTotalFunctionSpace = sizeof(ARE_Function_TableType)*_MaxTotalFunctionSize ;

	return _MaxSimultaneousNewFunctionSize ;
}


int BucketElimination::BEworkspace::ComputeMaxBucketFunctionWidth(void)
{
	_MaxBucketFunctionWidth = 0 ;
	for (int i = 0 ; i < _nBuckets ; i++) {
		BucketElimination::Bucket *b = _Buckets[i] ;
		ARE::Function & f = b->OutputFunction() ;
		if (f.N() > _MaxBucketFunctionWidth) _MaxBucketFunctionWidth = f.N() ;
		}
	return _MaxBucketFunctionWidth ;
}


int BucketElimination::BEworkspace::ComputeMaxNumVarsInBucket(void)
{
	_MaxNumVarsInBucket = 0 ;
	for (int i = 0 ; i < _nBuckets ; i++) {
		BucketElimination::Bucket *b = _Buckets[i] ;
		if (b->Width() > _MaxNumVarsInBucket) 
			_MaxNumVarsInBucket = b->Width() ;
		}
	return _MaxNumVarsInBucket ;
}


int BucketElimination::BEworkspace::ComputeNBucketsWithSingleChild(void)
{
	_nBucketsWithSingleChild = 0 ;
	for (int i = 0 ; i < _nBuckets ; i++) {
		BucketElimination::Bucket *b = _Buckets[i] ;
		if (1 == b->nChildBucketFunctions()) 
			_nBucketsWithSingleChild++ ;
		}
	return _nBucketsWithSingleChild ;
}


int BucketElimination::BEworkspace::ComputeMaxNumChildren(void)
{
	_MaxNumChildren = 0 ;
	for (int i = 0 ; i < _nBuckets ; i++) {
		BucketElimination::Bucket *b = _Buckets[i] ;
		if (b->nChildBucketFunctions() > _MaxNumChildren) 
			_MaxNumChildren = b->nChildBucketFunctions() ;
		}
	return _MaxNumChildren ;
}


int BucketElimination::BEworkspace::RunSimple(void)
{
	int i ;

	// mark all bucket functions as in-memory
	for (i = _nBuckets - 1 ; i >= 0 ; i--) {
		BucketElimination::Bucket *b = _Buckets[i] ;
		ARE::Function & f = b->OutputFunction() ;
		b->AllocateOutputFunctionBlockComputationResult(1000000, 0) ;
		f.ComputeTableSize() ;
		f.SetnTableBlocks(0) ;
		f.AllocateInMemoryAsSingleTableBlock() ;
		}

	// compute all bucket functions
	ARE::Function *MissingFunction ;
	__int64 MissingBlockIDX ;
	for (i = _nBuckets - 1 ; i >= 0 ; i--) {
		BucketElimination::Bucket *b = _Buckets[i] ;
		ARE::Function & f = b->OutputFunction() ;
		ARE::FunctionTableBlock *ftb = f.Table() ;
		if (NULL == ftb) 
			// this is error
			continue ;
		if (FN_COBINATION_TYPE_PROD == _FnCombinationType && VAR_ELIMINATION_TYPE_SUM == _VarEliminationType) 
			ftb->ComputeDataBEEM_ProdSum(MissingFunction, MissingBlockIDX) ;
		else if (FN_COBINATION_TYPE_PROD == _FnCombinationType && VAR_ELIMINATION_TYPE_MAX == _VarEliminationType) 
			ftb->ComputeDataBEEM_ProdMax(MissingFunction, MissingBlockIDX) ;
		else if (FN_COBINATION_TYPE_SUM == _FnCombinationType && VAR_ELIMINATION_TYPE_MAX == _VarEliminationType) 
			ftb->ComputeDataBEEM_SumMax(MissingFunction, MissingBlockIDX) ;
		else 
			// cannot combine/marginalize
			break ;
		}

	return 0 ;
}


int BucketElimination::BEworkspace::GenerateRandomBayesianNetworkStructure(int N, int K, int P, int C, int ProblemCharacteristic)
{
return 1 ;
/*
	if (NULL == _Problem) 
		return 1 ;
	_Problem->Destroy() ;

	if (0 != _Problem->GenerateRandomUniformBayesianNetworkStructure(N, K, P, C, ProblemCharacteristic)) 
		{ return 1 ; }

	if (0 != _Problem->ComputeGraph()) 
		{ return 1 ; }
	int nComponents = _Problem->nConnectedComponents() ;
	if (nComponents > 1) 
		{ return 1 ; }

	if (0 != _Problem->ComputeMinDegreeOrdering()) 
		{ return 1 ; }

	int width = _Problem->MinDegreeOrdering_InducedWidth() ;
//	if (width < MinWidth || width > MaxWidth) 
//		{ return 1 ; }

	if (0 != _Problem->TestVariableOrdering(_Problem->MinDegreeOrdering_VarList(), _Problem->MinDegreeOrdering_VarPos())) 
		{ return 1 ; }

	if (0 != CreateBuckets(_Problem->N(), _Problem->MinDegreeOrdering_VarList(), true)) 
		{ return 1 ; }

	__int64 spaceO = _Problem->ComputeFunctionSpace() ;
	__int64 spaceN = ComputeNewFunctionSpace() ;
	__int64 space = spaceO + spaceN ;
//	if (space < MinMemory || space > MaxMemory) 
//		{ return 1 ; }

	return 0 ;
*/
}


int BucketElimination::GenerateRandomBayesianNetworksWithGivenComplexity(int nProblems, int N, int K, int P, int C, int ProblemCharacteristic, __int64 MinSpace, __int64 MaxSpace)
{
return 1 ;
/*
	ARE::ARP p("test") ;
	// create BEEM workspace; this includes BE workspace.
	BEworkspace ws(NULL) ;
	ws.Initialize(p, NULL) ;

	time_t ttNow, ttBeginning ;
	time(&ttBeginning) ;
	time_t dMax = 3600 ;

	char s[256] ;
	int i ;
	for (i = 0 ; i < nProblems ; ) {
		time(&ttNow) ;
		time_t d = ttNow - ttBeginning ;
		if (d >= dMax) 
			break ;

		ws.DestroyBuckets() ;
		ws.GenerateRandomBayesianNetworkStructure(N, K, P, C, ProblemCharacteristic) ;

		if (p.nConnectedComponents() > 1) 
			continue ;

		int width = p.MinDegreeOrdering_InducedWidth() ;
		__int64 spaceO = p.ComputeFunctionSpace() ;
		__int64 spaceN = ws.ComputeNewFunctionSpace() ;
		__int64 space = spaceO + spaceN ;

		if (space < MinSpace || space > MaxSpace) 
			continue ;

		// generate nice name
		sprintf(s, "random-test-problem-%d-Space=%I64d", (int) ++i, space) ;
		p.SetName(s) ;

		// fill in original functions
		p.FillInFunctionTables() ;
		if (p.CheckFunctions()) 
			{ int error = 1 ; }

		// save UAI08 format
		p.SaveUAI08("c:\\UCI\\problems") ;
		}

	return i ;
*/
}

