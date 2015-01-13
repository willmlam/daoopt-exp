// CVO.cpp : Defines the entry point for the console application.
//

#include <stdlib.h>
#include <stdio.h>
#include <string>
#include <time.h>

#if defined WINDOWS || _WINDOWS
#include <process.h>    /* _beginthread, _endthread */
#endif // WINDOWS

#include "Utils/AVLtreeSimple.hxx"
#include "Utils/RandomProblemGenerator.hxx"
#include <Utils/MiscUtils.hxx>
#include <CVO/VariableOrderComputation.hxx>

#ifdef LINUX
pthread_mutex_t nRunsSumMutex = PTHREAD_MUTEX_INITIALIZER ;
pthread_mutex_t stopSignalMutex = PTHREAD_MUTEX_INITIALIZER ;
#endif

static void GetCurrentDTsec(char *strDT, time_t & ttNow)
{
	// get current date/time string so that we can log when stuff is happening
	if (0 == ttNow) time(&ttNow) ;
	struct tm *pTime = localtime(&ttNow) ;
	char *sDT = asctime(pTime) ; // from manual : The string result produced by asctime contains exactly 26 characters and has the form Wed Jan 02 02:03:55 1980\n\0.
	memcpy(strDT, sDT, 26) ;
	strDT[24] = 0 ;
}

static void GetCurrentDTmsec(char *strDT, INT64 & tNow)
{
	// get current date/time string so that we can log when stuff is happening
	if (0 == tNow) tNow = ARE::GetTimeInMilliseconds() ;
	sprintf(strDT, "%I64d", tNow) ;
}

int ARE::VarElimOrderComp::CVOcontext::NoteVarOrderComputationCompletion(int w_IDX, ARE::Graph & G)
{
	++_nRunsCompleted ;
	if ((StateSpaceSize==_ObjCode && (G._TotalVarElimComplexity < _BestOrder->_Complexity || (fabs(G._TotalVarElimComplexity - _BestOrder->_Complexity) < 0.01 && G._VarElimOrderWidth < _BestOrder->_Width))) ||  
		(Width==_ObjCode && (G._VarElimOrderWidth < _BestOrder->_Width || (G._VarElimOrderWidth == _BestOrder->_Width && G._TotalVarElimComplexity < _BestOrder->_Complexity)))) {
		INT64 tNow = ARE::GetTimeInMilliseconds() ;
		if (NULL != _fpLOG) {
			fprintf(_fpLOG, "\n%I64d worker %2d found better solution : width=%d complexity=%g space(#elements)=%g", tNow, (int) w_IDX, (int) G._VarElimOrderWidth, (double) G._TotalVarElimComplexity, (double) G._TotalNewFunctionStorageAsNumOfElements) ;
			fflush(_fpLOG) ;
			}

		ARE::VarElimOrderComp::ResultSnapShot & result_record = _Improvements[_nImprovements++] ;
		result_record._dt = tNow - _tStart ;
		result_record._width = G._VarElimOrderWidth ;
		result_record._complexity = G._TotalVarElimComplexity ;

		_BestOrder->_Width = G._VarElimOrderWidth ;
		_BestOrder->_nFillEdges = G._nFillEdges ;
		_BestOrder->_MaxSingleVarElimComplexity = G._MaxVarElimComplexity ;
		_BestOrder->_Complexity = G._TotalVarElimComplexity ;
		_BestOrder->_TotalNewFunctionStorageAsNumOfElements = G._TotalNewFunctionStorageAsNumOfElements ;
		for (int i = 0 ; i < _Problem->N() ; i++) 
			_BestOrder->_VarListInElimOrder[i] = (G._VarElimOrder)[i] ;
		}

	if (G._VarElimOrderWidth >= 0 && G._VarElimOrderWidth < 1024) {
		_Width2CountMap[G._VarElimOrderWidth]++ ;
		if (_Width2MaxComplexityMap[G._VarElimOrderWidth] > G._TotalVarElimComplexity) 
			_Width2MaxComplexityMap[G._VarElimOrderWidth] = G._TotalVarElimComplexity ;
		if (_Width2MaxComplexityMap[G._VarElimOrderWidth] < G._TotalVarElimComplexity) 
			_Width2MaxComplexityMap[G._VarElimOrderWidth] = G._TotalVarElimComplexity ;
		}

	return 0 ;
}


#if defined WINDOWS || _WINDOWS
typedef unsigned int (*pWorkerThreadFn)(void *X) ;
static unsigned int __stdcall WorkerThreadFn(void *X) 
#elif defined (LINUX)
typedef void *(*pWorkerThreadFn)(void *X) ;
static void *WorkerThreadFn(void *X)
#endif 
{
	ARE::VarElimOrderComp::Worker *w = (ARE::VarElimOrderComp::Worker *) X;
	ARE::VarElimOrderComp::CVOcontext & CVOcontext = *(w->_CVOcontext) ;
	ARE::ARP & Problem = *(CVOcontext._Problem) ;
	ARE::VarElimOrderComp::Order & best_order = *(CVOcontext._BestOrder) ;

	// initialize best order to something bad; we don't want anything worse than that.
	int bestWidth = 255 ;
	double bestComplexity = DBL_MAX ;
	{
	ARE::utils::AutoLock lock(CVOcontext._MasterGraphMutex) ;
	bestWidth = best_order._Width ;
	bestComplexity = best_order._Complexity ;
	}

	INT64 tNow = 0 ;
	int res ;

	if (NULL != CVOcontext._fpLOG) {
		tNow = ARE::GetTimeInMilliseconds() ;
		fprintf(CVOcontext._fpLOG, "\n%I64d worker %2d start; objCode=%d alg=%d ...", tNow, (int) w->_IDX, (int) CVOcontext._ObjCode, (int) CVOcontext._AlgCode) ;
		fflush(CVOcontext._fpLOG) ;
		}

	try {
		*(w->_G) = CVOcontext._MasterGraph ;
		if (! w->_G->_IsValid) 
			goto done ;
		}
	catch (...) {
		if (NULL != CVOcontext._fpLOG) {
			tNow = ARE::GetTimeInMilliseconds() ;
			fprintf(CVOcontext._fpLOG, "\n%I64d worker %2d prep0 exception ...", tNow, (int) w->_IDX) ;
			fflush(CVOcontext._fpLOG) ;
			}
		}

	int nCompleteRunsTodo ; nCompleteRunsTodo = 3 ;
	while (! w->_ThreadStop) {
		tNow = 0 ;
		++(w->_nRunsDone) ;
#if defined WINDOWS || _WINDOWS
		long v = InterlockedIncrement(&(CVOcontext._nRunsDone)) ;
#else
		pthread_mutex_lock(&nRunsSumMutex) ;
		long v = ++CVOcontext._nRunsDone ;
		pthread_mutex_unlock(&nRunsSumMutex) ;
#endif
		if (0 == (v % CVOcontext._LogIncrement)) {
			if (NULL != CVOcontext._fpLOG) {
				tNow = ARE::GetTimeInMilliseconds() ;
				fprintf(CVOcontext._fpLOG, "\n%I64d worker %2d will do run %d ...", tNow, (int) w->_IDX, (int) v) ;
				fflush(CVOcontext._fpLOG) ;
				}
			}
// DEBUGGG
//printf("\nworker %d starting; nRunsSum=%d ...", (int) w->_IDX, (int) v) ;
		try {
			bool earlyTerminationOk = nCompleteRunsTodo-- > 0 ? false : true ;
			int widthLimit = bestWidth ;
			if (CVOcontext._FindPracticalVariableOrder && widthLimit > CVOcontext._PracticalOrderLimit_W) 
				widthLimit = CVOcontext._PracticalOrderLimit_W ;
			int spaceLimit = bestComplexity ;
			if (CVOcontext._FindPracticalVariableOrder && spaceLimit > CVOcontext._PracticalOrderLimit_C) 
				spaceLimit = CVOcontext._PracticalOrderLimit_C ;
// 2014-03-21 KK : even if MinFill algorithm is used, run Simple(), not Simple_wMinFillOnly(), because Simple() will compute complexity also, as so we can try to minimize complexity too.
//			if (ARE::VarElimOrderComp::MinFill == CVOcontext._AlgCode) 
//				res = w->_G->ComputeVariableEliminationOrder_Simple_wMinFillOnly(widthLimit, earlyTerminationOk && CVOcontext._EarlyTerminationOfBasic_W, false, 1, CVOcontext._nRandomPick, CVOcontext._eRandomPick, w->_TempAdjVarSpace, TempAdjVarSpaceSize) ;
//			else 
				res = w->_G->ComputeVariableEliminationOrder_Simple(CVOcontext._AlgCode, widthLimit, earlyTerminationOk && CVOcontext._EarlyTerminationOfBasic_W, spaceLimit, earlyTerminationOk && CVOcontext._EarlyTerminationOfBasic_C, false, 1, CVOcontext._nRandomPick, CVOcontext._eRandomPick, w->_TempAdjVarSpaceSizeExtraArrayN, w->_TempAdjVarSpaceSizeExtraArray) ;
			}
		catch (...) {
			if (NULL != CVOcontext._fpLOG) {
				tNow = ARE::GetTimeInMilliseconds() ;
				fprintf(CVOcontext._fpLOG, "\n%I64d worker %2d exception ...", tNow, (int) w->_IDX) ;
				fflush(CVOcontext._fpLOG) ;
				goto done ;
				}
			}
// DEBUGGG
//		printf("\nworker %d finished res=%d; width=%d complexity=%I64d", (int) w->_IDX, (int) res, (int) w->_G->_VarElimOrderWidth, (INT64) w->_G->_TotalVarElimComplexity) ;
//		if (NULL != CVOcontext._fpLOG) {
//			fprintf(CVOcontext._fpLOG, "\nworker %d finished res=%d; width=%d complexity=%I64d", (int) w->_IDX, (int) res, (int) w->_G->_VarElimOrderWidth, (INT64) w->_G->_TotalVarElimComplexity) ;
//			fflush(CVOcontext._fpLOG) ;
//			}
		try {
			ARE::utils::AutoLock lock(CVOcontext._MasterGraphMutex) ;
//GetCurrentDTmsec(strDT, tNow) ;
//printf("\n%s worker %2d found width=%d complexity=%I64d space(#elements)=%I64d res=%d", strDT, (int) w->_IDX, (int) w->_G->_VarElimOrderWidth, (INT64) w->_G->_TotalVarElimComplexity, (INT64) w->_G->_TotalNewFunctionStorageAsNumOfElements, (int) res) ;
			if (0 == res) {
				CVOcontext.NoteVarOrderComputationCompletion(w->_IDX, *(w->_G)) ;
				}
			else {
				// DEBUGGG
//				printf("\nThread %d res=%d", w->_IDX, res) ;
				}
			bestWidth = best_order._Width ;
			bestComplexity = best_order._Complexity ;
			if (CVOcontext._nRunsDone >= CVOcontext._nRunsToDoMax) {
// DEBUGGG
//				printf("\nworker %d out of runs (%d >= %d) ...", (int) w->_IDX, (int) nRunsSum, (int) CVOcontext._nRunsToDoMax) ;
				goto done ;
				}
			if (CVOcontext._tToStop > 0) {
				if (0 == tNow) tNow = ARE::GetTimeInMilliseconds() ;
				if (tNow >= CVOcontext._tToStop) {
// DEBUGGG
//					printf("\nworker %d out of time ...", (int) w->_IDX) ;
					goto done ;
					}
				}

			}
		catch (...) {
			if (NULL != CVOcontext._fpLOG) {
				tNow = ARE::GetTimeInMilliseconds() ;
				fprintf(CVOcontext._fpLOG, "\n%I64d worker %d summary exception ...", tNow, (int) w->_IDX) ;
				fflush(CVOcontext._fpLOG) ;
				}
			}

		try {
			*(w->_G) = CVOcontext._MasterGraph ;
			}
		catch (...) {
			if (NULL != CVOcontext._fpLOG) {
				tNow = ARE::GetTimeInMilliseconds() ;
				fprintf(CVOcontext._fpLOG, "\n%I64d worker %d prep exception ...", tNow, (int) w->_IDX) ;
				fflush(CVOcontext._fpLOG) ;
				}
			}
		}

done :
	if (NULL != CVOcontext._fpLOG) {
		tNow = ARE::GetTimeInMilliseconds() ;
		fprintf(CVOcontext._fpLOG, "\n%I64d worker %d done ...", tNow, (int) w->_IDX) ;
		fflush(CVOcontext._fpLOG) ;
		}
	w->_ThreadHandle = 0 ;
	w->_ThreadIsRunning = false ;
#if defined WINDOWS || _WINDOWS
	_endthreadex(0) ;
	return 0  ;
#else
	return NULL ;
#endif
}


#if defined WINDOWS || _WINDOWS
typedef unsigned int (*pCVOThreadFn)(void *X) ;
static unsigned int __stdcall CVOThreadFn(void *X) 
#elif defined (LINUX)
typedef void *(*pCVOThreadFn)(void *X) ;
static void *CVOThreadFn(void *X)
#endif 
{
	ARE::VarElimOrderComp::CVOcontext *context = (ARE::VarElimOrderComp::CVOcontext *)(X) ;
	context->Reset() ;
	int nWorkers = context->_nThreads ;
	int nRunsToDoMax = context->_nRunsToDoMax ;
	ARE::ARP & p = *(context->_Problem) ;
	ARE::VarElimOrderComp::Order & BestOrder = *(context->_BestOrder) ;
	int & ret = context->_ret, i ;
	ARE::Graph & OriginalGraph = context->_OriginalGraph ;
	ARE::Graph & MasterGraph = context->_MasterGraph ;
	ARE::VarElimOrderComp::Order & best_order = *(context->_BestOrder) ;

	ARE::VarElimOrderComp::Worker *Workers = NULL ;


	int nRunning;
	int nWrunning;
	long stop_signalled;

	char strDT[64] ;
	INT64 tNow = 0 ; // ARE::GetTimeInMilliseconds() ;
	GetCurrentDTmsec(strDT, tNow) ;
	if (NULL != context->_fpLOG) {
		fprintf(context->_fpLOG, "\n%s CVO control thread; start preprocessing ...", strDT) ;
		fflush(context->_fpLOG) ;
		}

	if (NULL == BestOrder._VarListInElimOrder) {
		if (0 != BestOrder.Initialize(p.N())) {
			ret = 1000 ;
			if (NULL != context->_fpLOG) {
				fprintf(context->_fpLOG, "\n%I64d CVO control thread; bestorder init failed ...", tNow) ;
				fflush(context->_fpLOG) ;
				}
			goto done ;
			}
		}

	// create problem graph
	OriginalGraph.Create(p) ;
	if (! OriginalGraph._IsValid) {
		ret = 1001 ;
		goto done ;
		}

	// do all the easy eliminations; this will give us a starting point for large-scale randomized searches later.
	tNow = ARE::GetTimeInMilliseconds() ;
	if (NULL != context->_fpLOG) {
		fprintf(context->_fpLOG, "\n%I64d CVO control thread; eliminate easy vars from graph ...", tNow) ;
		fflush(context->_fpLOG) ;
		}
	MasterGraph = OriginalGraph ;
// 2014-03-21 KK : even if MinFill algorithm is used, run Simple(), not Simple_wMinFillOnly(), because Simple() will compute complexity also, as so ew can try to minimize complexity too.
//	if (ARE::VarElimOrderComp::MinFill == context->_AlgCode) 
//		i = MasterGraph.ComputeVariableEliminationOrder_Simple_wMinFillOnly(INT_MAX, false, false, 1, 1, 0.0, context->_TempAdjVarSpace, TempAdjVarSpaceSize) ;
//	else 
		i = MasterGraph.ComputeVariableEliminationOrder_Simple(0, INT_MAX, false, DBL_MAX, false, true, 1, 1, 0.0, context->_TempAdjVarSpaceSizeExtraArrayN, context->_TempAdjVarSpaceSizeExtraArray) ;
	MasterGraph.ReAllocateEdges() ;
	tNow = ARE::GetTimeInMilliseconds() ;
	if (NULL != context->_fpLOG) {
		fprintf(context->_fpLOG, "\n%I64d CVO control thread; %d vars eliminated, %d remaining ...", tNow, (int) MasterGraph._OrderLength, (int) MasterGraph._nRemainingNodes) ;
		fflush(context->_fpLOG) ;
		}

	context->_tStart = tNow ;
	if (context->_TimeLimitInMilliSeconds > 0) 
		context->_tToStop = context->_tStart + context->_TimeLimitInMilliSeconds ;
	if (NULL != context->_fpLOG) {
		fprintf(context->_fpLOG, "\n%I64d CVO control thread; preprocessing done; ready to run : nRunsToDo=%d/%d, tLimitmsec=%I64d (tAutoStop=%I64d) ...", tNow, context->_nRunsToDoMin, context->_nRunsToDoMax, context->_TimeLimitInMilliSeconds, context->_tToStop) ;
		fflush(context->_fpLOG) ;
		}

	// if strictly best order (whatever the width/complexity) is required, execute one run here to get some real bound on width/complexity.
	// when we launch multi-threaded search for good order, knowing a decent bound will help the threads right away.
	// when we want a practical order, width/complexity limits are (should be) quite low.
	best_order._Width = p.N() ;
	best_order._Complexity = DBL_MAX ;
	if (! context->_FindPracticalVariableOrder) {
		ARE::Graph g ;
		g = MasterGraph ;
		int widthLimit = context->_FindPracticalVariableOrder ? context->_PracticalOrderLimit_W : INT_MAX ;
		int spaceLimit = context->_FindPracticalVariableOrder ? context->_PracticalOrderLimit_C : DBL_MAX ;
		i = g.ComputeVariableEliminationOrder_Simple(0, widthLimit, context->_FindPracticalVariableOrder ? true : false, spaceLimit, context->_FindPracticalVariableOrder ? true : false, false, 10, 1, 0.0, context->_TempAdjVarSpaceSizeExtraArrayN, context->_TempAdjVarSpaceSizeExtraArray) ;
		++context->_nRunsDone ;
		tNow = 0 ;
		GetCurrentDTmsec(strDT, tNow) ;
		if (0 == i) {
			if (NULL != context->_fpLOG) {
				fprintf(context->_fpLOG, "\n%s Initial computation width=%d; MaxSingleVarElimComplexity=%g, TotalVarElimComplexity=%g, TotalNewFunctionStorageAsNumOfElements=%g", strDT, (int) g._VarElimOrderWidth, (double) g._MaxVarElimComplexity, (double) g._TotalVarElimComplexity, (double) g._TotalNewFunctionStorageAsNumOfElements) ;
				fflush(context->_fpLOG) ;
				}
			context->NoteVarOrderComputationCompletion(-1, g) ;
			}
		else {
			if (NULL != context->_fpLOG) {
				fprintf(context->_fpLOG, "\n%s Initial computation failed; res=%d; will continue ...", strDT, i) ;
				fflush(context->_fpLOG) ;
				}
			}
		}

	if (context->_nRunsDone >= context->_nRunsToDoMax) 
		goto done ;

	if (NULL != context->_fpLOG) {
		tNow = ARE::GetTimeInMilliseconds() ;
		fprintf(context->_fpLOG, "\n%I64d CVO control thread; nWorkers=%d", tNow, nWorkers) ;
		fflush(context->_fpLOG) ;
		}

	if (nWorkers < 1) {
		if (NULL != context->_fpLOG) {
			fprintf(context->_fpLOG, "\nNo threads; will quit ...") ;
			fflush(context->_fpLOG) ;
			goto done ;
			}
		}
//goto done ;

	// create workers; don't let them run yet
	Workers = new ARE::VarElimOrderComp::Worker[nWorkers] ;
	if (NULL == Workers) {
		ret = 1002 ;
		if (NULL != context->_fpLOG) {
			tNow = ARE::GetTimeInMilliseconds() ;
			fprintf(context->_fpLOG, "\n%I64d CVO control thread; failed to create workers, will quit ...", tNow) ;
			fflush(context->_fpLOG) ;
			}
		goto done ;
		}
	for (i = 0 ; i < nWorkers ; i++) {
		Workers[i]._IDX = i ;
		Workers[i]._G = new ARE::Graph ;
		Workers[i]._CVOcontext = context ;
		if (NULL == Workers[i]._G) {
			ret = 1003 ;
			goto done ;
			}
		}

	context->_LogIncrement = nRunsToDoMax/20 ;
	if (context->_LogIncrement < 1) context->_LogIncrement = 1 ;
  nWrunning = 0 ;
	for (i = 0 ; i < nWorkers ; i++) {
		Workers[i]._ThreadIsRunning = true ;
		Workers[i]._ThreadStop = false ;
#if defined WINDOWS || _WINDOWS
		Workers[i]._ThreadHandle = _beginthreadex(NULL, 0, WorkerThreadFn, &(Workers[i]), 0, NULL) ;
#else
		pthread_create(&(Workers[i]._ThreadHandle), NULL, WorkerThreadFn, &(Workers[i])) ; // TODO third argument
#endif 
		if (0 == Workers[i]._ThreadHandle) {
			Workers[i]._ThreadIsRunning = false ;
			Workers[i]._ThreadStop = true ;
			if (NULL != context->_fpLOG) {
				tNow = ARE::GetTimeInMilliseconds() ;
				fprintf(context->_fpLOG, "\n%I64d CVO control thread; FAILED to create thread %d", tNow, i) ;
				fflush(context->_fpLOG) ;
				}
			}
		else 
			nWrunning++ ;
		}
	if (0 == nWrunning) {
		if (NULL != context->_fpLOG) {
			tNow = ARE::GetTimeInMilliseconds() ;
			fprintf(context->_fpLOG, "\n%I64d CVO control thread; no worker threads, will quite ...", tNow) ;
			fflush(context->_fpLOG) ;
			}
		ret = 1004 ;
		goto done ;
		}

	// loop here waiting for stop signal of until all workers threads have exited
	nRunning = -1 ;
	stop_signalled = 0 ;
	while (true) {
		SLEEP(50) ;
#if defined WINDOWS || _WINDOWS
		stop_signalled = InterlockedCompareExchange(&(context->_StopAndExit), 1, 1) ;
#else
    pthread_mutex_lock(&stopSignalMutex);
    stop_signalled = context->_StopAndExit;
    pthread_mutex_unlock(&stopSignalMutex);
#endif
		if (0 != stop_signalled) {
			if (NULL != context->_fpLOG) {
				tNow = ARE::GetTimeInMilliseconds() ;
				fprintf(context->_fpLOG, "\n%I64d CVO control thread; stop_signalled ...", tNow) ;
				fflush(context->_fpLOG) ;
				}
			}
		nRunning = 0 ;
		for (i = 0 ; i < nWorkers ; i++) {
			if (0 == Workers[i]._ThreadHandle) continue ;
			if (0 != stop_signalled) {
				if (NULL != context->_fpLOG) {
					tNow = ARE::GetTimeInMilliseconds() ;
					fprintf(context->_fpLOG, "\n%I64d CVO control thread; stopping worker %d ...", tNow, i) ;
					fflush(context->_fpLOG) ;
					}
				Workers[i]._ThreadStop = true ;
				}
			++nRunning ;
			}
		if (nRunning <= 0 || 0 != stop_signalled) 
			break ;
		}
	// wait for threads to stop
	while (nRunning > 0) {
		SLEEP(50) ;
		nRunning = 0 ;
		for (i = 0 ; i < nWorkers ; i++) {
			if (0 == Workers[i]._ThreadHandle) continue ;
			++nRunning ;
			}
		}
	if (NULL != context->_fpLOG) {
		tNow = ARE::GetTimeInMilliseconds() ;
		fprintf(context->_fpLOG, "\n%I64d CVO control thread; all worker threads have closed ...", tNow) ;
		fflush(context->_fpLOG) ;
		}

/*
	if (TimeLimitInSeconds > 0) {
		Sleep(1000*TimeLimitInSeconds) ;
		for (i = 0 ; i < nWorkers ; i++) 
			Workers[i]._ThreadStop = true ;
		for (i = 0 ; i < nWorkers ; i++) {
			if (0 == Workers[i]._ThreadHandle) continue ;
			DWORD r = WaitForSingleObject((HANDLE) Workers[i]._ThreadHandle, 10) ;
			if (WAIT_TIMEOUT == r) 
				TerminateThread((HANDLE) Workers[i]._ThreadHandle, 0) ;
			CloseHandle((HANDLE) Workers[i]._ThreadHandle) ;
			Workers[i]._ThreadHandle = 0 ;
			}
		}
*/

	ret = 0 ;
done :
	tNow = context->_tEnd = ARE::GetTimeInMilliseconds() ;
	if (NULL != context->_fpLOG) {
		fprintf(context->_fpLOG, "\n%I64d CVO control thread; done; time=%I64dmsec ...", tNow, (int) (context->_tEnd - context->_tStart)) ;
		fflush(context->_fpLOG) ;
		}
//	for (i = 0 ; i < nWorkers ; i++) 
//		nRuns += Workers[i]._nRunsDone ;
	if (NULL != Workers) 
		delete [] Workers ;
	if (best_order._Width < p.N()) {
		// some ordering was found
		}

	// remove redundant fill edges
	if (false && best_order._Width > 0) {
		ARE::Graph g ;
		g = OriginalGraph ;
		g._OrderLength = p.N() ;
		for (i = 0 ; i < p.N() ; i++) {
			int v_i = best_order._VarListInElimOrder[i] ;
			g._VarType[i] = 0 ;
			g._PosOfVarInList[v_i] = i ;
			g._VarElimOrder[i] = v_i ;
			}
		g._nFillEdges = best_order._nFillEdges ;
		g.RemoveRedundantFillEdges() ;
		int nEdgesRemoved = (OriginalGraph._nEdges + best_order._nFillEdges) - g._nEdges ;
		fprintf(context->_fpLOG, "\n# of redundant fill edges = %d", nEdgesRemoved) ;
		fflush(context->_fpLOG) ;
		}

	context->_ThreadHandle = 0 ;
#if defined WINDOWS || _WINDOWS
	_endthreadex(0) ;
	return 0  ;
#else
	return NULL ;
#endif
}


int ARE::VarElimOrderComp::CVOcontext::CreateCVOthread(void)
{
#if defined WINDOWS || _WINDOWS
	_ThreadHandle = _beginthreadex(NULL, 0, CVOThreadFn, this, 0, NULL) ;
#else
	pthread_create(&_ThreadHandle, NULL, CVOThreadFn, this) ; // TODO third argument
#endif
	return 0 != _ThreadHandle ? 0 : 1 ;
}


int ARE::VarElimOrderComp::CVOcontext::StopCVOthread(void)
{
	if (0 == _ThreadHandle) {
		if (NULL != _fpLOG) {
			INT64 tNowLog = ARE::GetTimeInMilliseconds() ;
			fprintf(_fpLOG, "\n%I64d CVO_th : stop variable order computation; already stopped ...", tNowLog) ;
			fflush(_fpLOG) ;
			}
		return 0 ;
		}
#if defined WINDOWS || _WINDOWS
	InterlockedCompareExchange(&_StopAndExit, 1, 0) ;
#else
  pthread_mutex_lock(&stopSignalMutex);
  if (_StopAndExit == 0) {
    _StopAndExit = 1;
  }
  pthread_mutex_unlock(&stopSignalMutex);
#endif
	INT64 tStart = ARE::GetTimeInMilliseconds() ;
	if (NULL != _fpLOG) {
		fprintf(_fpLOG, "\n%I64d CVO_th : stop variable order computation; stop signalled, will wait ...", tStart) ;
		fflush(_fpLOG) ;
		}
	while (true) {
		SLEEP(50) ;
		if (0 == _ThreadHandle) 
			break ;
		INT64 tNow = ARE::GetTimeInMilliseconds() ;
		INT64 dt = tNow - tStart ;
		if (dt > 10000) {
			// we asked the thread to stop and waited for it to stop, but it won't stop, so kill the thread.
#if defined WINDOWS || _WINDOWS
			TerminateThread((HANDLE) _ThreadHandle, 0) ;
			CloseHandle((HANDLE) _ThreadHandle) ;
			_ThreadHandle = 0 ;
#else
			// TODO : handle linux
      pthread_exit(&_ThreadHandle);
#endif
			if (NULL != _fpLOG) {
				INT64 tNowLog = ARE::GetTimeInMilliseconds() ;
				fprintf(_fpLOG, "\n%I64d CVO_th : stop variable order computation, hard kill ...", tNowLog) ;
				fflush(_fpLOG) ;
				}
			break ;
			}
		}
	return 0 ;
}


int ARE::VarElimOrderComp::Compute(
	// IN
	const std::string & uaifile, 
	ARE::VarElimOrderComp::ObjectiveToMinimize objcode,
	ARE::VarElimOrderComp::NextVarPickCriteria algcode,
	int nthreads, 
	int nrunstodo, 
	INT64 TimeLimitInMilliSeconds, 
	int nRP, 
	double eRP, 
	bool PerformSingletonConsistencyChecking, 
	bool earlyterminationofbasic_W, 
	bool earlyterminationofbasic_C, 
	// OUT
	ARE::VarElimOrderComp::Order & BestOrder, 
	ARE::VarElimOrderComp::CVOcontext *Context
	)
{
	if (uaifile.length() < 1) 
		return 1 ;

	int ret = 1 ;

	char strDT[64] ;
	INT64 tNow = 0 ;

	INT64 tStart;
  INT64 tStopSignalled;
	int i;

	ARE::ARP p("") ;
	ARE::VarElimOrderComp::CVOcontext *cvocontext = Context ;
	ARE::VarElimOrderComp::CVOcontext *localCVOcontext = NULL ;
	if (NULL == cvocontext) {
		localCVOcontext = new ARE::VarElimOrderComp::CVOcontext ;
		if (NULL == localCVOcontext) 
			goto done ;
		cvocontext = localCVOcontext ;
		}
	GetCurrentDTmsec(strDT, tNow) ;
	printf("\n%s CVO : Start ...", strDT) ;
	if (NULL != cvocontext->_fpLOG) {
		fprintf(cvocontext->_fpLOG, "\n%s CVO : Start ...", strDT) ;
		fflush(cvocontext->_fpLOG) ;
		}

	cvocontext->_Problem = &p ;
	cvocontext->_AlgCode = algcode ;
	cvocontext->_ObjCode = objcode ;

	cvocontext->_nRunsToDoMax = nrunstodo ;
	if (cvocontext->_nRunsToDoMax < 1) 
		cvocontext->_nRunsToDoMax = 1 ;
	else if (cvocontext->_nRunsToDoMax > 1000000000) 
		cvocontext->_nRunsToDoMax = 1000000000 ;
	cvocontext->_nRunsToDoMin = cvocontext->_nRunsToDoMax ;

	if (TimeLimitInMilliSeconds < 0) 
		TimeLimitInMilliSeconds = 0 ;
	else if (TimeLimitInMilliSeconds > 86400000) 
		TimeLimitInMilliSeconds = 86400000 ;
	cvocontext->_TimeLimitInMilliSeconds = TimeLimitInMilliSeconds ;

	cvocontext->_eRandomPick = eRP ;
	cvocontext->_nRandomPick = nRP ;
	if (cvocontext->_nRandomPick < 1) 
		cvocontext->_nRandomPick = 1 ;

	cvocontext->_EarlyTerminationOfBasic_W = earlyterminationofbasic_W ;
	cvocontext->_EarlyTerminationOfBasic_C = earlyterminationofbasic_C ;
	cvocontext->_BestOrder = &BestOrder ;

	cvocontext->_nThreads = nthreads ;
	if (cvocontext->_nThreads > cvocontext->_nRunsToDoMax-1) 
		cvocontext->_nThreads = cvocontext->_nRunsToDoMax-1 ;
	else if (cvocontext->_nThreads < 1) 
		cvocontext->_nThreads = 1 ;
// DEBUGGG
//nThreads = 1 ;

	// at least one of nRunsToDo/TimeLimitInSeconds has to be given
	if (cvocontext->_nRunsToDoMin < 1 && cvocontext->_TimeLimitInMilliSeconds < 1) 
		cvocontext->_nRunsToDoMin = 1 ;

	// uaiFN may have dir in it; extract filename.
	i = uaifile.length() - 1 ;
	for (; i >= 0 ; i--) {
#ifdef LINUX
		if ('\\' == uaifile[i] || '/' == uaifile[i]) break ;
#else
    if ('\\' == uaifile[i] || '//' == uaifile[i]) break ;
#endif
		}
	{
	std::string fn(uaifile.substr(i+1)) ;
	p.SetName(fn) ;
	}

	if (0 != p.LoadFromFile(uaifile)) {
		ret = 2 ;
		printf("\nload failed ...") ;
		goto done ;
		}
	if (0 != p.PerformPostConstructionAnalysis()) {
		ret = 3 ;
		printf("\nPerformPostConstructionAnalysis failed ...") ;
		goto done ;
		}
	if (p.N() < 1) {
		ret = 4 ;
		printf("\nN=%d; will exit ...", p.N()) ;
		goto done ;
		}
/*
	if (p.N() > (1 << 16)) {
//		printf("\nN=%d; currently nVars limit is 64K since AVL tree (edgesadded) with key type long is used to store 2 variable indeces; will exit ...", p.N()) ;
		printf("\nN=%d; currently nVars limit is 64K since heap uses 2 bytes for cost and 2 bytes for var index, for total 4 bytes; change heap key to __in64; will exit ...", p.N()) ;
		goto done ;
		}
*/

	tNow = 0 ;
	GetCurrentDTmsec(strDT, tNow) ;
	printf("\n%s File loaded; start preprocessing ...", strDT) ;
	if (NULL != cvocontext->_fpLOG) {
		fprintf(cvocontext->_fpLOG, "\n%s File loaded; start preprocessing ...", strDT) ;
		fflush(cvocontext->_fpLOG) ;
		}

	// eliminate singleton-domain variables; do this before ordering is computed; this is easy and should be done by any algorithm processing the network
	printf("\n%s EliminateSingletonDomainVariables(): start, N=%d ...", strDT, (int) p.N()) ;
	if (NULL != cvocontext->_fpLOG) {
		fprintf(cvocontext->_fpLOG, "\n%s EliminateSingletonDomainVariables(): start, N=%d ...", strDT, (int) p.N()) ;
		fflush(cvocontext->_fpLOG) ;
		}
	p.EliminateSingletonDomainVariables() ;
	tNow = 0 ;
	GetCurrentDTmsec(strDT, tNow) ;
	printf("\n%s EliminateSingletonDomainVariables(): done; nSingletonDomainVariables=%d ...", strDT, (int) p.nSingletonDomainVariables()) ;
	if (NULL != cvocontext->_fpLOG) {
		fprintf(cvocontext->_fpLOG, "\n%s EliminateSingletonDomainVariables(): done; nSingletonDomainVariables=%d ...", strDT, (int) p.nSingletonDomainVariables()) ;
		fflush(cvocontext->_fpLOG) ;
		}

	// prune domains of variables by checking which values participate in no complete assignment with probability>0.
	// basically, we compute a minimal domain for each variable.
	if (PerformSingletonConsistencyChecking) {
		int nNewSingletonDomainVariablesAfterSingletonConsistency = 0 ;
		int res_sc = p.ComputeSingletonConsistency(nNewSingletonDomainVariablesAfterSingletonConsistency) ;
		tNow = 0 ;
		GetCurrentDTmsec(strDT, tNow) ;
		if (-1 == res_sc) {
			// domain of some variable is empty
			printf("\n%s ComputeSingletonConsistency(): domain of some variable is empty; will quit ...", strDT) ;
			if (NULL != cvocontext->_fpLOG) {
				fprintf(cvocontext->_fpLOG, "\n%s ComputeSingletonConsistency(): domain of some variable is empty; will quit ...", strDT) ;
				fflush(cvocontext->_fpLOG) ;
				}
			goto done ;
			}
		printf("\n%s ComputeSingletonConsistency(): nNewSingletonDomainVariablesAfterSingletonConsistency = %d", strDT, nNewSingletonDomainVariablesAfterSingletonConsistency) ;
		if (NULL != cvocontext->_fpLOG) {
			fprintf(cvocontext->_fpLOG, "\n%s ComputeSingletonConsistency(): nNewSingletonDomainVariablesAfterSingletonConsistency = %d", strDT, nNewSingletonDomainVariablesAfterSingletonConsistency) ;
			fflush(cvocontext->_fpLOG) ;
			}
		if (nNewSingletonDomainVariablesAfterSingletonConsistency > 0) 
			p.EliminateSingletonDomainVariables() ;
		}
	else {
		printf("\n%s Singleton-Consistency check not requested ...", strDT) ;
		if (NULL != cvocontext->_fpLOG) {
			fprintf(cvocontext->_fpLOG, "\n%s Singleton-Consistency check not requested ...", strDT) ;
			fflush(cvocontext->_fpLOG) ;
			}
		}

	tNow = 0 ;
	GetCurrentDTmsec(strDT, tNow) ;
	printf("\n%s Launching CVO thread ...", strDT) ;
	if (NULL != cvocontext->_fpLOG) {
		fprintf(cvocontext->_fpLOG, "\n%s Launching CVO thread ...", strDT) ;
		fflush(cvocontext->_fpLOG) ;
		}

	cvocontext->CreateCVOthread() ;
	if (0 == cvocontext->_ThreadHandle) {
		printf("\nFAILED to create cvo thread ...") ;
		if (NULL != cvocontext->_fpLOG) {
			fprintf(cvocontext->_fpLOG, "\nFAILED to create cvo thread ...") ;
			fflush(cvocontext->_fpLOG) ;
			}
		ret = 8 ;
		goto done ;
		}

	tStart = ARE::GetTimeInMilliseconds();
  tStopSignalled = 0 ;
	while (true) {
		SLEEP(100) ;
		if (0 == cvocontext->_ThreadHandle) 
			break ;
		INT64 tNow = ARE::GetTimeInMilliseconds() ;
		INT64 dt = tNow - tStart ;
		if (dt < cvocontext->_TimeLimitInMilliSeconds) 
			continue ;
		if (0 == tStopSignalled) {
			tStopSignalled = tNow ;
#if defined WINDOWS || _WINDOWS
      InterlockedCompareExchange(&(cvocontext->_StopAndExit), 1, 0) ;
#else
      pthread_mutex_lock(&stopSignalMutex);
      if (cvocontext->_StopAndExit == 0) {
        cvocontext->_StopAndExit = 1;
      }
      pthread_mutex_unlock(&stopSignalMutex);
#endif
			continue ;
			}
		dt = tNow - tStopSignalled ;
		if (dt > 10000) {
			// we asked the thread to stop and waited for it to stop, but it won't stop, so kill the thread.
#if defined WINDOWS || _WINDOWS
			TerminateThread((HANDLE) cvocontext->_ThreadHandle, 0) ;
			CloseHandle((HANDLE) cvocontext->_ThreadHandle) ;
			cvocontext->_ThreadHandle = 0 ;
#else
			// TODO
      pthread_exit(&cvocontext->_ThreadHandle);
#endif 
			break ;
			}
		}
	tNow = 0 ;
	GetCurrentDTmsec(strDT, tNow) ;
	printf("\n%s CVO : thread has closed ...", strDT) ;
	if (NULL != cvocontext->_fpLOG) {
		fprintf(cvocontext->_fpLOG, "\n%s CVO : thread has closed ...", strDT) ;
		fflush(cvocontext->_fpLOG) ;
		}

	ret = 0 ;
done :

	if (NULL != localCVOcontext) 
		delete localCVOcontext ;

	return ret ;
}

