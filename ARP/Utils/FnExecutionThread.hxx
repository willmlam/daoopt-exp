#ifndef FnExecutionThread_HXX_INCLUDED
#define FnExecutionThread_HXX_INCLUDED

#include <stdlib.h>
#include <Utils/MiscUtils.hxx>

#if defined WINDOWS || _WINDOWS
typedef unsigned int (__stdcall *pFnExecutionThreadFn)(void *X) ;
#elif defined (LINUX)
typedef void *(*pFnExecutionThreadFn)(void *X) ;
#endif 


class FnExecutionThread
{

public :
#if defined (LINUX)
  static pthread_mutex_t _stopSignalMutex ;
#endif
	// stop/exit computation thread
	long volatile _StopAndExit ;

	// ptr to the execution thread
	pFnExecutionThreadFn _Thread ;

#if defined WINDOWS || _WINDOWS
	uintptr_t _ThreadHandle ;
#elif defined (LINUX)
	pthread_t _ThreadHandle ;
#endif 

	INT64 _tStart, _tEnd, _tToStop, _RunTimeInMilliseconds ;

	int CreateThread(pFnExecutionThreadFn pThreadFn, void *X)
	{
		_tStart = ARE::GetTimeInMilliseconds() ;
#if defined WINDOWS || _WINDOWS
		_ThreadHandle = _beginthreadex(NULL, 0, 0 == pThreadFn ? _Thread : pThreadFn, X, 0, NULL) ;
#else
		pthread_create(&_ThreadHandle, NULL, pThreadFn, X) ; // TODO third argument
#endif
		return 0 != _ThreadHandle ? 0 : 1 ;
	}

	int StopThread(void)
	{
		if (0 == _ThreadHandle) 
			return 0 ;
#if defined WINDOWS || _WINDOWS
		InterlockedCompareExchange(&_StopAndExit, 1, 0) ;
#else
    pthread_mutex_lock(&_stopSignalMutex);
    if (_StopAndExit == 0) {
      _StopAndExit = 1;
    }
    pthread_mutex_unlock(&_stopSignalMutex);
#endif
		_tToStop = ARE::GetTimeInMilliseconds() ;
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
#else
				// TODO : handle linux
        pthread_exit(&_ThreadHandle);
#endif
				_ThreadHandle = 0 ;
				break ;
				}
			}
		_tEnd = ARE::GetTimeInMilliseconds() ;
		return 0 ;
	}

public :

	FnExecutionThread(void) :
		_StopAndExit(0), 
		_Thread(0), 
		_ThreadHandle(0), 
		_tStart(0), _tEnd(0), _tToStop(0)
	{
	}
	virtual ~FnExecutionThread(void)
	{
		StopThread() ;
	}
} ;

#endif // FnExecutionThread_HXX_INCLUDED
