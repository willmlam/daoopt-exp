#include <stdlib.h>
#include "time.h"
#include "math.h"

#if defined WINDOWS || _WINDOWS
#include "process.h"    /* _beginthread, _endthread */
#endif // WINDOWS

#include "Utils/Sort.hxx"

#include "Globals.hxx"
#include "Utils/MiscUtils.hxx"
#include "Utils/FnExecutionThread.hxx"

#if defined(LINUX)
pthread_mutex_t FnExecutionThread::_stopSignalMutex = PTHREAD_MUTEX_INITIALIZER ;
#endif