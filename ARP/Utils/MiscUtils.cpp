/*
	This module contains misc routines.

	Kalev Kask, March 2002.
*/

#include <stdlib.h>
#include "Utils/MiscUtils.hxx"

#if defined WINDOWS || _WINDOWS
#include <windows.h>
#endif

#include <cstring>
#include <string>
#include <time.h>
#include <errno.h>

#ifdef LINUX
#include <chrono>
#endif


INT64 ARE::GetTimeInMilliseconds(void)
{
#if defined WINDOWS || _WINDOWS
	FILETIME ft ;
	GetSystemTimeAsFileTime(&ft) ;
	__int64 t = ft.dwHighDateTime ;
	t <<= 32 ;
	t += ft.dwLowDateTime ;
	// 64-bit value representing the number of 100-nanosecond intervals since January 1, 1601 (UTC).
	// convert to milliseconds
	t /= 10000 ;
	return t ;
#else
	return std::chrono::duration_cast<std::chrono::milliseconds>(std::chrono::high_resolution_clock::now().time_since_epoch()).count();
#endif
}


int ARE::ExtractVarValuePairs(char *BUF, int L, std::list<std::pair<std::string,std::string>> & AssignmentList)
{
	AssignmentList.clear() ;
	if (NULL == BUF) 
		return 0 ;
	if (L < 0) 
		L = strlen(BUF) ;

	int i, j, n = 0 ;
	char *s = BUF ;
	std::string sN, sV ;
	for (i = 0 ; i < L ; i++) {
		if ('\n' != BUF[i]) continue ;
		int l = BUF + i - s ;
		if (l > 0) {
			// [s, s+l) is a string
			for (j = 0 ; j < l ; j++) 
				{ if ('=' == s[j]) break ; }
			if (j < l) {
				const char *sName = s ;
				const char *sValue = s + j + 1 ;
				int lN = j ;
				int lV = l - j - 1 ;
				if (lN > 0) sN.assign(sName, lN) ; else sN.erase() ;
				if (lV > 0) sV.assign(sValue, lV) ; else sV.erase() ;
				std::pair<std::string,std::string> assignment(sN,sV) ;
				AssignmentList.push_back(assignment) ;
				}
			}
		if ('\r' == BUF[i+1]) 
			i++ ;
		s = BUF + i + 1 ;
		}
	int l = BUF + i - s ;
	if (l > 0) {
		// [s, s+l) is a string
		for (j = 0 ; j < l ; j++) 
			{ if ('=' == s[j]) break ; }
		if (j < l) {
			const char *sName = s ;
			const char *sValue = s + j + 1 ;
			int lN = j ;
			int lV = l - j - 1 ;
			if (lN > 0) sN.assign(sName, lN) ; else sN.erase() ;
			if (lV > 0) sV.assign(sValue, lV) ; else sV.erase() ;
			std::pair<std::string,std::string> assignment(sN,sV) ;
			AssignmentList.push_back(assignment) ;
			}
		}

	return n ;
}


int ARE::ExtractParameterValue(/* IN */ std::string & Paramater, std::list<std::pair<std::string,std::string>> AssignmentList, /* OUT */ std::string & Value)
{
	Value.erase() ;
	for (std::list<std::pair<std::string,std::string>>::iterator i = AssignmentList.begin() ; i != AssignmentList.end() ; i++) {
		std::pair<std::string,std::string> & assignment = *i ;
		std::string & sN = assignment.first ;
		if (Paramater == sN) 
			{ Value = assignment.second ; return 0 ; }
		}

	return 1 ;
}


