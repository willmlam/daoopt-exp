/*
	This module contains misc utils declarations.

	Kalev Kask, March 2002.
*/

#ifndef MISCUTILS_HXX_INCLUDED
#define MISCUTILS_HXX_INCLUDED

#include <climits>
#include <list>
#include <string>

typedef int64_t INT64 ;

namespace ARE {

int64_t GetTimeInMilliseconds(void) ;

int ExtractVarValuePairs(/* IN */ char *BUF, int L, /* OUT */ std::list<std::pair<std::string,std::string>> & List) ;
int ExtractParameterValue(/* IN */ std::string & Paramater, std::list<std::pair<std::string,std::string>> AssignmentList, /* OUT */ std::string & Value) ;

}

#endif // SORT_HXX_INCLUDED
