#ifndef ARE_Globals_HXX_INCLUDED
#define ARE_Globals_HXX_INCLUDED

#include <stdlib.h>
#include <stdio.h>
#include <string>

#include "Utils/Mutex.h"

#ifdef LINUX
typedef int64_t __int64;
typedef unsigned int DWORD;
#include <cfloat>
#include <cstring>

#include <climits>
#include <limits>
#define _FPCLASS_NINF (- std::numeric_limits<double>::infinity())
#define _I64_MAX (std::numeric_limits<__int64>::max())

#include <cmath>
#endif

// macro to get system time; others.
#if defined LINUX
#include <unistd.h>
#define GETCLOCK() clock()
#define SLEEP(X) usleep(1000*X)
#define stricmp strcasecmp
#else
#define GETCLOCK() GetTickCount()
#define SLEEP(X) Sleep(X)
#endif

#define LOG_OF_SUM_OF_TWO_NUMBERS_GIVEN_AS_LOGS(sum,x,y) { if (x == y) sum = 0.30102999566398119521373889472449 + x ; else if (x > y) sum = x + log10(1.0 + pow(10.0, y - x)) ; else sum = y + log10(1.0 + pow(10.0, x - y)) ; }
#define LOG_OF_SUB_OF_TWO_NUMBERS_GIVEN_AS_LOGS(sub,x,y) { if (x == y) sub = - std::numeric_limits<double>::infinity() ; else if (x > y) sub = x + log10(1.0 - pow(10.0, y - x)) ; else sub = y + log10(1.0 - pow(10.0, x - y)) ; }

namespace ARE
{

// Variable Elimination Complexity = # operations to eliminate a variable (from a bucket/cluster e.g.).
// this typically is the product of domain sizes of all variables involved.
// set single-var elimination complexity limit at 10^15 = 1000 trillion
//#define InfiniteSingleVarElimComplexity 1000000000000000LL
//#define InfiniteVarElimComplexity       1000000000000000LL
#define InfiniteSingleVarElimComplexity_log 999.0
#define InfiniteVarElimComplexity_log       999.0
// note that _I64_MAX = 9,223,372,036,854,775,807
//          out limit =     1,000,000,000,000,000
// note that we want enough distance between MaxSingleVarElimComplexity and _I64_MAX so that can detect going over MaxSingleVarElimComplexity without overflowing the int64.

#define MAX_NUM_VARIABLES_PER_BUCKET	128
#define MAX_NUM_FUNCTIONS_PER_BUCKET	256
#define MAX_NUM_ARGUMENTS_PER_FUNCTION	128
#define MAX_NUM_VALUES_PER_VAR_DOMAIN	2048
#define MAX_NUM_FACTORS_PER_FN_DOMAIN	64
#define MAX_NUM_BUCKETS 65536
#define MAX_DEGREE_OF_GRAPH_NODE		1024
#define MAX_NUM_VARIABLES_PER_PROBLEM   1000000

#define ERRORCODE_generic									100
#define ERRORCODE_too_many_variables						101
#define ERRORCODE_too_many_functions						102
#define ERRORCODE_input_FTB_not_computed					103
#define ERRORCODE_input_FTB_computed_but_failed_to_fetch	104
#define ERRORCODE_fn_ptr_is_NULL							105
#define ERRORCODE_cannot_open_file							106
#define ERRORCODE_cannot_allocate_buffer_size_too_large		107
#define ERRORCODE_file_too_large							108
#define ERRORCODE_problem_type_unknown						109
#define ERRORCODE_EliminationComplexityWayTooLarge			110
#define ERRORCODE_EliminationWidthTooLarge					111
#define ERRORCODE_EliminationComplexityTooLarge				112
#define ERRORCODE_out_of_memory								113
#define ERRORCODE_out_of_memory_CVOedges					114
#define ERRORCODE_out_of_memory_CVOadjvarlist				115
#define ERRORCODE_memory_allocation_failure					116
#define ERRORCODE_NoVariablesLeftToPickFrom					117
#define ERRORCODE_DuplicateVarInFnArgumentList				118
#define ERRORCODE_BadVarInFnArgumentList					119
#define ERRORCODE_VarDegreeTooLarge							120
#define ERRORCODE_InvalidInputData							121
#define ERRORCODE_VarDomainSizeTooLarge						122

#define FN_COBINATION_TYPE_NONE		0
#define FN_COBINATION_TYPE_PROD		1
#define FN_COBINATION_TYPE_SUM		2

#define VAR_ELIMINATION_TYPE_NONE	0
#define VAR_ELIMINATION_TYPE_SUM	1
#define VAR_ELIMINATION_TYPE_MAX	2
#define VAR_ELIMINATION_TYPE_MIN	3

extern FILE *fpLOG ;
extern ARE::utils::RecursiveMutex LOGMutex ;

extern const double neg_infinity ;
extern       double pos_infinity ;

// returns 0 iff ok, error code otherwise.
int LoadVarOrderFromFile(
	// IN
	const std::string & fnVO, 
	int N,
	// OUT
	int *VarListBeforeEvidence, 
	int & Width
	) ;

int Initialize(void) ;

// remove from Set1 everything not in Set2; assume both Set1 and Set2 are sorted in increasing order. on return, Set1 is sorted.
inline void SetIntersection(int *Set1, int & N1, const int *Set2, int N2)
{
	int *s1 = Set1 ; const int *s2 = Set2 ;
	int *s1E = s1+N1 ; const int *s2E = s2+N2 ;
	while (s1 < s1E && s2 < s2E) {
		if (*s1 == *s2) 
			{ s1++ ; s2++ ; }
		else if (*s1 < *s2) {
			--s1E ;
			for (int *s = s1 ; s < s1E ; s++) *s = *(s+1) ;
			}
		else 
			{ s2++ ; }
		}
	N1 = s1 - Set1 ;
}

// remove from Set1 everything in Set2; assume both Set1 and Set2 are sorted in increasing order. on return, Set1 is sorted.
inline void SetMinus(int *Set1, int & N1, const int *Set2, int N2)
{
	int *s1 = Set1 ; const int *s2 = Set2 ;
	int *s1E = s1+N1 ; const int *s2E = s2+N2 ;
	while (s1 < s1E && s2 < s2E) {
		if (*s1 == *s2) {
			--s1E ;
			for (int *s = s1 ; s < s1E ; s++) *s = *(s+1) ;
			s2++ ;
			}
		else if (*s1 < *s2) 
			{ s1++ ; }
		else 
			{ s2++ ; }
		}
	N1 = s1E - Set1 ;
}

} // namespace ARE

// number factorization; done by division (not very efficient; ok for small numbers).
// output is sorted in increasing order.
char PrimeFactoringByTrialDivision(char N, char factors[128]) ;

#endif // ARE_Globals_HXX_INCLUDED
