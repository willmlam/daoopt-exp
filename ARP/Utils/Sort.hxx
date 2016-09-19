/*
	This module contains sorting function's declarations.

	Kalev Kask, March 2002.
*/

#ifndef SORT_HXX_INCLUDED
#define SORT_HXX_INCLUDED

#include <inttypes.h>

#ifdef LINUX
typedef int64_t INT64 ;
#else 
typedef __int64 INT64 ;
#endif

void QuickSortShort(short *input_key, uint32_t len, int32_t *input_data, 
	// stack; passed in, to make this function thread-safe
	int32_t left[32], int32_t right[32]) ;

void QuickSortDouble(double *input_key, uint32_t len, int32_t *input_data, 
	// stack; passed in, to make this function thread-safe
	int32_t left[32], int32_t right[32]) ;
void QuickSortDouble_Descending(double *input_key, uint32_t len, int32_t *input_data, 
	// stack; passed in, to make this function thread-safe
	int32_t left[32], int32_t right[32]) ;
int SortCheckDouble(double *key, uint32_t len) ;

void QuickSortLong(int32_t *input_key, uint32_t len, int32_t *input_data, 
	// stack; passed in, to make this function thread-safe
	int32_t left[32], int32_t right[32]) ;
void QuickSortLong_Descending(int32_t *input_key, uint32_t len, int32_t *input_data, 
	// stack; passed in, to make this function thread-safe
	int32_t left[32], int32_t right[32]) ;
void QuickSortLong2(int32_t *input_key, uint32_t len, 
	// stack; passed in, to make this function thread-safe
	int32_t left[32], int32_t right[32]) ;
int SortCheckint32_t(int32_t *key, uint32_t len) ;
void QuickSortLong_i64(int32_t *input_key, uint32_t len, int64_t *input_data,
	// stack; passed in, to make this function thread-safe
	int32_t left[32], int32_t right[32]) ;
void QuickSortLong_i64_Descending(int32_t *input_key, uint32_t len, int64_t *input_data, 
	// stack; passed in, to make this function thread-safe
	int32_t left[32], int32_t right[32]) ;

void QuickSorti64(int64_t *input_key, uint32_t len, int32_t *input_data,
	// stack; passed in, to make this function thread-safe
	int32_t left[32], int32_t right[32]) ;
void QuickSorti64_i64(int64_t *input_key, uint32_t len, int64_t *input_data,
	// stack; passed in, to make this function thread-safe
	int32_t left[32], int32_t right[32]) ;
int SortChecki64(int64_t *key, uint32_t len) ;

void QuickSortWchar_t(const wchar_t **input_key, uint32_t len, int32_t *input_data, 
	// stack; passed in, to make this function thread-safe
	int32_t left[32], int32_t right[32]) ;
void QuickSortWchar_t(const wchar_t **input_key, uint32_t len, 
	// stack; passed in, to make this function thread-safe
	int32_t left[32], int32_t right[32]) ;
void QuickSortWchar_t_Descending(const wchar_t **input_key, uint32_t len, int32_t *input_data, 
	// stack; passed in, to make this function thread-safe
	int32_t left[32], int32_t right[32]) ;

void QuickSortChar(const char **input_key, uint32_t len, int32_t *input_data, 
	// stack; passed in, to make this function thread-safe
	int32_t left[32], int32_t right[32]) ;

void QuickSortMem(const char **input_key, uint32_t keysize, uint32_t len, int32_t *input_data, 
	// stack; passed in, to make this function thread-safe
	int32_t left[32], int32_t right[32]) ;

// return true iff Obj1 is greater than Obj2 (i.e. Obj1 is after Obj2).
typedef bool (*QuickSortGreaterCompFn)(void *Obj1, void *Obj2) ;
void QuickSort(void *objarray[], uint32_t len, 
	// stack; passed in, to make this function thread-safe
	int32_t left[32], int32_t right[32],
	// fn to compare 2 objects
	QuickSortGreaterCompFn CompGreaterOperator) ;

#endif // SORT_HXX_INCLUDED
