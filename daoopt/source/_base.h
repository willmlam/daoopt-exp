/*
 * _base.h
 *
 *  Copyright (C) 2008-2012 Lars Otten
 *  This file is part of DAOOPT.
 *
 *  DAOOPT is free software: you can redistribute it and/or modify it
 *  under the terms of the GNU General Public License as published by
 *  the Free Software Foundation, either version 3 of the License, or
 *  (at your option) any later version.
 *
 *  DAOOPT is distributed in the hope that it will be useful, but
 *  WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU General Public License for more details.
 *
 *  You should have received a copy of the GNU General Public License
 *  along with DAOOPT.  If not, see <http://www.gnu.org/licenses/>.
 *  
 *  Created on: Oct 9, 2008
 *      Author: Lars Otten <lotten@ics.uci.edu>
 */


/* for easier debugging output and nicer code,
 * redefined for every file that includes _base.h */
#ifdef DIAG
#undef DIAG
#endif

#ifdef DEBUG
#define DIAG(X) { X ; }
#else
#define DIAG(X) ;
#endif

#ifndef _BASE_H_
#define _BASE_H_

/* some customizable definitions */
#include "DEFINES.h"
#include "float.h"

/* either LINUX or WINDOWS has to be defined -> should happen in the Makefile */
//#define LINUX
//#define WINDOWS
#if !defined(LINUX) && !defined(WINDOWS)
#error Either LINUX or WINDOWS identifiers have to be defined (for preprocessor)!
#endif

/* define if data files (solution and subproblem files) should be written in binary */
//#define BINARY_DATAFILES

/* define if gmp lib should be used for large numbers */
//#define USE_GMP


#ifdef NOTHREADS
#undef PARALLEL_DYNAMIC
#undef PARALLEL_STATIC
#endif

#ifdef PARALLEL_DYNAMIC

/* Boost thread libraries */
#include "boost/thread.hpp"
//#include "boost/thread/shared_mutex.hpp"

#define GETLOCK(X,Y) boost::mutex::scoped_lock Y ( X )
#define CONDWAIT(X,Y) ( X ).wait( Y )
#define NOTIFY(X) ( X ).notify_one()
#define NOTIFYALL(X) ( X ).notify_all()
#define INCREASE(X)  ++( X )
#define DECREASE(X) --( X )

/* static IO mutex for console output */
namespace daoopt {
static boost::mutex mtx_io;
}  // namespace daoopt

#else

#define GETLOCK(X,Y) ;
#define CONDWAIT(X,Y) ;
#define NOTIFY(X) ;
#define NOTIFYALL(X) ;
#define INCREASE(X) ;
#define DECREASE(X) ;

#endif


#ifdef PARALLEL_DYNAMIC
#define PAR_ONLY(X) X
#else
#define PAR_ONLY(X) ;
#endif

/* define if all computation should be done in log format (recommended) */
#define USE_LOG

#ifdef USE_LOG

#define OP_TIMES +
#define OP_TIMESEQ +=
#define OP_DIVIDE -
#define OP_DIVIDEEQ -=
#define OP_ROOT(X,N) X * (1.0 / N);
#define ELEM_ZERO (- std::numeric_limits<double>::infinity() )
#define ELEM_ONE 0.0

//#define OP_PLUS(X,Y) ( (X==Y) ? (X+log10(2.0)) : ( ( X > Y ) ? ( X + log10(1.0 + pow(10.0, Y-X ))) : ( Y + log10(1.0 +pow(10.0, X-Y ) )) ) )
#define OP_PLUS(X,Y) ( ( X > Y ) ? ( X + log10(1.0 + pow(10.0, Y-X ))) : ( Y + log10(1.0 +pow(10.0, X-Y ) )) )
//#define OP_PLUS(X,Y) ( log10( pow(10.0, X ) + pow(10.0, Y ) ) )
#define OP_MINUS(X,Y) ( ( X > Y ) ? ( X + log10(1.0 - pow(10.0, Y-X ))) : ( Y + log10(pow(10.0, X-Y ) - 1.0 )) )
//#define OP_MINUS(X,Y) ( ( X > Y ) ? ( X + log(1.0 - exp( Y-X ))) : ( Y + log(exp( X-Y ) - 1.0 )) )


#define ELEM_ENCODE(X) log10( X )
#define ELEM_DECODE(X) pow(10.0, X )
//#define ELEM_ENCODE(X) log( X )
//#define ELEM_DECODE(X) exp( X )

#define SCALE_LOG(X) ( X )
#define SCALE_NORM(X) pow(10.0, X )
//#define SCALE_NORM(X) exp( X )

#else

#define OP_TIMES *
#define OP_TIMESEQ *=
#define OP_DIVIDE /
#define OP_DIVIDEEQ /=
#define OP_ROOT(X,N) pow(X, 1.0/N)
#define ELEM_ZERO 0.0
#define ELEM_ONE 1.0

#define OP_PLUS(X,Y) (X+Y)
#define OP_MINUS(X,Y) (X-Y)

#define ELEM_ENCODE(X) ( X )
#define ELEM_DECODE(X) ( X )

#define SCALE_LOG(X) log10( X )
#define SCALE_NORM(X) ( X )

#endif


#define ELEM_NAN std::numeric_limits<double>::quiet_NaN()
#define ISNAN(x) ( x!=x )

#include <cassert>
#include <cmath>
#include <ctime>

/* for signal handling */
#include <csignal>

/* for limits */
#include <limits>

/* for typeinfo */
#include <typeinfo>

/* STL includes */
#include <iostream>
#include <iomanip>
#include <fstream>
#include <map>
#include <list>
#include <vector>
#include <stack>
#include <set>
#include <queue>
#include <string>
#include <utility>
#include <limits>
#include <sstream>
#include <algorithm>
#include <functional>
#include <locale>

#include <boost/foreach.hpp>
#include <boost/scoped_ptr.hpp>
#include <boost/scoped_array.hpp>
using boost::scoped_ptr;
using boost::scoped_array;


/* shorthand for convenience */
namespace daoopt {
typedef std::ostringstream oss;
}  // namespace daoopt

/* which hashtable to use? define only *one*  */
//#define HASH_CPP11
//#define HASH_BOOST
//#define HASH_TR1
//#define HASH_SGI
//#define HASH_GOOGLE_DENSE
//#define HASH_GOOGLE_SPARSE

/* type for storing contexts in binary */
namespace daoopt {
typedef std::vector<val_t> context_t;
}  // namespace daoopt

// google utilties
#include "base/integral_types.h"
/*
#include "base/logging.h"
#include "strings/split.h"
#include "strings/strcat.h"
#include "util/gtl/stl_util.h"
*/

#ifdef HASH_CPP11
#include <unordered_set>
#include <unordered_map>

namespace daoopt {
template<typename K>
using hash_set = std::unordered_set<K>;
template<typename K, typename V>
using hash_map = std::unordered_map<K, V>;
}  // namespace daoopt
#endif


#ifdef HASH_SGI
/* SGI hash set and map */
#include <backward/hash_set> // deprecated!
#include <backward/hash_map> // deprecated!

namespace daoopt {
using __gnu_cxx::hash_set;
using __gnu_cxx::hash_map;
}  // namespace daoopt
#endif

#ifdef HASH_TR1
/* TR1 hash sets and maps */
#include <tr1/unordered_set>
#include <tr1/unordered_map>

/* some renaming (crude hack) */
namespace daoopt {
template<typename K>
using hash_set = std::tr1::unordered_set<K>;
template<typename K, typename V>
using hash_map = std::tr1::unordered_map<K, V>;
}  // namespace daoopt
#endif

#ifdef HASH_BOOST
#include "boost/unordered_map.hpp"
#include "boost/unordered_set.hpp"

namespace daoopt {
template<typename K>
using hash_set = boost::unordered_set<K>;
template<typename K, typename V>
using hash_map = boost::unordered_map<K, V>;
}  // namespace daoopt
#endif

#ifdef HASH_GOOGLE_DENSE
#include "google/dense_hash_map"
#include "google/dense_hash_set"

namespace daoopt {
template<typename K>
using hash_set = google::dense_hash_set<K>;
template<typename K, typename V>
using hash_map = google::dense_hash_map<K, V>;
}  // namespace daoopt
#endif

#ifdef HASH_GOOGLE_SPARSE
#include "google/sparse_hash_map"
#include "google/sparse_hash_set"

namespace daoopt {
template<typename K>
using hash_set = google::sparse_hash_set<K>;
template<typename K, typename V>
using hash_map = google::sparse_hash_map<K, V>;
}  // namespace daoopt
#endif

#if defined HASH_SGI || defined HASH_GOOGLE_DENSE || defined HASH_GOOGLE_SPARSE
/* BEGIN FIX FOR STRINGS AND SGI HASH TABLE in __gnu_cxx */
namespace __gnu_cxx {
template<> struct hash<std::string> {
  size_t operator()(const std::string& x) const {
    return hash<const char*> ()(x.c_str());
  }
};
}  // namespace __gnu_cxx
/* END FIX */
#endif

/* might not be needed anymore? */
#ifdef __APPLE__
#include <sys/types.h>
#endif

/* for debugging */
#ifdef DEBUG
//#include "debug.h"
#endif

/* type for counting nodes, fixed precision for 32/64 bit machines. */
#include <stdint.h>
typedef uint64_t count_t;

/* LibGMP C++ interface */
#if defined PARALLEL_DYNAMIC || defined PARALLEL_STATIC
 #ifdef USE_GMP
 #include <gmpxx.h>
namespace daoopt {
typedef mpz_class bigint;
typedef mpf_class bigfloat;
typedef mpq_class bigfrac;
}  // namespace daoopt
 #else
namespace daoopt {
typedef unsigned long int bigint;
typedef double bigfloat;
typedef double bigfrac;
}  // namespace daoopt
 #endif
#endif

/* Boost random number library */
#include <boost/random/linear_congruential.hpp>
/* Boost lexical cast (for version string) */
#include <boost/lexical_cast.hpp>

namespace daoopt {
using namespace std;
}  // namespace daoopt

/*////////////////////////////////*/
/*////// MACRO DEFINITIONS ///////*/
/*////////////////////////////////*/

#define UNKNOWN -1
#define NOID -1
#define NONE -1

#define TRUE 1
#define FALSE 0

/* FUNCTION TYPES */
#define TYPE_BAYES 1
#define TYPE_WCSP 2

/* INPUT FORMATS */
#define INP_UAI 1
#define INP_ERGO 2

/* TASK TYPE */
#define TASK_MIN 1
#define TASK_MAX 2

/* PROBLEM TYPE */
#define PROB_MULT 1
#define PROB_ADD 2

/* SEARCH NODE TYPE */
#define NODE_AND 1
#define NODE_OR 2

namespace daoopt {
const int SUBPROB_WIDTH_INC = 0;
const int SUBPROB_WIDTH_DEC = 1;
const int SUBPROB_HEUR_INC = 2;
const int SUBPROB_HEUR_DEC = 3;
const string subprob_order[4]
  = {"width-inc","width-dec","heur-inc","heur-dec"};
}  // namespace daoopt

/*//////////////////////////////////////////////////////////////*/

// undefine min/max macros
#if defined(WINDOWS)
#ifdef min
#undef min
#endif
#ifdef max
#undef max
#endif
#endif

namespace daoopt {
/* static random number generator */
class rand {
private:
  static int state;
  static boost::minstd_rand _r;
public:

  static void seed(const int& s) { _r.seed(s); }
  static int seed() { return state; }
  static int max() { return _r.max(); }
  static int next() { return state=_r(); }
  static int next(const int& hi) {
    return static_cast<int>( (state=_r()) / (_r.max()+1.0) * hi );
  }
  // probability from a uniform distribution.
  static double next_unif() {
    return static_cast<double>( (state=_r()) / (_r.max()+1.0) );
  }

};
}  // namespace daoopt

/*//////////////////////////////*/

#if FALSE
/*
 * the following is based on code and explanations from
 * http://www.cygnus-software.com/papers/comparingfloats/Comparing%20floating%20point%20numbers.htm
 */

/* NOTE: Taking out the void* casts produces the gcc warning
 * "dereferencing type-punned pointer will break strict-aliasing rules"
 */

namespace daoopt {

/* floating point equality comparison (modulo floating point precision) */
inline bool fpEq(double A, double B, int64_t maxDist=2) {
  assert(sizeof(double) == sizeof(int64_t));
  if (A == B)
    return true;
  int64_t aInt = *(int64_t*)(void*)&A;
  if (aInt<0) aInt = 0x8000000000000000LL - aInt;
  int64_t bInt = *(int64_t*)(void*)&B;
  if (bInt<0) bInt = 0x8000000000000000LL - bInt;
  int64_t intDiff = std::abs(aInt - bInt);
  if (intDiff <= maxDist)
    return true;
  return false;
}

/* floating point "less than" (modulo fp. precision) */
inline bool fpLt(double A, double B, int64_t maxDist=2) {
  assert(sizeof(double) == sizeof(int64_t));
  if (A == B || A == INFINITY || B == -INFINITY)
    return false;
  if (A == -INFINITY || B == INFINITY)
    return true;
  int64_t aInt = *(int64_t*)(void*)&A;
  if (aInt<0) aInt = 0x8000000000000000LL - aInt;
  int64_t bInt = *(int64_t*)(void*)&B;
  if (bInt<0) bInt = 0x8000000000000000LL - bInt;
  int64_t intDiff = bInt - aInt;
  if (intDiff > maxDist)
    return true;
  return false;
}

/* floating point "less or equal than" (modulo fp. precision) */
inline bool fpLEq(double A, double B, int64_t maxDist=2) {
  assert(sizeof(double) == sizeof(int64_t));
  if (A == B || A == -INFINITY || B == INFINITY)
    return true;
  if (A == INFINITY || B == -INFINITY)
    return false;
  int64_t aInt = *(int64_t*)(void*)&A;
  if (aInt<0) aInt = 0x8000000000000000LL - aInt;
  int64_t bInt = *(int64_t*)(void*)&B;
  if (bInt<0) bInt = 0x8000000000000000LL - bInt;
  int64_t intDiff = aInt - bInt;
  if (intDiff > maxDist)
    return false;
  return true;
}

}  // namespace daoopt

#endif
/*//////////////////////////////////////////////////////////////*/

/* binary read and write of value X to/from stream S */
#ifdef BINARY_DATAFILES
#define BINWRITE(S,X) ( S ).write((char*)&( X ), sizeof( X ))
#define BINREAD(S,X) ( S ).read((char*)&( X ), sizeof( X ))
#define BINSKIP(S,T,N) ( S ).ignore( sizeof( T ) * ( N ) )
#else
#define BINWRITE(S,X) ( S ) << ( X ) << ' '
#define BINREAD(S,X)  ( S ) >> ( X )
#define BINSKIP(S,T,N) {T x; for (size_t _i=N; _i; --_i) (S)>>(x); }
#endif


#if FALSE
namespace daoopt {
/* encode doubles to 64 bit integers (and back) */
typedef int64_t int64bit;
inline std::string encodeDoubleAsInt(double d) {
  int64bit x = *(int64bit*)(void*)&d;
  std::ostringstream ss;
  ss << x;
  return ss.str();
}

inline double decodeDoubleFromString(std::string s) {
  int64bit x = 0;
  std::istringstream ss(s);
  ss >> x;
  return *(double*)(void*)&x;
}

}  // namespace daoopt

#endif /* false */

namespace daoopt {

inline std::string encodeDoubleAsInt(double d) {
  std::ostringstream ss;
  ss << d;
  return ss.str();
}

inline double decodeDoubleFromString(std::string s) {
  double x = -((double)INFINITY);
  std::istringstream ss(s);
  ss >> x;
  return x;
}

}  // namespace daoopt


/*///////////////////////////////////////////////////////////*/

namespace daoopt {
inline double mylog10(unsigned long a) {
  return log10(a);
}

#if defined PARALLEL_DYNAMIC || defined PARALLEL_STATIC
 #ifdef USE_GMP
 double mylog10(bigint a);
 #endif
#endif

}  // namespace daoopt


/*///////////////////////////////////////////////////////////*/

#if FALSE
/* hash function for pair<int,int> */
namespace __gnu_cxx {

  using std::size_t;

  template<> struct hash< std::pair<int,int> > {
    size_t operator () (const std::pair<int,int>& p) const {
      size_t h = 0;
      h ^= size_t(p.first) + 0x9e3779b9 + (h<<6) + (h>>2) ;
      h ^= size_t(p.second) + 0x9e3779b9 + (h<<6) + (h>>2) ;
      return h;
    }
  };

}  // namespace __gnu_cxx

#endif /* false */

#endif /* _BASE_H_ */



