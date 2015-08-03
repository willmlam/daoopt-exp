/////////////////////////////////////////////////////////////////////////////////////
// FactorSparse.h  --  class definition for matlab-compatible "sparse" factor class
//
// A few functions are defined only for MEX calls (construction & load from matlab)
// Most others can be used more generally.
// 
//////////////////////////////////////////////////////////////////////////////////////
//
// Written by Alex Ihler
// Copyright (C) 2014-15 Alexander Ihler; distributable under GPL -- see README.txt
//
//////////////////////////////////////////////////////////////////////////////////////
#ifndef __FACTOR_SPARSE_H
#define __FACTOR_SPARSE_H

#include <cmath>
#include <cstring>
#include <iostream>
#include <limits>
#include <algorithm>
#include <numeric>
#include <float.h>
#include <vector>
#include <map>

#include "mxUtil.h"
#include "enum.h"
#include "VarSet.h"
#include "subindex.h"
#include "vector.h"

// MEMORY define: keep track of global Factor memory usage (helpful in debugging other code)
//#define __FACTOR_H_MEMORY
#ifdef __FACTOR_H_MEMORY
  #define FACTOR_ADD_MEMORY(val)  { memused += val; mmax = std::max(mmax,memused); }
  #define FACTOR_SUB_MEMORY(val)  { memused -= val; }
#else
  #define FACTOR_ADD_MEMORY(val)
  #define FACTOR_SUB_MEMORY(val)
#endif

// !!! add bigUInt class
// !!! TODO: Resolve vindex vs size_t differences?


namespace mex {

class Factor : public virtual mxObject {
 public:
  typedef size_t BigInt;
  typedef double value;
  typedef std::map<BigInt, value> _tableType;
  typedef std::map<BigInt, value>::iterator _tableIter;
  typedef std::map<BigInt, value>::const_iterator _tableCIter;
  typedef VarSet::vindex vindex;   // variable identifiers (1...N)
  typedef VarSet::vsize  vsize;    // variable values (0...K-1)

  // Constructors //////////////////////////////////////////////////////////////////////////
  Factor(Factor const& f);                                 // copy ctor
  explicit Factor(const value s=1.0);                      // scalar constructor
  explicit Factor(VarSet const& vs, value s=1.0);          // constant factor over given vars
  //Factor(VarSet const& vs, value* T, bool littleEndian=true); // constructor from array of values
  Factor(VarSet const& vs, value* T);                      // constructor from array of values
  ~Factor();                                               // destructor

  // Assignments & copy constructors //////////////////////////////////////////////////////
  Factor& operator=(Factor const& rhs);                    // assignment (deep copy)
  void swap(Factor& F);                                    // object contents exchange
  virtual Factor* clone() const;                           // clone from pointer

  #ifdef __FACTOR_H_MEMORY
  static size_t memused;                     // helper f'ns / vars when memory tracking
  static size_t mmax;                        //  is enabled
  static size_t mem() { return memused; };
  #endif

  // MEX Class Wrapper Functions //////////////////////////////////////////////////////////
#ifdef MEX 
   bool        mxCheckValid(const mxArray*);   // check if matlab object is compatible with this object
   void        mxSet(mxArray*);     // associate with A by reference to data
   mxArray*    mxGet();             // get a pointer to the matlab object wrapper (creating if required)
   void        mxRelease();         // disassociate with a matlab object wrapper, if we have one
   void        mxDestroy();         // disassociate and delete matlab object
   void        mxSwap(Factor&);     // disassociate and delete matlab object
#endif
  void setDims();   // !!! protected?        // make variable & table dimensions consistent (matlab only)

  /////////////////////////////////////////////////////////////////////////////////////////
  // Accessor Functions ///////////////////////////////////////////////////////////////////
  /////////////////////////////////////////////////////////////////////////////////////////
  const size_t     nvar()      const;        // # of variables
  const VarSet&    vars()      const;        // variable IDs
  // Non-const version removed; unsafe? !!
  const VarSet&    variables() const;        // name difference (!!)
  const vsize*     dims()      const;        // variable dimensions
  const value*     table()     const;        // table of factor values
  size_t           nrStates()  const;        // get number of possible states
  size_t           numel()     const;        // get factor size

  // Boolean checks on factor properties //////////////////////////////////////////////////
  bool isempty()  const;                     // is factor empty (no values)?
  bool isnan()    const;                     // does factor contain NaNs?
  bool isfinite() const;                     // does factor contain inf/-inf?
  bool isscalar() const;                     // is factor a scalar?

  void cleanInf();                   // !!! hack?

  // Direct value accessor ////////////////////////////////////////////////////////////////
  value  operator[] (vsize v) const;                        // access table elements
  value& operator[] (vsize v)      ;                        // rewrite table elements

  value  get(vsize i)         const;                        // "safe" versions (check size, etc)
  void   set(vsize i, value v)     ;

  // Filling //////////////////////////////////////////////////////////////////////////////
  Factor& fill(value v);                                    // constant: F(x)=v for all x
  Factor& randomize();                                      // fill with (uniform) random values                
  Factor& fillUniform() { return fill(1.0/nrStates()); }    // uniform: fill with 1/N
  static Factor delta(const VarSet& vs, size_t idx); // dirac delta f'n: F(idx)=1, F(~idx)=0

  // Unary transformations (in-place); see outside class def'n for const versions /////////
  inline Factor& abs(void)   { transform(unOpAbs()); return *this; };       // F(x) = |F(x)|
  inline Factor& exp(void)   { transform(unOpExp()); return *this; };       // F(x) = exp[ F(x) ]
  inline Factor& log(void)   { transform(unOpLog()); return *this; };      // F(x) = log[ F(x) ]
  inline Factor& log0(void)  { transform(unOpLog0()); return *this; };       // F(x) = log[ F(x) ], or 0 if F(x)=0
  inline Factor& log2(void)  { transform(unOpLogL(std::log(2))); return *this; }; //= log2[ F(x) ]
  inline Factor& log10(void) { transform(unOpLog10()); return *this; };     // F(x) = log10[ F(x) ]
  // Normalize? (!!)

////////////////////////////////////////////////////////////////////////////////
// Basic factor operations (+,-,*,/)
////////////////////////////////////////////////////////////////////////////////
  Factor  operator+ (const Factor& B) const;    // add two factors, F+B
  Factor  operator- (const Factor& B) const;    // difference of two factors, F-B
  Factor  operator* (const Factor& B) const;    // product, F*B
  Factor  operator/ (const Factor& B) const;    // ratio, F/B
  Factor& operator+=(const Factor& B)      ;    // in-place versions of same (rewrite F)
  Factor& operator-=(const Factor& B)      ;
  Factor& operator*=(const Factor& B)      ;
  Factor& operator/=(const Factor& B)      ;

  Factor  operator+ (const value B) const  ;    // add factor + scalar
  Factor  operator- (const value B) const  ;    //  and so on
  Factor  operator* (const value B) const  ;
  Factor  operator/ (const value B) const  ;
  Factor  operator^ (const value B) const  ;
  Factor& operator+=(const value B)        ;    // in-place versions of same (rewrite F)
  Factor& operator-=(const value B)        ;
  Factor& operator*=(const value B)        ;
  Factor& operator/=(const value B)        ;
  Factor& operator^=(const value B)        ;

////////////////////////////////////////////////////////////////////////////////
// Comparative thresholding operators (useful? !!!)
////////////////////////////////////////////////////////////////////////////////
// operator> (const Factor& F)
// operator>=(const Factor& F)
// operator< (const Factor& F)
// operator<=(const Factor& F)

////////////////////////////////////////////////////////////////////////////////
// Partition function, entropy, and normalization
////////////////////////////////////////////////////////////////////////////////
  Factor& normalize();                           // normalize F (in-place) to sum to 1.0
  Factor  normalized()  const;                   // return a normalized version of F (copy)
  double logpartition() const;                   // compute log of the normalization constant
  value entropy(void) const ;                    // compute entropy (in nats)

////////////////////////////////////////////////////////////////////////////////
// Elimination operators (sum, max, min, ...)
////////////////////////////////////////////////////////////////////////////////
  Factor sum(VarSet const& sumOut) const;         // compute sum over a set of variables "sumOut"
  value  sum()                     const;         // sum over all variables
  Factor sumPower(VarSet const& sumOut,value pow)  const;  // power sum operation

  Factor logsumexp(const VarSet& sumOut) const;   // log-sum-exp (sum for log-factors)
  double logsumexp() const;
  // !!! lsePower?

  Factor max(VarSet const& sumOut) const;         // maximize over a set of variables "sumOut"
  value  max()                     const;

  Factor min(VarSet const& sumOut) const;         // minimize over a set of variables "sumOut"
  value  min()                     const;

  size_t argmax() const;                          // find the (or a) maximizer of F(x)
  size_t argmin() const;                          //   ""            minimizer of F(x)

  size_t argmax(const VarSet& vCond, vsize vState) const;    // argmax, conditioned on vCond=vState
  // !!! TODO: argmin( ... )

  Factor condition(const VarSet& vRem, vsize vState) const;  // factor with vRem=vState clamped
  Factor slice(const VarSet& vRem, vsize vState) const;      // "condition" by libDAI's name
  Factor embed(const VarSet& v) const;                       // expand scope of factor to "v"

  size_t sample() const;                                     // random draw proportional to F(x)
  vsize  draw() const;                                       // "sample" by libDAI's name

  Factor marginal(VarSet const& target) const;               // find the marginal of F over target
  Factor maxmarginal(VarSet const& target) const;            //   or the max-marginal or min-marginal
  Factor minmarginal(VarSet const& target) const;

  void marginalInto(VarSet const& target, Factor& F) const;  // In-place marginalization (TODO)

////////////////////////////////////////////////////////////////////////////////
// Misc other functions
////////////////////////////////////////////////////////////////////////////////

  MEX_ENUM( Distance , L1,L2,LInf,KL,HPM,MAS,OptGap );

  double distance(Factor const& F2, Distance type=Distance::L2) const;
  double norm(Distance type=Distance::L2) const;

  MEX_ENUM( Decomp , L2,L2_HPM,L2_MAS );

  std::vector<Factor> decompSum(std::vector<VarSet> vlist, Factor::Decomp method) const;
  std::vector<Factor> decompProd(std::vector<VarSet> vlist, Factor::Decomp method) const;

  ////////////////////////////////////////////////////////////////////////////////
  // File Input/Output formats
  ////////////////////////////////////////////////////////////////////////////////

  static vector<Factor> readUai10(std::istream& is);
  static void           writeUai10(std::ostream& os, const vector<Factor>&);
  static vector<Factor> readErgo(std::istream& is);
  static void           writeErgo(std::ostream& os, const vector<Factor>&);
  static vector<Factor> readWCSP(std::istream& is);
  static void           writeWCSP(std::ostream& os, const vector<Factor>&);

  friend std::ostream& operator<<( std::ostream& out, const mex::Factor& F);

  //////////////////////////////////////////////////////////////////////////////
  // Private object functions
  //////////////////////////////////////////////////////////////////////////////
protected:
  VarSet v_;      // variable list vector
  _tableType t_;  // lookup from (tuple)->(value)

#ifdef MEX 
  static inline bool isfinite(value v) { return mxIsFinite(v); }
  static inline bool isnan(value v)    { return mxIsNaN(v); }
  static inline value infty()          { return mxGetInf(); }
#else
  static inline bool isfinite(double v) { return (v <= DBL_MAX && v >= -DBL_MAX); }
  static inline bool isnan(value v)    { return (v!=v); }
  static inline value infty()          { return std::numeric_limits<value>::infinity(); }
#endif

  /////////////////////////////////////////////////////////////////////////////
  // Private members for implementation of various operations                //
  /////////////////////////////////////////////////////////////////////////////
  template <class xform>
  inline void transform(xform xf) { for (_tableIter i=t_.begin();i!=t_.end();++i) i->second = xf(i->second); };

  // Functors defined for unary operations (transformations) to the factor table
  struct unOpAbs   { value operator()(value a)  { return std::abs(a);  } };
  struct unOpExp   { value operator()(value a)  { return std::exp(a);  } };
  struct unOpLog0  { value operator()(value a)  { return (a!=0.0) ? std::log(a) : 0.0;  } };
  struct unOpLog   { value operator()(value a)  { return std::log(a);  } };
  struct unOpLogL  { value l; unOpLogL(value L):l(L) {}; value operator()(value a)  { return std::log(a)/l;  } };
  struct unOpLog10 { value operator()(value a)  { return std::log10(a);} };

  struct myMin { value operator()(value a,value b)  { return std::min(a,b);} };
  struct myMax { value operator()(value a,value b)  { return std::max(a,b);} };

  template<typename Function> Factor binaryOp( const Factor& B, Function Op) const;
  template<typename Function> Factor binaryOp( const value B, Function Op) const;
  template<typename Function> Factor& binaryOpIP( const Factor& B, Function Op);
  template<typename Function> Factor& binaryOpIP( const value B, Function Op);

  static size_t unionVals(const VarSet& v1, size_t r1, const VarSet& v2, size_t r2, const VarSet& vUnion);
  static size_t isectVals(const VarSet& v1, const VarSet& v12, size_t r1);
};

////////////////////////////////////////////////////////////////////////////////
// "Static" functions that operate on Factor class variables
////////////////////////////////////////////////////////////////////////////////

inline Factor  operator+ (const Factor::value B, const Factor& A) { return A+B; };
inline Factor  operator* (const Factor::value B, const Factor& A) { return A*B; };
inline Factor  operator- (const Factor::value B, const Factor& A) { Factor BF(A.vars(),B); return BF-=A; };
inline Factor  operator/ (const Factor::value B, const Factor& A) { Factor BF(A.vars(),B); return BF/=A; };

inline Factor abs(const Factor& A)   { Factor F=A; F.abs();   return F; };
inline Factor exp(const Factor& A)   { Factor F=A; F.exp();   return F; };
inline Factor log0(const Factor& A)  { Factor F=A; F.log0();  return F; };
inline Factor log(const Factor& A)   { Factor F=A; F.log();   return F; };
inline Factor log2(const Factor& A)  { Factor F=A; F.log2();  return F; };
inline Factor log10(const Factor& A) { Factor F=A; F.log10(); return F; };

template <class InputIterator>
inline Factor mean(InputIterator first, InputIterator last) { size_t N=0; Factor F(0.0); for (;first!=last;++first,++N) F+=*first; if (N) F/=N; return F; }
template <class InputIterator>
inline Factor geomean(InputIterator first, InputIterator last) { size_t N=0; Factor F(1.0); for (;first!=last;++first,++N) F*=*first; if (N) F^=1.0/N; return F; }


////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
// IMPLEMENTATION                                                             //
////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////



inline Factor::Factor(Factor const& f) : v_(f.v_), t_(f.t_) {                   // copy ctor
  FACTOR_ADD_MEMORY( t_.size()*sizeof(double) );
}
inline Factor::Factor(const value s) : v_(),   t_()  {             // scalar constructor
  t_[0]=s;
  FACTOR_ADD_MEMORY( t_.size()*sizeof(double) );
}
inline Factor::Factor(VarSet const& vs, value s) : v_(vs),t_() {      // constant factor over given vars
  if (s) for (size_t i=0, end=v_.nrStates();i<end;++i) t_[i]=s;     // TODO: sparse???
  FACTOR_ADD_MEMORY( t_.size()*sizeof(double) );
}
inline Factor::Factor(VarSet const& vs, value* T) : v_(vs),t_() {               // constructor from array of values
  if (T) for (size_t i=0, end=v_.nrStates();i<end;++i) t_[i]=T[i];          // TODO: sparse???
  FACTOR_ADD_MEMORY( t_.size()*sizeof(double) );
}
//TODO: Memory tracking: may grow (shrink?) with operations?
// Factor( from Vars and a vector of values)  !!
// Factor( from permuted vector of Var and vector of values) !!

inline Factor::~Factor() {                                                      // destructor
  FACTOR_SUB_MEMORY( t_.capacity()*sizeof(double) );
}

// Assignments & copy constructors //////////////////////////////////////////////////////

inline Factor& Factor::operator=(Factor const& rhs) {                   // assignment (deep copy)
  if (this!=&rhs) {
    FACTOR_SUB_MEMORY( t_.capacity()*sizeof(double) );
    #ifndef MEX
    { _tableType tmp; t_.swap(tmp); }                // force vector to release memory
    #endif
    v_ = rhs.v_; t_ = rhs.t_; setDims();                 // then reassign
    FACTOR_ADD_MEMORY( t_.capacity()*sizeof(double) );
  }
  return *this;
}

inline void Factor::swap(Factor& F) {                                  // object contents exchange
  if (&F!=this) { v_.swap(F.v_); t_.swap(F.t_); setDims(); F.setDims(); }
}

inline Factor* Factor::clone() const { Factor *F=new Factor(*this); return F; }  // clone from pointer

// MEX Class Wrapper Functions //////////////////////////////////////////////////////////
inline void Factor::setDims() {  // make variable & table dimensions consistent (matlab only)
  #ifdef MEX
  if (t_.mxAvail())
    if (v_.size()>1) mxSetDimensions(t_.mxGet(),v_.dims(),v_.size());
    else { mxSetM(t_.mxGet(),v_.size() ? v_.dims()[0] : 1); mxSetN(t_.mxGet(),1); }
  #endif
}

/////////////////////////////////////////////////////////////////////////////////////////
// Accessor Functions ///////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////
inline const size_t          Factor::nvar()      const { return v_.nvar();  };// # of variables
inline const VarSet&         Factor::vars()      const { return v_;   };      // variable IDs
// Non-const version removed; unsafe? !!
inline const VarSet&         Factor::variables() const { return v_;   };      // name difference (!!)
inline const Factor::vsize*  Factor::dims()      const { return v_.dims();};  // variable dimensions
inline const Factor::value*  Factor::table()     const { return NULL;   };  // TODO: throw exception
inline size_t                Factor::nrStates()  const { return v_.nrStates(); }  // # of states
inline size_t                Factor::numel()     const { return t_.size(); }; // # of values : name diff (!!)
  // TODO: differentiate between nrStates and # used / effective states... (nnz?)


// Boolean checks on factor properties //////////////////////////////////////////////////
inline bool Factor::isempty()  const { return t_.empty(); };             // empty factor
inline bool Factor::isnan()    const { bool b=false; for (_tableCIter i=t_.begin();i!=t_.end();++i) b|=isnan(i->second);    return b; };
inline bool Factor::isfinite() const { bool b=true;  for (_tableCIter i=t_.begin();i!=t_.end();++i) b&=isfinite(i->second); return b; };
inline bool Factor::isscalar() const { return (v_.nvar()==0 && !isempty()); };  // TODO: verify / do something else?

inline void Factor::cleanInf() { for (size_t i=0;i<nrStates();i++) if (!isfinite(t_[i])) t_[i]=0.0; }; // TODO:iterator

// Direct value accessor ////////////////////////////////////////////////////////////////
inline Factor::value  Factor::operator[] (vsize v) const { return t_.at(v); };   // access table elements
inline Factor::value& Factor::operator[] (vsize v)       { return t_[v]; };      // rewrite table elements
inline Factor::value  Factor::get(vsize i)         const { return t_.at(i); }    // "safe" versions
inline void           Factor::set(vsize i, value v)      { t_.at(i)=v; }

////////////////////////////////////////////////////////////////////////////////
// Inline definitions of some trivial functions
////////////////////////////////////////////////////////////////////////////////

// Filling
inline Factor& Factor::fill(value v) { for (_tableIter i=t_.begin();i!=t_.end();++i) i->second=v; return *this; } // TODO: "" ?
inline Factor& Factor::randomize()   { for (_tableIter i=t_.begin();i!=t_.end();++i) i->second=mex::randu(); return *this; } // TODO: "" ?

inline Factor  Factor::normalized()  const { Factor F=*this; F.normalize(); return F;  }
inline Factor& Factor::normalize() { double Z=sum(); if (Z!=0) *this/=Z; return *this; }

inline double Factor::logpartition() const { return std::log(sum()); };   // compute log of the normalization constant

inline Factor         Factor::sum(VarSet const& sumOut) const { VarSet t=v_ - sumOut; return marginal(t); }
inline Factor::value  Factor::sum()                     const { value v=0.0; for (_tableCIter i=t_.begin();i!=t_.end();++i) v+=i->second; return v; }

inline Factor         Factor::max(VarSet const& sumOut) const { VarSet t=v_ - sumOut; return maxmarginal(t); }
inline Factor::value  Factor::max()                     const { value v=-infty(); for (_tableCIter i=t_.begin();i!=t_.end();++i) v=(v>i->second)?v:i->second; return v; }

inline Factor         Factor::min(VarSet const& sumOut) const { VarSet t= v_ - sumOut; return minmarginal(t); }
inline Factor::value  Factor::min()                     const { value v=infty(); for (_tableCIter i=t_.begin();i!=t_.end();++i) v=(v<i->second)?v:i->second; return v; }

inline size_t Factor::argmax() const { return std::max_element(t_.begin(),t_.end())->first; }
inline size_t Factor::argmin() const { return std::min_element(t_.begin(),t_.end())->first; }

inline Factor Factor::slice(const VarSet& vRem, Factor::vsize vState) const { return condition(vRem,vState); } // !! name change

inline Factor::vsize Factor::draw() const { return sample(); }  // !! name difference

inline Factor Factor::embed(const VarSet& v) const { if (vars()==v) return *this; else return (*this)+Factor(v/vars(),0.0);}

} // end namespace mex

#undef FACTOR_ADD_MEMORY
#undef FACTOR_SUB_MEMORY

#endif
