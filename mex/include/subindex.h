#ifndef __MEX_SUBINDEX_H
#define __MEX_SUBINDEX_H

#include <iostream>

namespace mex {

/*****************************************************************************
Some basics on table-based factor representations:
Each entry in a table corresponds implicitly to some configuration, or tuple,
of the argument variables.  The order in which these tuples occur is critical.
Assuming the factor is F(x1,x2,x3,x4) and arguments are in "sorted" order,
we can imagine two standard orders we could use:
  Little Endian (Matlab)                   Big Endian (UAI)
  0000                                     0000
  1000                                     0001
  0100                                     0010
  1100                                     0011
  ...                                      ...
In other words, in LE format, the span of the first variable (x1) is small,
meaning that increments to x1's value correspond to small changes in the
accessed memory location, while the span of x4 is large (incrementing x4
corresponds to a much larger memory location).

We usually store our tables in Little Endian format, for Matlab compatibility.

*****************************************************************************
Basic operations:

The "subindex" class allows one to walk through a function with a large scope,
say F(x1,x2,x3), and use a subindex to simultaneously step through the 
corresponding elements of a function with smaller scope, e.g. G(x1,x3).

 Ex:  sub = [x1,x3]    =>  [0 0] [0 1] [0 0] [0 1] [1 0] [1 1] [1 0] [1 1]
     full = [x1,x2,x3] =>  [000] [001] [010] [011] [100] [101] [110] [111]

The "superindex" class allows one to walk through an implicitly specified
sub-function, or a "slice" of F:  F(x1, x2=X, x3)

The "permuteIndex" class allows walking through a function whose arguments
are specified in a different (unsorted) order, and/or has a different 
Endian-ness.

*****************************************************************************/



// Subindex : an iterator class for stepping through the entries of a small function "sub"
//   in the correct order to match linearly stepping through the entries of a larger function "full"
// Ex:  sub = [x1,x3]    =>  [0 0] [0 1] [0 0] [0 1] [1 0] [1 1] [1 0] [1 1]
//     full = [x1,x2,x3] =>  [000] [001] [010] [011] [100] [101] [110] [111]

// Notes:
// This is a "new" faster version that does not keep track of all variables' states, just
//   "runs" of variables that are jointly skipped or jointly not skipped
// The "sub-ptr" version uses pointers rather than indirect indexing, but does not appear to be faster.

class subindex {
  public:
   typedef VarSet::vsize vsize;

   vsize  idx_,end_;
   vsize  Nd;
   vsize  *state;           // vector of variable-indices (values) for the current position (1-based)
   vsize *dims;             // dimensions of each variable  (!!! delete in destructor)
   bool *skipped;           // variables not included in the sub-index
   vsize *add;              // how much to add when we increment
   vsize *subtract;         // how much to subtract when we wrap each variable
   
//
   subindex(const VarSet& full, const VarSet& sub) {
     assert( full >> sub );
     idx_=0; 
     end_=1;
     Nd=full.nvar(); 
     //state = new vsize[4*Nd]; dims=state+Nd; add=dims+Nd; subtract=add+Nd; 
     dims     = new vsize[Nd]; 
     state    = new vsize[Nd];
     add      = new vsize[Nd];
     subtract = new vsize[Nd];
     skipped  = new bool[Nd];
     // Compute reference index updates
     vsize i,j,k; 
     for (i=0;i<Nd;++i) dims[i]=full.dims()[i];        // copy off dimensions before modifying
     for (i=0,j=0,k=0;i<Nd;++i) {                      // i=index in full; j=index in sub; k=index in local vectors
       state[k]=1;                                     // initialize to state 1
       skipped[k]=( j>=sub.nvar() || sub[j]!=full[i] ); // does the subset include this variable?
       end_ *= dims[i];                                // track how large the full index count is
       if (i!=0 && skipped[k]==skipped[k-1]) dims[--k]*=dims[i]; // continuing to skip (or not) => combine into one
       else dims[k]=dims[i];                                     //  "variable"; otherwise transition to next
       if (k==0) add[k]=1;                             // how much does adding one to this var
       else add[k]=add[k-1]*(skipped[k-1]?1:dims[k-1]); //    add to our position?
       subtract[k]=add[k]*((skipped[k]?1:dims[k])-1);  // how much does wrapping back to 1 remove from pos?
       if (!skipped[k]) ++j;                           // var was found in sub => move to next var in sub
       ++k;
     }
     Nd = k;                                               // ended up with k "variables" after combining
     for (i=0;i<Nd;++i) if (skipped[i]) subtract[i]=add[i]=0; // allows unconditional add/subtract in ++
   }

   subindex(const subindex& S) {
     dims = S.dims; idx_=S.idx_; end_=S.end_;
     Nd = S.Nd;
     //state = new vsize[4*Nd]; dims=state+Nd; add=dims+Nd; subtract=add+Nd; 
     dims     = new vsize[Nd]; std::copy(S.dims,S.dims+Nd,dims);
     state    = new vsize[Nd]; std::copy(S.state,S.state+Nd,state);
     add      = new vsize[Nd]; std::copy(S.add,S.add+Nd,add);
     subtract = new vsize[Nd]; std::copy(S.subtract,S.subtract+Nd,subtract);
     skipped  = new bool[Nd];  std::copy(S.skipped,S.skipped+Nd,skipped);
   }

   ~subindex(void) {
     delete[] dims; delete[] state; delete[] skipped; delete[] add; delete[] subtract;
   }

   subindex& reset() {
     for (size_t i=0;i<Nd;++i) state[i]=1;
     idx_=0;
     return *this;
   };

   size_t end() const { return end_; }

#ifdef __SUB_PTR
   /// Prefix addition operator  (is using pointers directly faster?)
   subindex& operator++ (void) {
     vsize *State=state,*Dims=dims,*Add=add,*Sub=subtract,*End=state+Nd;
     bool *Skip=skipped;
     while (State!=End) {                        // for each variable                   
       if (*State==*Dims) {                      // if we reached the maximum, wrap around to 1
         *State=1;                               //   subtract wrap value from position
         //if (!*Skip) idx_ -= *Sub;               //   and continue to next variable in sequence  
         idx_ -= *Sub;                           //   and continue to next variable in sequence  
         ++State; ++Dims; ++Add; ++Sub; ++Skip;  
       } else {                                  // otherwise, increment variable index 
         ++(*State);                             //   add to our current position
         //if (!*Skip) idx_ += *Add;               //   and break (leave later vars the same)  
         idx_ += *Add;                           //   and break (leave later vars the same)  
         break;                                     
       }
     }
     return *this;
   }
#else 
   /// Prefix addition operator
   subindex& operator++ (void) {
     for (size_t i=0;i<Nd;++i) {                 // for each variable                   
       if (state[i]==dims[i]) {                  // if we reached the maximum, wrap around to 1
         state[i]=1;                             //   subtract wrap value from position
         //if (!skipped[i]) idx_ -= subtract[i]; //   and continue to next variable in sequence
         idx_ -= subtract[i];                    //   and continue to next variable in sequence
       } else {                                  // otherwise, increment variable index
         ++state[i];                             //   add to our current position
				 //if (!skipped[i]) idx_ += add[i];      //   and break (leave later vars the same)
         idx_ += add[i];                         //   and break (leave later vars the same)
         break;
       }
     }
     return *this;
   }
#endif
  /// Postfix addition : use prefix and copy
  subindex operator++ (int) { subindex S(*this); ++(*this); return S; }

  /// Conversion to index value
  operator size_t() const { return idx_; };

  /// Non-unit addition TODO:test
  subindex& operator+= (size_t C) {
    for (size_t i=0;i<Nd;++i) {
      size_t R = (state[i]+C)/dims[i];           // # rollovers: times we'll increment next var
      idx_ += (C-R)*add[i] - R*subtract[i];      // incrementing this var by C results in this
      state[i] = (state[i]-1+C)%dims[i]+1;       //   many adds & subtracts, and this state
      if (R==0) break;                           // if we incremented the next var, move on 
      C = R;                                     //   to that, with R increments on it
    }
    return *this;
  }

};


// Superindex : an iterator class for stepping through the entries of a large function "full"
//   in the correct order to match linearly stepping through the entries of a small function "sub"
// Ex: offset=[0 1 0] = 2
//      sub = [x1,x3]    =>  [0 0] [0 1] [1 0] [1 1] 
//     full = [x1,x2,x3] =>  [010] [011] [110] [111] 
// !!! TODO: shouldn't we pass the variables we want to slice on, rather than the ones to step over?
class superindex {
public:
  typedef VarSet::vsize vsize;

  vsize  idx_,end_;
  vsize  Ns;
  vsize  *state;           // vector of variable-indices (values) for the current position (1-based)
  const vsize *dims;       // dimensions of each variable
  vsize *add;              // how much to add when we increment
   

  superindex(const VarSet& full, const VarSet& sub, size_t offSub) {
    assert( full >> sub );
		size_t offset=0;
    idx_=offset; end_=full.nrStates();	// initially, span the whole table
    Ns = sub.nvar();
    const vsize* dimf = full.dims(); dims=sub.dims();
    state=new vsize[Ns]; 
    add  =new vsize[Ns]; 
    size_t i,j,d=1;
    for (i=0,j=0;i<full.nvar();++i) {
      if (j<sub.nvar() ? full[i]==sub[j] : false) {						// on variables in "sub", figure out how
        state[j]=1;											//  to span "full" with each step
        add[j]=d;
        j++;
      } else {													// on variables outside "sub", convert "offSub"
				size_t rem=(offSub % dimf[i]);  //  to a position in full table: 
				offSub -= rem; offSub/=dimf[i];	//  get the value of var vi in offSub
				offset += d*rem;                //  & locate that value in full var set
      }
      d*=dimf[i];
    }
		idx_ = offset;
    end_ = add[Ns-1]*dims[Ns-1]+offset;	// todo: fix (this is last+1, need last++)
  }

  superindex(const superindex& S) {
    dims = S.dims; idx_=S.idx_; end_=S.end_;
    Ns = S.Ns;
    state    = new vsize[Ns]; std::copy(S.state,S.state+Ns,state);
    add      = new vsize[Ns]; std::copy(S.add,S.add+Ns,add);
  }

  ~superindex(void) {
    delete[] state; delete[] add;
  }

   superindex& reset() {
     for (size_t i=0;i<Ns;++i) state[i]=1;
     idx_=0;
     return *this;
   };

   size_t end() { return end_; }

   /// Prefix addition operator
   superindex& operator++ (void) {
   for (size_t i=0;i<Ns;++i) {                  // for each variable                   
     if (state[i]==dims[i] && i<Ns-1) {         // if we reached the maximum, wrap around to 1
       state[i]=1;                              // 
       idx_ -= add[i]*(dims[i]-1);
     } else {                                   // otherwise, increment variable index
       ++state[i];                              //   add to our current position
       idx_+= add[i];
       break;
     }
   }
   return *this;
  }
  /// Postfix addition : use prefix and copy
  superindex operator++ (int) { superindex S(*this); ++(*this); return S; }

  /// Conversion to index value
  operator size_t() const { return idx_; };

};


class permuteIndex {
  public:
    // Construct permutation mapping from VarSet -> Order  (bigEndian: is 1st variable largest stride?)
    permuteIndex( const vector<Var>& order, bool bigEndian=false ) {
      _i = 0;
      VarSet _vs = VarSet(order.begin(),order.end(),order.size());   // compute implicit source order (VarSet)
      _pi.resize(order.size()); _dim.resize(order.size());
      for (size_t j=0;j<order.size();++j) _dim[j]=_vs[j].states();   // save dimensions in source order (VarSet)
      for (size_t j=0;j<order.size();++j) {                          // compute mapping from target order to source
        size_t jj = bigEndian ? order.size()-1-j : j;
        for (size_t k=0;k<order.size();++k)
          if (_vs[k] == order[j]) { _pi[jj]=k; break; }
      }
    }

    // get target index corresponding to current or specified source index
    operator      size_t()       { return convert(_i); }
    permuteIndex& set(size_t i)  { _i=i; return *this; };

    // convert a source index into a target index
    size_t convert(size_t i) {
      vector<size_t> I(_dim.size());
      size_t r=0, m=1;
      for (size_t v=0; v<_dim.size(); ++v) { I[v] = i%_dim[v]; i-=I[v]; i/=_dim[v]; }
      for (size_t j=0; j<_dim.size(); ++j) { r += m*I[_pi[j]]; m*=_dim[_pi[j]]; }
      return r;
    }

    // invert mapping from Order -> VarSet
    permuteIndex inverse() {
      permuteIndex inv(*this);
      for (size_t i=0;i<_pi.size();++i) { inv._pi[_pi[i]]=i; inv._dim[i]=_dim[_pi[i]]; }
      inv._i  = (size_t) *this;
      return inv;
    }

    // can be used as an iterator as well
    permuteIndex& operator++ (void) { ++_i; return *this; }
    permuteIndex  operator++ (int)  { permuteIndex r(*this); ++_i; return r; }
    permuteIndex& operator-- (void) { --_i; return *this; }
    permuteIndex  operator-- (int)  { permuteIndex r(*this); --_i; return r; }
  private:
    size_t _i;
    vector<size_t> _pi;
    vector<size_t> _dim;
};

} // end namespace mex
#endif
