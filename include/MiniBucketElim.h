/*
 * MiniBucket.h
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
 *  along with DAOOPT.  If not, see http://www.gnu.org/licenses/>.
 *  
 *  Created on: Nov 8, 2008
 *      Author: Lars Otten <lotten@ics.uci.edu>
 */

#ifndef MINIBUCKETELIM_H_
#define MINIBUCKETELIM_H_

#define m_augmented m_miniBucketFunctions.top().m_augmentedF
#define m_intermediate m_miniBucketFunctions.top().m_intermediateF

#include "Heuristic.h"
#include "Function.h"
#include "Problem.h"
#include "ProgramOptions.h"
#include "Pseudotree.h"
#include "utils.h"

#include "MiniBucket.h"


/* The overall minibucket elimination */
class MiniBucketFunctions;

class MiniBucketElim : public Heuristic {

  friend class MiniBucket;

protected:
  int m_ibound;                  // The ibound for this MB instance
  double m_globalUB;             // The global upper bound

/*
  // The augmented buckets that will store the minibucket functions (but not the original ones)
  vector<vector<Function*> > m_augmented;
  // Precompute and store, for each variable v, the relevant intermediate functions that are
  // generated in a pseudotree descendant and passed to an ancestor of v
  // (points to the same function objects as m_augmented)
  vector<vector<Function*> > m_intermediate;
  */

  stack<MiniBucketFunctions> m_miniBucketFunctions;

  vector<vector<int> > m_augmentedSource;
  vector<vector<int> > m_intermediateSource;

  bool m_momentMatching;
  bool m_dynamic;

protected:
  // Computes a dfs order of the pseudo tree, for building the bucket structure
  void findDfsOrder(vector<int>&) const;
  void findDfsOrder(vector<int>&, int var) const;

  // Compares the size of the scope of two functions
//  bool scopeIsLarger(Function*, Function*) const;

  // reset the data structures
  void reset() ;

public:

  // checks if the given i-bound would exceed the memlimit and lowers
  // it accordingly.
  size_t limitSize(size_t memlimit, const vector<val_t> * assignment);

  // builds the heuristic, limited to the relevant subproblem, if applicable.
  // if computeTables=false, only returns size estimate (no tables computed)
  size_t build(const vector<val_t>* assignment = NULL, bool computeTables = true);

  // builds the heuristic, restricted to the subtree rooted by the current assignment
  size_t buildSubproblem(int var, const map<int,val_t> &assignment, const vector<val_t> &vAssn, const vector<int> &elimOrder, bool computeTables = true);

  // returns the global upper bound
  double getGlobalUB() const { return m_globalUB; }

  // computes the heuristic for variable var given a (partial) assignment
  double getHeur(int var, const vector<val_t>& assignment);
  // computes heuristic values for all instantiations of var, given context assignment
  void getHeurAll(int var, const vector<val_t>& assignment, vector<double>& out);

  // reset the i-bound
  void setIbound(int ibound) { m_ibound = ibound; }
  // gets the i-bound
  int getIbound() const { return m_ibound; }

  // gets sum of tables sizes
  size_t getSize() const;

  // gets the width of the subproblem rooted at node i
  int getWidthSubproblem(int i) const;

  bool writeToFile(string fn) const;
  bool readFromFile(string fn);

  bool isAccurate();

public:
  MiniBucketElim(Problem* p, Pseudotree* pt, ProgramOptions* po, int ib);
  virtual ~MiniBucketElim();

};

// Class to store minibucket functions used for heuristics, along with its associated assignment
class MiniBucketFunctions {
  friend class MiniBucketElim;

  map<int,val_t> m_assignment;

  // Keep track of which variables are eliminated
  vector<int> m_elimOrder;
  // The augmented buckets that will store the minibucket functions (but not the original ones)
  vector<vector<Function*> > m_augmentedF;
  // Precompute and store, for each variable v, the relevant intermediate functions that are
  // generated in a pseudotree descendant and passed to an ancestor of v
  // (points to the same function objects as m_augmented)
  vector<vector<Function*> > m_intermediateF;


  // value is based on whether any partitioning was performed during function computation
  bool isAccurate;

public:
  MiniBucketFunctions() : isAccurate(false) {}
  MiniBucketFunctions(const map<int,val_t> &assignment, const vector<int> elimOrder) : m_assignment(assignment), m_elimOrder(elimOrder), isAccurate(false) {}
  MiniBucketFunctions(const vector<int> elimOrder) : m_elimOrder(elimOrder), isAccurate(false) {}
  ~MiniBucketFunctions() {
    for (vector<vector<Function*> >::iterator it=m_augmentedF.begin() ;it!=m_augmentedF.end(); ++it)
      for (vector<Function*>::iterator it2=it->begin(); it2!=it->end(); ++it2)
          delete (*it2);
  }

  const map<int,val_t> &getAssignment() const { return m_assignment; }

  // Check to see if the conditioning of these functions are compatible with the assignment 
  // to be evaluated. Used to see if the stack should be popped.
  bool isCompatible(const map<int,val_t> &assignment, const vector<int> &elimOrder) const {
      for (map<int,val_t>::const_iterator it=m_assignment.begin(); it!=m_assignment.end(); ++it)
          if (it->second != assignment.find(it->first)->second) {
              return false;
          }

      // m_elimOrder must be a superset of elimOrder
      unsigned int i=0,j=0;
      while(j<elimOrder.size()) {
          while(m_elimOrder[i] != elimOrder[j])
              if(++i >= m_elimOrder.size()) {
//                  cout << "var " << elimOrder[j] << " not found" << endl;
                  return false;
              }
          ++j;
      }
      return true;
  }

  void printAssignAndElim() const {
      for (map<int,val_t>::const_iterator it=m_assignment.begin(); it!=m_assignment.end(); ++it)
          cout << " " << it->first << " "<< int(it->second) << endl;
      cout << endl << "m_elimOrder: " << endl;
      for (unsigned int i=0; i<m_elimOrder.size(); ++i)
          cout << " " << m_elimOrder[i];
      cout << endl;
  }

  bool isEmpty() const {
      return m_augmentedF.size()==0 && m_intermediateF.size()==0;
  }
};

/* Inline definitions */

inline bool MiniBucketElim::isAccurate() {
  assert(m_pseudotree);
  return (m_pseudotree->getWidthCond() <= m_ibound);
}

inline MiniBucketElim::MiniBucketElim(Problem* p, Pseudotree* pt,
				      ProgramOptions* po, int ib) :
    Heuristic(p, pt, po), m_ibound(ib), m_globalUB(ELEM_ONE), m_momentMatching(po->match), m_dynamic(po->dynamic)
// , m_augmented(p->getN()), m_intermediate(p->getN())
  { }

inline MiniBucketElim::~MiniBucketElim() {
  // make sure to delete each function only once
  while(m_miniBucketFunctions.size()) m_miniBucketFunctions.pop();
}

inline bool scopeIsLarger(Function* p, Function* q) {
  assert(p && q);
  if (p->getArity() == q->getArity())
    return (p->getId() > q->getId());
  else
    return (p->getArity() > q->getArity());
}


#endif /* MINIBUCKETELIM_H_ */
