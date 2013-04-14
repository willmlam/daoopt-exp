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

#include "Heuristic.h"
#include "Function.h"
#include "Problem.h"
#include "ProgramOptions.h"
#include "Pseudotree.h"
#include "utils.h"

#include "MiniBucket.h"
#include "MBEHeuristicInstance.h"


/* The overall minibucket elimination */
//class MiniBucketFunctions;
class Scope;

enum DomRel { DOMINATES = 10, NEG_DOMINATES = -10, EQUAL = 1, INDETERMINATE = 0 };

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

  vector<stack<ConditionedMessages*> > m_cMessages;

  vector<vector<int> > m_mbCountSubtree; // cumulative on entire subtree
  vector<vector<int> > m_dupeCount; // number of mini buckets for each bucket

  // for each entry in this matrix M[i][j] =
  // 1 if i dominates j
  // -1 if j dominates i
  // 0 if i and j are equal
  // -10 if i and j are not comparable
  // condition: M[i][j] = -M[j][i] unless they are not comparable
  vector<vector<DomRel> > m_heuristicDominates;

  vector<int> m_mbCountAccurate;

  int m_currentGIter;            // Counter for managing granularity

  int m_maxDynHeur;

  int m_numHeuristics;
  
  int m_memlimit;

  vector<vector<int> > m_elimOrder;

  // Stores the root instance of the heuristic
  // also used by static setup
  MBEHeuristicInstance* m_rootHeurInstance;

  // Stores all conditioned heuristics
  // (Mostly a mechanism to properly free memory)
  vector<MBEHeuristicInstance*> m_heurCollection;

protected:
  // Computes a dfs order of the pseudo tree, for building the bucket structure
  void findDfsOrder(vector<int>&) const;
  void findDfsOrder(vector<int>&, int var) const;

  // Computes a bfs order of the pseudo tree
  void findBfsOrder(vector<int>&) const;

  // Compares the size of the scope of two functions
//  bool scopeIsLarger(Function*, Function*) const;

  // Remove messages from the augmented and intermediate sets that are stored by 
  // the given ConditionedMessages object
  void eraseMessages(ConditionedMessages *cm);

  // Inserts messages from the augmented and intermediate sets that are stored by
  // the given ConditionedMessages object
  void insertMessages(ConditionedMessages *cm, int bVar, vector<bool> &visited);

  // Calculate the number of extra variables induced by the minibucket heuristic 
  // computed at at variable u for evaluating a variable v
  int numberOfDuplicateVariables(int u, int v) {
      return m_mbCountSubtree[u][v] - m_mbCountAccurate[v];
  }

  // Calculate the number of variables not in the given context
  int getConditionedArity(const set<int> &scope, const set<int> &context) {
      int s = 0;
      for (typeof(scope.begin()) it = scope.begin(); it != scope.end(); ++it) {
          if (context.count(*it)==0) ++s;
      }
      return s;
  }

  bool meetsComputeConditions(int var, int varAncestor, int depth) {
      return m_numHeuristics < m_maxDynHeur &&
          m_options->dhDepth > depth &&
          depth % m_options->depthInterval == 0 &&
          (m_options->gNodes > 0 && m_currentGIter == 0) &&
          (numberOfDuplicateVariables(varAncestor,var) -
           numberOfDuplicateVariables(var,var)) >= m_options->dupeRed &&
          (m_options->strictDupeRed <= 0 || m_heuristicDominates[var][varAncestor] == DOMINATES) &&
          rand::next() < int(m_options->randDyn * rand::max());
  }

  bool meetsReuseCondition(int var, int varAncestor, int bucket) {
      switch(m_options->reuseLevel) {
          case 1:
              //return m_dupeCount[varAncestor][bucket] == m_dupeCount[var][bucket];
          case 2:
              return m_dupeCount[varAncestor][bucket] == 1;
          default:
              return false;
      }
  }

  void buildDominanceMatrix() {
      vector<int> ordering; 
      findBfsOrder(ordering);
      vector<int>::iterator it = ordering.begin();
      for (; it != ordering.end(); ++it) {
          PseudotreeNode *n = m_pseudotree->getNode(*it);
          m_heuristicDominates[*it][*it] = EQUAL;
          vector<bool> processed(ordering.size(),false);
          while ( (n = n->getParent()) ) {
              int v2 = n->getVar();
              // skip since variable was already processed via transitivity
              if (processed[v2]) continue;
                  
              // check dupe count on current variable only
              vector<int>::iterator vit = m_elimOrder[*it].begin();
              bool allEqualOrLess = true;
              int countStrictlyLess = 0;
              for (; vit != m_elimOrder[*it].end(); ++vit) {
                  if (*vit == m_pseudotree->getRoot()->getVar()) continue;
                  else if (m_dupeCount[*it][*vit] > m_dupeCount[v2][*vit]) {
                      allEqualOrLess = false;
                      break;
                  }
                  if (m_dupeCount[*it][*vit] < m_dupeCount[v2][*vit]) countStrictlyLess++;
              }
              if (allEqualOrLess) {
                  PseudotreeNode *n2 = m_pseudotree->getNode(v2);
                  if (countStrictlyLess >= m_options->strictDupeRed) {
                      m_heuristicDominates[*it][v2] = DOMINATES;
                      m_heuristicDominates[v2][*it] = NEG_DOMINATES;
                      while ( (n2 = n2->getParent()) ) {
                          int v3 = n2->getVar();
                          // if the relation is defined, it is transitively dominating
                          if (m_heuristicDominates[v2][v3] != INDETERMINATE) {
                              m_heuristicDominates[*it][v3] = DOMINATES;
                              m_heuristicDominates[v3][*it] = NEG_DOMINATES;
                              processed[v3] = true;
                          }
                      }
                  }
                  else {
                      m_heuristicDominates[*it][v2] = m_heuristicDominates[v2][*it] = EQUAL;
                      while ( (n2 = n2->getParent()) ) {
                          int v3 = n2->getVar();
                          // if the relation is defined, it is transitively the same relation
                          if (m_heuristicDominates[v2][v3] != INDETERMINATE) {
                              m_heuristicDominates[*it][v3] = m_heuristicDominates[v2][v3];
                              m_heuristicDominates[v3][*it] = m_heuristicDominates[v3][v2];
                              processed[v3] = true;
                          }
                      }
                  }
              }
          }
          processed[*it] = true;
      }
  }

  // reset the data structures
  void reset();

public:

  // checks if the given i-bound would exceed the memlimit and lowers
  // it accordingly.
  size_t limitSize(size_t memlimit, const vector<val_t> * assignment);

  // builds the heuristic, limited to the relevant subproblem, if applicable.
  // if computeTables=false, only returns size estimate (no tables computed)
  size_t build(const vector<val_t>* assignment = NULL, bool computeTables = true);

  // builds the heuristic, restricted to the subtree rooted by the current assignment
  size_t buildSubproblem(int var, const vector<val_t> &vAssn, MBEHeuristicInstance *ancHeur, MBEHeuristicInstance *curHeur, bool computeTables = true);

  // simulates building the subproblem heuristic, returning the number of buckets
  void simulateBuildSubproblem(int var, const vector<int> &elimOrder);

  // returns the global upper bound
  double getGlobalUB() const { return m_globalUB; }

  // computes the heuristic for variable var given a (partial) assignment
  double getHeur(int var, const vector<val_t>& assignment, SearchNode *n);
  // computes heuristic values for all instantiations of var, given context assignment
  void getHeurAll(int var, const vector<val_t>& assignment, SearchNode *n, vector<double>& out);


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

  int getNumHeuristics() const {
      return m_numHeuristics;
  }

  int getCurrentNumActive() const {
      return MBEHeuristicInstance::getCurrentNumActive();
  }

  int getMaxNumActive() const {
      return MBEHeuristicInstance::getMaxNumActive();
  }

public:
  MiniBucketElim(Problem* p, Pseudotree* pt, ProgramOptions* po, int ib);
  virtual ~MiniBucketElim();

};


class Scope {
    int m_id;
    set<int> m_scope;
public:
    Scope(int id, const set<int> &scope) : m_id(id) {
        m_scope.insert(scope.begin(),scope.end());
    };
    int getId() const { return m_id; }
    set<int> &getScope() { return m_scope; }
    int getArity() const { return m_scope.size(); }
    void erase(int i) {
        m_scope.erase(i);
    }
    ~Scope() {
        m_scope.clear();
    }
};


/*
// Class to store minibucket functions used for heuristics, along with its associated assignment
class MiniBucketFunctions {
  friend class MiniBucketElim;

  // Variable for heuristic
  int m_var;

  // Conditioning used in computing heuristic
  map<int,val_t> m_assignment;

  // Keep track of which variables are eliminated
  //vector<int> m_elimOrder;
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
  MiniBucketFunctions(int var) : m_var(var), isAccurate(false) {}
  MiniBucketFunctions(int var, const map<int,val_t> &assignment) : m_var(var), m_assignment(assignment), isAccurate(false) {}
  ~MiniBucketFunctions() {
    for (vector<vector<Function*> >::iterator it=m_augmentedF.begin() ;it!=m_augmentedF.end(); ++it)
      for (vector<Function*>::iterator it2=it->begin(); it2!=it->end(); ++it2)
          delete (*it2);
  }

  const map<int,val_t> &getAssignment() const { return m_assignment; }

  // Check to see if the conditioning of these functions are compatible with the assignment 
  // to be evaluated. Used to see if the stack should be popped.
  bool isCompatible(int var, const map<int,val_t> &assignment, Pseudotree *pt) const {
      // Always compatible if there is no associated assignment
      if (m_assignment.empty()) return true;

      // Check if m_var is an ancestor of var
      PseudotreeNode *n = pt->getNode(var);
      while (n->getVar() != m_var) {
          if (n != pt->getRoot())
              n = n->getParent();
          else
              return false;
      }
      map<int,val_t>::const_iterator it = m_assignment.begin();
      for (; it!=m_assignment.end() && it->second != assignment.find(it->first)->second; ++it);
      return it==m_assignment.end();
  }

  void printVarAndAssign() const {
      cout << "var: " << m_var << endl;
      for (map<int,val_t>::const_iterator it=m_assignment.begin(); it!=m_assignment.end(); ++it)
          cout << " " << it->first << " "<< int(it->second) << endl;
  }

  bool isEmpty() const {
      return m_augmentedF.size()==0 && m_intermediateF.size()==0;
  }
};
*/



/* Inline definitions */

inline bool MiniBucketElim::isAccurate() {
  assert(m_pseudotree);
  return (m_pseudotree->getWidthCond() <= m_ibound);
}

inline MiniBucketElim::MiniBucketElim(Problem* p, Pseudotree* pt,
				      ProgramOptions* po, int ib) :
    Heuristic(p, pt, po), m_ibound(ib), m_globalUB(ELEM_ONE), 
    m_cMessages(p->getN()), 
    m_mbCountSubtree(p->getN(),vector<int>(p->getN(), 0)),
    m_dupeCount(p->getN(),vector<int>(p->getN(), 0)),
    m_heuristicDominates(p->getN(),vector<DomRel>(p->getN(),INDETERMINATE)),
    m_mbCountAccurate(p->getN()),
    m_currentGIter(0), 
    m_maxDynHeur(po->maxDynHeur),
    m_numHeuristics(0)
  { 
      m_rootHeurInstance = new MBEHeuristicInstance(p->getN(), pt->getRoot()->getVar(), NULL);
      
      // If dynamic, precomupute all DFS elimination orders for each node
      // and precompute number of minibuckets used in each subproblem 
      // rooted by each node
      if (m_options->dynamic) {
          m_elimOrder.resize(p->getN());
          for (int i = 0 ; i < p->getN(); ++i) {
              findDfsOrder(m_elimOrder[i], i);
              m_mbCountAccurate[i] = m_elimOrder[i].size() - (i==p->getN()-1 ? 1 : 2);
              simulateBuildSubproblem(i, m_elimOrder[i]);
          }
          if (m_options->strictDupeRed > 0) buildDominanceMatrix();
      }
      else {
          m_elimOrder.resize(1);
          findDfsOrder(m_elimOrder[0]);
      }
  }

inline MiniBucketElim::~MiniBucketElim() {
  // make sure to delete each function only once
  //delete m_rootHeurInstance;
  /*
  for (unsigned i = 0; i < m_heurCollection.size(); ++i) {
      delete m_heurCollection[i];
  }
  */
}

inline bool scopeIsLarger(Function* p, Function* q) {
  assert(p && q);
  if (p->getArity() == q->getArity())
    return (p->getId() > q->getId());
  else
    return (p->getArity() > q->getArity());
}

inline bool scopeIsLargerS(Scope* p, Scope* q) {
  assert(p && q);
  if (p->getArity() == q->getArity())
    return (p->getId() > q->getId());
  else
    return (p->getArity() > q->getArity());
}

inline bool scopeIsLargerIF(const pair<int,Function*> &l, const pair<int,Function*> &r) {
  Function *p = l.second;
  Function *q = r.second;
  assert(p && q);
  if (l.first == r.first)
    return (p->getId() > q->getId());
  else
    return (l.first > r.first);
}

#endif /* MINIBUCKETELIM_H_ */
