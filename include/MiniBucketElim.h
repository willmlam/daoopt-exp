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


/* The overall minibucket elimination */
//class MiniBucketFunctions;
class ConditionedMessages;

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

  vector<vector<int> > m_mbCount;
  vector<int> m_mbCountAccurate;

  vector<int> m_mbCurrentHeur;

  // Sets of messages that change based on the conditioning
  vector<set<Function* > > m_augmented;
  vector<set<Function* > > m_intermediate;

  // Keeps track of whether outgoing messages are accurate
  vector<bool> m_accurateHeuristic;

  // Keeps track of whether incoming messages are accurate
  vector<bool> m_accurateHeuristicIn;

  bool m_momentMatching;
  bool m_dynamic;

  int m_gNodes;                  // compute the heuristic every g nodes
  int m_currentGIter;            // Counter for managing granularity

  int m_dhDepth;
  int m_depthInterval;
  int m_maxDupe;
  int m_dupeImp;

  int m_maxDynHeur;

  int m_buildSubCalled;

  vector<vector<int> > m_elimOrder;

protected:
  // Computes a dfs order of the pseudo tree, for building the bucket structure
  void findDfsOrder(vector<int>&) const;
  void findDfsOrder(vector<int>&, int var) const;

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
      return m_mbCount[u][v] - m_mbCountAccurate[v];
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
  size_t buildSubproblem(int var, const map<int,val_t> &assignment, const vector<val_t> &vAssn, const vector<int> &elimOrder, bool computeTables = true);

  // simulates building the subproblem heuristic, returning the number of buckets
  int simulateBuildSubproblem(int var, const vector<int> &elimOrder);

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

// Class to store the outgoing messages from a bucket
class ConditionedMessages {
    // The messages that are sent from this bucket
    vector<Function*> m_functions;

    // Contains the path that the corresponding message 
    vector<vector<int> *> m_paths;
    map<int,val_t> m_assignment;
    bool m_isAccurate;
    public:
    ConditionedMessages(const map<int,val_t> &assignment, bool isAccurate) 
        : m_assignment(assignment), m_isAccurate(isAccurate) {
    }

    const vector<Function*> &getFunctions() const {
        return m_functions;
    }

    const vector<vector<int> *> &getPaths() const {
        return m_paths;
    }

    // Return the bucket that message i goes into
    const int getTargetVar(int i) const {
        return m_paths[i]->back();
    }

    const map<int,val_t> &getAssignment() const {
        return m_assignment;
    }

    bool isAccurate() const {
        return m_isAccurate;
    }

    bool isConsistent(const map<int,val_t> &assignment) {
        map<int,val_t>::iterator it = m_assignment.begin();
        for(; it != m_assignment.end(); ++it) {
            const map<int,val_t>::const_iterator fit = assignment.find(it->first);
            if (fit == assignment.end() || fit->second != it->second) 
                return false;
        }
        return true;
    }

    void addFunction(Function *f, vector<int> *path) {
        m_functions.push_back(f);
        m_paths.push_back(path);
    }

    ~ConditionedMessages() {
        for (unsigned i = 0; i < m_functions.size(); ++i) {
            delete m_functions[i];
            delete m_paths[i];
        }
    }
};

class Scope {
    int m_id;
    set<int> m_scope;
public:
    Scope(int id, const set<int> &scope) : m_id(id), m_scope(scope) {
    };
    int getId() const { return m_id; }
    set<int> &getScope() { return m_scope; }
    int getArity() const { return m_scope.size(); }
    void erase(int i) {
        m_scope.erase(i);
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
    m_mbCount(p->getN(),vector<int>(p->getN(), 0)),
    m_mbCountAccurate(p->getN()),
    m_mbCurrentHeur(p->getN()),
    m_augmented(p->getN()), 
    m_intermediate(p->getN()), 
    m_accurateHeuristic(p->getN(), true),
    m_accurateHeuristicIn(p->getN(), true),
    m_momentMatching(po->match),
    m_dynamic(po->dynamic), 
    m_gNodes(po->gNodes), 
    m_currentGIter(0), 
    m_dhDepth(po->dhDepth),
    m_depthInterval(po->depthInterval),
    m_maxDupe(po->maxDupe),
    m_dupeImp(po->dupeImp),
    m_maxDynHeur(po->maxDynHeur),
    m_buildSubCalled(0)
  { 
      // If dynamic, precomupute all DFS elimination orders for each node
      // and precompute number of minibuckets used in each subproblem rooted by each node
      if (m_dynamic) {
          m_elimOrder.resize(p->getN());
          for (int i = 0 ; i < p->getN(); ++i) {
              findDfsOrder(m_elimOrder[i], i);
              m_mbCountAccurate[i] = m_elimOrder[i].size() - (i==p->getN()-1 ? 1 : 2);
              simulateBuildSubproblem(i, m_elimOrder[i]);
          }
      }
      else {
          m_elimOrder.resize(1);
          findDfsOrder(m_elimOrder[0]);
      }
  }

inline MiniBucketElim::~MiniBucketElim() {
  // make sure to delete each function only once
  for (int i = 0; i < m_problem->getN(); ++i) {
      while (!m_cMessages[i].empty()) {
          delete m_cMessages[i].top();
          m_cMessages[i].pop();
      }
  }
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

#endif /* MINIBUCKETELIM_H_ */
