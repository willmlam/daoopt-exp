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

#include "mex/mbe.h"
#include "FGLP.h"


/* The overall minibucket elimination */
//class MiniBucketFunctions;
class Scope;

enum class HeurRelation {EQUAL, BETTER, WORSE, _UNKNOWN};

ostream& operator <<(ostream &os, const HeurRelation &val);

/*

class ComputeCondition {
    friend MiniBucketElim;
    MiniBucketElim *m_mbe;
    bool met(int var, int varAncestor, int depth) = 0;
    ComputeCondition() : m_mbe(NULL) {}
    ComputeCondition(MiniBucketElim *mbe) : m_mbe(mbe) {}
};

class MaxPathCondition : public ComputeCondition {
    bool met(int var, int varAncestor, int depth) {
        return MBEHeuristicInstance::getCurrentNumActive() < m_options->maxPathHeur;
    }
};

class ExactFrontierCondition : public ComputeCondition {
};

class DepthCondition : public ComputeCondition {
};

class DepthIntervalCondition : public ComputeCondition {
};

class EveryNodesCondition : public ComputeCondtion {
};

class DuplicateVarsCondition : public ComputeCondition {
};

class StrictDuplicateVarsCondition : public ComputeCondition {
};

class RandComputeCondition : public ComputeCondition {
};
*/


class MiniBucketElim : public Heuristic {

  friend class MiniBucket;

protected:
  int m_ibound;                  // The ibound for this MB instance
  double m_globalUB;             // The global upper bound

  vector<double> tempOut;        // Used for temporary storage of heuristic values


/*
  // The augmented buckets that will store the minibucket functions (but not the original ones)
  vector<vector<Function*> > m_augmented;
  // Precompute and store, for each variable v, the relevant intermediate functions that are
  // generated in a pseudotree descendant and passed to an ancestor of v
  // (points to the same function objects as m_augmented)
  vector<vector<Function*> > m_intermediate;
  */

  vector<vector<int> > m_mbCountSubtree; // cumulative on entire subtree
  vector<vector<int> > m_dupeCount; // number of mini buckets for each bucket

  vector<int> m_mbCountAccurate;

  vector<int> m_subproblemWidth; // widths for each subproblem

  int m_currentGIter;            // Counter for managing granularity

  int m_maxDynHeur;

  int m_numHeuristics;
  
  int m_memlimit;

  vector<vector<int> > m_elimOrder;

  // Stores the root instance of the heuristic
  // also used by static setup
  MBEHeuristicInstance* m_rootHeurInstance;

  // Stores the root instance of FGLP if used
  FGLP *m_fglpRoot;

  // Stores all conditioned heuristics
  // (Mostly a mechanism to properly free memory)
  vector<MBEHeuristicInstance*> m_heurCollection;

  // Stores all number of times a variable is visited
  vector<unsigned int> m_countVarVisited;

  // Statistics: counts if the heuristic at this variable was better
  vector<int> m_heurBetter;
  
  // Estimation of the heuristic at this variable is expected to be
  // at [nodevar] : EQUAL, BETTER, WORST, _UNKNOWN
  vector<vector<HeurRelation>> m_heurBetterPre;

  // Stores the partitionings used for minibucket
  // [nodevar][var][partitionIdx] function set
  vector<vector<vector<set<int>>>> m_varPartitioning;

  vector<map<int,set<int>>> m_varPartitioningScopes;

  // Mechanism to make minibucket function ids unique
  int m_idMultiplier;

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
      for (auto it = scope.begin(); it != scope.end(); ++it) {
          if (context.count(*it)==0) ++s;
      }
      return s;
  }

  bool meetsComputeConditions(int var, int varAncestor, int depth) {
      /*
      return m_numHeuristics < m_options->maxDynHeur &&
          (m_options->maxPathHeur == -1 || 
           (m_options->maxPathHeur > 0 && MBEHeuristicInstance::getCurrentNumActive() < m_options->maxPathHeur) ||
           (m_options->computeExactFrontier && m_subproblemWidth[var] == m_options->ibound)) &&
          (depth > 0 && m_options->dhDepth >= depth) &&
          depth % m_options->depthInterval == 0 &&
          (m_options->gNodes > 0 && m_currentGIter == 0) &&
          (numberOfDuplicateVariables(varAncestor,var) -
           numberOfDuplicateVariables(var,var)) >= m_options->dupeRed &&
          (m_options->randDyn >= 1.0 || rand::next() < int(m_options->randDyn * rand::max()));
          */
      /*
      cout << "Subproblem width of " << var << " : " << m_subproblemWidth[var] << endl;
      cout << "Subibound : " << m_options->subibound << endl;
      cout << "Subibound distance : " << m_options->subwidthDistance << endl;
      */
      return (m_options->maxPathHeur == -1 || MBEHeuristicInstance::getCurrentNumActive() < m_options->maxPathHeur) && 
        (m_options->subwidthDistance == numeric_limits<int>::max() || 
         (m_subproblemWidth[var] - m_options->subibound <= m_options->subwidthDistance &&
         m_subproblemWidth[var] < m_subproblemWidth[m_pseudotree->getNode(var)->getParent()->getVar()]));
  }

  // Higher level does more reuse
  bool meetsReuseCondition(int var, int varAncestor, int bucket,
          MBEHeuristicInstance *curHeur) {
      bool conditionMet = false;
      switch(m_options->reuseLevel) {
          case 2: {
              // Check which are the new conditioning variables
              vector<int> vars;
              PseudotreeNode *n = m_pseudotree->getNode(var);
              do {
                  n = n->getParent();
                  vars.push_back(n->getVar());
              } while (n->getVar() != varAncestor);

              bool found = false;
              // Check if any exist in the original functions or incoming messages
              for (unsigned int i = 0; i < vars.size() && !found; ++i) {
                  const vector<Function*>& fnlist = m_pseudotree->getFunctions(bucket);
                  vector<Function*>::const_iterator itF = fnlist.begin();
                  for (; itF != fnlist.end() && !found; ++itF) {
                      if ((*itF)->hasInScope(vars[i]))
                          found = true;
                  }
                  const vector<Function*>& auglist = curHeur->getAugmented()[bucket];
                  itF = auglist.begin();
                  for (; itF != auglist.end() && !found; ++itF) {
                      if ((*itF)->hasInScope(vars[i]))
                          found = true;
                  }
              }
              conditionMet = conditionMet || !found;

          }
            
          case 1:
              conditionMet = conditionMet || (m_dupeCount[varAncestor][bucket] == 1);
          default:
              return conditionMet;
      }
  }

  // reset the data structures
  void reset();

public:

  // checks if the given i-bound would exceed the memlimit and lowers
  // it accordingly.
  size_t limitSize(size_t memlimit, const vector<val_t> * assignment);

  // checks if the given i-bound would exceed the memlimit and lowers
  // it accordingly. (for JGLP)
  size_t limitJGLPIBound(size_t memlimit, const vector<val_t> * assignment);

  // builds the heuristic, limited to the relevant subproblem, if applicable.
  // if computeTables=false, only returns size estimate (no tables computed)
  size_t build(const vector<val_t>* assignment = NULL, bool computeTables = true);

  // builds the heuristic, restricted to the subtree rooted by the current assignment
  size_t buildSubproblem(int var, const vector<val_t> &vAssn, MBEHeuristicInstance *ancHeur, MBEHeuristicInstance *curHeur, bool computeTables = true);

  // simulates building the subproblem heuristic, returning the number of buckets
  void simulateBuildSubproblem(int var, const vector<int> &elimOrder, int ibound);

  // returns the global upper bound
  double getGlobalUB() const { return m_globalUB; }

  // computes the heuristic for variable var given a (partial) assignment
  double getHeur(int var, const vector<val_t>& assignment, SearchNode *n);
  // computes heuristic values for all instantiations of var, given context assignment
  void getHeurAll(int var, const vector<val_t>& assignment, SearchNode *n, vector<double>& out);

  double getLabel(int var, const vector<val_t>& assignment, SearchNode *n);
  void getLabelAll(int var, const vector<val_t>& assignment, SearchNode *n, vector<double> &out);

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

  void printExtraStats() const {
      cout << "Heuristic stats" << endl;
      cout << "---------------" << endl;
      cout << "# Heuristics:    " << m_numHeuristics << endl;
      cout << "Max active:      " << MBEHeuristicInstance::getMaxNumActive() << endl;
      cout << "Total time:      " << m_heurCompTime << " seconds" << endl;
      cout << "Root time:       " << m_heurRootCompTime << " seconds" << endl;
      if (m_options->dynamic && m_numHeuristics > 1) {
          double dynTime = m_heurCompTime - m_heurRootCompTime;
          double avgTime = dynTime / (m_numHeuristics - 1);
          cout << "Dynamic time:    " << dynTime << " seconds" << endl;
          cout << "Avg. time (dyn): " << avgTime << " seconds" << endl;
      }
      cout << "Max memory:      " 
          << (MBEHeuristicInstance::getMaxMemory() / (1024*1024.0)) * sizeof(double) 
          << " MB" << endl;
      cout << "---------------" << endl;
      /*
      const auto &better = m_heuristic->getHeurBetter();
      const auto &varCount = m_heuristic->getVarTimesVisited();
      cout << "var,depth,width,#better,total,ratio" << endl;
      for (unsigned int i = 0; i < better.size(); ++i) {
          Pseudotree *temp = new Pseudotree(*m_pseudotree);
          temp->restrictSubproblem(i);
          cout << i << "," << m_pseudotree->getNode(i)->getDepth() << ","  << temp->getWidthCond() << "," << better[i] << "," << varCount[i] << "," << double(better[i])/varCount[i] << endl;
          delete temp;
      }
      */

  }

  int getNumHeuristics() const {
      return m_numHeuristics;
  }

  int getCurrentNumActive() const {
      return MBEHeuristicInstance::getCurrentNumActive();
  }

  int getMaxNumActive() const {
      return MBEHeuristicInstance::getMaxNumActive();
  }

  double getCurrentMemory() const {
      return (MBEHeuristicInstance::getCurrentMemory() / (1024*1024.0)) * sizeof(double);
  }

  double getMaxMemory() const {
      return (MBEHeuristicInstance::getMaxMemory() / (1024*1024.0)) * sizeof(double);
  }

  const vector<int> &getHeurBetter() const {
      return m_heurBetter;
  }

  const vector<unsigned int> &getVarTimesVisited() const {
      return m_countVarVisited;
  }

  // Preprocess problem using FGLP/JGLP
  bool doFGLP();
  bool doJGLP();

  // Copies and conditions the problem in the root.
  bool getNodeFGLPHeur(SearchNode *n, const vector<val_t> &assignment, vector<double> &out);

public:
  MiniBucketElim(Problem* p, Pseudotree* pt, ProgramOptions* po, int ib);
  virtual ~MiniBucketElim();

protected:
  mex::vector<mex::Factor> copyFactors(void);
  void rewriteFactors(const vector<mex::Factor>& factors);


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

    friend ostream& operator<<(ostream &os, const Scope &s) {
        os << "f" << s.m_id << ":" << s.m_scope;
        return os;
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
    m_mbCountSubtree(p->getN(),vector<int>(p->getN(), 0)),
    m_dupeCount(p->getN(),vector<int>(p->getN(), 0)),
    m_mbCountAccurate(p->getN()),
    m_subproblemWidth(p->getN(),-1),
    m_currentGIter(0), 
    m_numHeuristics(0),
    m_countVarVisited(p->getN()),
    m_heurBetter(p->getN(), 0),
    m_heurBetterPre(p->getN(), vector<HeurRelation>(p->getN(),HeurRelation::EQUAL))
    { 
    m_varPartitioning.resize(p->getN());
    for (auto &nodeVarVp : m_varPartitioning) {
        nodeVarVp.resize(p->getN());
    }
    m_varPartitioningScopes.resize(p->getN());
    m_idMultiplier = 1;
    while (m_idMultiplier < p->getN()) m_idMultiplier *= 10;

        m_rootHeurInstance = new MBEHeuristicInstance(p->getN(), pt->getRoot()->getVar(), NULL);
        m_fglpRoot = NULL;

        // If dynamic, precomupute all DFS elimination orders for each node
        // Also if using node FGLP
        // Otherwise, just compute one order for the entire problem
        if (m_options->dynamic || m_options->ndfglp > 0 || m_options->ndfglps > 0) {
            m_elimOrder.resize(p->getN());
            for (int i = 0 ; i < p->getN(); ++i) {
                findDfsOrder(m_elimOrder[i], i);
            }
        }
        else {
            m_elimOrder.resize(1);
            findDfsOrder(m_elimOrder[0]);
        }

        // Perform FGLP to preprocess problem original problem first
        if (m_options->mplp > 0 || m_options->mplps > 0) {
            if (doFGLP())
                m_pseudotree->addFunctionInfo(m_problem->getFunctions());
        }

        // If dynamic, precompute number of minibuckets used in each subproblem 
        // rooted by each node
        if (m_options->dynamic) {
            int i;
            for (i = 0 ; i < p->getN() - 1; ++i) {
                m_mbCountAccurate[i] = m_elimOrder[i].size() - (i==p->getN()-1 ? 1 : 2);
                simulateBuildSubproblem(i, m_elimOrder[i], m_options->subibound);
            }
            m_mbCountAccurate[i] = m_elimOrder[i].size() - (i==p->getN()-1 ? 1 : 2);
            simulateBuildSubproblem(i, m_elimOrder[i], m_ibound);
            /*
            for (int i = 0; i < p->getN(); ++i) {
                cout << i << ", " << m_pseudotree->getNode(i)->getDepth() << ", " << m_subproblemWidth[i] << endl;
            }
            */
            for (int i = 0; i < int(m_subproblemWidth.size()); ++i) {
                m_subproblemWidth[i] = m_pseudotree->computeSubproblemWidth(i);
            }

            /*
            for (int i = 0; i < p->getN(); ++i) {
                cout << pt->getNode(i)->getVar() << " " << pt->getNode(i)->getDepth() << " " << m_mbCountSubtree[pt->getRoot()->getVar()][i] << " " << m_mbCountSubtree[i][i] << endl;
            }
            */

            // DEBUG show partitioning for root
            const auto &rootVarPartition = m_varPartitioning.back();
            for (int i = 0; i < int(rootVarPartition.size()); ++i) {
                cout << i << " : | ";
                for (const auto &mbi : rootVarPartition[i]) {
                    cout << mbi << " | ";
                }
                cout << endl;
            }

            // for each nodevar
            for (int i=0; i<p->getN();++i) {
                auto &curList = m_heurBetterPre[i];
                // process each variable
                for (int j=0; j<p->getN();++j) {
                    if (curList[j] == HeurRelation::_UNKNOWN) continue; // skip update
                    HeurRelation rel = HeurRelation::_UNKNOWN;
                    const auto &curVarPartition = m_varPartitioning[i][j];

                    // check if identical (0)
                    if (curVarPartition.size() == rootVarPartition[j].size()) {
                        // check if the function contents are identical
                        bool allIdentical = true;
                        for (unsigned int k=0; 
                                allIdentical && k<curVarPartition.size(); ++k) {
                            if (curVarPartition[k] != rootVarPartition[j][k])  
                                allIdentical = false;
                        }
                        if (allIdentical) rel = HeurRelation::EQUAL;
                    }
                    // check if dynamic heuristic is better (1)
                    else if (curVarPartition.size() < rootVarPartition[j].size()) {
                        bool allHaveSuperset = true;
                        for (const auto &er : rootVarPartition[j]) {
                            bool exists = false;
                            for (const auto &ec : curVarPartition) {
                                if ( (exists = isSubset(er,ec)) ) break;
                            }
                            allHaveSuperset = allHaveSuperset && exists;
                            if (!allHaveSuperset) break;
                        }
                        if (allHaveSuperset) rel = HeurRelation::BETTER;
                    }
                    // check if dynamic heuristic is worse (-1)
                    else if (curVarPartition.size() > rootVarPartition[j].size()) {
                        bool allHaveSuperset = true;
                        for (const auto &ec : curVarPartition) {
                            bool exists = false;
                            for (const auto &er : rootVarPartition[j]) {
                                if ( (exists = isSubset(ec,er)) ) break;
                            }
                            allHaveSuperset = allHaveSuperset && exists;
                            if (!allHaveSuperset) break;
                        }
                        if (allHaveSuperset) rel = HeurRelation::WORSE;
                    }
                    // otherwise unknown, leave flag as default
                    
                    // record result and also propagate result to parents
                    PseudotreeNode *curNode = pt->getNode(j);
                    do {
                        int v = curNode->getVar();
                        if (rel==HeurRelation::_UNKNOWN) 
                            curList[v] = HeurRelation::_UNKNOWN;
                        else {
                            switch(curList[v]) {
                                case HeurRelation::EQUAL:
                                    curList[v] = rel;
                                    break;
                                case HeurRelation::BETTER:
                                    curList[v] = 
                                        (rel==HeurRelation::WORSE) ? 
                                        HeurRelation::_UNKNOWN : 
                                        HeurRelation::BETTER;
                                    break;
                                case HeurRelation::WORSE:
                                    curList[v] = 
                                        (rel==HeurRelation::BETTER) ? 
                                        HeurRelation::_UNKNOWN : 
                                        HeurRelation::WORSE;
                                    break;
                                default:
                                    break;
                            }
                        }
                    } while ( (curNode = curNode->getParent()) );
                }
            }

            for (int i=p->getN() - 1; i>=0;--i) {
                // process each variable in the reverse elimination ordering
                int nv = m_elimOrder.back()[i];
                cout << "nodevar: " << nv << endl;
                cout << "width: " << m_subproblemWidth[nv] << endl;
                cout << "depth: " << pt->getNode(nv)->getDepth() << endl;
                for (auto itV = m_elimOrder[nv].rbegin(); 
                        itV != m_elimOrder[nv].rend(); ++itV) {
                    cout << *itV << " " << m_heurBetterPre[nv][*itV] << endl;
                    cout << "depth: " << pt->getNode(*itV)->getDepth() << endl;

                    cout << "    |";
                    for(const auto &mbi : m_varPartitioning[nv][*itV]) {
                        cout << "{ ";
                        for (const auto &f : mbi) {
                            cout << "f" << f << ":" << m_varPartitioningScopes[nv][f] << " "; 
                        }
                        cout << "}";
                        cout << "|";
                    }
                    cout << endl;

                    cout << "    |";
                    for(const auto &mbi : m_varPartitioning.back()[*itV]) {
                        cout << "{ ";
                        for (const auto &f : mbi) {
                            cout << "f" << f << ":" << m_varPartitioningScopes.back()[f] << " "; 
                        }
                        cout << "}";
                        cout << "|";
                    }
                    cout << endl;

                }
                cout << endl << "======================" << endl;
            }
            
            
        }


}

inline MiniBucketElim::~MiniBucketElim() {
  // make sure to delete each function only once
  if (m_rootHeurInstance) delete m_rootHeurInstance;
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
