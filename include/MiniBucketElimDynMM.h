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

#ifndef MINIBUCKETELIMDYNMM_H_
#define MINIBUCKETELIMDYNMM_H_

//#define m_augmented m_miniBucketFunctions.top().m_augmentedF
//#define m_intermediate m_miniBucketFunctions.top().m_intermediateF

#include "Heuristic.h"
#include "Function.h"
#include "Problem.h"
#include "ProgramOptions.h"
#include "Pseudotree.h"
#include "utils.h"

#include "MiniBucket.h"

#undef DEBUG


/* The overall minibucket elimination */
class MiniBucketFunctionTree;

class MiniBucketElimDynMM : public Heuristic {

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

  MiniBucketFunctionTree *m_tree;

  bool m_initialized;

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

  // builds the minibucket tree structure:
  size_t initialize();

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
  MiniBucketElimDynMM(Problem* p, Pseudotree* pt, ProgramOptions* po, int ib);
  virtual ~MiniBucketElimDynMM();

};

// Class to store the minibucket functions, along with the partitioning
class MiniBucketFunctionTree {
    friend class MiniBucketElimDynMM;

    // These can store augmented or intermediate messages, used to keep track of messages and also their 
    // entire history on a stack
    class FunctionEdge {
        MiniBucket *m_mbSource; // The minibucket which the message came from
        MiniBucket *m_mbTarget; // The minibucket which the message is sent to
        int m_fIdx; // The index of the function in the minibucket

        stack<Function*> m_fStack;
        stack<map<int,val_t> > m_aStack;
    public:
        FunctionEdge(MiniBucket *source, MiniBucket *target, int fIdx) 
        : m_mbSource(source), m_mbTarget(target), m_fIdx(fIdx) {}
        Function *getFunction() const { return m_mbTarget->getFunctionRef(m_fIdx); }
        MiniBucket *getSource() const { return m_mbSource; }
        //void setSource(MiniBucket *mb) { m_mbSource = mb; }
        MiniBucket *getTarget() const { return m_mbTarget; }
        //void setTarget(MiniBucket *mb) { m_mbTarget = mb; }

        void pushFunction(Function *f, const map<int,val_t> &assn) {
            m_fStack.push(f);
            m_aStack.push(assn);
            m_mbTarget->getFunctionRef(m_fIdx) = m_fStack.top();
        }
        void popFunction() {
            assert(m_fStack.size() >= 2); // should never have to pop empty assignment
#ifdef DEBUG
            cout << "Popping function " << *(m_fStack.top()) << endl;
#endif
            delete m_fStack.top();
            m_fStack.pop();
            m_aStack.pop();
            m_mbTarget->getFunctionRef(m_fIdx) = m_fStack.top();
        }

        size_t getStackSize() {
            return m_fStack.size();
        }

        bool isInDummyState() {
            return m_fStack.empty();
        }

        bool isCompatible(const map<int,val_t> &assignment) const {
            if (m_aStack.empty()) return true;
#ifdef DEBUG
            cout << "current edge assign: " << m_aStack.top() << endl;
            cout << "stack size: " << m_aStack.size() << endl;
            cout << "assign to test: " << assignment << endl;
#endif
            for (map<int,val_t>::const_iterator it=m_aStack.top().begin(); 
                    it!=m_aStack.top().end(); ++it) {
                if (assignment.find(it->first) == assignment.end() ||
                        it->second != assignment.find(it->first)->second) {
                    return false;
                }
            }
            return true;
        }
        bool assignmentIsEqual(const map<int,val_t> &assignment) const {
            if (m_aStack.empty()) return false;
            return m_aStack.top() == assignment;
        }

        void backtrack(const map<int,val_t> &assignment) {
            while (!isCompatible(assignment)) {
                popFunction();
            }
        }

        virtual ~FunctionEdge() {
            m_mbTarget->getFunctionRef(m_fIdx) = NULL;
            while(!m_fStack.empty()) {
                assert(m_fStack.top() != NULL);
                delete m_fStack.top();
                m_fStack.pop();
            }
        }
    };

    // Stores the partitions of functions at each bucket
    vector<vector<MiniBucket *> > m_minibuckets;

    // Stores the separators between minibuckets for each variable
    vector<set<int> > m_separators;
    vector<size_t> m_sepSizes;

    // Preallocate memory to store max-marginals and average max-marginal
    vector<vector<Function *> > m_maxMarginals;
    vector<Function *> m_avgMaxMarginals;

    // Stores references to messages which are dynamic
    // vector index indicates the destination variable of the edge
    vector<vector<FunctionEdge *> > m_dynamicMessages;

    // Stores references to messages which are dynamic
    // vector index indicates the source variable of the edge
    // corresponds to message sent by a minibucket
    vector<vector<FunctionEdge *> > m_dynamicMessagesSource;

    // Stores the same edges but according to those which go to an ancestor
    // of the vector index variable
    vector<vector<FunctionEdge *> > m_intermediateMessages;

    Problem *m_problem;
    Pseudotree *m_pseudotree;
    
public:

    MiniBucketFunctionTree(Problem *p, Pseudotree *pt) : m_problem(p), m_pseudotree(pt) {
        m_minibuckets.resize(m_problem->getN());
        m_dynamicMessages.resize(m_problem->getN());
        m_dynamicMessagesSource.resize(m_problem->getN());
        m_intermediateMessages.resize(m_problem->getN());
        m_maxMarginals.resize(m_problem->getN());
        m_avgMaxMarginals.resize(m_problem->getN());
    }

    vector<MiniBucket*> &getVarMiniBuckets(int var) {
        return m_minibuckets[var];
    }

    const vector<FunctionEdge*> &getIncomingMessages(int var) {
        return m_dynamicMessages[var];
    }

    const vector<FunctionEdge*> &getIntermediateMessages(int var) {
        return m_intermediateMessages[var];
    }

    void addEdge(MiniBucket *source, MiniBucket *target, const vector<int> &intermediateVars, 
            int fIdx) {
        FunctionEdge *newEdge = new FunctionEdge(source, target, fIdx);
#ifdef DEBUG
        cout << "Function edge addr: " << newEdge;
        cout << ", fIdx: " << fIdx << ", addr: " << target->getFunctionRef(fIdx) << endl;
        cout << ", source: " << source->getVar() << ", target: " << target->getVar() << endl;
        cout << ", source addr: " << source << ", target addr: " << target << endl;
#endif
        assert(target->getFunctionRef(fIdx));
        m_dynamicMessages[target->getVar()].push_back(newEdge);
        m_dynamicMessagesSource[source->getVar()].push_back(newEdge);
        for (unsigned i = 0; i < intermediateVars.size(); ++i) {
            m_intermediateMessages[intermediateVars[i]].push_back(newEdge);
        }
    }

    // updates the message passed in by reprocessing the source bucket
    // and also any other messages that go into descendants of the variable 
    // that are affected by the update
    void firstOrderUpdate(FunctionEdge *message, const map<int,val_t> &assignment) {
        int source = message->getSource()->getVar();
        int target = message->getTarget()->getVar();
#ifdef DEBUG
        cout << "Performing first-order update on " << message << " from " << source << endl;
#endif
        queue<int> msgTargetQ;
        updateMessages(message->getSource()->getVar(), assignment);

        /*
        vector<FunctionEdge*>::iterator it = m_dynamicMessagesSource[source].begin();
        for (; it != m_dynamicMessagesSource[source].end(); ++it) {
            MiniBucket *targetBucket = (*it)->getTarget();
            if (*it == message || 
                    targetBucket->getJointScope().find(target) == targetBucket->getJointScope().end())
                continue;
            int newTarget = (*it)->getTarget()->getVar();
#ifdef DEBUG
            cout << "pushing variable " << newTarget << endl;
#endif
            msgTargetQ.push(newTarget);
        }
        */
        /*
        while (!msgTargetQ.empty()) {
            int currentSource = msgTargetQ.front();
            msgTargetQ.pop();
            updateMessages(currentSource, assignment);
            it = m_dynamicMessagesSource[currentSource].begin();
            for (; it != m_dynamicMessagesSource[currentSource].end(); ++it) {
                MiniBucket *targetBucket = (*it)->getTarget();
                if (targetBucket->getJointScope().find(target) == targetBucket->getJointScope().end()) 
                    continue;
                int newTarget = targetBucket->getVar();
                cout << "pushing variable " << newTarget << endl;
                msgTargetQ.push(newTarget);
            }
        }
        */
    }

    // Makes the edges coming into the variable consistent with the assignment
    void backtrack(int var, const map<int,val_t> &assignment) {
        vector<FunctionEdge*>::iterator it = m_dynamicMessages[var].begin();
        for (; it != m_dynamicMessages[var].end(); ++it)
            (*it)->backtrack(assignment);
    }

    bool assignmentIsEqual(int var, const map<int,val_t> &assignment) {
        vector<FunctionEdge*>::iterator it = m_dynamicMessagesSource[var].begin();
        for (; it != m_dynamicMessagesSource[var].end(); ++it) {
            if (!(*it)->assignmentIsEqual(assignment)) {
                return false;
            }
        }
        return true;
    }

    // Processes a bucket with a given assignment, updating its outgoing messages
    void updateMessages(int var, const map<int,val_t> &assignIn, bool forceUpdate=false) {

        // Make a copy and remove unnecessary conditioning
        map<int,val_t> assignment(assignIn);
        set<int>::iterator sit = m_separators[var].begin();
        for (; sit != m_separators[var].end(); ++sit) {
            assignment.erase(*sit);
        }
        // First check the bucket's current assignment by checking outgoing messages
        // no need to perform updates if all match
        if (assignmentIsEqual(var, assignment) && !forceUpdate) return;


#ifdef DEBUG
        cout << "Updating " << var << " with assignment: " << assignment << endl;
#endif

        // First ensure incoming messages have a consistent assignment
        backtrack(var, assignment);

        vector<FunctionEdge*> &outgoingMessages = m_dynamicMessagesSource[var];

        // no partitioning
        if (outgoingMessages.size() <= 1) {
            FunctionEdge* msg = *(outgoingMessages.begin());
            MiniBucket *mb = msg->getSource();
                
            // Compute for the first time
            if (msg->isInDummyState()) {
                delete msg->getFunction(); // delete dummy function
                msg->pushFunction(mb->eliminate(), assignment);
#ifdef DEBUG
                cout << "Message sent to " << msg->getTarget()->getVar() << " (address: " << msg->getTarget() << ")" << endl;
#endif
            }
            else if (forceUpdate) {
                msg->pushFunction(mb->eliminate(), assignment);
            }
#ifdef DEBUG
            else {
                cout << "No partitioning! No recomputation." << endl;
            }
#endif
        }
        else {

            vector<FunctionEdge*>::iterator itV;
#ifdef DEBUG
            cout << " ======" << endl;
            cout << "var: " << var << endl;
            cout << "assign: " << assignment << endl;
            cout << "orig sep: " << m_separators[var] << endl;
#endif

            /*
            size_t tempSepSize = m_sepSizes[var];
            set<int> tempSep = m_separators[var];

            // remove conditioned variables from separator
            for (set<int>::iterator it = tempSep.begin(); it != tempSep.end(); ) {
                set<int>::iterator toErase = it;
                ++it;
                if (assignment.find(*toErase) != assignment.end()) {
                    int erasedDomain = m_problem->getDomainSize(*toErase);
                    tempSep.erase(toErase);
                    tempSepSize /= erasedDomain;
                }
            }
            */

            // Compute conditioned max-marginals
            // For the minibucket of each outgoing message, compute the max-marginal
#ifdef DEBUG
            cout << " Conditioning: " << assignment << endl;
#endif
            vector<Function*> &maxMarginals = m_maxMarginals[var];
            int count = 0;
            for (itV = outgoingMessages.begin(); itV != outgoingMessages.end(); ++itV, ++count) {
                MiniBucket *mb = (*itV)->getSource();
                set<int> elimVar = setminus(mb->getJointScope(), m_separators[var]);
#ifdef DEBUG
                cout << "Computing a max-marginal." << endl;
                cout << "Joint scope: " << mb->getJointScope() << endl;
                cout << "Separator: "<< m_separators[var] << endl;
                cout << "Elim vars: "<< elimVar << endl;
#endif
                if (maxMarginals.size() != outgoingMessages.size()) {
                    Function *maxMarg = mb->conditionEliminate(true,assignment,elimVar);
                    maxMarginals.push_back(maxMarg);
                }
                else {
                    mb->conditionEliminateInPlace(true,assignment,elimVar,maxMarginals[count]->getTable());
                }
#ifdef DEBUG
                cout << maxMarginals[count]->getScopeSet().size() << endl;
#endif
            }

            // Compute average max-marginal
            double *avgMMTable = m_avgMaxMarginals[var]->getTable();
            for (unsigned i = 0; i < m_sepSizes[var]; ++i) {
                avgMMTable[i] = ELEM_ONE;
            }
            for (vector<Function*>::iterator itMM=maxMarginals.begin();
                    itMM!=maxMarginals.end(); ++itMM) {
                for (unsigned j = 0; j < m_sepSizes[var]; ++j) {
                    avgMMTable[j] OP_TIMESEQ (*itMM)->getTable()[j];
                }
            }
            for (unsigned i = 0; i < m_sepSizes[var]; ++i) {
                avgMMTable[i] = OP_ROOT(avgMMTable[i],maxMarginals.size());
            }
            /*
            Function *avgMaxMarginal = new FunctionBayes(var, m_problem, 
                    tempSep, avgMMTable, tempSepSize);
                    */

            assert(maxMarginals.size() == outgoingMessages.size());
            vector<Function*>::iterator itF;
            // Compute updated messages
            for (itV = outgoingMessages.begin(), itF = maxMarginals.begin(); itV != outgoingMessages.end(); ++itV, ++itF) {
                MiniBucket *mb = (*itV)->getSource();
                if ((*itV)->isInDummyState()) {
                    delete (*itV)->getFunction(); // delete dummy function
                }
                (*itV)->pushFunction(mb->conditionEliminateMM(true, assignment, 
                            *itF, m_avgMaxMarginals[var]), assignment);
            }
        }
    }

    void computeSeparators() {
        m_separators.clear();
        for (unsigned i = 0; i < m_minibuckets.size(); ++i) {
            if (m_minibuckets[i].empty()) {
                m_separators.push_back(set<int>());
                continue;
            }
            m_separators.push_back(m_minibuckets[i][0]->getJointScope());
            set<int> &currentSep = m_separators.back();
            for (int j = 1; j < int(m_minibuckets[i].size()); ++j) {
                set<int> newInter = intersection(currentSep, m_minibuckets[i][j]->getJointScope());
                currentSep = newInter;
            }
            m_sepSizes.push_back(1);
            for (set<int>::iterator sit=currentSep.begin(); sit!=currentSep.end(); ++sit) {
                m_sepSizes.back() *= m_problem->getDomainSize(*sit);
            }

            // ...also initialize average max marginal functions
            m_avgMaxMarginals[i] = new FunctionBayes(i, m_problem, currentSep, new double[m_sepSizes.back()], m_sepSizes.back());
        }
    }

    int getSumOfStackSizes() {
        int result = 0;
        for (unsigned i = 0; i < m_dynamicMessages.size(); ++i) {
            for (unsigned j = 0; j < m_dynamicMessages[i].size(); ++j) {
                result += m_dynamicMessages[i][j]->getStackSize();
            }
        }
        return result;
    }
    
public:
    virtual ~MiniBucketFunctionTree() {
        for (unsigned i = 0; i < m_dynamicMessages.size(); ++i) {
            for (unsigned j = 0; j < m_dynamicMessages[i].size(); ++j) {
                delete m_dynamicMessages[i][j];
            }
        }
        m_dynamicMessages.clear();
        for (unsigned i = 0; i < m_minibuckets.size(); ++i) {
            delete m_avgMaxMarginals[i];
            for (unsigned j = 0; j < m_minibuckets[i].size(); ++j) {
                delete m_minibuckets[i][j];
                delete m_maxMarginals[i][j];
            }
        }
    }
};

/* Inline definitions */

inline bool MiniBucketElimDynMM::isAccurate() {
  assert(m_pseudotree);
  return (m_pseudotree->getWidthCond() <= m_ibound);
}

inline MiniBucketElimDynMM::MiniBucketElimDynMM(Problem* p, Pseudotree* pt,
				      ProgramOptions* po, int ib) :
    Heuristic(p, pt, po), m_ibound(ib), m_globalUB(ELEM_ONE), m_initialized(false) { 
        m_tree = new MiniBucketFunctionTree(p,pt); 
    }

inline MiniBucketElimDynMM::~MiniBucketElimDynMM() {
  delete m_tree;
}

typedef pair<Function*,int> IndexedFunction;

inline bool scopeIsLargerIndexed(IndexedFunction p, IndexedFunction q) {
    assert(p.first && q.first);
  if (p.first->getArity() == q.first->getArity())
    return (p.first->getId() > q.first->getId());
  else
    return (p.first->getArity() > q.first->getArity());
}



#endif /* MINIBUCKETELIMDYNMM_H_ */
