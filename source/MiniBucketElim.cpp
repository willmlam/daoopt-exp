/*
 * MiniBucket.cpp
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
 *  Created on: Nov 8, 2008
 *      Author: Lars Otten <lotten@ics.uci.edu>
 */

#include "MiniBucketElim.h"

/* disables DEBUG output */
#undef DEBUG


#ifdef DEBUG
/* ostream operator for debugging */
ostream& operator <<(ostream& os, const list<Function*>& l) {
  list<Function*>::const_iterator it = l.begin();
  os << '[';
  while (it!=l.end()) {
    os << (**it);
    if (++it != l.end()) os << ',';
  }
  os << ']';
  return os;
}

ostream& operator <<(ostream& os, const vector<Function*>& l) {
  vector<Function*>::const_iterator it = l.begin();
  os << '[';
  while (it!=l.end()) {
    os << (**it);
    if (++it != l.end()) os << ',';
  }
  os << ']';
  return os;
}
#endif


/* computes the augmented part of the heuristic estimate */
double MiniBucketElim::getHeur(int var, const vector<val_t>& assignment) {

  assert( var >= 0 && var < m_problem->getN());

  // Rebuild heuristic with conditioning if dynamic
  if (m_dynamic) {
      const vector<int> &elimOrder = m_elimOrder[var]; // retrieve ordering

      // create map form of assignment
      map<int,val_t> mAssn;
      for (unsigned int i=0; i<assignment.size(); ++i)
      if (assignment[i] != -1 && 
          find(elimOrder.begin(), elimOrder.end(), i) == elimOrder.end()) {
          mAssn.insert(pair<int,val_t>(i, assignment[i]));
      }
      // while top of stack is not compatible pop functions
      while (!m_miniBucketFunctions.top()->isCompatible(var,mAssn,m_pseudotree))
          m_miniBucketFunctions.pop();
      // if the heuristic on top is not accurate, compute conditioned subproblem heuristics
#ifdef DEBUG
      cout << "stack size: " << m_miniBucketFunctions.size() << endl;
#endif
      if (m_miniBucketFunctions.empty() || !m_miniBucketFunctions.top()->isAccurate) {
#ifdef DEBUG
          cout << "Ancestor heuristic is not accurate!" << endl;
#endif
          m_miniBucketFunctions.push(new MiniBucketFunctions(var,mAssn));
          buildSubproblem(var, mAssn, assignment, elimOrder);
      }
#ifdef DEBUG
      else {
          cout << "Ancestor heuristic is accurate!" << endl;
          m_miniBucketFunctions.top()->printVarAndAssign();
          cout << endl;
      }
#endif
  }

  double h = ELEM_ONE;

  // go over augmented and intermediate lists and combine all values
  vector<Function*>::const_iterator itF = m_augmented[var].begin();
  for (; itF!=m_augmented[var].end(); ++itF) {
    h OP_TIMESEQ (*itF)->getValue(assignment);
  }

  itF = m_intermediate[var].begin();
  for (; itF!=m_intermediate[var].end(); ++itF) {
    h OP_TIMESEQ (*itF)->getValue(assignment);
  }

  return h;
}


void MiniBucketElim::getHeurAll(int var, const vector<val_t>& assignment, vector<double>& out) {
    /*
    cout << "var :" << var << endl;
    cout << "Width subproblem: " << getWidthSubproblem(var) << endl;
    */

  // Rebuild heuristic with conditioning if dynamic
  if (m_dynamic) {
      const vector<int> &elimOrder = m_elimOrder[var]; // will hold dfs order

/*
      cout << "elimOrder:" << endl;
      for (unsigned int i=0; i<elimOrder.size(); ++i)
          cout << " " << elimOrder[i];
      cout << endl;
      */
      // create map form of assignment

      map<int,val_t> mAssn;
      for (unsigned int i=0; i<assignment.size(); ++i) {
          if (assignment[i] != -1 && 
                  find(elimOrder.begin(), elimOrder.end(), i) == elimOrder.end()) {
              mAssn.insert(pair<int,val_t>(i, assignment[i]));
          }
      }

      // DEBUG PRINT mAssn
/*
      
      cout << "DEBUG PRINT mAssn" << endl;
      for(map<int,val_t>::iterator it=mAssn.begin(); it!=mAssn.end(); ++it) {
          cout << " " << it->first << " " << int(it->second) << endl;
      }
      cout << endl;

      cout << "DEBUG PRINT assignment" << endl;
      for (unsigned int i=0; i<assignment.size();++i)
          if (assignment[i] != -1) cout << " " << i << " " << int(assignment[i]) << endl;
      cout << endl;

*/
      // ======
//      cout << "stack size before: " << m_miniBucketFunctions.size() << endl;
      // while top of stack is not compatible pop functions
      assert(m_miniBucketFunctions.size());
      while (m_miniBucketFunctions.size() && !m_miniBucketFunctions.top()->isCompatible(var,mAssn,m_pseudotree)) {
//          cout << "not compatible" << endl;
          assert(m_dhDepth >=0);
/*
          cout << m_miniBucketFunctions.top()->getAssignment().size() << endl;
          m_miniBucketFunctions.top()->printAssignAndElim();
          cout << endl;
          */
          delete m_miniBucketFunctions.top();
          m_miniBucketFunctions.pop();
          assert(m_miniBucketFunctions.size());
      }
//      cout << "stack size after: " << m_miniBucketFunctions.size() << endl;
      // if the heuristic on top is not accurate, compute conditioned subproblem heuristics
      cout << "var: " << var << endl;
      int currentDepth = m_pseudotree->getNode(var)->getDepth();
      cout << "depth: " << currentDepth << endl;
      if (m_miniBucketFunctions.empty() || 
              (!m_miniBucketFunctions.top()->isAccurate && m_dhDepth > currentDepth && currentDepth % m_depthInterval == 0)) {
//          cout << "Ancestor heuristic is not accurate!" << endl;
//          cout << "Rebuilding to evaluate..." << endl;
          /*
          for(map<int,val_t>::iterator it=mAssn.begin(); it!=mAssn.end(); ++it) {
              cout << " " << it->first << " " << int(it->second) << endl;
          }
          cout << endl;
          cout << m_miniBucketFunctions.top().getAssignment().size() << endl;
          m_miniBucketFunctions.top().printAssignAndElim();
          cout << endl;
          */
          if (m_gNodes > 0 && m_currentGIter == 0) {
              m_miniBucketFunctions.push(new MiniBucketFunctions(var,mAssn));
              buildSubproblem(var, mAssn, assignment, elimOrder);
          }
          m_currentGIter = (m_currentGIter + 1) % m_gNodes;
      }
      else {
//          cout << "Ancestor heuristic is accurate!" << endl;
//          cout << "Using heuristic with this assignment..." << endl;
          /*
          cout << m_miniBucketFunctions.top().getAssignment().size() << endl;
          m_miniBucketFunctions.top().printAssignAndElim();
          cout << "to evaluate..." << endl;
          for(map<int,val_t>::iterator it=mAssn.begin(); it!=mAssn.end(); ++it) {
              cout << " " << it->first << " " << int(it->second) << endl;
          }
          cout << endl;
          cout << endl;
          */
      }
  }

#ifdef DEBUG
  cout << "Using this heuristic: " << endl;
  m_miniBucketFunctions.top()->printVarAndAssign();
#endif

  out.clear();
  out.resize(m_problem->getDomainSize(var), ELEM_ONE);
  vector<double> funVals;
  vector<Function*>::const_iterator itF;
  for (itF = m_augmented[var].begin(); itF!=m_augmented[var].end(); ++itF) {
    (*itF)->getValues(assignment, var, funVals);
    for (size_t i=0; i<out.size(); ++i)
      out[i] OP_TIMESEQ funVals[i];
  }
  for (itF = m_intermediate[var].begin(); itF!=m_intermediate[var].end(); ++itF) {
    (*itF)->getValues(assignment, var, funVals);
    for (size_t i=0; i<out.size(); ++i)
      out[i] OP_TIMESEQ funVals[i];
  }
}


void MiniBucketElim::reset() {

  while(m_miniBucketFunctions.size()) {
      delete m_miniBucketFunctions.top();
      m_miniBucketFunctions.pop();
  }
//  m_miniBucketFunctions.push(MiniBucketFunctions());

/*
  vector<vector<Function*> > empty;
  m_augmented.swap(empty);

  vector<vector<Function*> > empty2;
  m_intermediate.swap(empty2);
  */

}


size_t MiniBucketElim::build(const vector<val_t> * assignment, bool computeTables) {

#ifdef DEBUG
  cout << "$ Building MBE(" << m_ibound << ")" << endl;
#endif

  this->reset();

  vector<int> elimOrder; // will hold dfs order
  findDfsOrder(elimOrder); // computes dfs ordering of relevant subtree
  m_miniBucketFunctions.push(new MiniBucketFunctions(m_pseudotree->getRoot()->getVar()));

  m_miniBucketFunctions.top()->isAccurate = true;

  m_augmented.resize(m_problem->getN());
  m_intermediate.resize(m_problem->getN());

  // keep track of total memory consumption
  size_t memSize = 0;

  // ITERATES OVER BUCKETS, FROM LEAVES TO ROOT
  for (vector<int>::reverse_iterator itV=elimOrder.rbegin(); itV!=elimOrder.rend(); ++itV) {

#ifdef DEBUG
    cout << "$ Bucket for variable " << *itV << endl;
#endif

    // collect relevant functions in funs
    vector<Function*> funs;
    const vector<Function*>& fnlist = m_pseudotree->getFunctions(*itV);
    funs.insert(funs.end(), fnlist.begin(), fnlist.end());
    funs.insert(funs.end(), m_augmented[*itV].begin(), m_augmented[*itV].end());
#ifdef DEBUG
    for (vector<Function*>::iterator itF=funs.begin(); itF!=funs.end(); ++itF)
      cout << ' ' << (**itF);
    cout << endl;
#endif

    // compute global upper bound for root (dummy) bucket
    if (*itV == elimOrder[0]) {// variable is dummy root variable
      if (computeTables && assignment) { // compute upper bound if assignment is given
        m_globalUB = ELEM_ONE;
        for (vector<Function*>::iterator itF=funs.begin(); itF!=funs.end(); ++itF)
          m_globalUB OP_TIMESEQ (*itF)->getValue(*assignment);
        cout << "    MBE-ALL  = " << SCALE_LOG(m_globalUB) << " (" << SCALE_NORM(m_globalUB) << ")" << endl;
        m_globalUB OP_DIVIDEEQ m_problem->globalConstInfo();  // for backwards compatibility of output
        cout << "    MBE-ROOT = " << SCALE_LOG(m_globalUB) << " (" << SCALE_NORM(m_globalUB) << ")" << endl;
      }
      continue; // skip the dummy variable's bucket
    }

    // sort functions by decreasing scope size
    sort(funs.begin(), funs.end(), scopeIsLarger);

    // partition functions into minibuckets
    vector<MiniBucket> minibuckets;
//    vector<Function*>::iterator itF; bool placed;
    for (vector<Function*>::iterator itF = funs.begin(); itF!=funs.end(); ++itF) {
//    while (funs.size()) {
      bool placed = false;
      for (vector<MiniBucket>::iterator itB=minibuckets.begin();
            !placed && itB!=minibuckets.end(); ++itB)
      {
        if (itB->allowsFunction(*itF)) { // checks if function fits into bucket
          itB->addFunction(*itF);
          placed = true;
        }
      }
      if (!placed) { // no fit, need to create new bucket
        MiniBucket mb(*itV,m_ibound,m_problem);
        mb.addFunction(*itF);
        minibuckets.push_back(mb);
      }
//      funs.pop_front();
      //funs.erase(itF);
    }


    // minibuckets for current bucket are now ready, process each
    // and place resulting function


    // Compute max-marginals for each bucket

    vector<Function*> maxMarginals;
    double *avgMMTable;
    Function *avgMaxMarginal;
    
    if (m_momentMatching && minibuckets.size() > 1) {
        // Find intersection of scopes (the scope of all max-marginals)
        vector<MiniBucket>::iterator itB=minibuckets.begin();
        set<int> intersectScope(itB->getJointScope());
        itB++;
        for (; itB!=minibuckets.end(); ++itB)
        {
            set<int> newInter = intersection(intersectScope, itB->getJointScope());
            intersectScope = newInter;
        }
        for (vector<MiniBucket>::iterator itB=minibuckets.begin();
                itB!=minibuckets.end(); ++itB)
        {
            set<int> elimVars = setminus(itB->getJointScope(), intersectScope);
#ifdef DEBUG
            Function *eF = itB->eliminate(computeTables, elimVars);
            cout << "Max-marginal values: " << endl;
            for (int i=0; i < eF->getTableSize(); ++i) {
                cout << ' ' << ELEM_DECODE(eF->getTable()[i]) << endl;
            }
            cout << endl;
            //        cin.get();
#endif
            maxMarginals.push_back(itB->eliminate(computeTables, elimVars));
        }

        // Find average max-marginal (geometric mean)
        size_t tablesize = 1;
        for (set<int>::iterator sit=intersectScope.begin(); 
                sit!=intersectScope.end(); ++sit) {
            tablesize *= m_problem->getDomainSize(*sit);
        }

        avgMMTable = new double[tablesize];
#ifdef DEBUG
        cout << "Tablesize: " << tablesize << ", buckets: " << minibuckets.size() << endl;
#endif
        for (unsigned int i = 0; i < tablesize; ++i) avgMMTable[i] = ELEM_ONE;
        for (vector<Function*>::iterator itMM=maxMarginals.begin();
                itMM!=maxMarginals.end(); ++itMM) {
            for (unsigned int i = 0; i < tablesize; ++i) {
                avgMMTable[i] OP_TIMESEQ (*itMM)->getTable()[i];
            }
        }
        for (unsigned int i = 0; i < tablesize; ++i) {
            avgMMTable[i] = OP_ROOT(avgMMTable[i],minibuckets.size());
        }
#ifdef DEBUG
        cout << "Avg max-marginal: " << endl;
        for (unsigned int i = 0; i < tablesize; ++i) {
            cout << ' ' << avgMMTable[i] << endl;
        }
#endif
    
        int ammid = 0;
        avgMaxMarginal = new FunctionBayes(ammid, m_problem, intersectScope, avgMMTable, tablesize);
    }

    int bucketIdx = 0;

    if (m_miniBucketFunctions.top()->isAccurate && minibuckets.size() > 1) 
        m_miniBucketFunctions.top()->isAccurate = false;

    for (vector<MiniBucket>::iterator itB=minibuckets.begin();
          itB!=minibuckets.end(); ++itB, ++bucketIdx)
    {

      // Replace this to generate moment-matched version if #minibuckets > 1
      Function* newf;
      if (!m_momentMatching || minibuckets.size() <= 1)
          newf = itB->eliminate(computeTables); // process the minibucket
      else {
          newf = itB->eliminateMM(computeTables,
                  maxMarginals[bucketIdx],avgMaxMarginal); // process the minibucket
      }

      const set<int>& newscope = newf->getScopeSet();
      memSize += newf->getTableSize();
      // go up in tree to find target bucket
      PseudotreeNode* n = m_pseudotree->getNode(*itV)->getParent();
      while (newscope.find(n->getVar()) == newscope.end() && n != m_pseudotree->getRoot() ) {
        m_intermediate[n->getVar()].push_back(newf);
        n = n->getParent();
      }
      // matching bucket found OR root of pseudo tree reached
      m_augmented[n->getVar()].push_back(newf);
    }
    // all minibuckets processed and resulting functions placed

    // free up memory used by max-marginals
    if (m_momentMatching && minibuckets.size() > 1) {
      for (unsigned int i = 0; i < maxMarginals.size(); ++i)
        delete maxMarginals[i];
      maxMarginals.clear();
      if(avgMaxMarginal) delete avgMaxMarginal;
    }
  }

#ifdef DEBUG
  // output augmented and intermediate buckets
  if (computeTables)
    for (int i=0; i<m_problem->getN(); ++i) {
      cout << "$ AUG" << i << ": " << m_augmented[i] << " + " << m_intermediate[i] << endl;
    }
#endif

  // clean up for estimation mode
  if (!computeTables) {
    for (vector<vector<Function*> >::iterator itA = m_augmented.begin(); itA!=m_augmented.end(); ++itA)
      for (vector<Function*>::iterator itB = itA->begin(); itB!=itA->end(); ++itB)
        delete *itB;
    m_augmented.clear();
  }

  return memSize;
}

size_t MiniBucketElim::buildSubproblem(int var, const map<int,val_t> &assignment, const vector<val_t> &vAssn, const vector<int> &elimOrder, 
bool computeTables) {


#ifdef DEBUG
  cout << "$ Building MBE(" << m_ibound << ")" << endl;
#endif

//  this->reset();
  assert(m_miniBucketFunctions.top()->isEmpty());
  m_miniBucketFunctions.top()->isAccurate = true;


  m_augmented.resize(m_problem->getN());
  m_intermediate.resize(m_problem->getN());

  // keep track of total memory consumption
  size_t memSize = 0;

  // ITERATES OVER BUCKETS, FROM LEAVES TO ROOT
  for (vector<int>::const_reverse_iterator itV=elimOrder.rbegin(); itV!=elimOrder.rend(); ++itV) {

#ifdef DEBUG
    cout << "$ Bucket for variable " << *itV << endl;
#endif

    // collect relevant functions in funs
    vector<Function*> funs;
    const vector<Function*>& fnlist = m_pseudotree->getFunctions(*itV);
    vector<Function*> condfnlist;
    for(vector<Function*>::const_iterator itF=fnlist.begin(); itF!=fnlist.end(); ++itF) {
        condfnlist.push_back((*itF)->substitute(assignment));
    }
    funs.insert(funs.end(), condfnlist.begin(), condfnlist.end());
    funs.insert(funs.end(), m_augmented[*itV].begin(), m_augmented[*itV].end());
#ifdef DEBUG
    for (vector<Function*>::iterator itF=funs.begin(); itF!=funs.end(); ++itF)
      cout << ' ' << (**itF);
    cout << endl;
#endif

// test
    if (funs.size() == 0) continue;
// ===


    // compute global upper bound for root (dummy) bucket
    if (*itV == elimOrder[0]) {// variable is dummy root variable
      if (computeTables) { // compute upper bound if assignment is given
        m_globalUB = ELEM_ONE;
        for (vector<Function*>::iterator itF=funs.begin(); itF!=funs.end(); ++itF)
          m_globalUB OP_TIMESEQ (*itF)->getValue(vAssn);
//        cout << "    MBE-ALL  = " << SCALE_LOG(m_globalUB) << " (" << SCALE_NORM(m_globalUB) << ")" << endl;
        m_globalUB OP_DIVIDEEQ m_problem->globalConstInfo();  // for backwards compatibility of output
//        cout << "    MBE-ROOT = " << SCALE_LOG(m_globalUB) << " (" << SCALE_NORM(m_globalUB) << ")" << endl;
      }
      continue; // skip the dummy variable's bucket
    }

    // sort functions by decreasing scope size
    sort(funs.begin(), funs.end(), scopeIsLarger);

    // partition functions into minibuckets
    vector<MiniBucket> minibuckets;
//    vector<Function*>::iterator itF; bool placed;
    for (vector<Function*>::iterator itF = funs.begin(); itF!=funs.end(); ++itF) {
//    while (funs.size()) {
      bool placed = false;
      for (vector<MiniBucket>::iterator itB=minibuckets.begin();
            !placed && itB!=minibuckets.end(); ++itB)
      {
        if (itB->allowsFunction(*itF)) { // checks if function fits into bucket
          itB->addFunction(*itF);
          placed = true;
        }
      }
      if (!placed) { // no fit, need to create new bucket
        MiniBucket mb(*itV,m_ibound,m_problem);
        mb.addFunction(*itF);
        minibuckets.push_back(mb);
      }
//      funs.pop_front();
      //funs.erase(itF);
    }


    // minibuckets for current bucket are now ready, process each
    // and place resulting function

    // Compute max-marginals for each bucket

    vector<Function*> maxMarginals;
    double *avgMMTable; 
    Function *avgMaxMarginal;
    
    if (m_momentMatching && minibuckets.size() > 1) {
        // Find intersection of scopes (the scope of all max-marginals)
        vector<MiniBucket>::iterator itB=minibuckets.begin();
        set<int> intersectScope(itB->getJointScope());
        itB++;
        for (; itB!=minibuckets.end(); ++itB)
        {
            set<int> newInter = intersection(intersectScope, itB->getJointScope());
            intersectScope = newInter;
        }
#ifdef DEBUG
        cout << "Intersection: " << endl;
        for (set<int>::iterator it = intersectScope.begin(); it!=intersectScope.end(); ++it) {
            cout << ' ' << *it;
        }
        cout << endl;
#endif
        for (vector<MiniBucket>::iterator itB=minibuckets.begin();
                itB!=minibuckets.end(); ++itB)
        {
            set<int> elimVars = setminus(itB->getJointScope(), intersectScope);
#ifdef DEBUG
            Function *eF = itB->eliminate(computeTables, elimVars);
            cout << "Max-marginal values: " << endl;
            for (int i=0; i < eF->getTableSize(); ++i) {
                cout << ' ' << ELEM_DECODE(eF->getTable()[i]) << endl;
            }
            cout << endl;
            //        cin.get();
#endif
            maxMarginals.push_back(itB->eliminate(computeTables, elimVars));
        }

        // Find average max-marginal (geometric mean)
        size_t tablesize = 1;
        for (set<int>::iterator sit=intersectScope.begin(); 
                sit!=intersectScope.end(); ++sit) {
            tablesize *= m_problem->getDomainSize(*sit);
        }

        avgMMTable = new double[tablesize];
#ifdef DEBUG
        cout << "Tablesize: " << tablesize << ", buckets: " << minibuckets.size() << endl;
#endif
        for (unsigned int i = 0; i < tablesize; ++i) avgMMTable[i] = 0;
        for (vector<Function*>::iterator itMM=maxMarginals.begin();
                itMM!=maxMarginals.end(); ++itMM) {
            for (unsigned int i = 0; i < tablesize; ++i) {
                avgMMTable[i] OP_TIMESEQ (*itMM)->getTable()[i];
            }
        }
        for (unsigned int i = 0; i < tablesize; ++i) {
            avgMMTable[i] = OP_ROOT(avgMMTable[i],minibuckets.size());
        }
#ifdef DEBUG
        cout << "Avg max-marginal: " << endl;
        for (unsigned int i = 0; i < tablesize; ++i) {
            cout << ' ' << avgMMTable[i] << endl;
        }
#endif
    
        int ammid = 0;
        avgMaxMarginal = new FunctionBayes(ammid, m_problem, intersectScope, avgMMTable, tablesize);
    }

    int bucketIdx = 0;

    if (minibuckets.size() > 1) 
        m_miniBucketFunctions.top()->isAccurate = false;

    for (vector<MiniBucket>::iterator itB=minibuckets.begin();
          itB!=minibuckets.end(); ++itB, ++bucketIdx)
    {

      // Replace this to generate moment-matched version if #minibuckets > 1
      Function* newf;
      if (!m_momentMatching || minibuckets.size() <= 1)
          newf = itB->eliminate(computeTables); // process the minibucket
      else {
          newf = itB->eliminateMM(computeTables,
                  maxMarginals[bucketIdx],avgMaxMarginal); // process the minibucket
      }

      const set<int>& newscope = newf->getScopeSet();
      memSize += newf->getTableSize();
      // go up in tree to find target bucket
      PseudotreeNode* n = m_pseudotree->getNode(*itV)->getParent();
      while (newscope.find(n->getVar()) == newscope.end() && n != m_pseudotree->getRoot() ) {
        m_intermediate[n->getVar()].push_back(newf);
        n = n->getParent();
      }
      // matching bucket found OR root of pseudo tree reached
      m_augmented[n->getVar()].push_back(newf);
    }
    // all minibuckets processed and resulting functions placed

    // free up memory used by max-marginals
    if (m_momentMatching && minibuckets.size() > 1) {
      for (unsigned int i = 0; i < maxMarginals.size(); ++i)
        delete maxMarginals[i];
      maxMarginals.clear();
      if(avgMaxMarginal) delete avgMaxMarginal;
    }
  }

#ifdef DEBUG
  // output augmented and intermediate buckets
  if (computeTables)
    for (int i=0; i<m_problem->getN(); ++i) {
      cout << "$ AUG" << i << ": " << m_augmented[i] << " + " << m_intermediate[i] << endl;
    }
#endif

  // clean up for estimation mode
  if (!computeTables) {
    for (vector<vector<Function*> >::iterator itA = m_augmented.begin(); itA!=m_augmented.end(); ++itA)
      for (vector<Function*>::iterator itB = itA->begin(); itB!=itA->end(); ++itB)
        delete *itB;
    m_augmented.clear();
    m_augmented.clear();
  }

  return memSize;
}


/* finds a dfs order of the pseudotree (or the locally restricted subtree)
 * and writes it into the argument vector */
void MiniBucketElim::findDfsOrder(vector<int>& order) const {
  order.clear();
  stack<PseudotreeNode*> dfs;
  dfs.push(m_pseudotree->getRoot());
  PseudotreeNode* n = NULL;
  while (!dfs.empty()) {
    n = dfs.top();
    dfs.pop();
    order.push_back(n->getVar());
    for (vector<PseudotreeNode*>::const_iterator it=n->getChildren().begin();
          it!=n->getChildren().end(); ++it) {
      dfs.push(*it);
    }
  }
}

/* finds a dfs order of the pseudotree (or the locally restricted subtree)
 * and writes it into the argument vector 
 * this version takes in the variable which is the root of the restricted subtree */
void MiniBucketElim::findDfsOrder(vector<int>& order, int var) const {
  order.clear();
  stack<PseudotreeNode*> dfs;
  if (m_pseudotree->getRoot()->getVar() != var) {
      order.push_back(m_pseudotree->getRoot()->getVar());
  }
  dfs.push(m_pseudotree->getNode(var));
  PseudotreeNode* n = NULL;
  while (!dfs.empty()) {
    n = dfs.top();
    dfs.pop();
    order.push_back(n->getVar());
    for (vector<PseudotreeNode*>::const_iterator it=n->getChildren().begin();
          it!=n->getChildren().end(); ++it) {
      dfs.push(*it);
    }
  }
}


size_t MiniBucketElim::limitSize(size_t memlimit, const vector<val_t> * assignment) {

  // convert to bits
  memlimit *= 1024 *1024 / sizeof(double);

  int ibound = m_options->ibound;

  cout << "Adjusting mini bucket i-bound..." << endl;
  this->setIbound(ibound);
  size_t mem = this->build(assignment, false);
  cout << " i=" << ibound << " -> " << ((mem / (1024*1024.0)) * sizeof(double) )
       << " MBytes" << endl;

  while (mem > memlimit && ibound > 1) {
    this->setIbound(--ibound);
    mem = this->build(assignment, false);
    cout << " i=" << ibound << " -> " << ((mem / (1024*1024.0)) * sizeof(double) )
         << " MBytes" << endl;
  }

  m_options->ibound = ibound;
  return mem;
}


size_t MiniBucketElim::getSize() const {
  size_t S = 0;
  for (vector<vector<Function*> >::const_iterator it=m_augmented.begin(); it!= m_augmented.end(); ++it) {
    for (vector<Function*>::const_iterator itF=it->begin(); itF!=it->end(); ++itF)
      S += (*itF)->getTableSize();
  }
  return S;
}

int MiniBucketElim::getWidthSubproblem(int i) const {
    const vector<int>& condset = m_pseudotree->getNode(i)->getFullContextVec();
    stack<PseudotreeNode*> stck; 
    stck.push(m_pseudotree->getNode(i));
    int width = -1;
    while(stck.size()) {
        PseudotreeNode *n = stck.top();
        stck.pop();
        int x = setminusSize(n->getFullContextVec(), condset);
        width = max(width,x);
        for (vector<PseudotreeNode*>::const_iterator it=n->getChildren().begin(); it!=n->getChildren().end(); ++it) {
            stck.push(*it);
        }
    }
    return width;
}


/*
 * mini bucket file format (all data in binary):
 * - size_t: no. of variables
 * - int: i-bound
 * - double: global upper bound
 * for every variable:
 *   - size_t: number of functions in bucket structure
 *   for every such function:
 *     - int: function ID
 *     - size_t: scope size
 *     for every scope variable:
 *       - int: variable index
 *     - size_t: table size
 *     for every table entry:
 *       - double: CPT entry
 * for every variable:
 *   - size_t: number of intermediate function pointers
 *   for every function pointer:
 *     - size_t: function index (implicit order from above)
 */

bool MiniBucketElim::writeToFile(string fn) const {

  ogzstream out(fn.c_str());
  if ( ! out ) {
    cerr << "Error writing mini buckets to file " << fn << endl;
    return false;
  }

  // used later
  int x = NONE;
  size_t y = NONE;

  // number of variables
  size_t sz = m_augmented.size();
  out.write((char*)&( sz ), sizeof( sz ));

  // i-bound
  out.write((char*)&( m_ibound ), sizeof( m_ibound ));

  // global UB
  out.write((char*)&( m_globalUB ), sizeof( m_globalUB ));

  map<const Function*,size_t> funcMap;

  // over m_augmented
  for (size_t i=0; i<sz; ++i) {
    size_t sz2 = m_augmented[i].size();
    out.write((char*)&( sz2 ), sizeof( sz2 ));

    vector<Function*>::const_iterator itF = m_augmented[i].begin();
    for (size_t j=0; j<sz2; ++j, ++itF) {
      const Function* f = *itF;
      funcMap.insert(make_pair(f,funcMap.size()));

      // function ID
      int id = f->getId();
      out.write((char*)&( id ), sizeof( id ));


      // scope
      size_t sz3 = f->getScopeVec().size();
      out.write((char*)&( sz3 ), sizeof( sz3 ));
 // scope size
      for (vector<int>::const_iterator it=f->getScopeVec().begin(); it!=f->getScopeVec().end(); ++it) {
        x = *it;
        out.write((char*)&( x ), sizeof( x ));
 // vars from scope
      }

      // table size
      sz3 = f->getTableSize();
      out.write((char*)&( sz3 ), sizeof( sz3 ));

      // table
      out.write((char*) ( f->getTable() ), sizeof( double ) * sz3);

    }


  }

  // over m_intermediate
  for (size_t i=0; i<sz; ++i) {
    size_t sz2 = m_intermediate[i].size();
    out.write((char*)&( sz2 ), sizeof( sz2 ));

    vector<Function*>::const_iterator itF = m_intermediate[i].begin();
    for (size_t j=0; j<sz2; ++j, ++itF) {
      y = funcMap.find(*itF)->second;
      out.write((char*) &( y ), sizeof(y));
    }

  }

  out.close();

  return true;

}


bool MiniBucketElim::readFromFile(string fn) {

  ifstream inTemp(fn.c_str());
  inTemp.close();
  if (inTemp.fail()) { // file not existent yet
    return false;
  }

  igzstream in(fn.c_str());

  this->reset();

  // used later
  int x = NONE;
  size_t y = NONE;
  vector<Function*> allFuncs;

  // no. of variables
  size_t sz;
  in.read((char*) &( sz ), sizeof( sz ));

  if (sz != (size_t) m_problem->getN()) {
    cerr << "Number of variables in mini bucket file doesn't match" << endl;
    return false;
  }

  m_augmented.resize(sz);
  m_intermediate.resize(sz);

  // i-bound
  int ibound;
  in.read((char*) &(ibound), sizeof(ibound));
  m_ibound = ibound;

  // global UB
  double ub;
  in.read((char*) &(ub), sizeof( ub ));
  m_globalUB = ub;

  // over variables for m_augmented
  for (size_t i=0; i<sz; ++i) {
    size_t sz2;
    in.read((char*) &(sz2), sizeof(sz2));

    // over functions
    for (size_t j=0; j<sz2; ++j) {
      int id;
      in.read((char*) &( id ), sizeof(id));

      // scope
      size_t sz3;
      in.read((char*) &(sz3), sizeof(sz3));
      set<int> scope;
      for (size_t k=0; k<sz3; ++k) {
        in.read((char*) &( x ), sizeof(x));
        scope.insert(x);
      }

      // table size and table
      in.read((char*) &( sz3 ), sizeof(sz3));
      double* T = new double[sz3];
      in.read((char*) ( T ), sizeof(double)*sz3);

      // create function and store it
      Function* f = new FunctionBayes(id,m_problem,scope,T,sz3);
      m_augmented[i].push_back(f);
      allFuncs.push_back(f);
    }
  }

  for (size_t i=0; i<sz; ++i) {
    // no. of function pointers
    size_t sz2;
    in.read((char*) &(sz2), sizeof(sz2));

    for (size_t j=0; j<sz2; ++j) {
      // function index
      in.read((char*) &(y), sizeof(y));
      m_intermediate[i].push_back(allFuncs.at(y));
    }
  }

  in.close();
  cout << "Read mini bucket with i-bound " << ibound << " from file " << fn << endl;
  return true;
}


