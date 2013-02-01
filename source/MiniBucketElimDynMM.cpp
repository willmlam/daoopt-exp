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

#include "MiniBucketElimDynMM.h"

/* disables DEBUG output */
#undef DEBUG


#ifdef DEBUGGA
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
double MiniBucketElimDynMM::getHeur(int var, const vector<val_t>& assignment) {

  assert( var >= 0 && var < m_problem->getN());


  map<int,val_t> assign;
  const vector <int> &relVars = m_pseudotree->getNode(var)->getFullContextVec();
  for (unsigned i = 0; i < relVars.size(); ++i) {
      assign[relVars[i]] = assignment[relVars[i]];
  }

  double h = ELEM_ONE;

  // use minibucketfunction tree for messages
  const vector<MiniBucketFunctionTree::FunctionEdge* > &messages = m_tree->getIncomingMessages(var);

  vector<MiniBucketFunctionTree::FunctionEdge*>::const_iterator it;

//  cout << "Assignment: " << assign << endl;

  for (it = messages.begin(); it != messages.end(); ++it) {
      //m_tree->firstOrderUpdate(*it,assign);
      h OP_TIMESEQ (*it)->getFunction()->getValue(assignment);
  }

  const vector<MiniBucketFunctionTree::FunctionEdge* > &ancMessages = m_tree->getIntermediateMessages(var);
  for (it = ancMessages.begin(); it != ancMessages.end(); ++it) {
      //m_tree->firstOrderUpdate(*it,assign);
      h OP_TIMESEQ (*it)->getFunction()->getValue(assignment);
  }
#ifdef DEBUG
  cout << "number of messages stored: " << m_tree->getSumOfStackSizes();
  cout << ", depth: " << m_pseudotree->getNode(var)->getDepth() << endl;
#endif
  return h;
}


void MiniBucketElimDynMM::getHeurAll(int var, const vector<val_t>& assignment, vector<double>& out) {

  out.clear();
  out.resize(m_problem->getDomainSize(var), ELEM_ONE);
  /*
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
  */
}


void MiniBucketElimDynMM::reset() {

    //m_tree->reset();
//  m_miniBucketFunctions.push(MiniBucketFunctions());

/*
  vector<vector<Function*> > empty;
  m_augmented.swap(empty);

  vector<vector<Function*> > empty2;
  m_intermediate.swap(empty2);
  */

}

// initializes the mini bucket tree structure
size_t MiniBucketElimDynMM::initialize() {
    this->reset();
    vector<int> elimOrder; 
    findDfsOrder(elimOrder);

    size_t memSize = 0;

    // Temporary storage for messages
    vector<vector<Function*> > augmented;
    // Temporary storage indicating the source of the message
    vector<vector<MiniBucket*> > source;
    // Storage storing all of the variables on the pseudo tree path
    vector<vector<vector<int> > > intermediateVars;
    augmented.resize(m_problem->getN());
    source.resize(m_problem->getN());
    intermediateVars.resize(m_problem->getN());
    

    // iterate over buckets, from leaves to root
    for (vector<int>::reverse_iterator itV=elimOrder.rbegin(); itV!=elimOrder.rend(); ++itV) {
        //if (*itV == elimOrder[0]) continue; // Skip root
#ifdef DEBUG
        cout << "$ Bucket for variable " << *itV << endl;
#endif
        vector<IndexedFunction> funs;
        const vector<Function*>& fnlist = m_pseudotree->getFunctions(*itV);
        for (unsigned i = 0; i < fnlist.size(); ++i) {
            funs.push_back(IndexedFunction(fnlist[i],funs.size()));
        }

        // The first index that deals with messages
        int offset = funs.size();
        for (unsigned i = 0; i < augmented[*itV].size(); ++i) {
            funs.push_back(IndexedFunction(augmented[*itV][i],funs.size()));
        }


        // sort functions by decreasing scope size
        sort(funs.begin(), funs.end(), scopeIsLargerIndexed);

        // partition functions into minibuckets
        vector<MiniBucket*> &minibuckets = m_tree->getVarMiniBuckets(*itV);
        //    vector<Function*>::iterator itF; bool placed;
        
        for (vector<IndexedFunction>::iterator itF = funs.begin(); itF!=funs.end(); ++itF) {
            //cout << itF->first << " " << itF->second << " " << itF->first->getTable() <<endl;
            bool placed = false;
            for (vector<MiniBucket*>::iterator itB=minibuckets.begin();
                    !placed && itB!=minibuckets.end(); ++itB)
            {
                if ((*itB)->allowsFunction(itF->first)) { // checks if function fits into bucket
                    int idx = (*itB)->addFunction(itF->first);
                    if (itF->second >= offset) {
                        //cout << "Above is a message" << endl;
                        assert((*itB)->getFunctionRef(idx)->getTable() == NULL);
                        m_tree->addEdge(source[*itV][itF->second-offset],*itB,
                                intermediateVars[*itV][itF->second-offset],idx);
                        //cout << "Adding edge with Function addr: " << (*itB)->getFunctionRef(idx) << endl;
                    }
                    placed = true;
                }
            }
            if (!placed) { // no fit, need to create new bucket
                minibuckets.push_back(new MiniBucket(*itV,m_ibound,m_problem));
                int idx = minibuckets.back()->addFunction(itF->first);
                if (itF->second >= offset) {
                    //cout << "Above is a message" << endl;
                    assert(minibuckets.back()->getFunctionRef(idx)->getTable() == NULL);
                    m_tree->addEdge(source[*itV][itF->second-offset],minibuckets.back(),
                            intermediateVars[*itV][itF->second-offset],idx); 
                    //cout << "Adding edge with Function addr: " 
                        //<< minibuckets.back()->getFunctionRef(idx) << endl;
                }
            }
        }

        if (*itV == elimOrder[0]) continue;
        // Ready to send messages
        for (vector<MiniBucket*>::iterator itB=minibuckets.begin();
                itB!=minibuckets.end(); ++itB)
        {
            // process the minibucket to get a dummy function with scope information
            Function *newf = (*itB)->eliminate(false); 

            const set<int>& newscope = newf->getScopeSet();
            memSize += newf->getTableSize();
            // go up in tree to find target bucket
            PseudotreeNode* n = m_pseudotree->getNode(*itV)->getParent();
            //PseudotreeNode* child = m_pseudotree->getNode(*itV);
            vector<int> pathVars;
            while (newscope.find(n->getVar()) == newscope.end() && n != m_pseudotree->getRoot() ) {
                pathVars.push_back(n->getVar());
                n = n->getParent();
            }
            // matching bucket found OR root of pseudo tree reached
            
            augmented[n->getVar()].push_back(newf);
            source[n->getVar()].push_back(*itB);
            intermediateVars[n->getVar()].push_back(pathVars);
        }
        // all minibuckets processed and resulting functions placed
    }

    /*
#ifdef DEBUG
    // output augmented and intermediate buckets
    if (computeTables)
    for (int i=0; i<m_problem->getN(); ++i) {
    cout << "$ AUG" << i << ": " << m_augmented[i] << " + " << m_intermediate[i] << endl;
    }
#endif
*/

    // clean up for estimation mode
    
    m_tree->computeSeparators(); 
    m_initialized = true;
    return memSize;
}


size_t MiniBucketElimDynMM::build(const vector<val_t> * assignment, bool computeTables) {

#ifdef DEBUG
  cout << "$ Building MBE(" << m_ibound << ")" << endl;
#endif
  if (!m_initialized) initialize();


  vector<int> elimOrder; // will hold dfs order
  findDfsOrder(elimOrder); // computes dfs ordering of relevant subtree

  // keep track of total memory consumption
  size_t memSize = 0;

  map<int,val_t> assign;

  // ITERATES OVER BUCKETS, FROM LEAVES TO ROOT
  for (vector<int>::reverse_iterator itV=elimOrder.rbegin(); itV!=elimOrder.rend(); ++itV) {
    // compute global upper bound for root (dummy) bucket
    if (*itV == elimOrder[0]) {// variable is dummy root variable
        // use minibucketfunction tree for messages
        m_globalUB = ELEM_ONE;
        const vector<MiniBucketFunctionTree::FunctionEdge* > &messages = 
            m_tree->getIncomingMessages(*itV);
        for (unsigned i=0; i<messages.size(); ++i) {
            m_globalUB OP_TIMESEQ messages[i]->getFunction()->getValue(*assignment);
        }
        cout << "    MBE-ALL  = " << SCALE_LOG(m_globalUB) << " (" << SCALE_NORM(m_globalUB) << ")" << endl;
        m_globalUB OP_DIVIDEEQ m_problem->globalConstInfo();  // for backwards compatibility of output
        cout << "    MBE-ROOT = " << SCALE_LOG(m_globalUB) << " (" << SCALE_NORM(m_globalUB) << ")" << endl;

      continue; // skip the dummy variable's bucket
    }
    //cout << "Updating " << *itV << endl;
    m_tree->updateMessages(*itV, assign);

  }

/*
#ifdef DEBUG
  // output augmented and intermediate buckets
  if (computeTables)
    for (int i=0; i<m_problem->getN(); ++i) {
      cout << "$ AUG" << i << ": " << m_augmented[i] << " + " << m_intermediate[i] << endl;
    }
#endif
    */

  // clean up for estimation mode
  return memSize;
}

/* finds a dfs order of the pseudotree (or the locally restricted subtree)
 * and writes it into the argument vector */
void MiniBucketElimDynMM::findDfsOrder(vector<int>& order) const {
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
void MiniBucketElimDynMM::findDfsOrder(vector<int>& order, int var) const {
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

size_t MiniBucketElimDynMM::limitSize(size_t memlimit, const vector<val_t> * assignment) {

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


// TO DO LATER
size_t MiniBucketElimDynMM::getSize() const {
  size_t S = 0;
  /*
  for (vector<vector<Function*> >::const_iterator it=m_augmented.begin(); it!= m_augmented.end(); ++it) {
    for (vector<Function*>::const_iterator itF=it->begin(); itF!=it->end(); ++itF)
      S += (*itF)->getTableSize();
  }
  */
  return S;
}

int MiniBucketElimDynMM::getWidthSubproblem(int i) const {
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


bool MiniBucketElimDynMM::writeToFile(string fn) const {
  return false;
}


bool MiniBucketElimDynMM::readFromFile(string fn) {
  return false;
}


