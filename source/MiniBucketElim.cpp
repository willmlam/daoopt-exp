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


//#ifdef DEBUG
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

ostream& operator <<(ostream& os, const set<Function*>& l) {
  set<Function*>::const_iterator it = l.begin();
  os << '[';
  while (it!=l.end()) {
    os << (**it);
    if (++it != l.end()) os << ',';
  }
  os << ']';
  return os;
}
//#endif


/* computes the augmented part of the heuristic estimate */
// TO DO: fix for new data structure
double MiniBucketElim::getHeur(int var, const vector<val_t>& assignment, SearchNode *n) {

  assert( var >= 0 && var < m_problem->getN());

  double h = ELEM_ONE;

  /*
  // go over augmented and intermediate lists and combine all values
  set<Function*>::const_iterator itF = m_augmented[var].begin();
  for (; itF!=m_augmented[var].end(); ++itF) {
    h OP_TIMESEQ (*itF)->getValue(assignment);
  }

  itF = m_intermediate[var].begin();
  for (; itF!=m_intermediate[var].end(); ++itF) {
    h OP_TIMESEQ (*itF)->getValue(assignment);
  }
  */

  return h;
}


void MiniBucketElim::getHeurAll(int var, const vector<val_t>& assignment, SearchNode *n, vector<double>& out) {
    /*
    cout << "var :" << var << endl;
    cout << "Width subproblem: " << getWidthSubproblem(var) << endl;
    */
    m_countVarVisited[var]++;

    SearchNodeOR *sNode = dynamic_cast<SearchNodeOR*>(n);
    assert(sNode);
    // Handle initial root search node
    if (sNode->getHeurInstance() == NULL) {
        sNode->setHeurInstance(m_rootHeurInstance);
        //m_rootHeurInstance->setOwner(sNode);
    }
    int currentDepth = m_pseudotree->getNode(var)->getDepth();

    /*
      map<int,val_t> mAssn;
      const vector<int> &relVars = m_pseudotree->getNode(var)->getFullContextVec();
      for (unsigned i = 0; i < relVars.size(); ++i) {
          mAssn[relVars[i]] = assignment[relVars[i]];
      }
      cout << "Assignment: " << mAssn << endl;
      cout << "Assignment: " << assignment << endl;
      */

  // Rebuild heuristic with conditioning if dynamic
  if (m_options->dynamic) {
      /*
      if (sNode->isHeuristicLocked()) {
          cout << "Heuristic is locked." << endl;
      }
      */
      if (m_pseudotree->getRoot()->getVar() != var &&
              !(sNode->isHeuristicLocked()) &&
              !(m_options->reuseLevel != 0 && sNode->getHeurInstance()->getAccurateHeurIn()[var]) &&
              meetsComputeConditions(var, sNode->getHeurInstance()->getVar(), currentDepth)) {
          MBEHeuristicInstance *newHeur = 
              new MBEHeuristicInstance(m_problem->getN(),var,sNode->getHeurInstance());
          newHeur->setDepth(currentDepth);
          m_heurCollection.push_back(newHeur);
          buildSubproblem(var, 
                  assignment, 
                  sNode->getHeurInstance(), 
                  newHeur);
          sNode->setHeurInstance(newHeur);
          newHeur->setOwner(sNode);
          /*
             cout << "Compiled heuristic at (var,depth): " 
             << "(" << var << "," << currentDepth << ")" << endl;
             */
      }
      m_currentGIter = (m_currentGIter + 1) % m_options->gNodes;
      assert(sNode->getHeurInstance()->getDepth() <= currentDepth);
  }

  /*
  const vector<vector<Function*> >&augmented = sNode->getHeurInstance()->getAugmented();
  const vector<vector<Function*> >&intermediate = sNode->getHeurInstance()->getIntermediate();
  */
  MBEHeuristicInstance *curHeur = sNode->getHeurInstance();
  /*
  if (curHeur->getParent()) 
      curHeur = curHeur->getParent();
      */
  const vector<vector<Function*> >&augmented = curHeur->getAugmented();
  const vector<vector<Function*> >&intermediate = curHeur->getIntermediate();

  /*
  if (sNode->getHeurInstance() == m_rootHeurInstance) {
      cout << "using root heuristic" << endl;
  }
  else {
      cout << sNode->getHeurInstance() << endl;
      //cout << mAssn << endl;
      cout << "using conditioned heuristic" << endl;
  }
  */


  out.clear();
  out.resize(m_problem->getDomainSize(var), ELEM_ONE);


  vector<double> funVals;
  vector<Function*>::const_iterator itF;
  for (itF = augmented[var].begin(); itF!=augmented[var].end(); ++itF) {
    (*itF)->getValues(assignment, var, funVals);
    for (size_t i=0; i<out.size(); ++i)
      out[i] OP_TIMESEQ funVals[i];
  }
  for (itF = intermediate[var].begin(); itF!=intermediate[var].end(); ++itF) {
    (*itF)->getValues(assignment, var, funVals);
    for (size_t i=0; i<out.size(); ++i)
      out[i] OP_TIMESEQ funVals[i];
  }

  // Check bound if the parent was used instead if this was a newly computed heuristic
  MBEHeuristicInstance *parentHeur = curHeur->getParent();
  if (!sNode->isHeuristicLocked() && parentHeur) {
      vector<double> tempOut;
      tempOut.resize(m_problem->getDomainSize(var), ELEM_ONE);
      const vector<vector<Function*> >&augmentedAnc = parentHeur->getAugmented();
      const vector<vector<Function*> >&intermediateAnc = parentHeur->getIntermediate();
      for (itF = augmentedAnc[var].begin(); itF!=augmentedAnc[var].end(); ++itF) {
          (*itF)->getValues(assignment, var, funVals);
          for (size_t i=0; i<tempOut.size(); ++i)
              tempOut[i] OP_TIMESEQ funVals[i];
      }
      for (itF = intermediateAnc[var].begin(); itF!=intermediateAnc[var].end(); ++itF) {
          (*itF)->getValues(assignment, var, funVals);
          for (size_t i=0; i<tempOut.size(); ++i)
              tempOut[i] OP_TIMESEQ funVals[i];
      }
      vector<double> diff;
      bool printResults = false;
      diff.resize(m_problem->getDomainSize(var));
      for (size_t i=0; i<diff.size(); ++i) {
          if (isinf(out[i]) && isinf(tempOut[i])) 
              diff[i] = 0;
          else 
              diff[i] = out[i] - tempOut[i];
          if (diff[i] != 0) printResults = true;
      }

      double maxDecrease = diff[0];
      for (size_t i=1; i<diff.size(); ++i) {
          maxDecrease = min(maxDecrease, diff[i]);
      }
      curHeur->setMostRecentDecrease(maxDecrease);
      if (false && printResults && m_pseudotree->getNode(var)->getDepth() <= 30) {
          cout << "var,depth: " << var << "," << m_pseudotree->getNode(var)->getDepth() << endl;
          cout << "Maximum decrease: " << maxDecrease << endl;
          cout << endl;
      }
  }

  if (m_options->useSimpleHeurSelection) {
      vector<double> tempOut;
      tempOut.resize(m_problem->getDomainSize(var), ELEM_ONE);
      // modified to look at all parent heuristics as well
      // guarantees heuristic computed dominates
      while ( (curHeur = curHeur->getParent()) ) {
          const vector<vector<Function*> >&augmentedAnc = curHeur->getAugmented();
          const vector<vector<Function*> >&intermediateAnc = curHeur->getIntermediate();
          for (itF = augmentedAnc[var].begin(); itF!=augmentedAnc[var].end(); ++itF) {
              (*itF)->getValues(assignment, var, funVals);
              for (size_t i=0; i<tempOut.size(); ++i)
                  tempOut[i] OP_TIMESEQ funVals[i];
          }
          for (itF = intermediateAnc[var].begin(); itF!=intermediateAnc[var].end(); ++itF) {
              (*itF)->getValues(assignment, var, funVals);
              for (size_t i=0; i<tempOut.size(); ++i)
                  tempOut[i] OP_TIMESEQ funVals[i];
          }
          for (size_t i=0; i<out.size(); ++i) {
              out[i] = min(out[i],tempOut[i]);
              tempOut[i] = ELEM_ONE;
          }
      }
  }



  /*
  const vector<vector<Function*> >&augmentedR = m_rootHeurInstance->getAugmented();
  const vector<vector<Function*> >&intermediateR = m_rootHeurInstance->getIntermediate();

  vector<double> outRoot;
  outRoot.resize(m_problem->getDomainSize(var), ELEM_ONE);

  for (itF = augmentedR[var].begin(); itF!=augmentedR[var].end(); ++itF) {
    (*itF)->getValues(assignment, var, funVals);
    for:(size_t i=0; i<outRoot.size(); ++i)
      outRoot[i] OP_TIMESEQ funVals[i];
  }
  for (itF = intermediateR[var].begin(); itF!=intermediateR[var].end(); ++itF) {
    (*itF)->getValues(assignment, var, funVals);
    for (size_t i=0; i<out.size(); ++i)
      outRoot[i] OP_TIMESEQ funVals[i];
  }
  if (true || sNode->getHeurInstance() != m_rootHeurInstance) {
      cout << "var: " << var << endl;
      cout << mAssn << endl;
      //cout << sNode->getHeurInstance()->getAssignment() << endl;
      cout << out << endl;
  }
  */

}


void MiniBucketElim::reset() {
    if (m_rootHeurInstance) delete m_rootHeurInstance;
    m_rootHeurInstance = new MBEHeuristicInstance(m_problem->getN(), m_pseudotree->getRoot()->getVar(), NULL);

//  m_miniBucketFunctions.push(MiniBucketFunctions());

/*
  vector<vector<Function*> > empty;
  m_augmented.swap(empty);

  vector<vector<Function*> > empty2;
  m_intermediate.swap(empty2);
  */

}


// Computes the root heuristic instance
size_t MiniBucketElim::build(const vector<val_t> * assignment, bool computeTables) {

#ifdef DEBUG
  cout << "$ Building MBE(" << m_ibound << ")" << endl;
#endif

  this->reset();

  time_t heurCompStart, heurCompEnd;
  time(&heurCompStart);

  const vector<int> &elimOrder = 
      (m_options->dynamic ? 
       m_elimOrder[m_pseudotree->getRoot()->getVar()] :
       m_elimOrder[0]);

  //m_miniBucketFunctions.push(new MiniBucketFunctions(m_pseudotree->getRoot()->getVar()));

  //m_miniBucketFunctions.top()->isAccurate = true;

  // keep track of total memory consumption
  size_t memSize = 0;

  vector<bool> visited(m_problem->getN(), false);

  vector<vector<Function*> > &augmented = m_rootHeurInstance->getAugmented();
  vector<ConditionedMessages*> &computedMessages = m_rootHeurInstance->getComputedMessages();
  vector<bool> &accurateHeurIn = m_rootHeurInstance->getAccurateHeurIn();

  // ITERATES OVER BUCKETS, FROM LEAVES TO ROOT
  for (vector<int>::const_reverse_iterator itV=elimOrder.rbegin(); itV!=elimOrder.rend(); ++itV) {

#ifdef DEBUG
    cout << "$ Bucket for variable " << *itV << endl;
#endif

    // collect relevant functions in funs
    vector<Function*> funs;
    const vector<Function*>& fnlist = m_pseudotree->getFunctions(*itV);
    funs.insert(funs.end(), fnlist.begin(), fnlist.end());
    funs.insert(funs.end(), augmented[*itV].begin(), augmented[*itV].end());
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
        m_rootHeurInstance->setUpperBound(m_globalUB);
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
    Function *avgMaxMarginal = NULL;

    // messages are accurate if no partitioning and incoming messages were also accurate
    bool accurateHeur = accurateHeurIn[*itV] && minibuckets.size() == 1;
    
    if (computeTables && m_options->match && minibuckets.size() > 1) {
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
            maxMarginals.push_back(itB->eliminate(computeTables, elimVars));
        }

        // Find average max-marginal (geometric mean)
        size_t tablesize = 1;
        for (set<int>::iterator sit=intersectScope.begin(); 
                sit!=intersectScope.end(); ++sit) {
            tablesize *= m_problem->getDomainSize(*sit);
        }

        avgMMTable = new double[tablesize];
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
    
        int ammid = 0;
        avgMaxMarginal = new FunctionBayes(ammid, m_problem, intersectScope, avgMMTable, tablesize);
    }

    int bucketIdx = 0;

    map<int,val_t> emptyAssignment;
    computedMessages[*itV] = new ConditionedMessages(m_rootHeurInstance, accurateHeur);

    for (vector<MiniBucket>::iterator itB=minibuckets.begin();
          itB!=minibuckets.end(); ++itB, ++bucketIdx)
    {

      // Replace this to generate moment-matched version if #minibuckets > 1
      Function* newf;
      if (!computeTables || !m_options->match || minibuckets.size() <= 1)
          newf = itB->eliminate(computeTables); // process the minibucket
      else {
          newf = itB->eliminateMM(computeTables,
                  maxMarginals[bucketIdx],avgMaxMarginal); // process the minibucket
      }

      /*
      for (int i = 0; i < newf->getTableSize(); ++i) {
        cout << " " << newf->getTable()[i] << endl;
      }
      cout << endl;
      */

      const set<int>& newscope = newf->getScopeSet();
      memSize += newf->getTableSize();

      vector<int> *path = new vector<int>();
      // go up in tree to find target bucket
      PseudotreeNode* n = m_pseudotree->getNode(*itV)->getParent();
      while (newscope.find(n->getVar()) == newscope.end() && n != m_pseudotree->getRoot() ) {
        path->push_back(n->getVar());
        /*
        m_accurateHeuristicIn[n->getVar()] = 
            m_accurateHeuristicIn[n->getVar()] && m_accurateHeuristic[*itV];
        */
        n = n->getParent();
      }
      // matching bucket found OR root of pseudo tree reached
      path->push_back(n->getVar());
      m_rootHeurInstance->getNumValues()[n->getVar()] += newf->getTableSize();
      m_rootHeurInstance->getNumZeroValues()[n->getVar()] += newf->getTableSize() - newf->getTightness();
      /*
      m_accurateHeuristicIn[n->getVar()] = 
          m_accurateHeuristicIn[n->getVar()] && m_accurateHeuristic[*itV];
      */
      if (*itV == 190) {
          cout << newf;
      }
      computedMessages[*itV]->addFunction(newf, path);
    }
    m_rootHeurInstance->populateMessages(*itV,visited);
    // all minibuckets processed and resulting functions placed

    // free up memory used by max-marginals
    if (m_options->match && minibuckets.size() > 1) {
      for (unsigned int i = 0; i < maxMarginals.size(); ++i)
        delete maxMarginals[i];
      maxMarginals.clear();
      if(avgMaxMarginal) delete avgMaxMarginal;
    }
  }

#ifdef DEBUG
  // output augmented and intermediate buckets
  vector<vector<Function*> > &intermediate = m_rootHeurInstance->getIntermediate();
  if (computeTables)
    for (int i=0; i<m_problem->getN(); ++i) {
      cout << "$ AUG" << i << ": " << augmented[i] << " + " << intermediate[i] << endl;
    }
#endif

  // clean up for estimation mode
  /*
  if (!computeTables) {
    for (vector<vector<Function*> >::iterator itA = m_augmented.begin(); itA!=m_augmented.end(); ++itA)
      for (vector<Function*>::iterator itB = itA->begin(); itB!=itA->end(); ++itB)
        delete *itB;
    m_augmented.clear();
  }
  */

  // adaptively set maxDynHeur
  if (m_options->dynamic && m_options->memlimit > 0) {
      //cout << memSize << endl;
      int maxDynHeur = m_options->memlimit / (memSize / (1024*1024.0)) * sizeof(double);
      if (maxDynHeur < m_maxDynHeur) {
          cout << "Cannot fit " << m_maxDynHeur << " heuristics in memory." << endl;
          cout << "Adjusting maximum number of heuristics to: " << maxDynHeur << endl;
          m_maxDynHeur = maxDynHeur;
      }
  }

  time(&heurCompEnd);
  m_heurCompTime += difftime(heurCompEnd,heurCompStart);
  m_numHeuristics++;

  /*
  cout << "Root heuristic computed. " << endl;
  for (int i = 0; i < m_problem->getN(); ++i) {
      cout << i << "," << m_pseudotree->getNode(i)->getDepth() << "," << m_rootHeurInstance->getNumValues()[i] << "," << m_rootHeurInstance->getNumZeroValues()[i] << "," << m_rootHeurInstance->getConstrainedness(i) << endl;
  }
  cout << endl;
  */
  return memSize;
}

size_t MiniBucketElim::buildSubproblem(int var, const vector<val_t> &vAssn, 
        MBEHeuristicInstance *ancHeur,
        MBEHeuristicInstance *curHeur, bool computeTables) {
    time_t heurCompStart, heurCompEnd;
    time(&heurCompStart);

    /*
    cout << "buildsub called: " << m_numHeuristics << endl;
    cout << "var,depth: " << var << ", " << m_pseudotree->getNode(var)->getDepth() << endl;
    */

#ifdef DEBUG
  cout << "$ Building MBE(" << m_ibound << ")" << endl;
#endif


  // should only be here if heuristic is not accurate or if not reusing messages
  assert(m_options->reuseLevel == 0 || !ancHeur->getAccurateHeurIn()[var]); 

  // should never be recomputing the root heuristic
  assert(var != m_pseudotree->getRoot()->getVar());


  map<int,val_t> assignment;
  const vector<int> &relVars = m_pseudotree->getNode(var)->getFullContextVec();
  const set<int> &context = m_pseudotree->getNode(var)->getFullContext();
  for (unsigned i = 0; i < relVars.size(); ++i) {
      assignment[relVars[i]] = vAssn[relVars[i]];
  }
  const vector<int> &elimOrder = m_elimOrder[var];// will hold dfs order
  /*
  cout << elimOrder << endl;
  //cout << "Assignment: " << assignment << endl;
  cout << "BuildSubproblem" << endl;
  cout << "var: " << var << endl;
  //cout << "assign: " << assignment << endl;
  */

//  this->reset();

  // keep track of total memory consumption
  size_t memSize = 0;

  vector<bool> visited(m_problem->getN(), false);
  vector<bool> forceUpdate(m_problem->getN(), false);

  vector<vector<Function*> > &augmented = curHeur->getAugmented();
  vector<ConditionedMessages*> &computedMessages = curHeur->getComputedMessages();
  vector<bool> &accurateHeurIn = curHeur->getAccurateHeurIn();



  // ITERATES OVER BUCKETS, FROM LEAVES TO ROOT
  for (vector<int>::const_reverse_iterator itV=elimOrder.rbegin(); itV!=elimOrder.rend(); ++itV) {

    if (*itV == var || *itV == elimOrder[0]) {
        /*
        cout << "will skip var "  << var << endl;
        cout << "assign " << assignment << endl;
        */
        continue;
    }
#ifdef DEBUG
    cout << "$ Bucket for variable " << *itV << endl;
#endif

    // reuse if child messages were not updated and conditions met
    if (!forceUpdate[*itV] && meetsReuseCondition(var,ancHeur->getVar(),*itV)) {
            computedMessages[*itV] = ancHeur->getComputedMessages()[*itV];   
            curHeur->populateMessages(*itV,visited);
            const vector<Function*> &funs = computedMessages[*itV]->getFunctions();
            const vector<vector <int> *> &paths = computedMessages[*itV]->getPaths();
            for (unsigned i = 0 ; i < funs.size(); ++i) {
                    int destVar = paths[i]->back();
                    curHeur->getNumValues()[destVar] += 
                        funs[i]->getTableSize();
                    curHeur->getNumValues()[destVar] += 
                        funs[i]->getTableSize() - funs[i]->getTightness();
            }

            continue;
    }

    // otherwise start recomputing bucket

    // collect relevant functions in funs
    vector<pair<int, Function*> > funs;
    const vector<Function*>& fnlist = m_pseudotree->getFunctions(*itV);
    /*
    for(vector<Function*>::const_iterator itF=fnlist.begin(); itF!=fnlist.end(); ++itF) {
        funs.push_back((*itF)->substitute(assignment));
    }
    */
    int fcount = 0;
    for(vector<Function*>::const_iterator itF=fnlist.begin(); itF!=fnlist.end(); ++itF, ++fcount) {
        funs.push_back(pair<int,Function*>(
                    getConditionedArity((*itF)->getScopeSet(), 
                        context), 
                    *itF));
    }
    //funs.insert(funs.end(),fnlist.begin(),fnlist.end());
    /*
    for(set<Function*>::const_iterator itF=m_augmented[*itV].begin(); 
            itF!=m_augmented[*itV].end(); ++itF) {
        funs.push_back((*itF)->substitute(assignment));
    }
    */
    for(vector<Function*>::const_iterator itF=augmented[*itV].begin(); 
            itF!=augmented[*itV].end(); ++itF, ++fcount) {
        funs.push_back(pair<int,Function*>(
                    getConditionedArity((*itF)->getScopeSet(), 
                        context),
                    *itF));
    }
    //funs.insert(funs.end(),m_augmented[*itV].begin(),m_augmented[*itV].end());
    /*
    // Check later: better to condition augmented?
    funs.insert(funs.end(), m_augmented[*itV].begin(), m_augmented[*itV].end());
    */
#ifdef DEBUGGGA
    if (!ancHeur->getComputedMessages()[*itV]->isAccurate()) {
    for (vector<pair<int,Function*> >::iterator itF=funs.begin(); itF!=funs.end(); ++itF) {
        cout << ' ' << itF->first << endl;
      cout << ' ' << (*(itF->second)) << endl;
    }
    cout << endl;
    }
    /*
    cout << assignment << endl;
    for (vector<pair<int,Function*> >::iterator itF=funs.begin(); itF!=funs.end(); ++itF)
      cout << ' ' << (*(itF->second)->substitute(assignment));
    cout << endl;
    */
#endif

// test
 /*
    if (funs.size() == 0) {
        cout << "empty bucket" << endl;
        continue;
    }
    */
// ===


    // compute global upper bound for root (dummy) bucket
    if (*itV == elimOrder[0]) {// variable is dummy root variable
      if (computeTables) { // compute upper bound if assignment is given
        double subproblemUB = ELEM_ONE;
        for (vector<pair<int,Function*> >::iterator itF=funs.begin(); itF!=funs.end(); ++itF)
          subproblemUB OP_TIMESEQ (itF->second)->getValue(vAssn);
//        cout << "    MBE-ALL  = " << SCALE_LOG(m_globalUB) << " (" << SCALE_NORM(m_globalUB) << ")" << endl;
        subproblemUB OP_DIVIDEEQ m_problem->globalConstInfo();  // for backwards compatibility of output
//        cout << "    MBE-ROOT = " << SCALE_LOG(m_globalUB) << " (" << SCALE_NORM(m_globalUB) << ")" << endl;
        curHeur->setUpperBound(subproblemUB);
      }
      continue; // skip the dummy variable's bucket
    }


    // sort functions by decreasing scope size
    sort(funs.begin(), funs.end(), scopeIsLargerIF);

    // partition functions into minibuckets
    vector<MiniBucket> minibuckets;
//    vector<Function*>::iterator itF; bool placed;
    for (vector<pair<int,Function*> >::iterator itF = funs.begin(); itF!=funs.end(); ++itF) {
//    while (funs.size()) {
      bool placed = false;
      int bcount = 0;
      for (vector<MiniBucket>::iterator itB=minibuckets.begin();
            !placed && itB!=minibuckets.end(); ++itB, ++bcount)
      {
          //cout << "Trying bucket " << bcount << endl;
        if (itB->allowsFunction(itF->second,context)) { // checks if function fits into bucket
          itB->addFunction(itF->second);
          //cout << itF->first << " " << itF->second->getArity() << endl;
          /*
          if (!ancHeur->getComputedMessages()[*itV]->isAccurate()) {
              cout << itF->first << endl;
              cout << *(itF->second) << endl;
              cout << "(placed in bucket " << bcount << ") " << endl;
          }
          */
          placed = true;
        }
      }
      if (!placed) { // no fit, need to create new bucket
        MiniBucket mb(*itV,m_ibound,m_problem);
        mb.addFunction(itF->second);
          //cout << itF->first << " " << itF->second->getArity() << endl;
          /*
          if (!ancHeur->getComputedMessages()[*itV]->isAccurate()) {
              cout << itF->first << endl;
              cout << *(itF->second) << endl;
              cout << "(placed in bucket " << bcount << ") " << endl;
          }
          */
          
        minibuckets.push_back(mb);
      }
    }

    // minibuckets for current bucket are now ready, process each
    // and place resulting function
    
    // messages are accurate if no partitioning and incoming messages were also accurate
    bool accurateHeur = accurateHeurIn[*itV] && minibuckets.size() == 1;

    // Compute max-marginals for each bucket

    vector<Function*> maxMarginals;
    double *avgMMTable; 
    Function *avgMaxMarginal = NULL;
    
    if (m_options->match && minibuckets.size() > 1) {
        // Find intersection of scopes (the scope of all max-marginals)
        vector<MiniBucket>::iterator itB=minibuckets.begin();
        set<int> intersectScope = setminus(itB->getJointScope(), context);
        itB++;
        for (; itB!=minibuckets.end(); ++itB)
        {
            intersectScope = intersection(intersectScope, 
                    setminus(itB->getJointScope(), context));
            //intersectScope = newInter;
        }
        for (vector<MiniBucket>::iterator itB=minibuckets.begin();
                itB!=minibuckets.end(); ++itB)
        {
            set<int> elimVars = setminus(setminus(itB->getJointScope(),context), intersectScope);
            maxMarginals.push_back(itB->conditionEliminate(computeTables, assignment, elimVars));
        }

        // Find average max-marginal (geometric mean)
        size_t tablesize = 1;
        for (set<int>::iterator sit=intersectScope.begin(); 
                sit!=intersectScope.end(); ++sit) {
            tablesize *= m_problem->getDomainSize(*sit);
        }

        avgMMTable = new double[tablesize];
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
    
        int ammid = 0;
        avgMaxMarginal = new FunctionBayes(ammid, m_problem, intersectScope, avgMMTable, tablesize);
    }

    int bucketIdx = 0;

    computedMessages[*itV] = new ConditionedMessages(curHeur, accurateHeur);

    for (vector<MiniBucket>::iterator itB=minibuckets.begin();
          itB!=minibuckets.end(); ++itB, ++bucketIdx)
    {

      // Replace this to generate moment-matched version if #minibuckets > 1
      Function* newf;
      if (!m_options->match || minibuckets.size() <= 1)
          newf = itB->conditionEliminate(computeTables,assignment); // process the minibucket
      else {
          newf = itB->conditionEliminateMM(computeTables,assignment,
                  maxMarginals[bucketIdx],avgMaxMarginal); // process the minibucket
      }

      const set<int>& newscope = newf->getScopeSet();
      memSize += newf->getTableSize();

      /*
      */
      /*
      if (!ancHeur->getComputedMessages()[*itV]->isAccurate()) {
          cout << assignment << endl;
          cout << " " << *newf << endl;
          for (unsigned i = 0; i < newf->getTableSize(); ++i) {
              cout << "  " << newf->getTable()[i] << endl;
          }
      }
      */


      // go up in tree to find target bucket
      vector<int> *path = new vector<int>();

      PseudotreeNode* n = m_pseudotree->getNode(*itV)->getParent();
      while (newscope.find(n->getVar()) == newscope.end() && n != m_pseudotree->getRoot() && n->getVar() != var) {
        path->push_back(n->getVar());
        n = n->getParent();
      }
      // matching bucket found OR root of pseudo tree reached
      path->push_back(n->getVar());
      curHeur->getNumValues()[n->getVar()] += newf->getTableSize();
      curHeur->getNumZeroValues()[n->getVar()] += newf->getTableSize() - newf->getTightness();
      forceUpdate[n->getVar()] = true;
      computedMessages[*itV]->addFunction(newf, path);
    }
    curHeur->populateMessages(*itV,visited);
    // all minibuckets processed and resulting functions placed
    //

    // free up memory used by max-marginals
    if (m_options->match && minibuckets.size() > 1) {
      for (unsigned int i = 0; i < maxMarginals.size(); ++i)
        delete maxMarginals[i];
      maxMarginals.clear();
      if(avgMaxMarginal) delete avgMaxMarginal;
    }

  }

#ifdef DEBUG
  vector<vector<Function*> > &intermediate = curHeur->getIntermediate();
  // output augmented and intermediate buckets
  if (computeTables)
    for (int i=0; i<m_problem->getN(); ++i) {
      cout << "$ AUG" << i << ": " << augmented[i] << " + " << intermediate[i] << endl;
    }
#endif

  /*
  // clean up for estimation mode
  if (!computeTables) {
    for (vector<vector<Function*> >::iterator itA = m_augmented.begin(); itA!=m_augmented.end(); ++itA)
      for (vector<Function*>::iterator itB = itA->begin(); itB!=itA->end(); ++itB)
        delete *itB;
    m_augmented.clear();
    m_augmented.clear();
  }
  */
  time(&heurCompEnd);
  m_heurCompTime += difftime(heurCompEnd, heurCompStart);
  m_numHeuristics++;

  /*
  cout << "Heuristic computed. " << endl;
  cout << "var,depth: " << var << "," << m_pseudotree->getNode(var)->getDepth() << endl;
  for (int i = 0; i < m_problem->getN(); ++i) {
      cout << i << "," << m_pseudotree->getNode(i)->getDepth() << "," << curHeur->getNumValues()[i] << "," << curHeur->getNumZeroValues()[i] << "," << curHeur->getConstrainedness(i) << endl;
  }
  cout << endl;
  */

  return memSize;
}

void MiniBucketElim::simulateBuildSubproblem(int var, const vector<int> &elimOrder) {

  //if (m_pseudotree->getRoot()->getVar() == var) return 0;

#ifdef DEBUGSIM
  cout << "$ Building MBE(" << m_ibound << ")" << endl;
  cout << "simulating at " << var << endl;
#endif


  /*
  cout << "BuildSubproblem" << endl;
  cout << "var: " << var << endl;
  cout << "assign: " << assignment << endl;
  */
  //assert(!m_accurateHeuristicIn[var]); // should only be here if heuristic is not accurate

//  this->reset();


  //vector<bool> visited(m_problem->getN(), false);
  
  // sets of scopes for each variable
  vector<vector<Scope*> > augmentedSim(m_problem->getN());
  set<int> subproblemVars(elimOrder.begin(),elimOrder.end());
  subproblemVars.erase(var);

  // ITERATES OVER BUCKETS, FROM LEAVES TO ROOT
  for (vector<int>::const_reverse_iterator itV=elimOrder.rbegin(); itV!=elimOrder.rend(); ++itV) {

    if (*itV == var) {
        /*
        cout << "will skip var "  << var << endl;
        cout << "assign " << assignment << endl;
        */
        continue;
    }
#ifdef DEBUGSIM
    cout << "$ Bucket for variable " << *itV << endl;
#endif

    // collect relevant scopes in funs
    vector<Scope*> funs;
    const vector<Function*>& fnlist = m_pseudotree->getFunctions(*itV);
    for(vector<Function*>::const_iterator itF=fnlist.begin(); itF!=fnlist.end(); ++itF) {
        Scope *s = new Scope((*itF)->getId(), 
                intersection(subproblemVars, (*itF)->getScopeSet()));
        funs.push_back(s);
    }
    for(vector<Scope*>::const_iterator itF=augmentedSim[*itV].begin(); 
            itF!=augmentedSim[*itV].end(); ++itF) {
        Scope *s = new Scope((*itF)->getId(), 
                intersection(subproblemVars, (*itF)->getScope()));
        funs.push_back(s); 
    }

    // account for "empty" minibuckets
    if (funs.size() == 0) {
        m_mbCountSubtree[var][m_pseudotree->getNode(*itV)->getParent()->getVar()] += 
            m_mbCountSubtree[var][*itV];
        //cntMiniBuckets++;
        continue;
    }
// ===


    // compute global upper bound for root (dummy) bucket
    if (*itV == elimOrder[0]) {// variable is dummy root variable
      continue; // skip the dummy variable's bucket
    }


    // sort functions by decreasing scope size
    sort(funs.begin(), funs.end(), scopeIsLargerS);

    // partition functions into minibuckets
    //vector<MiniBucket> minibuckets;
    vector<Scope*> minibuckets;
//    vector<Function*>::iterator itF; bool placed;

    set<int> intersectionScope;
    for (vector<Scope*>::iterator itF = funs.begin(); itF!=funs.end(); ++itF) {
      intersectionScope.insert((*itF)->getScope().begin(), (*itF)->getScope().end());
      bool placed = false;
      for (vector<Scope*>::iterator itB=minibuckets.begin();
            !placed && itB!=minibuckets.end(); ++itB)
      {
        set<int> &a=(*itB)->getScope(), b=(*itF)->getScope();
        set<int>::iterator ita=a.begin(), itb=b.begin();
        int s=0,d=0;

        while (ita != a.end() && itb != b.end()) {
            d = *ita - *itb;
            if(d>0) {
                ++s; ++itb;
            } else if (d<0) {
                ++s; ++ita;
            } else {
                ++s; ++ita; ++itb;
            }
        }
        while (ita != a.end()) {
            ++s; ++ita;
        }
        while (itb != b.end()) {
            ++s; ++itb;
        }
        if (s == (int) a.size() || s <= m_ibound+1) {// checks if function fits into bucket
            a.insert(b.begin(),b.end()); 
            placed = true;
        }
      }
      if (!placed) { // no fit, need to create new bucket
        minibuckets.push_back(new Scope(-(*itV), (*itF)->getScope()));
      }
    }
    m_subproblemWidth[var] = max(m_subproblemWidth[var], int(intersectionScope.size() - 1));
    for (unsigned i = 0; i < funs.size(); ++i) {
        delete funs[i];
    }
    funs.clear();

    m_dupeCount[var][*itV] = minibuckets.size();
    

    for (vector<Scope*>::iterator itB=minibuckets.begin();
          itB!=minibuckets.end(); ++itB)
    {

      (*itB)->erase(*itV);
        
      const set<int> &newscope = (*itB)->getScope();

      // go up in tree to find target bucket

      PseudotreeNode* n = m_pseudotree->getNode(*itV)->getParent();
      while (newscope.find(n->getVar()) == newscope.end() && n != m_pseudotree->getRoot() && 
              n->getVar() != var) {
        n = n->getParent();
      }
      // matching bucket found OR root of pseudo tree reached
      augmentedSim[n->getVar()].push_back(*itB);

    }
    m_mbCountSubtree[var][m_pseudotree->getNode(*itV)->getParent()->getVar()] += 
        (m_mbCountSubtree[var][*itV] + minibuckets.size());
    // all minibuckets processed and resulting functions placed
    minibuckets.clear();

  }
  for (int i = 0; i < m_problem->getN(); ++i) {
      for (int j = 0; j < int(augmentedSim[i].size()); ++j) {
          delete augmentedSim[i][j];
      }
  }

  //return m_mbCountSubtree[var][var];
}

// reparameterize problem using mplp on the original factors
bool MiniBucketElim::doFGLP() {
  bool changedFunctions = false;

  if (m_options!=NULL && (m_options->mplp > 0 || m_options->mplps > 0)) {
    mex::mplp _mplp( copyFactors() );  // copyFactors()
    _mplp.setProperties("Schedule=Fixed,Update=Var,StopIter=100,StopObj=-1,StopMsg=-1,StopTime=-1");
    _mplp.init();
    char opt[50];
    if (m_options->mplp > 0)  { sprintf(opt,"StopIter=%d",m_options->mplp); _mplp.setProperties(opt); }
    if (m_options->mplps > 0) { sprintf(opt,"StopTime=%f",m_options->mplps); _mplp.setProperties(opt); }

    _mplp.run();
    rewriteFactors( _mplp.beliefs() );
    //_mbe.setModel( _mplp.beliefs() );

    changedFunctions = true;
  }

  return changedFunctions;
}

bool MiniBucketElim::doJGLP() {
    return false;
}

// Copy DaoOpt Function class into mex::Factor class structures
mex::vector<mex::Factor> MiniBucketElim::copyFactors( void ) {
  mex::vector<mex::Factor> fs(m_problem->getC());
  for (int i=0;i<m_problem->getC(); ++i) fs[i] = m_problem->getFunctions()[i]->asFactor().exp();
  return fs;
}

// Mini-bucket may have re-parameterized the original functions; if so, replace them
void MiniBucketElim::rewriteFactors( const vector<mex::Factor>& factors) {
  vector<Function*> newFunctions(factors.size()); // to hold replacement, reparameterized functions
  for (size_t f=0;f<factors.size();++f) {         // allocate memory, copy variables into std::set
    double* tablePtr = new double[ factors[f].nrStates() ];
    std::set<int> scope;
    for (mex::VarSet::const_iterator v=factors[f].vars().begin(); v!=factors[f].vars().end(); ++v)
      scope.insert(v->label());
    newFunctions[f] = new FunctionBayes(f,m_problem,scope,tablePtr,factors[f].nrStates());
    newFunctions[f]->fromFactor( log(factors[f]) );    // write in log factor functions
  }
  m_problem->replaceFunctions( newFunctions );                // replace them in the problem definition
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

/* finds a bfs order of the pseudotree */
void MiniBucketElim::findBfsOrder(vector<int>& order) const {
  order.clear();
  queue<PseudotreeNode*> bfs;
  bfs.push(m_pseudotree->getRoot());
  PseudotreeNode* n = NULL;
  while (!bfs.empty()) {
    n = bfs.front();
    bfs.pop(); 
    order.push_back(n->getVar());
    for (vector<PseudotreeNode*>::const_iterator it=n->getChildren().begin();
          it!=n->getChildren().end(); ++it) {
      bfs.push(*it);
    }
  }
}

size_t MiniBucketElim::limitSize(size_t memlimit, const vector<val_t> * assignment) {
    if (m_options->dynamic) return 0;
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
  const vector<vector<Function*> > &augmented = m_rootHeurInstance->getAugmented();
  for (vector<vector<Function*> >::const_iterator it=augmented.begin(); it!= augmented.end(); ++it) {
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
  vector <vector<Function*> > &augmented = m_rootHeurInstance->getAugmented();
  size_t sz = augmented.size();
  out.write((char*)&( sz ), sizeof( sz ));

  // i-bound
  out.write((char*)&( m_ibound ), sizeof( m_ibound ));

  // global UB
  out.write((char*)&( m_globalUB ), sizeof( m_globalUB ));

  map<const Function*,size_t> funcMap;

  // over m_augmented
  for (size_t i=0; i<sz; ++i) {
    size_t sz2 = augmented[i].size();
    out.write((char*)&( sz2 ), sizeof( sz2 ));

    vector<Function*>::const_iterator itF = augmented[i].begin();
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
  vector <vector<Function*> > &intermediate = m_rootHeurInstance->getIntermediate();

  // over m_intermediate
  for (size_t i=0; i<sz; ++i) {
    size_t sz2 = intermediate[i].size();
    out.write((char*)&( sz2 ), sizeof( sz2 ));

    vector<Function*>::const_iterator itF = intermediate[i].begin();
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

  vector <vector<Function*> > &augmented = m_rootHeurInstance->getAugmented();
  vector <vector<Function*> > &intermediate = m_rootHeurInstance->getIntermediate();

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
      augmented[i].push_back(f);
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
      intermediate[i].push_back(allFuncs.at(y));
    }
  }

  in.close();
  cout << "Read mini bucket with i-bound " << ibound << " from file " << fn << endl;
  return true;
}


