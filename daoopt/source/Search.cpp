/*
 * Search.cpp
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
 *  Created on: Sep 22, 2009
 *      Author: Lars Otten <lotten@ics.uci.edu>
 */

#undef DEBUG

#include "Search.h"
#include "ProgramOptions.h"

#include <chrono>
#include <gflags/gflags.h>
#include <iomanip>

//#define DECOMPOSE_H_INTO_INDEPENDENT_SUBPROBLEMS
// Note: We will later refactor AOBB specific parts into BranchAndBound.
// i.e. caching

using namespace std::chrono;
using namespace minisat;


// NONE, UNIT, FULL (not yet implemented).
DECLARE_string(cp_type);

namespace daoopt {

extern high_resolution_clock::time_point _time_start; // from Main.cpp

Search::Search(Problem* prob, Pseudotree* pt, SearchSpace* s, Heuristic* h,
    BoundPropagator* prop, ProgramOptions* po) :
    m_problem(prob), m_pseudotree(pt), m_space(s), m_heuristic(h),
    m_prop(prop),
    m_options(po), m_foundFirstPartialSolution(false)
#ifdef PARALLEL_DYNAMIC
  , m_nextSubprob(NULL)
#endif
{
  // initialize the array for counting nodes per level
  m_nodeProfile.resize(m_pseudotree->getHeight()+1, 0);
  // initialize the array for counting leaves per level
  m_leafProfile.resize(m_pseudotree->getHeight()+1, 0);

  // initialize the local assignment vector for BaB
  m_assignment.resize(m_problem->getN(),NONE);

  // Preallocate space for expansion vector. 128 should be plenty.
  m_expand.reserve(128);
}


SearchNode* Search::initSearch() {
  DoConstraintPropagation();
  if (!m_space->root) {
    // Create root OR node (dummy variable)
    PseudotreeNode* ptroot = m_pseudotree->getRoot();
    SearchNode* node = new SearchNodeOR(NULL, ptroot->getVar(), -1);
    m_space->root = node;
    return node;
  } else {
    if (m_space->root->getChildCountAct())  // TODO: not needed anymore?
      return m_space->root->getChildren()[0];
    else
      return m_space->root;
  }
}


#ifndef NO_HEURISTIC
void Search::finalizeHeuristic() {
  assert(m_space && m_space->getTrueRoot());
  assignCostsOR(m_space->getTrueRoot());
}
#endif


bool Search::doProcess(SearchNode* node) {
  assert(node);
  if (node->getType() == NODE_AND) {
    m_space->stats.numProcAND += 1;
    m_space->stats.numProcANDVar[node->getVar()] += 1;
    // 0-labeled nodes should not get generated in the first place
    assert(node->getLabel() != ELEM_ZERO);  // if this fires, something is wrong!
    DIAG( ostringstream ss; ss << "[" << m_space->stats.numProcAND << "] " << *node << " (l=" << node->getLabel() << ")\n"; myprint(ss.str()) );
    int var = node->getVar();
    int val = node->getVal();
    m_assignment[var] = val; // record assignment
  } else { // NODE_OR
    m_space->stats.numProcOR += 1;
    m_space->stats.numProcORVar[node->getVar()] += 1;
    DIAG( ostringstream ss; ss << "[" << m_space->stats.numProcOR << "] " << *node << "\n"; myprint(ss.str()) );
  }
  return false; // default
}


bool Search::doCaching(SearchNode* node) {

  assert(node);
  int var = node->getVar();
  PseudotreeNode* ptnode = m_pseudotree->getNode(var);

  if (node->getType() == NODE_AND) { // AND node -> reset associated adaptive cache tables

    const list<int>& resetList = ptnode->getCacheReset();
    for (list<int>::const_iterator it=resetList.begin(); it!=resetList.end(); ++it)
      m_space->cache->reset(*it);

  } else { // OR node, try actual caching
if (!ptnode->getParent() || !ptnode->getParent()->getParent())
      return false;  // pseudo tree root or one of its direct children

    if (ptnode->getFullContextVec().size() <= ptnode->getParent()->getFullContextVec().size()) {
      // add cache context information
      addCacheContext(node,ptnode->getCacheContextVec());
      //DIAG( myprint( str("    Context set: ") + ptnode->getCacheContext() + "\n" ) );
      // try to get value from cache
#ifndef NO_ASSIGNMENT
      const pair<const double, const vector<val_t> >& entry =
          m_space->cache->read(var, node->getCacheInst(), node->getCacheContext());
      if (!ISNAN(entry.first)) {
        node->setValue( entry.first ); // set value
        node->setOptAssig( entry.second ); // set assignment
#else
      const double entry = m_space->cache->read(var, node->getCacheInst(), node->getCacheContext());
      if (!ISNAN(entry)) {
        node->setValue( entry ); // set value
#endif
        node->setLeaf(); // mark as leaf
#ifdef DEBUG
        ostringstream ss;
        ss << "-Read " << *node << " with value " << node->getValue()
#ifndef NO_ASSIGNMENT
        << " and opt. solution " << node->getOptAssig()
#endif
        << endl;
        myprint(ss.str());
#endif
        return true;
	  } else {
        node->setCachable(); // mark for caching later
      }
    }

  } // if on node type

  return false; // default, no caching applied

} // Search::doCaching


bool Search::doPruning(SearchNode* node) {

#ifdef NO_HEURISTIC
  return false;
#endif

  assert(node);
  int var = node->getVar();
  PseudotreeNode* ptnode = m_pseudotree->getNode(var);
  int depth = ptnode->getDepth();

  if (canBePruned(node)) {
    DIAG( myprint("\t !pruning \n") );
    node->setLeaf();
    m_space->stats.numPruned += 1;
    node->setPruned();
    if (node->getType() == NODE_AND) {
      // count 1 leaf AND node
      if (depth >= 0) m_leafProfile.at(depth) += 1;
#if defined PARALLEL_DYNAMIC
      node->setSubLeaves(1);
#endif
      node->setValue(ELEM_ZERO);
    } else { // NODE_OR
      if ( ISNAN(node->getValue()) ) // value could be set by LDS
        node->setValue(ELEM_ZERO);
      // assume all AND children would have been created and pruned
      if (depth >= 0) m_leafProfile.at(depth) += m_problem->getDomainSize(var);
#if defined PARALLEL_DYNAMIC
      node->addSubLeaves(m_problem->getDomainSize(var));
#endif
    }
    return true;
  }

  return false; // default false

} // Search::doPruning



SearchNode* Search::nextLeaf() {

  SearchNode* node = this->nextNode();
  while (node) {
    // dependant on derived class how nodes are processed
    if (doCompleteProcessing(node)) {
      return node;
    }

    if (m_options->prop_heuristic &&
        propHeuristic(node)) { // early termination?
      return nullptr;
    }

    // Then move to next node.
    node = this->nextNode();
    high_resolution_clock::time_point time_now = high_resolution_clock::now();
    double time_elapsed = duration_cast<duration<double>>(time_now - _time_start).count();
    if (time_elapsed > m_options->maxTime) {
      cout << "Timed out at " << time_elapsed << " seconds." << endl;
      cout << "Stats at timeout: " << endl;
      cout << "================= " << endl;
      cout << "OR nodes:           " << m_space->stats.numExpOR << endl;
      cout << "AND nodes:          " << m_space->stats.numExpAND << endl;
      cout << "OR processed:       " << m_space->stats.numProcOR << endl;
      cout << "AND processed:      " << m_space->stats.numProcAND << endl;
      cout << "Leaf nodes:         " << m_space->stats.numLeaf << endl;
      cout << "Pruned nodes:       " << m_space->stats.numPruned << endl;
      cout << "Deadend nodes:      " << m_space->stats.numDead << endl;
      cout << "Deadend nodes (CP): " << m_space->stats.numDeadCP << endl;
      exit(0);
    }
  }
  return NULL;
}


bool Search::canBePruned(SearchNode* n) {
  DIAG(oss ss; ss << std::setprecision(20) << "\tcanBePruned(" << *n << ")" << " h=" << n->getHeur() << endl; myprint(ss.str());)


  n->PruningGap() = DBL_MAX ;

  // never prune the root node (is there a better solution maybe?)
  if (n->getDepth() < 0) return false;

  // use constraint propagation
  if (n->getType() == NODE_AND) {
    if (FLAGS_cp_type == "FULL") {
      vector<int> sat_vars;

      SearchNode* cur = n;
      SearchNode* parent = cur->getParent();
      for (;;) {
        int vvar = cur->getVar();
        int vval = cur->getVal();
        sat_vars.push_back(var2sat_[vvar][vval]);

        cur = parent->getParent();
        if (cur) break;
        parent = cur->getParent();
      }

      if (!DoSATPropagate(sat_vars)) {
        ++m_space->stats.numDeadCP;
        return true;
      }
    }
    else if (FLAGS_cp_type == "UNIT") {
      zchaff_solver_.dlevel()++;
      if (!DoCPLookahead(n->getVar(), n->getVal(), n->changes(),
          m_pseudotree->getNode(n->getVar())->getSubprobVars())) {
        ++m_space->stats.numDeadCP;
        return true;
      }
    }
  }

      // heuristic is upper bound, prune if zero
  if (n->getHeur() == ELEM_ZERO)
	  { n->PruningGap() = ELEM_ZERO ; return true; }

  double curPSTVal = n->getHeur();  // includes label in case of AND node
  SearchNode* curOR = (n->getType() == NODE_OR) ? n : n->getParent();

  double pg = curPSTVal - curOR->getValue() ;
  n->PruningGap() = pg ;
  if (pg <= 0.0)  // simple pruning case
	  return true;

  SearchNode* curAND = NULL;

  while (curOR->getParent()) {  // climb all the way up to root node, if we have to

    curAND = curOR->getParent();  // climb up one, update values
    curPSTVal OP_TIMESEQ curAND->getLabel();  // collect AND node label
    DIAG ( ostringstream ss; ss << *curAND << endl; myprint(ss.str()); )
    curPSTVal OP_TIMESEQ curAND->getSubSolved();  // incorporate already solved sibling OR nodes

    // incorporate new not-yet-solved sibling OR nodes through their heuristic
    NodeP* children = curAND->getChildren();
    for (size_t i = 0; i < curAND->getChildCountFull(); ++i) {
      if (!children[i] || children[i] == curOR) continue;
      else curPSTVal OP_TIMESEQ children[i]->getHeur();
    }
    curOR = curAND->getParent();

    DIAG( ostringstream ss; ss << std::setprecision(20) << "\t ?PST root: " << *curOR << " pst=" << curPSTVal << " v=" << curOR->getValue() << endl; myprint(ss.str()); )

    //if ( fpLt(curPSTVal, curOR->getValue()) ) {
	pg = curPSTVal - curOR->getValue() ;
	if (pg < n->PruningGap()) n->PruningGap() = pg ;
    if ( pg <= 0.0 /*|| fabs(curPSTVal - curOR->getValue()) < 1e-10*/ ) {
        for (SearchNode* nn = (n->getType() == NODE_OR) ? n : n->getParent();
            nn != curOR; nn = nn->getParent()->getParent())
            nn->setNotOpt();  // mark possibly not optimally solved subproblems
      return true;  // pruning is possible!
    }
  }
  return false;  // default, no pruning possible

} // Search::canBePruned


void Search::syncAssignment(const SearchNode* node) {
  // only accept OR nodes
  assert(node && node->getType()==NODE_OR);

  while (node->getParent()) {
    node = node->getParent(); // AND node
    m_assignment.at(node->getVar()) = node->getVal();
    node = node->getParent(); // OR node
  }
} // Search::syncAssignment


bool Search::generateChildrenAND(SearchNode* n, vector<SearchNode*>& chi) {

  assert(n && n->getType() == NODE_AND);

  if (n->getChildren()) {  // node previously expanded
    if (n->getChildCountAct() == 0) {
      n->clearChildren();  // previously expanded, but all children deleted
    } else {
      for (size_t i = 0; i < n->getChildCountFull(); ++i) {
        if (n->getChildren()[i])
          chi.push_back(n->getChildren()[i]);
      }
      return false;
    }
  }

  m_space->stats.numExpAND += 1;
  m_space->stats.numANDVar[n->getVar()] += 1;

  int var = n->getVar();
  PseudotreeNode* ptnode = m_pseudotree->getNode(var);
  int depth = n->getDepth();

#if defined PARALLEL_DYNAMIC || defined PARALLEL_STATIC
  n->setSubCount(1);
#endif
  if (depth>=0) m_nodeProfile[depth] +=1; // ignores dummy node

  // we have stored heuristic from each child
  const PseudotreeNode *node = m_pseudotree->getNode(var) ;
  int nChildren = node->getChildren().size() ;

  // create new OR children (going in reverse due to reversal on stack)
  int idxChild = nChildren - 1 ;
  for (vector<PseudotreeNode*>::const_reverse_iterator it=ptnode->getChildren().rbegin(); it!=ptnode->getChildren().rend(); ++it, --idxChild) {
	// _heurValueForEachIndSubproblem[idxChild] is the portion of this node's h that we got from that child.
    int vChild = (*it)->getVar();
    SearchNodeOR* c = new SearchNodeOR(n, vChild, depth+1);
    chi.push_back(c);
#ifdef DECOMPOSE_H_INTO_INDEPENDENT_SUBPROBLEMS
	c->heurValueOfParentFromThisNode() = ((SearchNodeAND*)n)->getHeurValueForEachIndSubproblem(idxChild) ;
#endif // DECOMPOSE_H_INTO_INDEPENDENT_SUBPROBLEMS
#ifndef NO_HEURISTIC
    // Compute and set heuristic estimate, includes child labels
    if (assignCostsOR(c) == ELEM_ZERO) {  // dead end, clean up and exit
      for (vector<NodeP>::iterator it=chi.begin(); it!=chi.end(); ++it)
        delete (*it);
      chi.clear();
      n->setLeaf();
      n->setValue(ELEM_ZERO);
      m_space->stats.numLeaf += 1;
      if (depth!=-1) m_leafProfile[depth] += 1;
#if defined PARALLEL_DYNAMIC
      n->setSubLeaves(1);
#endif
      return true;
    }
#endif

#ifndef NO_HEURISTIC
    // store initial lower bound on subproblem (needed for master init.)
#if defined PARALLEL_DYNAMIC //|| defined PARALLEL_STATIC
    c->setInitialBound(lowerBound(c));
#endif
#endif

  } // for loop over new OR children

  if (chi.empty()) {
    n->setLeaf(); // terminal node
    n->setValue(ELEM_ONE);
    m_space->stats.numLeaf += 1;
    if (depth!=-1) m_leafProfile[depth] += 1;
#if defined PARALLEL_DYNAMIC
    n->setSubLeaves(1);
#endif
    return true; // no children
  }

  // order subproblems by their heuristic (if applicable)
  // (inverse due to stack reversal)
  if (m_space->options->subprobOrder == SUBPROB_HEUR_INC)
    sort(chi.rbegin(), chi.rend(), SearchNode::heurLess);
  else if (m_space->options->subprobOrder == SUBPROB_HEUR_DEC)
    sort(chi.begin(), chi.end(), SearchNode::heurLess);

  n->addChildren(chi);

  return false; // default

} // Search::generateChildrenAND


bool Search::generateChildrenOR(SearchNode* n, vector<SearchNode*>& chi) {

  assert (n && n->getType() == NODE_OR);

  if (n->getChildren()) {  // node previously expanded
    if (n->getChildCountAct() == 0) {
      n->clearChildren();  // previously expanded, but all children deleted
    } else {
      for (size_t i = 0; i < n->getChildCountFull(); ++i) {
        if (n->getChildren()[i])
          chi.push_back(n->getChildren()[i]);
      }
      return false;
    }
  }

  m_space->stats.numExpOR += 1;
  m_space->stats.numORVar[n->getVar()] += 1;

  int var = n->getVar();
  int depth = n->getDepth();
//  PseudotreeNode* ptnode = m_pseudotree->getNode(var);

#ifndef NO_HEURISTIC
  // retrieve precomputed labels and heuristic values
  double* heur = n->getHeurCache();
  double* ordering_heur = n->getOrderingHeurCache();
#endif

  val_t vDomain = m_problem->getDomainSize(var) ;
  int subprobH_ = 2*vDomain ;
  const PseudotreeNode *node = m_pseudotree->getNode(var) ;
  int nChildren = node->getChildren().size() ;

  for (val_t i = vDomain - 1 ; i >= 0 ; --i) {
#ifdef NO_HEURISTIC
    // compute label value for new child node
    m_assignment[var] = i;


    double d = ELEM_ONE;
    const list<Function*>& funs = m_pseudotree->getFunctions(var);
    for (list<Function*>::const_iterator it = funs.begin(); it != funs.end(); ++it)
      d OP_TIMESEQ (*it)->getValue(m_assignment);
    if (d == ELEM_ZERO) {
      m_space->nodesLeaf += 1;
      m_leafProfile[depth] += 1;
#if defined PARALLEL_DYNAMIC
      n->addSubLeaves(1);
#endif
      continue; // label=0 -> skip
    }

    SearchNodeAND* c = new SearchNodeAND(n, i, d);
#else
    // early pruning if heuristic is zero (since it's an upper bound)
    if (heur[2*i+1] == ELEM_ZERO) { // 2*i=heuristic, 2*i+1=label
      m_space->stats.numDead += 1;
      m_leafProfile[depth] += 1;
#if defined PARALLEL_DYNAMIC
      n->addSubLeaves(1);
#endif
      continue;
    }

    // Found inconsistent by constriant propagation
    if (!current_domains_.empty() && !current_domains_[var][i]) {
      continue;
    }
    SearchNodeAND* c = new SearchNodeAND(n, i, heur[2*i+1]); // uses cached label
    // set cached heur. value
    c->setHeur( heur[2*i] );
    c->setOrderingHeur( ordering_heur[i] );
#ifdef DECOMPOSE_H_INTO_INDEPENDENT_SUBPROBLEMS
	c->setHeurValueForEachIndSubproblem(nChildren > 0 ? heur + (subprobH_ + i*nChildren) : NULL) ;
#endif // DECOMPOSE_H_INTO_INDEPENDENT_SUBPROBLEMS
#endif
    chi.push_back(c);
  }

#ifndef NO_HEURISTIC
//  n->clearHeurCache();
#endif

  if (chi.empty()) {  // should never happen
    n->setLeaf();
    n->setValue(ELEM_ZERO);
    return true; // no children
  }

#ifndef NO_HEURISTIC
  // sort new nodes by increasing heuristic value
  sort(chi.begin(), chi.end(), SearchNode::heurLess);
#endif
#ifdef DEBUG
// =================
    for (int i = 0;i < chi.size(); ++i) {
        cout << setprecision(20);
        cout << *chi[i] << " " << chi[i]->getHeur() << " " << chi[i]->getLabel() << endl;
    }
// =================
#endif
  n->addChildren(chi);

  return false; // default

} // Search::generateChildrenOR


/* define the following to enable fetching of function values in bulk */
#define GET_VALUE_BULK
double Search::assignCostsOR(SearchNode* n)
{
  int v = n->getVar();
  int old_value = m_assignment[v];
  const PseudotreeNode *node = m_pseudotree->getNode(v) ;
  int nChildren = node->getChildren().size() ;
  int vDomain = m_problem->getDomainSize(v);
  // in dv, [0 ... 2*vDomain) are label+heuristic/label values; [2*vDomain ... 2*vDomain + vDomain*nChildren) are heuristic values, for each value of this var and each child.
#ifdef DECOMPOSE_H_INTO_INDEPENDENT_SUBPROBLEMS
  double* dv = new double[vDomain*2 + vDomain*nChildren];
#else
  double* dv = new double[vDomain*2];
#endif // DECOMPOSE_H_INTO_INDEPENDENT_SUBPROBLEMS
  double* ordering_cache = new double[vDomain];

  // Set first so computed values can be used right away
  n->setHeurCache(dv);
  n->setOrderingHeurCache(ordering_cache);

  for (int i = 0; i < vDomain; ++i) dv[2*i+1] = ELEM_ONE;
  for (int i = 0; i < vDomain; ++i) ordering_cache[i] = 0;
  double h = ELEM_ZERO; // the new OR nodes h value
  int argmax_h = 0; // Initialize the first as max since they may all be zero.
  m_heuristic->noteOrNodeExpansionBeginning(v, m_assignment, n) ;

#ifdef GET_VALUE_BULK
  // Need to get correct costs for function after shifting
  // need to request from heuristic class instead
  vector<double> labelAll ;
  labelAll.clear();
  labelAll.resize(vDomain, ELEM_ONE);
  m_heuristic->getLabelAll(v, m_assignment, n, labelAll);

  m_costTmp.clear();
  m_costTmp.resize(vDomain, ELEM_ONE);
  m_heuristic->getHeurAll(v, m_assignment, n, m_costTmp);
  for (int i=0; i<vDomain; ++i) {
    dv[2*i] = m_costTmp[i];
    dv[2*i+1] = labelAll[i];
  }

  for (int i=0; i<vDomain; ++i) {
    m_assignment[v] = i;
    dv[2*i] = dv[2*i+1] OP_TIMES dv[2*i];
    if (dv[2*i] > h) {
      argmax_h = i;
      h = dv[2*i]; // keep max. for OR node heuristic
    }
    ordering_cache[i] = m_heuristic->getOrderingHeur(v, m_assignment, n);
  }
#else
  double label, heuristic;
  double ordering_heuristic;
  int j ;
#ifdef DECOMPOSE_H_INTO_INDEPENDENT_SUBPROBLEMS
  std::vector<double> subprobH ;
  int idxSubPH_ = 2*vDomain ;
#endif // DECOMPOSE_H_INTO_INDEPENDENT_SUBPROBLEMS
  for (val_t i = 0; i < vDomain; ++i) {
    m_assignment[v] = i;

    // get label from based on heuristic class instead
    label = m_heuristic->getLabel(v,m_assignment,n);

	// compute heuristic value
#ifdef DECOMPOSE_H_INTO_INDEPENDENT_SUBPROBLEMS
    heuristic = m_heuristic->getHeurPerIndSubproblem(v, m_assignment, n, label, subprobH);
#else
    heuristic = m_heuristic->getHeur(v, m_assignment, n);
#endif // DECOMPOSE_H_INTO_INDEPENDENT_SUBPROBLEMS


    // store label and heuristic into cache table
    dv[2*i+1] = label; // label
    dv[2*i] = heuristic OP_TIMES label; // heuristic
    if (dv[2*i] > h) {
      argmax_h = i;
      h = dv[2*i]; // keep max. for OR node heuristic
    }
    ordering_cache[i] = m_heuristic->getOrderingHeur(v, m_assignment, n);

#ifdef DECOMPOSE_H_INTO_INDEPENDENT_SUBPROBLEMS
    // store h of each subproblem
    int idxSubPH_End = idxSubPH_ + nChildren ; j = 0 ;
    for (; idxSubPH_ < idxSubPH_End ; idxSubPH_++, j++)
      dv[idxSubPH_] = subprobH[j] ;
#endif // DECOMPOSE_H_INTO_INDEPENDENT_SUBPROBLEMS
  }
#endif

  n->setHeur(h);
  m_assignment[v] = old_value;
  n->setOrderingHeur(ordering_cache[argmax_h]);


  return h;

} // Search::assignCostsOR


#ifndef NO_CACHING
void Search::addCacheContext(SearchNode* node, const vector<int>& ctxt) const {

  context_t sig;
  sig.reserve(ctxt.size());
  for (vector<int>::const_iterator itC=ctxt.begin(); itC!=ctxt.end(); ++itC) {
    sig.push_back(m_assignment[*itC]);
  }

  node->setCacheContext(sig);
#ifdef PARALLEL_DYNAMIC
  node->setCacheInst( m_space->cache->getInstCounter(node->getVar()) );
#endif

}
#endif /* NO_CACHING */



#if defined PARALLEL_DYNAMIC || defined PARALLEL_STATIC
void Search::addSubprobContext(SearchNode* node, const vector<int>& ctxt) const {

  context_t sig;
  sig.reserve(ctxt.size());
  for (vector<int>::const_iterator itC=ctxt.begin(); itC!=ctxt.end(); ++itC) {
    sig.push_back(m_assignment[*itC]);
  }
  node->setSubprobContext(sig);
}
#endif /* PARALLEL_DYNAMIC */


double Search::lowerBound(const SearchNode* node) const {

  // assume OR node
  assert(node->getType() == NODE_OR);

  double maxBound = ELEM_ZERO;
  maxBound = max(maxBound, node->getValue()); // to deal with NaN
  double curBound = ELEM_ZERO;
  double pstVal = ELEM_ONE;

  vector<double> pst;
  node->getPST(pst); // will store PST into pst

  // climb up to search space root node
  for (vector<double>::iterator it=pst.begin(); it!=pst.end(); ) {

    // iterator is at AND position
    pstVal OP_TIMESEQ (*it); // add current label

    ++it; // move to OR value
    curBound = (*it) OP_DIVIDE pstVal; // current bound
    maxBound = max(maxBound,curBound); // order is important (for NaN)

    ++it; // move to next AND label
  }

  return maxBound;

}


bool Search::loadInitialBound(string file) {

  // See if file can be opened
  ifstream inTest(file.c_str());
  inTest.close();
  if (inTest.fail()) {
    cerr << "ERROR reading from file " << file << endl;
    return false;
  }

  // Read the actual value from file
  igzstream infile(file.c_str(), ios::in | ios::binary);
  if (infile) {
    double bound;
    BINREAD(infile, bound);

#ifndef NO_ASSIGNMENT

    // no. of and and or nodes (not relevant here) // TODO breaks old versions
    count_t noOr, noAnd;
    BINREAD(infile, noOr);
    BINREAD(infile, noAnd);

    // read optimal assignment from file
    int noVars;
    BINREAD(infile, noVars);

    if (noVars != m_problem->getNOrg()) {
      cerr << "ERROR reading SLS solution, number of vars doesn't match" << endl;
      return false;
    }

    // read full tuple
    vector<int> tuple; tuple.reserve(noVars);
    for (int i=0; i<noVars; ++i) {
      int x;
      BINREAD(infile,x);
      tuple.push_back(x);
    }

    // filter evidence and unary domain variables that were removed during preprocessing
    vector<val_t> reduced;
    for (int i=0; i<noVars; ++i) {
      if (!m_problem->isEliminated(i)) {
        reduced.push_back( (val_t) tuple.at(i));
      }
    }

    // dummy node = 0
    reduced.push_back((val_t) 0);

    if ((int) reduced.size() != m_problem->getN()) {
      cerr << "ERROR reading SLS solution, reduced problem size doesn't match" << endl;
      return false;
    }
#endif

    // store solution/bound into search space and problem instance
#ifdef NO_ASSIGNMENT
    this->updateSolution(bound);
    m_problem->updateSolution(getCurOptValue(), NULL, true);
#else
    this->updateSolution(bound, reduced);
    m_problem->updateSolution(getCurOptValue(), getCurOptTuple(), NULL, true);
#endif


  }
  infile.close();
  return true;

}


/* restricts the search space to a subproblem as specified:
 * @rootVar: the subproblem root variable
 * @assig: variable assignment, possibly partial, needed for context instantiation
 * @pst: the parent partial solution tree, interpreted *top-down* !
 */
int Search::restrictSubproblem(int rootVar, const vector<val_t>& assig, const vector<double>& pst) {

  // adjust the pseudo tree, returns original depth of new root node
  int depth = m_pseudotree->restrictSubproblem(rootVar);

  // resize node count vectors for subproblem
  m_nodeProfile.clear();
  m_nodeProfile.resize(m_pseudotree->getHeightCond()+1);
  m_nodeProfile.clear();
  m_leafProfile.resize(m_pseudotree->getHeightCond()+1);
  // reset search stats
  m_space->stats = SearchStats();

  // set context assignment
  const vector<int>& context = m_pseudotree->getNode(rootVar)->getFullContextVec();
  for (vector<int>::const_iterator it = context.begin(); it!=context.end(); ++it) {
    m_assignment[*it] = assig[*it];
  }

  // generate structure of bogus nodes to hold partial solution tree copied
  // from master search space
  SearchNode *next = NULL, *node = NULL;
  delete m_space->root;  // delete previous search space root
  int pstSize = pst.size() / 2;
  int dummyVar = m_problem->getN() - 1;

  // build structure top-down, first entry in pst array is assumed to correspond to
  // highest OR value, last entry is the lowest AND label
  for (int i=0; i<pstSize; ++i) {
    // dummy OR node with solution bound from master PST
    next = new SearchNodeOR(node, dummyVar, -1) ;
    next->setValue(pst.at(2*i));
    DIAG(cout << "- Created  OR dummy with value " << pst.at(2*i) << endl;)
    if (i > 0) node->setChild(next);
    else m_space->root = next;
    node = next;

    // dummy AND node with label from master PST
    next = new SearchNodeAND(node, 0, pst.at(2*i+1)) ;
    DIAG(cout << " -Created AND dummy with label " << pst.at(2*i+1) << endl;)
    node->setChild(next);
    node = next;
  }

  // create the OR node for the actual (non-dummy) subproblem root variable
  next = new SearchNodeOR(node, rootVar, 0);
  node->setChild(next);
  m_space->subproblemLocal = next;
  // empty existing queue/stack/etc. and add new node
  this->reset(next);

  return depth;
}


bool Search::updateSolution(double d
#ifndef NO_ASSIGNMENT
    ,const vector<val_t>& tuple
#endif
  ) const {
  assert(m_space && m_space->root);
  if (ISNAN(d))
    return false;
  double curValue = m_space->root->getValue();
  if (!ISNAN(curValue) && d <= curValue)
    return false;
  m_space->root->setValue(d);
#ifndef NO_ASSIGNMENT
  m_space->root->setOptAssig(tuple);
#endif
  return true;
}


bool Search::restrictSubproblem(string file) {
  assert(!file.empty());

  size_t id = NONE;
  size_t i = file.rfind(":");
  if (i != string::npos) {
    string t = file.substr(i+1);
    id = atoi(t.c_str());
    file = file.substr(0,i);
  }

  { // Check if file access works
    ifstream inTemp(file.c_str());
    inTemp.close();

    if (inTemp.fail()) {
      cerr << "Error opening subproblem specification " << file << '.' << endl;
      return false;
    }
  }

  igzstream fs(file.c_str(), ios::in | ios::binary);

  int rootVar = UNKNOWN;
  int x = UNKNOWN;
  val_t y = UNKNOWN;

  if (id!= (size_t) NONE) {
    count_t count = NONE;
    count_t z = NONE;
    BINREAD(fs, count);
    BINREAD(fs, z);
    while (id != z) {
      BINREAD(fs, rootVar);
      BINREAD(fs, x); // context size
      BINSKIP(fs,int,x); // skip context assignment
      BINREAD(fs, x); // PST size
      x = (x<0)? -2*x : 2*x;
      BINSKIP(fs,double,x);

      if (fs.eof()) {
        cerr << "ERROR reading subproblem, id not found" << endl;
        return false;
      }

      BINREAD(fs,z); // next id
    }
  }

  //////////////////////////////////////
  // FIRST: read context of subproblem

  // root variable of subproblem
  BINREAD(fs, rootVar);
  if (rootVar<0 || rootVar>= m_problem->getN()) {
    cerr << "ERROR reading subproblem specification, variable index " << rootVar << " out of range" << endl;
    return false;
  }

  cout << "Restricting to subproblem with root node " << rootVar << endl;

  // context size
  BINREAD(fs, x);
  const vector<int>& context = m_pseudotree->getNode(rootVar)->getFullContextVec();
  if (x != (int) context.size()) {
    cerr << "ERROR reading subproblem specification, context size doesn't match." << endl;
    return false;
  }

  // Read context assignment into temporary vector
  vector<val_t> assignment(m_problem->getN(), UNKNOWN);
  int z = UNKNOWN;
  cout << "Subproblem context:";
  for (vector<int>::const_iterator it = context.begin(); it!=context.end(); ++it) {
    BINREAD(fs, z); // files always contain ints, convert to val_t
    y = (val_t) z;
    if (y<0 || y>=m_problem->getDomainSize(*it)) {
      cerr <<endl<< "ERROR reading subproblem specification, variable value " << (int)y << " not in domain." << endl;
      return false;
    }
    cout << ' ' << (*it) << "->" << (int) y;
    assignment[*it] = y;
  }
  cout << endl;

  ///////////////////////////////////////////////////////////////////////////
  // SECOND: read information about partial solution tree in parent problem

  // read encoded PST values for advanced pruning, and save into temporary vector.
  // the assumed order later is top-down; however, if pstSize from file is negative,
  // it's bottom-up in the file, so in that case we need to reverse
  double valueOR, labelAND;
  bool reverse = false;

  int pstSize = UNKNOWN;
  BINREAD(fs,pstSize);
  if (pstSize < 0) {
    reverse = true; // negative value indicates order needs to be reversed
    pstSize *= -1;
  }
  cout << "Reading parent partial solution tree of size " << pstSize << endl;

  vector<double> pstVals(pstSize*2, ELEM_NAN);
  if (reverse) { // fill array in reverse
    int i=pstSize-1;
    while (i >= 0) {
      BINREAD(fs, labelAND);
      BINREAD(fs, valueOR);
      pstVals.at(2*i+1) = labelAND;
      pstVals.at(2*i) = valueOR;
      i -= 1;
    }
  } else { // no reverse necessary
    int i=0;
    while (i<pstSize) {
      BINREAD(fs, valueOR);
      BINREAD(fs, labelAND);
      pstVals.at(2*i) = valueOR;
      pstVals.at(2*i+1) = labelAND;
      i += 1;
    }
  }
  fs.close();

  // call function to actually condition this search instance
  int depth = this->restrictSubproblem(rootVar, assignment, pstVals);
  cout << "Restricted to subproblem with root node " << rootVar << " at depth " << depth  << endl;

  return true; // success
}

bool Search::propHeuristic(SearchNode* node) {
  double or_delta = ELEM_ONE;
  // Propagate heuristics here.
  SearchNode* cur = node;
  SearchNode* prev;
  bool prop = true;
  uint32 prop_iteration = 0;
  do {
    //      cout << "Prop iteration: " << ++prop_iteration << endl;
    if (cur->getType() == NODE_AND) {
      // if this is where the propagation starts, we have freshly
      // generated OR nodes with new heuristics.
      if (prop_iteration == 0) {
        double d = cur->getSubSolved();
        assert(d == ELEM_ONE); // OR nodes are fresh!
        DIAG( ostringstream ss;
            ss << "value stored (subsolved only): " << d << endl;
            myprint(ss.str()); )
          NodeP* children = cur->getChildren();
        uint32 count_missing_children = 0;
        for (uint32 i = 0; i < cur->getChildCountFull(); ++i) {
          if (children[i]) {
            d OP_TIMESEQ children[i]->getHeur();
            //            cout << setprecision(20) << children[i]->getHeur() << endl;
          } else {
            ++count_missing_children;
          }
        }

        d OP_TIMESEQ cur->getLabel();

        double old_h = cur->getHeur();
        DIAG( ostringstream ss;
            ss << setprecision(20);
            ss << *cur <<" [AND] old_h: " << old_h << ", new_h: " << d << endl;
            myprint(ss.str()); )
          if (d - old_h > 1e-16) {
            if (d - old_h > 1e-12) {
              cout << "subsolved: " << cur->getSubSolved() << endl;
              cout << "missing children: " << count_missing_children << endl;
              cout << "warning: bad AND h would have been propagated. gap:  "
                << fabs(d - old_h) << endl;
              { ostringstream ss;
                ss << setprecision(20);
                ss << *cur <<" [AND] old_h: " << old_h
                  << ", new_h: " << d
                  << ", this label: " << cur->getLabel()
                  << endl;
                myprint(ss.str()); }
              //          cout << "gap: " << d - old_h << endl;
              //          cin.get();
            }
          } else {
            cur->setHeur(d);
          }
      } else { // otherwise we only need to propagate the differences
        double old_h = cur->getHeur();
        if (or_delta < ELEM_ONE) {
          double new_h = cur->getHeur() OP_TIMES or_delta;
          DIAG( ostringstream ss;
              ss << setprecision(20);
              ss << *cur << " [AND] old_h: " << old_h
              << ", new_h: " << new_h << endl; )

            cur->setHeur(new_h);
        } else {
          prop = false;
        }
      }
    } else { // OR node
      double old_h = cur->getHeur();
      double max_h = -std::numeric_limits<double>::infinity();
      NodeP* children = cur->getChildren();
      for (uint32 i = 0; i < cur->getChildCountFull(); ++i) {
        if (children[i]) {
          double new_h =
            //                children[i]->getLabel() OP_TIMES children[i]->getHeur();
            children[i]->getHeur();
          if (new_h > max_h) {
            max_h = new_h;
          }
        }
      }

      double delta = max_h OP_DIVIDE old_h;
      DIAG( ostringstream ss;
          ss << setprecision(20);
          ss << *cur << " [OR] old_h: " << old_h << ", new_h: " << max_h << endl;
          myprint(ss.str()); )
        assert(delta <= ELEM_ONE || fabs(delta) < 1e-12);
      if (max_h < old_h) {
        cur->setHeur(max_h);
        or_delta = delta;
      } else {
        prop = false;
      }

    }
    if (prop) {
      prev = cur;
      cur = cur->getParent();
      prop_iteration++;
    } else {
      break;
    }
  } while (cur);

  if (prop || !cur) {
    m_problem->updateUpperBound(prev->getHeur(), &m_space->stats, true);
  }

  // Check if LB >= UB. Early termination.
  if (m_problem->getSolutionCost() >= m_problem->getUpperBound()) {
    return true;
  }
  return false;
}

void Search::DoConstraintPropagation() {
  if (FLAGS_cp_type == "FULL") {
    vec<Lit> lits;

    // create SAT problem
    var2sat_.resize(m_problem->getN());
    for (int var = 0; var < m_problem->getN(); ++var) {
      var2sat_[var].resize(m_problem->getDomainSize(var));
      for (int val = 0; val < m_problem->getDomainSize(var); ++val) {
        Var vid = minisat_solver_.newVar();
        var2sat_[var][val] = vid;
      }
    }

    vector<val_t> assignment;
    assignment.resize(m_problem->getN(), UNKNOWN);

    // add nogoods
    for (const Function* f : m_problem->getFunctions()) {
      const vector<int>& scope_vars = f->getScopeVec();
      int scope_size = scope_vars.size();
      if (scope_size == 0) continue;
      vector<int> values(scope_size, 0);
      values.back() = -1;

      int i;
      for (;;) {
        for (i = scope_size - 1; i >= 0; --i) {
          int a = scope_vars[i];
          int last = m_problem->getDomainSize(a) - 1;
          if (values[i] < last) break;
          values[i] = 0;
        }

        if (i < 0) break;
        ++values[i];

        for (int a = 0; a < scope_size; ++a) {
          assignment[scope_vars[a]] = values[a];
        }

        double val = f->getValue(assignment);
        if (val == ELEM_ZERO) {
          lits.clear();
          for (int a = 0; a < scope_size; ++a) {
            int vid = var2sat_[scope_vars[a]][values[a]];
            lits.push(~mkLit(vid));
          }
          minisat_solver_.addClause_(lits);
        }
      }
      assignment.resize(m_problem->getN(), UNKNOWN);
    }

    // at-most-one clauses
    for (int var = 0; var < m_problem->getN(); ++var) {
      int domain_size = m_problem->getDomainSize(var);
      for (int v1 = 0; v1 < domain_size - 1; ++v1) {
        for (int v2 = v1 + 1; v2 < domain_size; ++v2) {
          int v1sat = var2sat_[var][v1];
          int v2sat = var2sat_[var][v2];

          lits.clear();
          lits.push(~mkLit(v1sat));
          lits.push(~mkLit(v2sat));
          minisat_solver_.addClause_(lits);
        }
      }
    }

    // at-least-one clauses
    for (int var = 0; var < m_problem->getN(); ++var) {
      lits.clear();
      for (int val = 0; val < m_problem->getDomainSize(var); ++val) {
        lits.push(mkLit(var2sat_[var][val]));
      }
      minisat_solver_.addClause_(lits);
    }

    minisat_solver_.simplify();
    cout << "Created SAT (minisat) instance with " << minisat_solver_.nVars()
         << " variables and " << minisat_solver_.nClauses() << " clauses"
         << endl;
    bool ok = minisat_solver_.solve();
    if (ok) {
      cout << "Initial SAT instance is satisfiable." << endl;
    } else {
      cout << "Initial SAT instance is not satisfiable." << endl;
    }
  }
  else if (FLAGS_cp_type == "UNIT") {
    int var_id = 1;
    int n_sat_vars = 0;
    for (int i = 0; i < m_problem->getN(); ++i) {
      n_sat_vars += m_problem->getDomainSize(i);
    }

    // Create SAT problem
    var2sat_.resize(m_problem->getN());
    sat2var_.resize(n_sat_vars + 1, std::make_pair(-1, -1));
    zchaff_solver_.set_variable_number(n_sat_vars);
    for (int i = 0; i < m_problem->getN(); ++i) {
      int domain_size = m_problem->getDomainSize(i);
      var2sat_[i].resize(domain_size, 0);
      for (int k = 0; k < domain_size; ++k, ++var_id) {
        var2sat_[i][k] = var_id;
        sat2var_[var_id] = std::make_pair(i, k);
      }
    }

    zchaff_solver_.set_var_maps(var2sat_, sat2var_);

    // Create clauses for 0 probability tuples
    int num_cl = 0;
    for (const Function* f : m_problem->getFunctions()) {
      int arity = f->getArity();
      if (arity == 0) continue;
      int* scope_vars = new int[arity];
      int l = 0;
      for (int v : f->getScopeSet()) {
        scope_vars[l++] = v;
      }

      int* values = new int[arity];
      for (int i = 0 ; i < arity - 1; ++i) {
        values[i] = 0;
      }
      values[arity - 1] = -1;

      // Enumerate instantiations
      vector<val_t> assignment(m_problem->getN(), -1);
      int i;
      for (;;) {
        for (i = arity - 1; i >= 0; --i) {
          int a = scope_vars[i];
          int last = m_problem->getDomainSize(a) - 1;
          if (values[i] < last) break;
          values[i] = 0;
        }
        if (i < 0) break;
        ++values[i];

        for (int a = 0; a < arity; ++a) {
          assignment[scope_vars[a]] = values[a];
        }

        // Nogood clause if assignment is 0
        if (f->getValue(assignment) == ELEM_ZERO) {
          zchaff_solver_.add_nogood_clause(arity, scope_vars, values);
          ++num_cl;
        }
      }

      // Reset assignment
      assignment.resize(m_problem->getN(), -1);
      delete[] values;
      delete[] scope_vars;
    }

    // Add at-most-one clauses for domain values
    // (assignment exclusivity constraint)
    for (int i = 0; i < m_problem->getN(); ++i) {
      int domain_size = m_problem->getDomainSize(i);
      for (int v1 = 0; v1 < domain_size - 1; ++v1) {
        for (int v2 = v1 + 1; v2 < domain_size; ++v2) {
          zchaff_solver_.add_atmostone_clause(i, v1, v2);
          ++num_cl;
        }
      }
    }

    // Add at-least-one clauses to ensure something is assigned
    for (int i = 0; i < m_problem->getN(); ++i) {
      int domain_size = m_problem->getDomainSize(i);
      zchaff_solver_.add_atleastone_clause(i, domain_size);
      ++num_cl;
    }

    zchaff_solver_.init_solve();
    zchaff_solver_.preprocess();

    cout << "Created SAT (zchaff) instance with "
         << zchaff_solver_.num_variables() << " variables and "
         << zchaff_solver_.num_clauses() << " clauses" << endl;

    current_domains_.resize(m_problem->getN());
    for (int i = 0; i < m_problem->getN(); ++i) {
      current_domains_[i].resize(m_problem->getDomainSize(i), true);
    }

    list<pair<int, int>> changes;
    zchaff_solver_.update_domains(current_domains_, changes);

    m_prop->setSatSolver(&zchaff_solver_);
  } else {
    m_prop->setSatSolver(nullptr);
    assert(current_domains_.empty());
  }
}

bool Search::DoSATPropagate(const vector<int>& vars) {
  assert(FLAGS_cp_type == "FULL");

  vec<Lit> lits;
  for (int v : vars) {
    lits.push(mkLit(v));
  }

  return minisat_solver_.solve(lits);
}

bool Search::DoCPLookahead(int var, int val, list<pair<int, int>>& changes,
                           const vector<int>& subtree) {
  assert(FLAGS_cp_type == "UNIT");

  // Unit resolution
  int res = zchaff_solver_.propagate(var, val);
  zchaff_solver_.update_domains(current_domains_, changes);
  if (res == zchaff::CONFLICT) {
    return false;
  }

  // Future empty domain in a subproblem?
  for (int desc : subtree) {
    if (var == desc) {
      continue;
    }

    vector<bool>& domain = current_domains_[desc];
    size_t sz = domain.size();
    size_t pruned = 0;
    for (size_t j = 0; j < sz; ++j) {
      if (!domain[j]) {
        ++pruned;
      }
    }
    if (pruned == sz) {
      return false;
    }
  }
  return true;
}

}  // namespace daoopt
