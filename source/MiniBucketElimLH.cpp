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

#include <cstddef>
#include <fstream>

#include "MiniBucketElimLH.h"

/* disables DEBUG output */
#undef DEBUG

/*
extern int64 nd1SpecialCalls;
extern int64 nd1GeneralCalls;
*/

//#define DEBUG_BUCKET_ERROR

#ifndef OUR_OWN_nInfinity
#define OUR_OWN_nInfinity (-std::numeric_limits<double>::infinity())
#endif // OUR_OWN_nInfinity
#ifndef OUR_OWN_pInfinity
#define OUR_OWN_pInfinity (std::numeric_limits<double>::infinity())
#endif // OUR_OWN_pInfinity

namespace daoopt {

static int computeMBOutFnArgsVectorPtrMap(int elim_var,
                                          vector<Function *> &Functions,
                                          set<int> &scope, int &n,
                                          val_t *&tuple,
                                          vector<vector<val_t *>> &idxMap) {
  idxMap.clear();
  idxMap.resize(Functions.size());

  int i;
  uint64 j;

  vector<Function *>::const_iterator fit;
  set<int>::const_iterator sit;

  scope.clear();
  for (fit = Functions.begin(); fit != Functions.end(); ++fit) {
    scope.insert((*fit)->getScopeVec().begin(), (*fit)->getScopeVec().end());
  }
  scope.erase(elim_var);

  if (NULL == tuple) {
    n = scope.size();          // new function arity
    tuple = new val_t[n + 1];  // this keeps track of the tuple assignment
  }
  val_t *elimVal = tuple + n;  // n+1 is elimVar

  // holds iterators .begin() and .end() for all function scopes
  vector<pair<set<int>::const_iterator, set<int>::const_iterator>> iterators;
  iterators.reserve(Functions.size());
  for (fit = Functions.begin(); fit != Functions.end(); ++fit) {
    // store begin() and end() for each function scope
    iterators.push_back(
        make_pair((*fit)->getScopeSet().begin(), (*fit)->getScopeSet().end()));
  }

  // collect pointers to tuple values
  bool bucketVarPassed = false;
  for (i = 0, sit = scope.begin(); i < n; ++i, ++sit) {
    if (!bucketVarPassed && *sit > elim_var) {  // just went past bucketVar
      for (j = 0; j < Functions.size(); ++j) {
        idxMap[j].push_back(elimVal);
        ++(iterators[j]
               .first);  // skip bucketVar in the original function scope
      }
      bucketVarPassed = true;
    }

    for (j = 0; j < Functions.size(); ++j) {
      //      cout << "  f" << funs[j]->getId() << ' ' << *sit << " == " <<
      // *(iterators[j].first) << endl;
      if (iterators[j].first != iterators[j].second  // scope iterator != end()
          && *sit == *(iterators[j].first)) {        // value found
        idxMap[j].push_back(&tuple[i]);
        ++(iterators[j].first);
      }
    }
  }

  if (!bucketVarPassed && elim_var < INT_MAX) {  // bucketVar has highest index
    for (j = 0; j < Functions.size(); ++j) {
      idxMap[j].push_back(elimVal);
    }
  }

  return 0;
}

void MiniBucketElimLH::reset(void) {
  MiniBucketElim::reset();
  _Stats.reset();
}

void MiniBucketElimLH::Delete() {
  for (MiniBucketElimLHError& be_error : _ErrorHelper) {
    be_error.Delete();
  }
  deleteLocalErrorFNs();

  for (vector<MiniBucket>& mb : _MiniBuckets) {
    mb.clear();
  }
  for (set<int>& scope : _BucketScopes) {
    scope.clear();
  }
  for (vector<Function*>& b_fns : _BucketFunctions) {
    b_fns.clear();
  }
  _MiniBuckets.clear();
  _BucketScopes.clear();
  _BucketFunctions.clear();

  for (auto& aug_fns : m_augmented) {
    for (Function* fn : aug_fns) {
      delete fn;
    }
  }
  m_augmented.clear();
  m_intermediate.clear();
}

size_t MiniBucketElimLH::build(const std::vector<val_t> *assignment,
                               bool computeTables) {

#if defined DEBUG || _DEBUG
  cout << "$ Building MBEX(" << m_ibound << ")" << endl;
#endif

  this->reset();
  _Stats._iBound = m_ibound;

  vector<int> elimOrder;    // will hold dfs order
  findDfsOrder(elimOrder);  // computes dfs ordering of relevant subtree

  // to define m_augmented/m_intermediate, let v be a variable being eliminated,
  // j be the mini-bucket of v being processed,
  // and u be the closest (to the bucket of v) bucket that contains a variable
  // in fn_v_j, where fn_v_j is the fn produced
  // by minibucket j of v.
  // m_augmented is a list of functions generated by the minibucket algorithm
  // that are placed in the bucket of u.
  // m_intermediate is a list of functions generated by the minibucket algorithm
  // that are placed in bucket between (u,v).
  m_augmented.resize(m_problem->getN());
  m_intermediate.resize(m_problem->getN());

  // keep track of total memory consumption
  _Stats._MemorySize = 0;

  _BucketScopes.resize(m_problem->getN());
  _BucketFunctions.resize(m_problem->getN());
  //	_BucketHFnScopes.resize(m_problem->getN()) ;
  _MiniBuckets.resize(m_problem->getN());
  _ErrorHelper.resize(m_problem->getN());
  _BucketErrorQuality.resize(m_problem->getN());
  for (int i = m_problem->getN() - 1; i >= 0; i--) _BucketErrorQuality[i] = -1;
  _distToClosestDescendantWithMBs.resize(m_problem->getN());
  for (int i = m_problem->getN() - 1; i >= 0; i--)
    _distToClosestDescendantWithMBs[i] = INT_MAX;
  _distToClosestDescendantWithLE.resize(m_problem->getN());
  for (int i = m_problem->getN() - 1; i >= 0; i--)
    _distToClosestDescendantWithLE[i] = INT_MAX;
  deleteLocalErrorFNs();

  // ITERATES OVER BUCKETS, FROM LEAVES TO ROOT
  for (vector<int>::reverse_iterator itV = elimOrder.rbegin();
       itV != elimOrder.rend(); ++itV) {
    int v = *itV;  // this is the variable being eliminated

    vector<Function *> &funs =
        _BucketFunctions[v];  // this is the list of all (original +
                              // MB-generated) in the bucket of v
    funs.clear();
    set<int> &jointScope = _BucketScopes[v];
    jointScope.clear();

#if defined DEBUG || _DEBUG
    cout << "$ Bucket for variable " << *itV << endl;
#endif
    if (NULL != m_options ? NULL != m_options->_fpLogFile : false) {
      fprintf(m_options->_fpLogFile, "\nBucket for variable %d(domainsize=%d)",
              (int)*itV, (int)m_problem->getDomainSize(v));
    }

    // collect relevant functions in funs
    /*		set<int> & jointHScope = _BucketHFnScopes[v] ;
          jointHScope.clear() ;
    // add h fn functions
    funs.insert(funs.end(), m_augmented[v].begin(), m_augmented[v].end()) ;
    funs.insert(funs.end(), m_intermediate[v].begin(), m_intermediate[v].end())
    ;
    // compute width
    int nHFNs = funs.size() ;
    for (vector<Function*>::iterator itF = funs.begin() ; itF != funs.end() ;
    ++itF)
    jointHScope.insert((*itF)->getScopeVec().begin(),
    (*itF)->getScopeVec().end()) ;
    if (_Stats._HPseudoWidth < jointHScope.size())
    _Stats._HPseudoWidth = jointHScope.size() ;*/
    // add original functions
    const vector<Function *>& fnlist = m_pseudotree->getFunctions(v);
    funs.insert(funs.end(), fnlist.begin(), fnlist.end());
    // compute width
    //		int nOriginalFNs = funs.size() ;
    for (vector<Function*>::iterator itF = funs.begin(); itF != funs.end();
         ++itF)
      jointScope.insert((*itF)->getScopeVec().begin(),
                        (*itF)->getScopeVec().end());
    double stats_Width = jointScope.size();
    if (_Stats._Width < int(jointScope.size()))
      _Stats._Width = jointScope.size();
    // add MB-generated functions
    funs.insert(funs.end(), m_augmented[v].begin(), m_augmented[v].end());
    // compute pseudo-width
    //		int nTotalFNs = funs.size() ;
    for (vector<Function *>::iterator itF = m_augmented[v].begin();
         itF != m_augmented[v].end(); ++itF)
      jointScope.insert((*itF)->getScopeVec().begin(),
                        (*itF)->getScopeVec().end());
    double stats_PseudoWidth = jointScope.size();
    if (_Stats._PseudoWidth < int(jointScope.size()))
      _Stats._PseudoWidth = jointScope.size();
#if defined DEBUG || _DEBUG
    for (vector<Function *>::iterator itF = funs.begin(); itF != funs.end();
         ++itF)
      cout << ' ' << (**itF);
    cout << endl;
#endif
    size_t tablesize = 1;
    if (NULL != m_options ? NULL != m_options->_fpLogFile : false) {
      for (const int scope_var : jointScope) {
        if (v == scope_var) {
          continue;
        }
        tablesize *= m_problem->getDomainSize(scope_var);
      }
      fprintf(m_options->_fpLogFile,
              " Width=%d PseudoWidth=%d OUTtablesize=%lld", (int)stats_Width,
              (int)stats_PseudoWidth, (int64)tablesize);
    }

    // compute global upper bound for root (dummy) bucket
    if (v == elimOrder[0]) {  // variable is dummy root variable
      if (computeTables &&
          assignment) {  // compute upper bound if assignment is given
        m_globalUB = ELEM_ONE;
        for (vector<Function *>::iterator itF = funs.begin(); itF != funs.end();
             ++itF)
          m_globalUB OP_TIMESEQ(*itF)->getValue(*assignment);
        cout << "    MBE-ALL  = " << SCALE_LOG(m_globalUB) << " ("
             << SCALE_NORM(m_globalUB) << ")" << endl;
        m_globalUB OP_DIVIDEEQ m_problem->globalConstInfo();  // for backwards
                                                              // compatibility
                                                              // of output
        cout << "    MBE-ROOT = " << SCALE_LOG(m_globalUB) << " ("
             << SCALE_NORM(m_globalUB) << ")" << endl;
      }
      continue;  // skip the dummy variable's bucket
    }

    // sort functions by decreasing scope size
    sort(funs.begin(), funs.end(), scopeIsLarger);

    // partition functions into minibuckets
    vector<MiniBucket>& minibuckets = _MiniBuckets[v];
    minibuckets.clear();
    //		vector<Function*>::iterator itF; bool placed ;
    for (vector<Function*>::iterator itF = funs.begin(); itF != funs.end();
         ++itF) {
      //			while (funs.size()) {
      bool placed = false;
      for (vector<MiniBucket>::iterator itB = minibuckets.begin();
           !placed && itB != minibuckets.end(); ++itB) {
        if (itB->allowsFunction(*itF)) {  // checks if function fits into bucket
          itB->addFunction(*itF);
          placed = true;
        }
      }
      if (!placed) {  // no fit, need to create new bucket
        MiniBucket mb(v, m_ibound, m_problem);
        mb.addFunction(*itF);
        minibuckets.push_back(mb);
      }
      //			funs.pop_front() ;
      //			funs.erase(itF) ;
    }

    // Moment-matching step, performed only if we have partitioning.
    vector<Function *> max_marginals;
    std::unique_ptr<Function> average_mm_function;
    if (computeTables && m_options->match && minibuckets.size() > 1) {

      set<int> scope_intersection;
      bool first_mini_bucket = true;
      for (const MiniBucket &mini_bucket : minibuckets) {
        if (first_mini_bucket) {
          scope_intersection = mini_bucket.getJointScope();
          first_mini_bucket = false;
        } else {
          scope_intersection =
              intersection(scope_intersection, mini_bucket.getJointScope());
        }
      }
      for (MiniBucket &mini_bucket : minibuckets) {
        set<int> elim_vars =
            setminus(mini_bucket.getJointScope(), scope_intersection);
        max_marginals.push_back(
            mini_bucket.eliminate(computeTables, elim_vars));
      }

      // Find average max-marginals (geometric mean)
      size_t table_size = 1;
      for (const int var : scope_intersection) {
        table_size *= m_problem->getDomainSize(var);
      }

      double *average_mm_table = new double[table_size];
      for (size_t i = 0; i < table_size; ++i) {
        average_mm_table[i] = ELEM_ONE;
      }
      for (const Function *max_marginal : max_marginals) {
        for (size_t i = 0; i < table_size; ++i) {
          average_mm_table[i] OP_TIMESEQ max_marginal->getTable()[i];
        }
      }
      for (size_t i = 0; i < table_size; ++i) {
        average_mm_table[i] = OP_ROOT(average_mm_table[i], minibuckets.size());
      }
      int dummy_id = 0;
      average_mm_function.reset(
          new FunctionBayes(dummy_id, m_problem, scope_intersection,
                            average_mm_table, table_size));
    }

    int nMBs = minibuckets.size();
    if (_Stats._MaxNumMBs < nMBs) _Stats._MaxNumMBs = nMBs;
    if (nMBs > 1) ++(_Stats._NumBucketsWithMoreThan1MB);

    if (NULL != m_options ? NULL != m_options->_fpLogFile : false) {
      fprintf(m_options->_fpLogFile, " nMBs=%d", (int)minibuckets.size());
    }

    // minibuckets for current bucket are now ready, process each and place
    // resulting function
    uint32 bucket_idx = 0;
    for (vector<MiniBucket>::iterator itB = minibuckets.begin();
         itB != minibuckets.end(); ++itB, ++bucket_idx) {
      Function *newf;
      if (!computeTables || !m_options->match || minibuckets.size() <= 1) {
        newf = itB->eliminate(computeTables);  // process the minibucket
      } else {
        newf = itB->eliminateMM(computeTables, max_marginals[bucket_idx],
                                average_mm_function.get());
      }
      const set<int> &newscope = newf->getScopeSet();
      _Stats._MemorySize += newf->getTableSize();
      // go up in tree to find target bucket
      PseudotreeNode *n = m_pseudotree->getNode(v)->getParent();
      while (newscope.find(n->getVar()) == newscope.end() &&
             n != m_pseudotree->getRoot()) {
        m_intermediate[n->getVar()].push_back(newf);
        n = n->getParent();
      }
      // matching bucket found OR root of pseudo tree reached
      m_augmented[n->getVar()].push_back(newf);
    }
    // all minibuckets processed and resulting functions placed
  }

#if defined DEBUG || _DEBUG
  // output augmented and intermediate buckets
  if (computeTables) {
    /*		for (int i = 0; i < m_problem->getN() ; i++) {
          cout << "$ AUG" << i << ": " << m_augmented[i] << " + " <<
       m_intermediate[i] << endl ;
          }*/
  }
#endif

  // for each bucket, compute closest distance to descendant with more than 1
  // MB.
  for (vector<int>::reverse_iterator itV = elimOrder.rbegin();
       itV != elimOrder.rend(); ++itV) {
    int v = *itV;  // this is the variable being eliminated
    vector<MiniBucket> &minibuckets = _MiniBuckets[v];
    PseudotreeNode *parent_n = m_pseudotree->getNode(v)->getParent();
    if (NULL == parent_n) continue;
    int parent_v = parent_n->getVar(), d2parent_v = INT_MAX;
    if (minibuckets.size() > 1)
      d2parent_v = 0;
    else
      d2parent_v = _distToClosestDescendantWithMBs[v];
    if (d2parent_v >= INT_MAX) continue;
    if (++d2parent_v < _distToClosestDescendantWithMBs[parent_v])
      _distToClosestDescendantWithMBs[parent_v] = d2parent_v;
  }

  // clean up for estimation mode
  if (!computeTables) {
    Delete();
  } else {
    double total_memory_limit =
        m_options->lookahead_local_error_all_tables_total_limit;
    double table_memory_limit =
        m_options->lookahead_local_error_single_table_limit;
    computeLocalErrorTables(true, total_memory_limit, table_memory_limit);

    // for each bucket, compute closest distance to descendant with >0 bucket
    // error.
    for (vector<int>::reverse_iterator itV = elimOrder.rbegin();
         itV != elimOrder.rend(); ++itV) {
      int v = *itV;  // this is the variable being eliminated
      PseudotreeNode *parent_n = m_pseudotree->getNode(v)->getParent();
      if (NULL == parent_n) continue;
      int parent_v = parent_n->getVar(), d2parent_v = INT_MAX;
      if (0 != _BucketErrorQuality[v])
        d2parent_v = 0;
      else
        d2parent_v = _distToClosestDescendantWithLE[v];
      if (d2parent_v >= INT_MAX) continue;
      if (++d2parent_v < _distToClosestDescendantWithLE[parent_v])
        _distToClosestDescendantWithLE[parent_v] = d2parent_v;
    }
    
    // for each bucket, compute closest distance to descendant with >0 bucket
    // error
    for (auto itV = elimOrder.rbegin(); itV != elimOrder.rend(); ++itV) {
      int v = *itV;
      _ErrorHelper[v].Initialize(*this, v, m_options->lookahead_depth);
    }
  }

  cout << StrCat("Pseudowidth: ", _Stats._PseudoWidth - 1) << endl;;
  double minibucket_mem_mb = _Stats._MemorySize *
      sizeof(double) / (1024.0 * 1024);
  cout << StrCat("Minibucket Memory (MB): ", minibucket_mem_mb) << endl;
  cout << StrCat("Local Error Memory (MB): ", _Stats._LEMemorySizeMB) << endl;
  cout << StrCat("Total Heuristic Memory (MB): ",
                 minibucket_mem_mb + _Stats._LEMemorySizeMB) << endl;
  return _Stats._MemorySize;
}

double MiniBucketElimLH::getHeur(int var, std::vector<val_t> & assignment, SearchNode *search_node) 
{
#if defined DEBUG || _DEBUG
	int nSet = 0;
	for (int ii = 0; ii < assignment.size(); ii++) {
		if (assignment[ii] >= 0) nSet++;
  }
#endif

	// normally (wo lookahead) we would go over augmented[var] and
	// intermediate[var] lists and combine all values.
	// i.e. we would process bucket[var].
	// however, for us here, we would have to look at one level below that in the
	// bucket tree.

	double h = MiniBucketElim::getHeur(var, assignment, search_node);
	if (m_options->lookahead_depth <= 0) return h;

  /*
	const int local_lookaheadDepth = m_options->lookahead_depth - 1;
	double DH = 0.0, dh;
	if (0 == local_lookaheadDepth) {  // if looking ahead just 1 step, logic is a little easier
		if (1 != _distToClosestDescendantWithLE[var]) 
			return h;  // there is no child with non-0 BucketError
		const PseudotreeNode *n = m_pseudotree->getNode(var);
		const vector<PseudotreeNode *> &children = n->getChildren();
		//		const int nChildren = children.size() ;
		for (vector<PseudotreeNode *>::const_iterator it = children.begin(); it != children.end(); ++it) {
			int child = (*it)->getVar();
			if (0 == _BucketErrorQuality[child])
				continue;                             // most likely, child has nMBs=1.
			dh = getLocalError(child, assignment);  // note : dh should be >= 0.0
			DH += dh;
			}
// DEBUGGG
//    if (DH < 0.0) printf("\nERROR 111");
		} 
	else {
// DEBUGGGG
exit(999);
		if (_distToClosestDescendantWithLE[var] > m_options->lookahead_depth)
			return h;  // there is no child with non-0 BucketError within the range of
		// lookahead depth
		const PseudotreeNode *n = m_pseudotree->getNode(var);
		const vector<PseudotreeNode *> &children = n->getChildren();
		//		const int nChildren = children.size() ;
		for (vector<PseudotreeNode *>::const_iterator it = children.begin();it != children.end(); ++it) {
			int child = (*it)->getVar();
			dh = getHeuristicError(child, assignment,local_lookaheadDepth);  // note : dh should be >= 0.0
			DH += dh;
			}
		}
	return h - DH;
  */
  if (_distToClosestDescendantWithLE[var] > m_options->lookahead_depth) {
    return h;
  }
//  ++nd1GeneralCalls;
  double DH = _ErrorHelper[var].Error(assignment);
  return h - DH;
}

void MiniBucketElimLH::getHeurAll(int var, vector<val_t> &assignment,
                                  SearchNode *search_node,
                                  vector<double> &out) {
  MiniBucketElim::getHeurAll(var, assignment, search_node, out);
  if (m_options->lookahead_depth <= 0) return;

  /*
  const int local_lookaheadDepth = m_options->lookahead_depth - 1;
  if (0 == local_lookaheadDepth) {  // if looking ahead just 1 step, logic is a
                                    // little easier
    if (1 != _distToClosestDescendantWithLE[var])
      return;  // there is no child with non-0 BucketError
    double var_original_value = assignment[var];
    PseudotreeNode *n = m_pseudotree->getNode(var);
    const vector<PseudotreeNode *> &children = n->getChildren();
    //		int nChildren = children.size() ;
    for (val_t i = m_problem->getDomainSize(var) - 1; i >= 0; i--) {
      assignment[var] = i;
      double DH = 0.0, dh;  // ELEM_ONE ;
      for (vector<PseudotreeNode *>::const_iterator it = children.begin();
           it != children.end(); ++it) {
        int child = (*it)->getVar();
        if (0 == _BucketErrorQuality[child])
          continue;  // most likely, child has nMBs=1.
        dh = getLocalError(child, assignment);  // note : dh should be >= 0.0
		DH += dh;
      }
      out[i] -= DH;
    }
    assignment[var] = var_original_value;
  } else {
  */
  if (_distToClosestDescendantWithLE[var] > m_options->lookahead_depth) {
    return;  // there is no child with non-0 BucketError within the range of
            // lookahead depth
  }
//  ++nd1GeneralCalls;
  _ErrorHelper[var].Error(assignment, out);
}

double MiniBucketElimLH::getLocalError(int var, vector<val_t> & assignment)
{
	// check if number of MBs is 1; if yes, bucket error is 0.
	vector<MiniBucket>& minibuckets = _MiniBuckets[var];
	if (minibuckets.size() <= 1 || 0 == _BucketErrorQuality[var]) return 0.0;

#ifndef DEBUG_BUCKET_ERROR
	if (NULL != _BucketErrorFunctions[var]) {
		double dh = _BucketErrorFunctions[var]->getValue(assignment);
		return dh;
  }
#endif // DEBUG_BUCKET_ERROR

  int var_domain_size = m_problem->getDomainSize(var);

  // combine MB output FNs
  double tableentryMB = ELEM_ONE;
  for (const MiniBucket& mb : minibuckets) {
    Function *fMB = mb.output_fn();
    if(!fMB) {
      continue;
    }
    tableentryMB OP_TIMESEQ fMB->getValue(assignment);
  }
  if (OUR_OWN_nInfinity == tableentryMB) {
    return 0.0; // it must be that tableentryB is also -inf and error is 0.
  }

  // compute bucket output for the given assignment; this is the exact bucket
  // value, given its functions and assignment.
  vector<double> funVals;
  vector<double> sumVals;
  sumVals.resize(m_problem->getDomainSize(var), ELEM_ONE);
	vector<Function*>& funs_B = _BucketFunctions[var];
  for (Function* fn : funs_B) {
    fn->getValues(assignment, var, funVals);
    for (int i = 0; i < var_domain_size; ++i) {
      sumVals[i] OP_TIMESEQ funVals[i];
    }
  }
  double tableentryB = sumVals.front();
  for (int i = 1; i < var_domain_size; ++i) {
    tableentryB = max(tableentryB, sumVals[i]);
  }
  /*
	// enumerate over all bucket var values
	double var_original_value = assignment[var];
	double tableentryB = ELEM_ZERO;
  for (int i = 0; i < m_problem->getDomainSize(var); ++i) {
		assignment[var] = i;
		// combine all bucket FNs
		double zB = ELEM_ONE;
    for (Function* fn : funs_B) {
			zB OP_TIMESEQ fn->getValue(assignment);
    }
		tableentryB = max(tableentryB, zB);
  }
	assignment[var] = var_original_value;
  */

#ifdef DEBUG_BUCKET_ERROR
  if (NULL != _BucketErrorFunctions[var]) {
    double dh = _BucketErrorFunctions[var]->getValue(assignment);
    double dh_ = (tableentryMB <= tableentryB) ? 0.0 : (tableentryMB - tableentryB) ;
    if (fabs(dh - dh_) > 1.0e-6) {
//      cout << StrCat("var: ", var, " : ", fabs(dh - dh_)) << endl;
      exit(998) ; // error : online/offline versions are different
    }
  }
#endif // DEBUG_BUCKET_ERROR

	// error must be >= 0 since MB approximation is an upper bound on tableentryB,
	// i.e. Bucket fn value is <= combination of MB output fn values.
	// note, in log-space, tableentryB could be -infinity, but in that case
	// tableentryMB should also be -infinity.
	if (tableentryMB <= tableentryB) return 0.0 ; // '<' is an error, '=' is ok.

/*
#ifdef USE_LOG
//	double e = log10(pow(10.0, tableentryMB) - pow(10.0, tableentryB)) ;
double e = tableentryMB + log10(1 - pow(10.0, tableentryB - tableentryMB)) ;
#else
double e = (tableentryMB - tableentryB) ;
#endif
*/
	return tableentryMB - tableentryB;
}

void MiniBucketElimLH::getLocalError(int parent, int var,
                                     vector<val_t>& assignment,
                                     std::vector<double>& out) {
  // check if number of MBs is 1; if yes, bucket error is 0.
  vector<MiniBucket>& minibuckets = _MiniBuckets[var];
  if (minibuckets.size() <= 1 || 0 == _BucketErrorQuality[var]) {
    return;
  }

  vector<double> fun_vals;
  int j, k, parent_domain_size = m_problem->getDomainSize(parent);

  if (NULL != _BucketErrorFunctions[var]) {
    _BucketErrorFunctions[var]->getValues(assignment, parent, fun_vals);
    for (k = 0; k < parent_domain_size; ++k) {
      out[k] OP_DIVIDEEQ fun_vals[k];
    }
    return;
  }

  // we could try to compute directly : out - MB + B, but that would have 
  // two problems:
  // 1) we could not catch the rate cases of floating-point errors where MB<B,
  // 2) when MB/B are -inf, the math goes crazy.

  // combine MB output FNs for each value of parent
  vector<double> mb_vals;
  mb_vals.resize((vector<double>::size_type) parent_domain_size, ELEM_ONE);

  for (const MiniBucket& mb : minibuckets) {
    Function* fMB = mb.output_fn();
    if (!fMB) {
      continue;
    }
    fMB->getValues(assignment, parent, fun_vals);
    for (k = 0; k < parent_domain_size; ++k) {
      mb_vals[k] OP_TIMESEQ fun_vals[k];
    }
  }

  // compute bucket output for the given assignment; this is the exact bucket 
  // value, given its functions and assignment.
  // note here we have two unassigned variables : parent and child (var), we 
  // have to enumerate over one of them.
  vector<Function*>& funs_B = _BucketFunctions[var];
  double parent_original_value = assignment[parent];
  int var_domain_size = m_problem->getDomainSize(var);
  vector<double> sum_vals;
  sum_vals.resize((vector<double>::size_type)var_domain_size);
  for (k = 0; k < parent_domain_size; ++k) {
    assignment[parent] = k;
    for (j = 0; j < var_domain_size; ++j) {
      sum_vals[j] = ELEM_ONE;
    }
    for (Function* fn : funs_B) {
      fn->getValues(assignment, var, fun_vals);
      for (j = 0; j < var_domain_size; ++j) {
        sum_vals[j] OP_TIMESEQ fun_vals[j];
      }
    }
    double tableentryB = sum_vals[0];
    for (j = 1; j < var_domain_size; ++j) {
      tableentryB = max(tableentryB, sum_vals[j]);
    }

    // process wrt MB value
    if (mb_vals[k] <= tableentryB) {
      continue;
    }
    out[k] OP_DIVIDEEQ (mb_vals[k] - tableentryB);
  }
  assignment[parent] = parent_original_value;
}


int MiniBucketElimLH::computeLocalErrorTable(
    int var, bool build_table, double TableMemoryLimitAsNumElementsLog,
    double &TableSizeLog, double &avgError, double &avgExact,
    Function *&errorFn) {
  // assumption : each minibucket has its output fn computed
  // algorithm : for each entry in the output table (over 'jointscope'-var)
  // compute the max-product/sum of MB_output_FNs and max-product/sum of all fns
  // in the bucket, and then the different between the two. this should be >= 0.

  // safety check
  if (TableMemoryLimitAsNumElementsLog > 9.0)
    TableMemoryLimitAsNumElementsLog = 9.0;

  avgError = avgExact = DBL_MAX;
  errorFn = NULL;
  TableSizeLog = OUR_OWN_nInfinity;

  int i, k;

  //	int N = m_problem->getN() ;
  vector<MiniBucket> &minibuckets = _MiniBuckets[var];
  set<int> &jointScope = _BucketScopes[var];

  if (minibuckets.size() <= 1) {
    if (NULL != m_options ? NULL != m_options->_fpLogFile : false) {
      fprintf(m_options->_fpLogFile,
              "\n   Computing localError for var=%d, nMBs = 1, avg error = 0",
              (int)var);
    }
    _BucketErrorQuality[var] = 0;
    avgError = 0.0;
    TableSizeLog = OUR_OWN_nInfinity;
    return 0;
  }

  // make scope of bucket output fn
  set<int> scope = jointScope;
  // remove elimVar from the new scope
  scope.erase(var);
  int scope_size = scope.size();  // new function arity

  // compute table size; do it is log space since sometimes the table size can
  // overflow int64.
  TableSizeLog = 0.0;
  set<int>::const_iterator sit;
  for (sit = scope.begin(); sit != scope.end(); ++sit) {
    double dsLog = log10(m_problem->getDomainSize(*sit));
    if (dsLog < 0.0) {
      _BucketErrorQuality[var] = 0;
      TableSizeLog = OUR_OWN_nInfinity;
      avgError = 0.0;
      return 0;
    }
    TableSizeLog += dsLog;
  }
  if (TableSizeLog < 0.0) {  // this should not happen; we check just in case.
    _BucketErrorQuality[var] = 0;
    avgError = 0.0;
    TableSizeLog = OUR_OWN_nInfinity;
    return 0;
  }
  if (TableSizeLog > TableMemoryLimitAsNumElementsLog) {  // table would be too
                                                          // large; will not fit
                                                          // in memory.
    //		if (NULL != m_options ? NULL != m_options->_fpLogFile : false) {
    //			fprintf(m_options->_fpLogFile, "\n   Computing localError for
    //var=%d, nMBs = %d, avg error = NA, tablesize = %g entries", (int) var,
    //(int) minibuckets.size(), (double) TableSizeLog) ;
    //			}
    // leave _BucketErrorQuality[var] as -1 (or whatever it is).
    return 1;
  }

  if (_Stats._EnforcedPseudoWidth < scope_size) {
    _Stats._EnforcedPseudoWidth = scope_size;
  }

  // construct vector of scope domains; compute new table size and collect
  // domain sizes
  vector<val_t> domains;
  domains.reserve(scope_size);
  int64 TableSize = 1;
  for (sit = scope.begin(); sit != scope.end(); ++sit) {
    TableSize *= ((int64)m_problem->getDomainSize(*sit));
    domains.push_back(m_problem->getDomainSize(*sit));
  }
  size_t bucket_var_domain_size = m_problem->getDomainSize(var);
  if (TableSize <= 0) {  // this should not happen; we check just in case.
    _BucketErrorQuality[var] = 0;
    avgError = 0.0;
    TableSizeLog = OUR_OWN_nInfinity;
    return 0;
  }

  // collect all FNs
  vector<Function *> &funs_B = _BucketFunctions[var];
  vector<Function *> funs_MB;
  for (const MiniBucket &mini_bucket : minibuckets) {
    Function *fMB = mini_bucket.output_fn();
    if (NULL == fMB) continue;
    funs_MB.push_back(fMB);
  }

  set<int> scopeB;
  int n = 0;
  val_t *tuple = NULL;
  vector<vector<val_t*>> idxMapB;
  //	int resB = computeMBOutFnArgsVectorPtrMap(var, funs_B, scopeB, n, tuple,
  //idxMapB) ;
  computeMBOutFnArgsVectorPtrMap(var, funs_B, scopeB, n, tuple, idxMapB);
  set<int> scopeMB;
  vector<vector<val_t*>> idxMapMB;
  // int resMB = computeMBOutFnArgsVectorPtrMap(INT_MAX, funs_MB, scopeMB, n,
  // tuple, idxMapMB) ;
  computeMBOutFnArgsVectorPtrMap(INT_MAX, funs_MB, scopeMB, n, tuple, idxMapMB);
  // size of tuple is n+1

  // safety check : scopes must be the same
  if (scopeB.size() != scopeMB.size()) { // scope of all bucket FNs - var should
                                         // be the same as (the union of) scopes
                                         // of MB output FNs.
    return 1;
  }
  if (scopeB.size() != scope.size()) {
    return 1;
  }
  if (int(scope.size()) != n) {
    return 1;
  }

  // we will iterate over all combinations of new fn scope value stored here
  for (i = n - 2; i >= 0; --i) {
    tuple[i] = 0;
  }
  tuple[n - 1] = -1;

  double *newTable = NULL;
  if (build_table) {
    newTable = new double[TableSize];
  }

#if defined DEBUG || _DEBUG
  size_t nBadErrorValues = 0;
#endif

  if (NULL != m_options ? NULL != m_options->_fpLogFile : false) {
    fprintf(
        m_options->_fpLogFile,
        "\n   MiniBucketElimLH::computeLocalErrorTable var=%d tablesize=%lld",
        (int)var, (int64)TableSize);
  }
  printf("\n   MiniBucketElimLH::computeLocalErrorTable var=%d tablesize=%lld",
         (int)var, (int64)TableSize);

  // enumerate all new fn scope combinations
  double e, numErrorItems = 0.0;
  avgError = avgExact = 0.0;
  for (int64 j = 0; j < TableSize; j++) {
    // find next combination
    for (i = n - 1; i >= 0; i--) {
      if (++tuple[i] < domains[i]) {
        break;
      }
      tuple[i] = 0;
    }

    // enumerate over all bucket var values; combine all bucket FNs
    double tableentryB = ELEM_ZERO, zB;
    for (tuple[n] = 0; tuple[n] < int(bucket_var_domain_size); tuple[n]++) {
      zB = ELEM_ONE;
      for (k = 0; k < int(funs_B.size()); ++k)
        zB OP_TIMESEQ funs_B[k]->getValuePtr(idxMapB[k]);
      tableentryB = max(tableentryB, zB);
    }

    // combine MB output FNs
    double tableentryMB = ELEM_ONE;
    for (k = 0; k < int(funs_MB.size()); k++) {
      tableentryMB OP_TIMESEQ funs_MB[k]->getValuePtr(idxMapMB[k]);
    }

    /*		// note that it should be that tableentryMB >= tableentryB, regardless
    of whether log/normal space.
    // note that tableentryB==-infinity could happen (in case of log scale).
    if (OUR_OWN_nInfinity == tableentryMB || OUR_OWN_nInfinity == tableentryB) {
    // if OUR_OWN_mInfinity == tableentryMB then it also must be that
    OUR_OWN_mInfinity == tableentryB; in that case e=infinity.
    // if OUR_OWN_mInfinity == tableentryB then e=infinity.
    e = OUR_OWN_pInfinity ;
    }
    else */
    if (tableentryMB <= tableentryB) {  // '<' is an error, '=' is ok.
      e = 0.0;
#if defined DEBUG || _DEBUG
      if (tableentryMB < tableentryB) nBadErrorValues++;
#endif
    } else  // note tableentryB may be -infinity.
      e = tableentryMB - tableentryB;

// at this point, it should be that e >= 0

#if defined DEBUG_BUCKER_ERROR
    // TESTING : that tables return same value(s); we use getValue() fn of the
    // tables; getValue() takes a complete assignment (that may be partial) as
    // input.
    // in actual exe, this piece of code can be commented out, since its purpose
    // is debugging.
    vector<val_t> assignment;
    assignment.resize(N);
    for (j = 0; j < N; ++j) assignment[j] = -1;

    // set up (partial) complete assignment
    for (sit = scopeB.begin(), k = 0; sit != scopeB.end(); ++sit, ++k)
      assignment[*sit] = tuple[k];
    // combine MB output FNs
    double tableentryMB_ = ELEM_ONE;
    for (const Function* fn : funs_MB) {
      tableentryMB_ OP_TIMESEQ fn->getValue(assignment);
    }
    // enumerate over all bucket var values; combine all bucket FNs
    double tableentryB_ = ELEM_ZERO;
    for (tuple[n] = 0; tuple[n] < bucket_var_domain_size; tuple[n]++) {
      assignment[var] = tuple[n];
      zB = ELEM_ONE;
      for (const Function* fn : funs_B) {
        zB OP_TIMESEQ fn->getValue(assignment);
      }
      tableentryB_ = max(tableentryB_, zB);
    }
    double e_ = tableentryMB_ - tableentryB_;
#endif

    numErrorItems += 1.0;
    avgError += e;
    //		if (tableentryB >= DBL_MIN) // don't add -infinity, which may happen
    //in case of log-space representation
    avgExact += tableentryB;
    if (NULL != newTable) newTable[j] = e;
  }

  if (numErrorItems > 0.0) {
    avgError = avgError / numErrorItems;
    avgExact = avgExact / numErrorItems;
  } else if (build_table)
    build_table = false;

  // if avgError is 0, don't build table. note avgError could be infinity (in
  // case of log space representation).
  if (avgError <= DBL_MIN) {
    _BucketErrorQuality[var] = 0;
    if (build_table) build_table = false;
  } else {
    _BucketErrorQuality[var] = 1;
  }

  if (build_table && NULL != newTable) {
    errorFn = new FunctionBayes(-var, m_problem, scope, newTable, TableSize);
    if (NULL != errorFn)
      newTable = NULL;  // table belongs to errorFn
    else if (NULL != errorFn) {
      delete errorFn;
      errorFn = NULL;
    }
  }

  if (NULL != m_options ? NULL != m_options->_fpLogFile : false) {
    double rel_error =
        avgExact > 0.0 ? fabs(100.0 * avgError / avgExact) : DBL_MAX;
    fprintf(m_options->_fpLogFile,
            "\n   Computing localError for var=%d, nMBs = %d, avg error = %g, "
            "avg exact = %g, tablesize = %lld entries",
            (int)var, (int)minibuckets.size(), (double)avgError,
            (double)avgExact, (int64)TableSize);
    if (rel_error < 1.0e+100)
      fprintf(m_options->_fpLogFile, ", rel avg error = %g%c",
              (double)rel_error, '%');
  }

  _Stats._LEMemorySizeMB += TableSize * sizeof(double) / (1024.0 * 1024);

  if (NULL != newTable) delete[] newTable;
  delete[] tuple;

  return 0;
}

int MiniBucketElimLH::computeLocalErrorTables(
    bool build_tables, double TotalMemoryLimitAsNumElementsLog,
    double TableMemoryLimitAsNumElementsLog) {
  vector<int> elimOrder;    // will hold dfs order
  findDfsOrder(elimOrder);  // computes dfs ordering of relevant subtree

  if (NULL != m_options ? NULL != m_options->_fpLogFile : false) {
    fprintf(m_options->_fpLogFile,
            "\n\nWILL COMPUTE LOCAL ERROR for each bucket ... error=(B value - "
            "(combine_all_MB : MB value))");
    fflush(m_options->_fpLogFile);
  }
  printf(
      "\n\nWILL COMPUTE LOCAL ERROR for each bucket ... error=(B value - "
      "(combine_all_MB : MB value))");

  deleteLocalErrorFNs();
  _BucketErrorFunctions.resize(m_problem->getN());
  _BuckerErrorFnTableSizes_Total = _BuckerErrorFnTableSizes_Precomputed =
      _BuckerErrorFnTableSizes_Ignored =
          -10.0;  // don't make -infinity; math will blow up.
  _nBucketsWithNonZeroBuckerError = _nBucketsWithMoreThan1MB = 0;

  // compute total of all bucket error table sizes
  for (vector<int>::reverse_iterator itV = elimOrder.rbegin();
       itV != elimOrder.rend(); ++itV) {
    int v = *itV;  // this is the variable being eliminated
    vector<MiniBucket> &minibuckets = _MiniBuckets[v];
    if (minibuckets.size() <= 1) continue;
    ++_nBucketsWithMoreThan1MB;
    set<int> &jointScope = _BucketScopes[v];
    double TableSize = 0.0;
    for (set<int>::const_iterator sit = jointScope.begin();
         sit != jointScope.end(); ++sit) {
      int u = *sit;
      if (u == v) continue;
      double dsLog = log10(m_problem->getDomainSize(u));
      if (dsLog < 0) {
        TableSize = -1.0;
        break;
      }
      TableSize += dsLog;
    }
    if (TableSize >= 0.0)
      _BuckerErrorFnTableSizes_Total =
          _BuckerErrorFnTableSizes_Total +
          log10(1.0 + pow(10.0, TableSize - _BuckerErrorFnTableSizes_Total));
  }
  if (NULL != m_options ? NULL != m_options->_fpLogFile : false)
    fprintf(m_options->_fpLogFile,
            "\n   BuckerErrorFnTableSizes total = %g, total_memory_limit = %g, "
            "total_memory_limit = %g",
            (double)_BuckerErrorFnTableSizes_Total,
            (double)TotalMemoryLimitAsNumElementsLog,
            (double)TableMemoryLimitAsNumElementsLog);
  printf(
      "\n   BuckerErrorFnTableSizes total = %g, total_memory_limit = %g, "
      "total_memory_limit = %g",
      (double)_BuckerErrorFnTableSizes_Total,
      (double)TotalMemoryLimitAsNumElementsLog,
      (double)TableMemoryLimitAsNumElementsLog);

  for (vector<int>::reverse_iterator itV = elimOrder.rbegin();
       itV != elimOrder.rend(); ++itV) {
    int v = *itV;  // this is the variable being eliminated
    double avgError, E;
    Function *errorFn = NULL;
    double tableSize = -1.0;
    double table_space_left =
        TotalMemoryLimitAsNumElementsLog > _BuckerErrorFnTableSizes_Precomputed
            ? TotalMemoryLimitAsNumElementsLog +
                  log10(1 - pow(10.0, _BuckerErrorFnTableSizes_Precomputed -
                                          TotalMemoryLimitAsNumElementsLog))
            : OUR_OWN_nInfinity;
    double table_size_actual_limit =
        table_space_left < TableMemoryLimitAsNumElementsLog
            ? table_space_left
            : TableMemoryLimitAsNumElementsLog;
    _BucketErrorQuality[v] = -1;
    // int res = computeLocalErrorTable(v, build_tables,
    // table_size_actual_limit, tableSize, avgError, E, errorFn) ;
    computeLocalErrorTable(v, build_tables, table_size_actual_limit, tableSize,
                           avgError, E, errorFn);
    _BucketErrorFunctions[v] = errorFn;
    if (_BucketErrorQuality[v] > 0) _nBucketsWithNonZeroBuckerError++;
    if (NULL != errorFn ? NULL != errorFn->getTable() : false) {
      //			if (NULL != m_options->_fpLogFile)
      //				fprintf(m_options->_fpLogFile, "\n
      //MiniBucketElimLH::localError var=%d, avgError=%g avgExact=%g", v, e, E)
      //;
      double *table = errorFn->getTable();
      if (tableSize > 0)
        //				_BuckerErrorFnTableSizes_Precomputed +=
        //tableSize ;
        _BuckerErrorFnTableSizes_Precomputed =
            _BuckerErrorFnTableSizes_Precomputed +
            log10(1.0 +
                  pow(10.0, tableSize - _BuckerErrorFnTableSizes_Precomputed));
      // compute standard deviation
      double stdDev = -DBL_MAX;
      size_t ts = errorFn->getTableSize();
      double sum = 0.0;
      for (size_t ti = 0; ti < ts; ti++) {
        double d = avgError - table[ti];
        double ds = d * d;
        sum += ds;
      }
#if defined DEBUG || _DEBUG
      if (isnan(sum)) {
        int bad = 1;
      }
#endif
      double variance = sum / ts;
      stdDev = sqrt(variance);
      if (NULL != m_options ? NULL != m_options->_fpLogFile : false)
        fprintf(m_options->_fpLogFile, ", stdDev = %g", (double)stdDev);
    } else if (tableSize > 0) {
      //			_BuckerErrorFnTableSizes_Ignored += tableSize ;
      _BuckerErrorFnTableSizes_Ignored =
          _BuckerErrorFnTableSizes_Ignored +
          log10(1.0 + pow(10.0, tableSize - _BuckerErrorFnTableSizes_Ignored));
      // printf("\n>>> %lld %lld", (int64) _BuckerErrorFnTableSizes_Ignored,
      // (int64) tableSize) ;
      if (_BucketErrorQuality[v] > 0) {
        // this would be bad; something is wrong. maybe table was too large.
      }
    }
  }
  if (NULL != m_options ? NULL != m_options->_fpLogFile : false) {
    fprintf(m_options->_fpLogFile,
            "\n   BuckerErrorFnTableSizes (precomputed/ignored/total) = "
            "%g/%g/%g entries",
            (double)_BuckerErrorFnTableSizes_Precomputed,
            (double)_BuckerErrorFnTableSizes_Ignored,
            (double)_BuckerErrorFnTableSizes_Total);
    fprintf(
        m_options->_fpLogFile,
        "\n   nBucketsWithNonZeroBuckerError (nMB>1/total) = %lld (%lld/%lld)",
        (int64)_nBucketsWithNonZeroBuckerError, (int64)_nBucketsWithMoreThan1MB,
        (int64)m_problem->getN());
    fprintf(m_options->_fpLogFile, "\n");
  }
  printf(
      "\n   BuckerErrorFnTableSizes (precomputed/ignored/total) = %g/%g/%g "
      "entries",
      (double)_BuckerErrorFnTableSizes_Precomputed,
      (double)_BuckerErrorFnTableSizes_Ignored,
      (double)_BuckerErrorFnTableSizes_Total);
  printf("\n   nBucketsWithNonZeroBuckerError (nMB>1/total) = %lld (%lld/%lld)",
         (int64)_nBucketsWithNonZeroBuckerError,
         (int64)_nBucketsWithMoreThan1MB, (int64)m_problem->getN());
  printf("\n");

  return 0;
}

}  // namespace daoopt
