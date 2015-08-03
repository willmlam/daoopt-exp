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

#include <vector>
#include <cstddef>
#include <fstream>
#include <chrono>

using namespace std::chrono;

#include "MiniBucketElimLH.h"

// 2015-07-17 KK : IGNORE_COPY_SUBTREE is applicable in a case when a LH subtree is a (proper/strict) subset of a LH subtree higher in a bucket tree.
// in this case, between the two nodes in the bucket/search tree, there is no information gain - the second (lower) LH subtree captures 
// exactly the same errors are the higher (earlier) LH subtree.
// #define IGNORE_COPY_SUBTREE

/* disables DEBUG output */
#undef DEBUG

/*
extern int64 nd1SpecialCalls;
extern int64 nd1GeneralCalls;
*/

//#define DEBUG_BUCKET_ERROR
//#define NO_LH_PREPROCESSING

#ifndef OUR_OWN_nInfinity
#define OUR_OWN_nInfinity (-std::numeric_limits<double>::infinity())
#endif // OUR_OWN_nInfinity
#ifndef OUR_OWN_pInfinity
#define OUR_OWN_pInfinity (std::numeric_limits<double>::infinity())
#endif // OUR_OWN_pInfinity

//#define USE_H_RESIDUAL

int daoopt::MBLHSubtreeNode::ComputeidxVarMapping(std::vector<val_t> & assignment)
{
	if (0 == _RelevantFunctions.size()) 
		{ _idxVarMappingComputed = true ; return 0 ; }
	if (NULL == _OutputFunction) 
		return 1 ;
	_idxVarMapping.clear() ;
	_idxVarMapping.reserve(_RelevantFunctions.size()) ;
	_idxVarMapping.resize(_RelevantFunctions.size()) ;
	val_t *assignment_data = assignment.data() ;
	// compute ptr2value mapping for each function
	int i = 0 ;
	for (vector<Function*>::const_iterator fit = _RelevantFunctions.begin(); fit != _RelevantFunctions.end(); ++fit, i++) {
		Function *f = *fit ;
		std::vector<val_t*> & ptr2value = _idxVarMapping[i] ;
		const std::vector<int> & scope = f->getScopeVec() ;
		ptr2value.clear() ;
		ptr2value.reserve(scope.size()) ;
		for (int j = 0 ; j < scope.size() ; j++)
			ptr2value.push_back(assignment_data + scope[j]) ;
		}
	// compute ptr2value mapping for output function scope
	_idxOutputFunctionScopeMapping.clear() ;
	_idxOutputFunctionScopeMapping.reserve(_OutputFunctionScope.size()) ;
	const std::vector<int> & scopeOF = _OutputFunction->getScopeVec() ;
	for (int j = 0 ; j < scopeOF.size() ; j++)
		_idxOutputFunctionScopeMapping.push_back(assignment_data + scopeOF[j]) ;
	_idxVarMappingComputed = true ;
	return 0 ;
}

int daoopt::MBLHSubtreeNode::PrepOutputFunction(std::set<int> & InstantiatedVariables)
{
	daoopt::Problem *problem = _H->getProblem() ;

	// compute output fn scope
	_OutputFunctionScope.clear() ;
	set<int> node_vars ; // note we cannot use _H->BucketScope(v) since here we may have more vars since the subtree does exact computation, not MB computation
	for (vector<Function*>::iterator itF = _RelevantFunctions.begin(); itF != _RelevantFunctions.end(); ++itF)
		node_vars.insert((*itF)->getScopeVec().begin(), (*itF)->getScopeVec().end()) ;
	node_vars.erase(_v) ;
	for (set<int>::const_iterator itI = node_vars.begin() ; itI != node_vars.end() ; ++itI) {
		int u = *itI ;
		const bool is_in = InstantiatedVariables.find(u) != InstantiatedVariables.end() ;
		if (! is_in) 
			_OutputFunctionScope.insert(u) ;
		}
	int n = _OutputFunctionScope.size() ;

//	std::sort(_OutputFunctionScope.begin(), _OutputFunctionScope.end()) ; // make sure it is sorted
	_OutputFunctionScopeDomains.reserve(n) ;
	_OutputFunctionSize = 1 ;
	// we assume here that when we walk through the set using begin/end, we get an increasing list of elements.
	for (set<int>::const_iterator sit = _OutputFunctionScope.begin(); sit!=_OutputFunctionScope.end(); ++sit) {
		_OutputFunctionSize *= problem->getDomainSize(*sit) ;
		_OutputFunctionScopeDomains.push_back(problem->getDomainSize(*sit)) ;
		}
	double *table = new double[_OutputFunctionSize] ;
	if (NULL == table) 
		return 1 ;
	_OutputFunction = new FunctionBayes(-1000000-_v,problem,_OutputFunctionScope,table,_OutputFunctionSize) ;
	if (NULL == _OutputFunction) 
		{ delete table ; return 1 ; }
	return 0 ;
}

bool daoopt::MBLHSubtreeNode::IsSubtreeMBEfunction(Function & f, std::stack<MBLHSubtreeNode *> & dfsHelper)
{
	while (dfsHelper.size() > 0) dfsHelper.pop() ;
	std::vector<std::vector<MiniBucket>> & MBlist = _H->MiniBuckets() ;
	dfsHelper.push(this) ;
	while (dfsHelper.size() > 0) {
		MBLHSubtreeNode *N = dfsHelper.top() ; dfsHelper.pop() ;
		for (MBLHSubtreeNode *n : N->_Children) {
			int u = n->_v ;
			std::vector<MiniBucket> & minibuckets = MBlist[u] ;
			for (const MiniBucket & mb : minibuckets) {
				Function *fnMB = mb.output_fn() ;
				if (&f == fnMB) return true ;
				}
			dfsHelper.push(n) ;
			}
		}
	return false ;
}

int daoopt::MBLHSubtreeNode::Delete(void) 
{
	_H = NULL ;
	_RootVar = -1 ;
	_v = -1 ;
	_Parent = NULL ;
	_idxWRTparent = -1 ;
	_k = -1 ;
	_depth2go = -1 ;
	for (vector<MBLHSubtreeNode *>::iterator it = _Children.begin(); it!=_Children.end(); ++it) 
		delete *it ;
	_Children.clear() ;
	_RelevantFunctions.clear() ;
	_OutputFunctionScope.clear() ;
	_OutputFunctionScope.clear() ;
	_OutputFunctionScopeDomains.clear() ;
	_OutputFunctionSize = -1 ;
	if (NULL != _OutputFunction) { delete _OutputFunction ; _OutputFunction = NULL ; }
	_idxVarMappingComputed = false ;
	_idxVarMapping.clear() ;
	_idxOutputFunctionScopeMapping.clear() ;
	return 0 ;
}

int daoopt::MBLHSubtreeNode::ComputeOutputFunction(std::vector<val_t> & assignment)
{
	if (NULL == _OutputFunction) 
		return 0 ;
	if (! _idxVarMappingComputed) {
		int res = ComputeidxVarMapping(assignment) ;
		if (0 != res) 
			return res ;
		}

	// if the node has no instantiated variables, we need to compute it just once
	if (0 == _nVarsInstantiated) {
		if (_nTimesComputed > 0) 
			return 0 ;
		}

	++_nTimesComputed ;

	// enumerate over all combinations of _OutputFunctionScope[] values; for each eliminate _v by max
	int i, k, n = _idxOutputFunctionScopeMapping.size() ;
	if (n > 0) {
		for (i = n-2 ; i >= 0 ; i--)  
			*(_idxOutputFunctionScopeMapping[i]) = 0 ;
		*(_idxOutputFunctionScopeMapping[n-1]) = -1 ;
		}
	double *output_table = _OutputFunction->getTable(), z ;
	for (size_t j = 0 ; j < _OutputFunctionSize ; j++) {
		// move to next combination of output fn scope value combination
		for (i = n - 1 ; i >= 0 ; i--) {
			if (++(*(_idxOutputFunctionScopeMapping[i])) < _OutputFunctionScopeDomains[i]) break ;
			*(_idxOutputFunctionScopeMapping[i]) = 0 ;
			}
		output_table[j] = ELEM_ZERO ;
		// eliminate var of this node
		for (k = 0 ; k < _k ; k++) {
			assignment[_v] = k ;
			z = ELEM_ONE ;
			for (i = _RelevantFunctions.size()-1 ; i >= 0 ; i--) 
				z OP_TIMESEQ _RelevantFunctions[i]->getValuePtr(_idxVarMapping[i]) ;
			if (z > output_table[j]) output_table[j] = z ;
			}
		}
	return 0 ;
}

bool daoopt::MBLHSubtree::IsSubtreeVariable(int v)
{
	for (vector<MBLHSubtreeNode *>::iterator itRB = _SubtreeNodes.begin(); itRB!=_SubtreeNodes.end(); ++itRB) 
		{ if (v == (*itRB)->_v) return true ; }
	return false ;
}

int daoopt::MBLHSubtree::GetAncestors(bool IncludeRootVar, std::set<int> & Ancestors)
{
	Ancestors.clear() ;
	daoopt::Pseudotree *pseudotree = _H->getPseudotree() ;
	int v = _RootVar ;
	do {
		if (IncludeRootVar ? true : v != _RootVar) Ancestors.insert(v) ;
		const PseudotreeNode *n = pseudotree->getNode(v) ;
		if (NULL == n) break ;
		const PseudotreeNode *p = n->getParent() ;
		if (NULL == p) break ;
		v = p->getVar() ;
		} while (v >= 0) ;
	return 0 ;
}

int daoopt::MBLHSubtree::ComputeOutputFnIndependence(std::set<int> & InstantiatedVariables)
{
	daoopt::Pseudotree *pseudotree = _H->getPseudotree() ;
	std::vector<std::vector<MiniBucket>> & MBlist = _H->MiniBuckets() ;

	_nSubtreeNodesIndependentOfContext = 0 ;
	// note : traverse the subtree bottom up, so that when a node is processed, all children are processed.
	// this is important, because it sometimes may be that this node is independent in that it contains no context variables, 
	// but some child node may contain (i.e. child node is not independent), then this node is also not independent.
	// to determine that, we need to process all children first, so that by the time we process this node we know the status of all children.
	for (vector<MBLHSubtreeNode *>::reverse_iterator itRB = _SubtreeNodes.rbegin(); itRB!=_SubtreeNodes.rend(); ++itRB) {
		MBLHSubtreeNode *n = *itRB ;
		n->_nVarsInstantiated = 0 ;
		int v = n->_v ;
		const PseudotreeNode *bucket = pseudotree->getNode(v) ;
		if (NULL == bucket) continue ;
		set<int> node_vars ; // note we cannot use _H->BucketScope(v) since here we may have more vars since the subtree does exact computation, not MB computation
		for (vector<Function*>::iterator itF = n->_RelevantFunctions.begin(); itF != n->_RelevantFunctions.end(); ++itF)
			node_vars.insert((*itF)->getScopeVec().begin(), (*itF)->getScopeVec().end()) ;
		for (set<int>::const_iterator itI = node_vars.begin() ; itI != node_vars.end() ; ++itI) {
			int u = *itI ;
			const bool is_in = InstantiatedVariables.find(u) != InstantiatedVariables.end() ;
			if (is_in) 
				(n->_nVarsInstantiated)++ ;
			}
		// test how many of output fn scope are instantiated; should be 0.
		int nOFinstantiated = 0 ;
		for (set<int>::const_iterator itI = n->_OutputFunctionScope.begin() ; itI != n->_OutputFunctionScope.end() ; ++itI) {
			int u = *itI ;
			const bool is_in = InstantiatedVariables.find(u) != InstantiatedVariables.end() ;
			if (is_in) 
				nOFinstantiated++ ;
			}
		// note : node_vars = (n->_nVarsInstantiated) + _OutputFunctionScope + 1
		int nTest = 1 + n->_nVarsInstantiated + n->_OutputFunctionScope.size() ;
		assert(node_vars.size() == nTest) ;
		// we have a problem right now : it could be that at a subtree leaf, we have some instantiated variables, but the output fn of this node has no instantiated variables, so its parent node may think that it has no instantiated variables, and compute it only the very first time, and not compute during later calls, leading to incorrect values.
		// fix this, if a child of a node has nInstantiated>0, then mark the node as having instantiated variables.
		if (0 == n->_nVarsInstantiated) {
			for (MBLHSubtreeNode *child : n->_Children) {
				if (child->_nVarsInstantiated > 0) 
					{ n->_nVarsInstantiated = INT_MAX ; break ; }
				}
			}
		if (0 == n->_nVarsInstantiated) 
			_nSubtreeNodesIndependentOfContext++ ;
		}
	// test : if a node has no instantiated variables, all children also must have instantiated variables
	for (vector<MBLHSubtreeNode *>::iterator itRB = _SubtreeNodes.begin(); itRB!=_SubtreeNodes.end(); ++itRB) {
		MBLHSubtreeNode *n = *itRB ;
		if (n->_nVarsInstantiated > 0) continue ;
		for (MBLHSubtreeNode *child : n->_Children) {
			if (child->_nVarsInstantiated > 0) {
				int bug = 1 ;
				}
			}
		}
	return 0 ;
}

int daoopt::MBLHSubtree::Delete(void) 
{
	_RootNode.Delete() ;
	_SubtreeNodes.clear() ;
	_IsCopyOfEarlierSubtree = NULL ;
	return 0 ;
}

int daoopt::MBLHSubtree::PrepOutputFunctions(std::set<int> & InstantiatedVariables)
{
	// prep output functions of all nodes; do in reverse order (since nodelist is DFS sorted this should be ok); add output fn to parent.
	for (vector<MBLHSubtreeNode *>::reverse_iterator itRB = _SubtreeNodes.rbegin(); itRB!=_SubtreeNodes.rend(); ++itRB) {
		MBLHSubtreeNode *n = *itRB ;
		n->PrepOutputFunction(InstantiatedVariables) ;
		if (NULL == n->_OutputFunction) 
			{ Delete() ; return 1 ; }
		if (NULL != n->_Parent) 
			((n->_Parent)->_RelevantFunctions).push_back(n->_OutputFunction) ;
		}
	return 0 ;
}

int daoopt::MBLHSubtree::FillInSubtreeNodeRelevantFunctions(void)
// fill in _RelevantFunctions of each LHsubtreenode
{
	daoopt::Pseudotree *pseudotree = _H->getPseudotree() ;
	std::stack<MBLHSubtreeNode *> dfsHelper ;
	for (vector<MBLHSubtreeNode *>::iterator itRB = _SubtreeNodes.begin(); itRB!=_SubtreeNodes.end(); ++itRB) {
		MBLHSubtreeNode *n = *itRB ;
		int u = n->_v ;
		std::vector<Function *> & funs = n->_RelevantFunctions ;
		// add OF_n
		const vector<Function *> & OFlist = pseudotree->getFunctions(u) ;
		funs.insert(funs.end(), OFlist.begin(), OFlist.end()) ;
		// add AF_n^{>d}; i.e. those functions in AF_n that came from a bucket below the LHsubtree.
		vector<Function *> & AFlist = _H->m_augmented[u] ;
		for (vector<Function*>::iterator itF = AFlist.begin() ; itF != AFlist.end() ; ++itF) {
			Function *f = *itF ;
			// note : we could just check what the originating bucket of 'fn' is; MBE-generated name of 'fn' should say that, but we are not sure the names are reliable (i.e. a MBE-generated by var 'v' fn name was supposed to be -v, but what if v=0???).
			if (! n->IsSubtreeMBEfunction(*f, dfsHelper)) 
				funs.push_back(f) ;
			}
		}
	return 0 ;
}

int daoopt::MBLHSubtree::FillInRootNodeRelevantFunctions(void)
// fill in _IntermediateSubtreeFunctions - IF/AF (in bucket of _RootVar) that came from outside of the LH-subteee
{
	daoopt::Pseudotree *pseudotree = _H->getPseudotree() ;
	std::stack<MBLHSubtreeNode *> dfsHelper ;
	std::vector<Function *> & rootvarfuns = _RootNode._RelevantFunctions ;
	vector<Function *> & rootvarAFlist = _H->m_augmented[_RootVar] ;
	for (vector<Function*>::iterator itF = rootvarAFlist.begin() ; itF != rootvarAFlist.end() ; ++itF) {
		Function *f = *itF ;
		if (! _RootNode.IsSubtreeMBEfunction(*f, dfsHelper)) 
			rootvarfuns.push_back(f) ;
		}
	vector<Function *> & rootvarIFlist = _H->m_intermediate[_RootVar] ;
	for (vector<Function*>::iterator itF = rootvarIFlist.begin() ; itF != rootvarIFlist.end() ; ++itF) {
		Function *f = *itF ;
		if (! _RootNode.IsSubtreeMBEfunction(*f, dfsHelper)) 
			rootvarfuns.push_back(f) ;
		}
	return 0 ;
}

int daoopt::MBLHSubtree::ComputeSubtree(void)
{
	daoopt::Problem *problem = _H->getProblem() ;
	daoopt::Pseudotree *pseudotree = _H->getPseudotree() ;
	daoopt::PseudotreeNode *nRootVar = pseudotree->getNode(_RootVar) ;

	// allocate memory for subtree nodes; assuming avg branching factor (e.g. 5) is problematic, since the memory would grow exponentially with d.
	int memory = _depth*2 ; if (memory > problem->getN()) memory = problem->getN() ; // assume two nodes at each level
	_SubtreeNodes.clear() ;
	_SubtreeNodes.reserve(memory) ; // try to allocate some memory to make it faster

	// build LH subtree of nodes
//	std::deque<MiniBucketElimLHobsoleteErrorNode *> bfsLHsubtreeTraversal ;
	std::stack<MBLHSubtreeNode *> dfsLHsubtreeTraversal ; // using stack will turn it into a DFS processing of the tree
	try { dfsLHsubtreeTraversal.push(&_RootNode) ; } catch (...) { return 1 ; }
	while (dfsLHsubtreeTraversal.size() > 0) {
		MBLHSubtreeNode *N = dfsLHsubtreeTraversal.top() ; dfsLHsubtreeTraversal.pop() ;
		const PseudotreeNode *N_ = pseudotree->getNode(N->_v) ;
		const vector<PseudotreeNode *> & children = N_->getChildren() ;
		N->_Children.reserve(children.size()) ; // reserver space for all children (some will not be added as children) so that reallocation is not needed
		int depth2go = N->_depth2go - 1 ; // remaining depth from each child
		for (vector<PseudotreeNode*>::const_iterator itC = children.begin() ; itC != children.end(); ++itC) {
			int child = (*itC)->getVar() ;
#ifndef USE_FULL_LOOKAHEAD_SUBTREE
			if (_H->_BucketErrorQuality[child] <= 1 ? _H->_distToClosestDescendantWithLE[child] > depth2go : false) continue ; // _BucketErrorQuality <= 1 means it is not proven that there is substantial bucket error
#endif // USE_FULL_LOOKAHEAD_SUBTREE
			MBLHSubtreeNode *n = NULL ; try { n = new MBLHSubtreeNode ; } catch (...) { return 1 ; } if (NULL == n) return 1 ;
			n->_H = _H ; n->_RootVar = _RootVar ; n->_v = child ; n->_Parent = N ; n->_k = problem->getDomainSize(child) ; n->_depth2go = depth2go ;
			try { _SubtreeNodes.push_back(n) ; } catch (...) { delete n ; return 1 ; }
			try { (N->_Children).push_back(n) ; } catch (...) { delete n ; return 1 ; }
			if (depth2go > 0) { try { dfsLHsubtreeTraversal.push(n) ; } catch (...) { delete n ; return 1 ; }}
			n->_idxWRTparent = N->_Children.size() - 1 ;
/*			// scope of output function of n is scope of its parent + its parent's variable
			MBLHSubtreeNode *p = n->_Parent ;
			int vParent = NULL != p ? p->_v : -1 ;
			if (vParent >= 0 && vParent != _RootVar) {
				n->_OutputFunctionScope = p->_OutputFunctionScope ;
				(n->_OutputFunctionScope).insert(vParent) ;
				}*/
			}
		}

	// compute signature for each subtree node
	_RootNode.SerializeSignature() ;
	for (vector<MBLHSubtreeNode *>::iterator itST = _SubtreeNodes.begin(); itST != _SubtreeNodes.end(); ++itST) 
		(*itST)->SerializeSignature() ;

	// check if this subtree is a subset of the subtree of some ancestor
	std::vector<MBLHSubtree> & lhArray = _H->LH() ;
	daoopt::PseudotreeNode *p = nRootVar->getParent() ;
	_IsCopyOfEarlierSubtree = NULL ;
	int dUp = 1 ;
	while (NULL != p) {
		// check if subtree of p contains signature of this subtree
		int v = p->getVar() ;
		MBLHSubtree & lhV = lhArray[v] ;
		MBLHSubtreeNode *n = lhV.GetSubtreeNode(_RootVar) ;
		if (NULL == n) 
			// subtree of v does not contain node '_RootVar'.
			break ;
		if (_RootNode._Signature != n->_Signature) 
			// signature of '_RootVar' in n does not match signature of this subtree; we cannot use n to compute the value of this subtree; quit.
			break ;
		// we can use n to compute the value of this subtree, since n and this subtree match!!!
		_IsCopyOfEarlierSubtree = n ;
		// go up on level and check again
		p = p->getParent() ;
		++dUp ;
		}
	if (NULL != _IsCopyOfEarlierSubtree) {
    if (_H->m_options->_fpLogFile) {
      fprintf(_H->m_options->_fpLogFile, "\nMATCH : subtree of rootvar %d is a subset of subtree of %d at distance=%d", _RootVar, _IsCopyOfEarlierSubtree->_RootVar, dUp) ;
      fflush(_H->m_options->_fpLogFile) ;
    }
  }

	return 0 ;
}

int daoopt::MBLHSubtree::Initialize(MiniBucketElimLH & H, int v, int depth)
{
	Delete() ;

	_H = &H ; _RootVar = v ; _depth = depth ; _RootNode._H = &H ; _RootNode._v = v ; _RootNode._k = (H.getProblem())->getDomainSize(v) ; _RootNode._depth2go = depth ;

	// always check if there is any point in computing error
	if (H._distToClosestDescendantWithLE[v] > depth) 
		return 0 ;

	// create subtree
	if (0 != ComputeSubtree()) 
		return 1 ;

	if (NULL == _IsCopyOfEarlierSubtree) {
		// fill in _RelevantFunctions of each LHsubtreenode
		if (0 != FillInSubtreeNodeRelevantFunctions()) 
			return 1 ;

		// compute instantiated variables for this subtree
		std::set<int> instantiated_variables ;
		if (0 != GetAncestors(false, instantiated_variables)) 
			return 1 ;

		// for each LH subtree node, bottom-up, compute output fn scope and create/allocate space for the output fn
		if (0 != PrepOutputFunctions(instantiated_variables)) 
			return 1 ;

		// figure out which subtree nodes are independent of the context; they can be computed during search just once.
		if (0 != ComputeOutputFnIndependence(instantiated_variables)) 
			return 1 ;
		}
	else {
		// since this node is a copy of a subtree of ancestor LH subtree, add the output functions of the children of the original node (whose copy this subtree is) to the list of relevant functions of this subtree
		for (MBLHSubtreeNode *child : _IsCopyOfEarlierSubtree->_Children) {
			if (NULL != child->_OutputFunction) 
				_RootNode._RelevantFunctions.push_back(child->_OutputFunction) ;
			}
		}

	// fill in _IntermediateSubtreeFunctions - IF/AF (in bucket of _RootVar) that came from outside of the LH-subteee
	if (0 != FillInRootNodeRelevantFunctions()) 
		return 1 ;

	return 0 ;
}

void daoopt::MBLHSubtree::ComputeHeuristic(std::vector<val_t> & assignment)
{
	// compute subtree output functions, bottom-up, only if this LH subtree is not a copy of a (subtree of) an earlier (higher in the bucket tree) LH subtree
	if (NULL == _IsCopyOfEarlierSubtree) {
		++_RootNode._nTimesComputed ;

		// compute output functions of all nodes, back-to-front; since _SubtreeNodes[] is in order where u<v means v is a descendant of u, this is ok.
		for (int i = _SubtreeNodes.size()-1 ; i >= 0 ; i--) {
			MBLHSubtreeNode *n = _SubtreeNodes[i] ;
			n->ComputeOutputFunction(assignment) ;
			}
		}
}

double daoopt::MBLHSubtree::GetHeuristic(std::vector<val_t> & assignment)
{
	if (0 == _SubtreeNodes.size()) 
		return _H->MiniBucketElim::getHeur(_RootVar, assignment, NULL) ; // most likely, all the children(buckets) up to depth d have exactly 1 mini-bucket.

#ifdef IGNORE_COPY_SUBTREE
	if (NULL != _IsCopyOfEarlierSubtree) {
		return _H->MiniBucketElim::getHeur(_RootVar, assignment, NULL) ; // most likely, all the children(buckets) up to depth d have exactly 1 mini-bucket.
		}
#endif // IGNORE_COPY_SUBTREE

	// add up all relevant functions of the root var. they contain all functions of the root_var bucket that came from outside the LH subtree
	// as well as output functions of all child nodes in the look-ahead subtree.
	double h = ELEM_ONE ;
	vector<Function*>::const_iterator itF_end = _RootNode._RelevantFunctions.end() ;
	for (vector<Function*>::const_iterator itF = _RootNode._RelevantFunctions.begin() ; itF != itF_end ; ++itF) {
		daoopt::Function *f = *itF ;
		h OP_TIMESEQ f->getValue(assignment) ;
		if (OUR_OWN_nInfinity == h) 
			return OUR_OWN_nInfinity ; // no point in continuing with the computation
		}

	return h ;
}

int daoopt::SetupLookaheadStructure(daoopt::MiniBucketElimLH & H, int Depth)
{
	daoopt::Problem *problem = H.getProblem() ;
	daoopt::Pseudotree *pseudotree = H.getPseudotree() ;
	daoopt::ProgramOptions *options = H.getProgramOptions() ;
	std::vector<MBLHSubtree> & lhArray = H.LH() ;

	// compute order of variables, so that a parent is before its children wrt the bucket tree; this can be used to traverse the bucket tree top-down or bottom-up.
	std::vector<int> btOrder ;
	btOrder.reserve(problem->getN()) ;
	const PseudotreeNode *ptRoot = pseudotree->getRoot() ;
	btOrder.push_back(ptRoot->getVar()) ;
	int i = 0 ; // how many nodes in btOrder have been processed
	while (i < btOrder.size()) {
		int v = btOrder[i++] ;
		const PseudotreeNode *n = pseudotree->getNode(v) ;
		const vector<PseudotreeNode *> & children = n->getChildren() ;
		for (vector<PseudotreeNode*>::const_iterator itC = children.begin() ; itC != children.end(); ++itC) {
			int child = (*itC)->getVar() ;
			btOrder.push_back(child) ;
			}
		}

	// create LH subtree for all nodes; process a parent before its children so that a child can use parents data.
	for (auto itV = btOrder.begin(); itV != btOrder.end(); ++itV) {
		int v = *itV;
		int depth = Depth ;
		int lhInitRes = lhArray[v].Initialize(H, v, Depth) ;
		if (0 != lhInitRes) {
			printf("\n\nERROR : lookahead init failed; v=%d", v) ;
      if (options->_fpLogFile) {
        fprintf(options->_fpLogFile, "\n\nERROR : lookahead init failed; v=%d", v) ;
        fflush(options->_fpLogFile) ;
      }
			exit(1) ;
			}
		}

	return 0 ;
}

