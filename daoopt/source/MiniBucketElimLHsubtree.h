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
 *  along with DAOOPT.  If not, see <http://www.gnu.org/licenses/>.
 *  
 *  Created on: Nov 8, 2008
 *      Author: Lars Otten <lotten@ics.uci.edu>
 */

#ifndef MiniBucketElimLHsubtree_H_
#define MiniBucketElimLHsubtree_H_

#include "MiniBucketElim.h"

#undef DEBUG

#ifndef OUR_OWN_nInfinity
#define OUR_OWN_nInfinity (-std::numeric_limits<double>::infinity())
#endif // OUR_OWN_nInfinity

namespace daoopt {

class MiniBucketElimLH ;

// this class represents a node in the look-ahead subtree in the mini-bucket tree.
class MBLHSubtreeNode
{
public :
	// MBE heuristic with LH
	MiniBucketElimLH *_H ;
	// root variable of the LH subtree that this node belongs to
	int _RootVar ;
	// variable of this node.
	int _v ;
	// parent of this LHsubtree node
	MBLHSubtreeNode *_Parent ;
	// idx of this node as child of its parent; this is wrt the childlen list of the parent.
	int _idxWRTparent ;
	// domain size of this variable.
	int _k ;
	// depth still to go from this node; as number of edges. this is used when depth-based LH subtree is built.
	int _depth2go ;
	// children of this node in the bucket tree, that have non-0 bucket error.
	std::vector<MBLHSubtreeNode *> _Children ;
	// signature of this node, serialized as a string. signature is a concatenated string of IDs of all nodes in the subtree rooted at this node, separated by ';'. 
	std::string _Signature ;
	// for nodes (except rootnode) _RelevantFunctions[] contains 3 kinds of functions, 
	// 1) original problem specification functions (places in the bucket of this var(node)),
	// 2) MB-generated functions that contain this var and come from outside (below) the subtree that this node is part of,
	// 3) output functions of the children of this node (i.e. child nodes of this node).
	// these 3 kinds of functions go into the computation of the output function of this node; this mimics the bucket-tree computation of a bucket.
	// for rootnode, we don't have 1) and instead have 
	// 4) MB-generated functions that don't contain this var and come from outside (below) the subtree that this node is part of (so called m_intermediate[thisvar] function in the language of daoopt library).
	std::vector<Function *> _RelevantFunctions ;
	// helper data used when OutputFunction is computed :
	bool _idxVarMappingComputed ;
	std::vector<std::vector<val_t*>> _idxVarMapping ; // for each _RelevantFunctions[i], mapping to ptr2value for each argument. mapping is wrt "std::vector<val_t> & assignment" used during search.
	std::vector<val_t*> _idxOutputFunctionScopeMapping ; // ptr2value for all output function arguments
	// output function scope; ordered so that if u and v are in the scope, and u is before v, u is an ancestor of v.
	std::set<int> _OutputFunctionScope ;
	// output function scope domain sizes
	std::vector<val_t> _OutputFunctionScopeDomains ;
	// output function table size
	size_t _OutputFunctionSize ;
	// output function; computed during look-ahead during search; then computed, we assume all ancestor variables (of the root variable of the LH subtree) + the root variable are instantiated.
	daoopt::Function *_OutputFunction ;
	// number of times this node is computed
	int64 _nTimesComputed ;
	// number of variables in this node that are instantiated, i.e. in the context.
	// if this number is 0, then this node needs to be computed just once. also, if it is 0, then all its descendants must have 0 too.
	int _nVarsInstantiated ;

  // marks whether the subtree node's output function is valid for the current
  // context
  bool _IsValidForCurrentContext;

public :
	// signature is a set of IDs of all nodes in the subtree rooted at this node, separated by ';'. 
	// note we will not clear/erase the input Signature array, so that it can be used in recursive calls.
	inline int ComputeSignature(std::set<int> & Signature)
	{
		Signature.insert(_v) ;
		for (int i = 0 ; i < _Children.size() ; i++) {
			MBLHSubtreeNode *c = _Children[i] ;
			c->ComputeSignature(Signature) ;
			}
		return 0 ;
	}
	inline int SerializeSignature(void)
	{
		std::set<int> sig ;
		ComputeSignature(sig) ;
		_Signature.erase() ;
		for (set<int>::const_iterator itI = sig.begin() ; itI != sig.end() ; ++itI) {
			int id = *itI ;
			if (_Signature.length() > 0) _Signature += ';' ;
			_Signature += std::to_string(id) ;
			}
		return 0 ;
	}
	int ComputeOutputFunction(std::vector<val_t> & assignment) ;
	int ComputeidxVarMapping(std::vector<val_t> & assignment) ;
	int PrepOutputFunction(std::set<int> & InstantiatedVariables) ;
	// return true iff the given function was generated by a minibucket in a node (bucket) in a subtree rooted at this node.
	bool IsSubtreeMBEfunction(Function & f, std::stack<MBLHSubtreeNode *> & dfsHelper) ;
	int Delete(void) ;

	inline MBLHSubtreeNode(void) : _H(NULL), _RootVar(-1), _v(-1), _Parent(NULL), _idxWRTparent(-1), _k(-1), _depth2go(-1), _OutputFunctionSize(-1), _OutputFunction(NULL), _idxVarMappingComputed(false), _nTimesComputed(0), _nVarsInstantiated(-1), _IsValidForCurrentContext(false)
	{
	}
	inline ~MBLHSubtreeNode(void)
	{
		Delete() ;
	}
} ;

class MBLHSubtree
{
public :
	// MBE heuristic with LH
	MiniBucketElimLH *_H ;
	// root variable of the LH subtree
	int _RootVar ;
	// depth to go from root; used when depth-based LH is done; this is the depth of the LH. if -1 then depth-based LH is not used.
	int _depth ;
	// limit on the size of the LH subtree. if -1, then size-based LH is not used.
	int _sizelimit ;
public :
	// node of the root variable
	// RelevantFunctions of the root variable are all augmented/intermediate functions, of the root variable bucket, that came from below the subtree.
	// here we also have output functions of all child nodes, so that in order to compute look-ahead heuristic for the root var, all we need to do is combine relevant functions of the root node.
	MBLHSubtreeNode _RootNode ;
	// this LH subtree may be a copy of a subtree of a LH subtree computed earlier (i.e. at an ancestor) variable. 
	MBLHSubtreeNode *_IsCopyOfEarlierSubtree ;
	// variables/buckets in the (minimal) subtree. they are here for reference. each node belongs to its parent (and should be deleted by the parent).
	// note that RootVar(Node) does not belong to the SubtreeNodes list.
	// note : we assume that _SubtreeNodes entries are stored in DFS order.
	// 2015-10-18 KK : actually, we assume that in _SubtreeNodes[], parent is before its children. not full DFS order; e.g. BFS order is fine too.
	std::vector<MBLHSubtreeNode *> _SubtreeNodes ;
	// number of subtree nodes that have no context variables.
	int _nSubtreeNodesIndependentOfContext ;

  // this is used to mark if lookahead was performed for the current context,
  // so this subtree would have a compatible instantiation with its descendants.
  bool _ComputedForCurrentContext;
public :
	// return true iff if the given variable is the variable of a node in _SubtreeNodes.
	bool IsSubtreeVariable(int v) ;
	// compute number of subtree nodes evaluated
	inline int64 nSubtreeNodesEvaluated(void)
	{
		int64 nComputations = _RootNode._nTimesComputed ;
		for (int i = _SubtreeNodes.size()-1 ; i >= 0 ; i--) {
			MBLHSubtreeNode *n = _SubtreeNodes[i] ;
			nComputations += n->_nTimesComputed ;
			}
		return nComputations ;
	}
	inline MBLHSubtreeNode *GetSubtreeNode(int var)
	{
		for (vector<MBLHSubtreeNode *>::iterator it = _SubtreeNodes.begin(); it != _SubtreeNodes.end(); ++it) {
			MBLHSubtreeNode *n = *it ;
			if (var == n->_v) 
				return n ;
			}
		return NULL ;
	}
	// a set of variables, from the parent of the RootVar to the root of the bucket tree.
	// we will use this to check if a subtree node depends on the ancestors of not; that determines if the subtree nodes can be evaluated just once (if it does not depend on ancestors).
	int GetAncestors(bool IncludeRootVar, std::set<int> & Ancestors) ;
public :
	// compute look-ahead heuristic for the given variable assignment, wrt root variable's subtree.
	void ComputeHeuristic(std::vector<val_t> & assignment) ;

  // compute the look-ahead heuristic for the given variable assignment, but 
  // only for the part of the subtree rooted by sub_root.
  void ComputeHeuristicSubset(std::vector<val_t>& assignment,
      MBLHSubtreeNode* sub_root);
	// get look-ahead heuristic for the given variable assignment, wrt root variable's subtree.
	double GetHeuristic(std::vector<val_t> & assignment) ;
public :
	// fill in RelevantFunctions array
	int FillInSubtreeNodeRelevantFunctions(void) ;
	int FillInRootNodeRelevantFunctions(void) ;
	// create space for the output function
	int PrepOutputFunctions(std::set<int> & InstantiatedVariables) ;
	// compute which output nodes are independent of the context (instantiated variables = ancestors of the root variable).
	int ComputeOutputFnIndependence(std::set<int> & InstantiatedVariables) ;
	// compute the subtree, i.e. the array of subtree nodes.
	int ComputeSubtree(void) ;
	// set up a lookahead subtree for the given bucket tree variable.
	int Initialize(MiniBucketElimLH & H, int v, int depth, int sizelimit) ;

  // note that all of the subtree nodes are no longer valid for the current
  // context
  void Invalidate();
public :
	int Delete(void) ;

	inline MBLHSubtree(void) : _H(NULL), _RootVar(-1), _depth(-1), _sizelimit(-1), _IsCopyOfEarlierSubtree(NULL), _nSubtreeNodesIndependentOfContext(0)	{
	}
	inline ~MBLHSubtree(void)
	{
		Delete() ;
	}
} ;

int SetupLookaheadStructure(MiniBucketElimLH & H, int Depth, int SizeLimit) ;

// here we mark (as buckets with error), after sorting buckets by BE, certain percentage of buckets with largest bucket error.
int SetupLookaheadStructure_FractionOfLargestAbsErrorNodesOnly(MiniBucketElimLH & H, int Depth, int SizeLimit) ;

}  // namespace daoopt

#endif // MiniBucketElimLHsubtree_H_
