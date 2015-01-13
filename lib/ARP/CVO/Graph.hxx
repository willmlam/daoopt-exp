#ifndef ARE_Graph_HXX_INCLUDED
#define ARE_Graph_HXX_INCLUDED

#include <stdlib.h>
#include <stdio.h>
#include <string>
#include <set>

#include "Utils/Heap.hxx"
#include "Utils/AVLtreeSimple.hxx"
#include "Utils/AVLtree.hxx"
#include "Utils/MersenneTwister.h"

#include "Problem/Problem.hxx"

#define ARE_GRAPH_VAR_ORDER_COMP_MAX_NUM_EDGES_ADDED 100000

namespace ARE
{

class ARP ;

class AdjVar
{
public :
	int _V ;
	int _IterationEdgeAdded ;
	AdjVar *_NextAdjVar ;
public :
	AdjVar(void) : _V(-1), _IterationEdgeAdded(-1), _NextAdjVar(NULL) { }
} ;

class Node
{
public :
	int _Degree ;
	double _LogK ; // log of the domain size of the variable
	int _MinFillScore ;
	double _EliminationScore ; // log of the product of domain sizes of variables in N_G[v]
	AdjVar *_Neighbors ;
public :
	Node(int Degree = 0) : _Degree(Degree), _LogK(0.0), _MinFillScore(-1), _EliminationScore(0.0), _Neighbors(NULL) { }
} ;

class Graph
{
public :
	ARP *_Problem ;
public :
	int _nNodes ;
	Node *_Nodes ; // we assume that the adj-node-list for each node is sorted in increasing order of node indeces
public :
	int _nEdges ; // this is the actual number of edges in the graph; size of _StaticAdjVarTotalList is 2*_nEdges.
	AdjVar *_StaticAdjVarTotalList ;
public :
	// this function is used when the initial edge set has changed (e.g. because some variables have been eliminated), 
	// and we want to consolidate/cleanup the edge set.
	int ReAllocateEdges(void) ;
	// this function is used when the initial edge set has changed and number of edges has changed, 
	// and we want to consolidate/cleanup the edge set.
	// note that NewNumEdges may be larger than current number of edges.
	int ReAllocateEdges(int NewNumEdges) ;
	// this function adds an edge (both ways) between two variables, using the provided AdjVar objects; if uv/vu data member _V is -1, then this obj was not used (probably because already adjacent before).
	int AddEdge(int u, int v, AdjVar & uv, AdjVar & vu) ;
	// this function adds an edge u->v, using the provided AdjVar object
	int AddEdge(int u, int v, AdjVar & uv) ;
	// this function removes an edge between two variables, returning the AdjVar objects that were used for the edges
	int RemoveEdge(int u, int v, AdjVar * & uv, AdjVar * & vu) ;
public :
	int _VarElimOrderWidth ;
	double _MaxVarElimComplexity ; // log og
	double _TotalVarElimComplexity ; // log of
	double _TotalNewFunctionStorageAsNumOfElements ; // log of
	int _nFillEdges ;
public :
	inline double ComputeEliminationComplexity(int v) 
	{
		double score = _Nodes[v]._LogK ;
		for (AdjVar *av = _Nodes[v]._Neighbors ; NULL != av ; av = av->_NextAdjVar) 
			score += _Nodes[av->_V]._LogK ;
		return score ;
/*
		__int64 score = _Problem->K(v) ;
		for (AdjVar *av = _Nodes[v]._Neighbors ; NULL != av ; av = av->_NextAdjVar) {
			score *= _Problem->K(av->_V) ;
			if (score >= InfiniteSingleVarElimComplexity) 
				{ score = InfiniteSingleVarElimComplexity ; return score ; }
			}
		return score ;
*/
	}
	inline void AdjustScoresForArcAddition(int u, int v, int CurrentIteration)
	{
		/* TODO improvement :
			we will move nodes whose score is decremented to 0, into _MFSchangelist[] list.
			scores whose score is incremented are in general list and stay there, so there is no need to add them to the avlVars2CheckScore list.
		*/

		// do a comparative scan of neighbor lists of u and v
		AdjVar *av_u = _Nodes[u]._Neighbors ;
		AdjVar *av_v = _Nodes[v]._Neighbors ;
move_on2 :
		if (NULL == av_v) { // all nodes left in av_u are adjacent to u but not to v; add 1 (to u) for each neighbor.
			for (; NULL != av_u ; av_u = av_u->_NextAdjVar) {
				if (v != av_u->_V) {
					_Nodes[u]._MinFillScore++ ;
					if (0 == _MFShaschanged[u]) 
						{ _MFShaschanged[u] = 1 ; _MFSchangelist[_nMFSchanges++] = u ; }
					}
				}
			return ;
			}
		if (u == av_v->_V) {
			av_v = av_v->_NextAdjVar ;
			goto move_on2 ;
			}
		if (NULL == av_u) { // all nodes left in av_v are adjacent to v but not to u; add 1 (to u) for each neighbor.
			for (; NULL != av_v ; av_v = av_v->_NextAdjVar) {
				if (u != av_v->_V) {
					_Nodes[v]._MinFillScore++ ;
					if (0 == _MFShaschanged[v]) 
						{ _MFShaschanged[v] = 1 ; _MFSchangelist[_nMFSchanges++] = v ; }
					}
				}
			return ;
			}
		if (v == av_u->_V) {
			av_u = av_u->_NextAdjVar ;
			goto move_on2 ;
			}
		if (av_u->_V < av_v->_V) { // av_u->_V is adjacent to u but not to v.
			if (v != av_u->_V) {
				_Nodes[u]._MinFillScore++ ;
				if (0 == _MFShaschanged[u]) 
					{ _MFShaschanged[u] = 1 ; _MFSchangelist[_nMFSchanges++] = u ; }
				}
			av_u = av_u->_NextAdjVar ;
			goto move_on2 ;
			}
		if (av_u->_V > av_v->_V) { // av_v->_V is adjacent to v but not to u.
			if (u != av_v->_V) {
				_Nodes[v]._MinFillScore++ ;
				if (0 == _MFShaschanged[v]) 
					{ _MFShaschanged[v] = 1 ; _MFSchangelist[_nMFSchanges++] = v ; }
				}
			av_v = av_v->_NextAdjVar ;
			goto move_on2 ;
			}
		// now it must be that av_u->_V == av_v->_V; if this common-adjacent-node existing before this iteration, subtract 1.
		// Note that here w = av_u->_V == av_v->_V may on adjacent to X ("on the circle") or not ("outside the circle").
		// in case that is it on the circle, we need to make sure that edges (u,w) and (v,w) existed before.
		if (u < v && av_u->_IterationEdgeAdded < CurrentIteration && av_v->_IterationEdgeAdded < CurrentIteration) {
			int w = av_u->_V ;
			--_Nodes[w]._MinFillScore ; 
			if (0 == _MFShaschanged[w]) 
				{ _MFShaschanged[w] = 1 ; _MFSchangelist[_nMFSchanges++] = w ; }
			}
		av_u = av_u->_NextAdjVar ;
		av_v = av_v->_NextAdjVar ;
		goto move_on2 ;
	}
public :
	// 0=ordered, 1=Trivial, 2=MinFillScore0, 3=General
	char *_VarType ;
	int *_PosOfVarInList ;
	// some variables should be ignored; they should go at the end of the elimination list
	int _nIgnoreVariables ;
	int _IgnoreVariables[1] ; // 2014-04-21 KK : for now, we have at most 1 ignore variable.
	// list of ordered variables; in elimination order.
	int _OrderLength ;
	int *_VarElimOrder ;
	// temporary space for holding trivial nodes (degree(X) <= 1)
	int _nTrivialNodes ;
	int *_TrivialNodesList ;
	// temporary holding space for MinFillScore=0 nodes. During OrderingCreation, we rely on the fact that once a node's MinFill score becomes 0, it stays zero.
	// note that if node is trivial, it should not be in MinFillScore=0 list.
	int _nMinFillScore0Nodes ;
	int *_MinFill0ScoreList ;
	// temporary holding space for remining nodes.
	int _nRemainingNodes ;
	int *_RemainingNodesList ;
public :
	// temporary space for storing edges to add during each variable ordering computation iteration
	short _EdgeU[ARE_GRAPH_VAR_ORDER_COMP_MAX_NUM_EDGES_ADDED] ;
	short _EdgeV[ARE_GRAPH_VAR_ORDER_COMP_MAX_NUM_EDGES_ADDED] ;
	// we need to track vars whose MFS changes.
	char *_MFShaschanged ; // a boolean for each var, whether its MFS has changed of not
	int *_MFSchangelist ; // a list of vars whose MFS has changed
	int _nMFSchanges ;
public :
	inline bool IsIgnoreVariable(int X)
	{
		for (int iii = 0 ; iii < _nIgnoreVariables ; iii++) {
			if (X == _IgnoreVariables[iii]) 
				return true ;
			}
		return false ;
	}
	inline void RemoveVarFromList(int X)
	{
		int i = _PosOfVarInList[X] ;
		// it must be that i>=0
		_PosOfVarInList[X] = -1 ;
		if (3 == _VarType[X]) {
			if (--_nRemainingNodes != i) {
				_RemainingNodesList[i] = _RemainingNodesList[_nRemainingNodes] ;
				_PosOfVarInList[_RemainingNodesList[i]] = i ;
				}
			}
		else if (2 == _VarType[X]) {
			if (--_nMinFillScore0Nodes != i) {
				_MinFill0ScoreList[i] = _MinFill0ScoreList[_nMinFillScore0Nodes] ;
				_PosOfVarInList[_MinFill0ScoreList[i]] = i ;
				}
			}
		else if (1 == _VarType[X]) {
			if (--_nTrivialNodes != i) {
				_TrivialNodesList[i] = _TrivialNodesList[_nTrivialNodes] ;
				_PosOfVarInList[_TrivialNodesList[i]] = i ;
				}
			}
		// else X is ordered; we will not remove these variables.
	}
	inline void ProcessPostEliminationNodeListLocation(int u) 
	{
		// check if u has become special; this function is called, after X is eliminated, for all nodes u whose MinFill score has changed

		if (2 == _VarType[u]) { // u is in MilFillScore0 list; check if it can be moved up to trivial list
			// check if u has become trivial
			if (_Nodes[u]._Degree <= 1) {
// DEBUGGG
//printf("\nMinFill : round %d picked var=%d; making u=%d trivial, VarToMinFill0ScoreListMap[u]=%d ...", (int) nOrdered, (int) X, (int) u, (int) _VarToMinFill0ScoreListMap[u]) ;
				RemoveVarFromList(u) ;
				_PosOfVarInList[u] = _nTrivialNodes ;
				_TrivialNodesList[_nTrivialNodes++] = u ;
				_VarType[u] = 1 ;
				}
			}
		else if (3 == _VarType[u]) { // u is in general list; check if it can be moved up to MinFillScore=0 or trivial list
			// check if u has become trivial
			if (_Nodes[u]._Degree <= 1) {
				RemoveVarFromList(u) ;
// DEBUGGG
//printf("\nMinFill : round %d picked var=%d; making u=%d trivial, VarToMinFill0ScoreListMap[u]=%d ...", (int) nOrdered, (int) X, (int) u, (int) _VarToMinFill0ScoreListMap[u]) ;
				_PosOfVarInList[u] = _nTrivialNodes ;
				_TrivialNodesList[_nTrivialNodes++] = u ;
				_VarType[u] = 1 ;
				}
			// check if u has become MinFillScore=0
			else if (_Nodes[u]._MinFillScore <= 0) {
				RemoveVarFromList(u) ;
				_PosOfVarInList[u] = _nMinFillScore0Nodes ;
				_MinFill0ScoreList[_nMinFillScore0Nodes++] = u ;
				_VarType[u] = 2 ;
				}
			}
	}
public :
	bool _IsValid ;
protected :
	MTRand _RNG ;
public :
	int ComputeVariableEliminationOrder_Simple(
		char CostFunction, // 0=MinFill, 1=MinDegree, 2=MinComplexity
		// width/complexity of the best know order; used to cut off search when the elimination order we found is not very good
		int WidthLimit, 
		bool EarlyTermination_W, 
		double TotalComplexityLimit, // log of the complexity limit
		bool EarlyTermination_C, 
		// if true, we will quit after all Trivial/MinFillScore0 nodes have been used up. This is typically used to generate a starting point for large-scale randomized searches.
		bool QuitAfterEasyIsDone, 
		// width that is considered trivial; whenever a variale with this width is found, it can be eliminated right away, even if it is not the variable with the smallest width
		int EasyWidth, 
		// when no easy variables (e.g. degree <=1, MinFillScore=0) are to pick, we pick randomly among a few best nodes.
		int nRandomPick, 
		double eRandomPick, 
		// temp AdjVar space; size of each block is TempAdjVarSpaceSize.
		int & TempAdjVarSpaceSizeExtraArrayN, AdjVar *TempAdjVarSpaceSizeExtraArray[]
		) ;
	int ComputeVariableEliminationOrder_Simple_wMinFillOnly(
		// width/complexity of the best know order; used to cut off search when the elimination order we found is not very good
		int WidthLimit, 
		bool EarlyTermination_W, 
		// if true, we will quit after all Trivial/MinFillScore0 nodes have been used up. This is typically used to generate a starting point for large-scale randomized searches.
		bool QuitAfterEasyIsDone, 
		// width that is considered trivial; whenever a variale with this width is found, it can be eliminated right away, even if it is not the variable with the smallest width
		int EasyWidth, 
		// when no easy variables (e.g. degree <=1, MinFillScore=0) are to pick, we pick randomly among a few best nodes.
		int nRandomPick, 
		double eRandomPick, 
		// temp AdjVar space; size of each block is TempAdjVarSpaceSize.
		int & TempAdjVarSpaceSizeExtraArrayN, AdjVar *TempAdjVarSpaceSizeExtraArray[]
		) ;
public :
	int RemoveRedundantFillEdges(void) ;
public :
	int operator=(const Graph & G) ;
	int Test(int MaxWidthAcceptableForSingleVariableElimination) ;
public :
	Graph(ARP *Problem = NULL) ;
	~Graph(void) ;
	int Destroy(void) ;

	// create a graph from a problem
	int Create(ARP & Problem) ;

//	// create a graph from the given set of fn signatures. 
//	// we assume each signature is correct : var indeces range [0, nNodes), there are no repetitions.
//	int Create(int nNodes, std::vector<std::set<int>> & fn_signatures) ;
} ;

} // namespace ARE

#endif // ARE_Graph_HXX_INCLUDED
