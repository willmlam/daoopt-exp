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

#ifndef MiniBucketElimLH_H_
#define MiniBucketElimLH_H_

#include "MiniBucketElim.h"
#include "MiniBucketElimLHsubtree.h"

#undef DEBUG

#define OUR_OWN_nInfinity (-std::numeric_limits<double>::infinity())

//#define USE_FULL_LOOKAHEAD_SUBTREE

namespace daoopt {

class MiniBucketElimLH;

// old/obsolete LH code
class MiniBucketElimLHErrorNode;
class MiniBucketElimLHError;
class MiniBucketElimLHheuristicNode;
class MiniBucketElimLHheuristic;

// new/current LH code
class MBLHSubtreeNode ;
class MBLHSubtree ;

class MiniBucketElimLHStatistics
{
public :
	int _iBound ;
	int _Width ; // num of vars when only original functions placed in the bucket are considered
	int _PseudoWidth ; // num of vars when all functions (original and MB generated) placed in the bucket are considered (this is the scope of bucket fn).
	int _EnforcedPseudoWidth; // "pseudowidth" based on the local error functions that were actually computed.
//	int _HPseudoWidth ; // width when only heuristic function h is done
	size_t _MemorySize ;
	double _LEMemorySizeMB;
	int _MaxNumMBs ;
	int _NumBucketsWithMoreThan1MB ;
	vector<size_t> _NumNodesLookahead; // actual number of lookahead nodes during search by variable.
	double _LookaheadTotalTime;
public :
	void reset(void)
	{
		_iBound = 0 ;
		_Width = 0 ;
		_PseudoWidth = 0 ;
		_EnforcedPseudoWidth = 0;
//		_HPseudoWidth = 0 ;
		_MemorySize = 0 ;
		_LEMemorySizeMB = 0 ;
		_MaxNumMBs = 0 ;
		_NumBucketsWithMoreThan1MB = 0 ;
	}
public :
	MiniBucketElimLHStatistics(void) :
		_iBound(-1), 
		_Width(-1), 
		_PseudoWidth(-1), 
		_EnforcedPseudoWidth(-1),
//		_HPseudoWidth(-1), 
		_MemorySize(0), 
		_LEMemorySizeMB(0),
		_MaxNumMBs(0), 
		_NumBucketsWithMoreThan1MB(0)
	{
	}
} ;

/* enhanced minibucket elimination */
class MiniBucketElimLH : public MiniBucketElim
{
	friend class MiniBucketElimLHError ;
	friend class MiniBucketElimLHErrorNode ;
	friend class MiniBucketElimLHheuristic ;
	friend class MiniBucketElimLHheuristicNode ;
	friend class MBLHSubtree ;
	friend class MBLHSubtreeNode ;

 protected:

	MiniBucketElimLHStatistics _Stats ;

	// in the bucket tree, distance of the variable from the root, as the number of edges needed to travel to reach the bucket. 0 for the root.
	std::vector<int> _Depth ;
	int _MaxDepth ;

	// joint scope of each bucket, including the bucket's var
	std::vector<std::set<int>> _BucketScopes ;
	// original + augmented function of each bucket; this is the set of functions that minibucket algorithm works with.
	std::vector<std::vector<Function *>> _BucketFunctions ;

	std::vector<Function*> _BucketErrorFunctions ;
	double _BuckerErrorFnTableSizes_Total ; // log
	double _BuckerErrorFnTableSizes_Precomputed ; // log
	double _BuckerErrorFnTableSizes_Ignored ; // log
	int64 _nBucketsWithNonZeroBuckerError ;
	int64 _nBucketsWithMoreThan1MB ;

	// for each variable, its bucket error indicator :
	//	   -99  = marked by user not to be included in the LH subtree
	//		-1  = unknown
	//		 0  = bucket error is 0
	//		 1  = bucket error is >0 (but perhaps trivial)
	//		 2+ = bucket error is >0 and non-trivial
	//		99  = marked by default (by user) to be included in the LH subtree
	// note that bucket error is MBvalue - B(exact)value. given that we solve maximization problem, this should be always >=0, since MBvalue>=Bvalue.
	std::vector<signed char> _BucketErrorQuality ;

	// number of times a LH call was made on each variable
	std::vector<int64> _nLHcalls ;

	// avg abs of bucket error. when initialized this is -infinity (unknown). actual value should be >=0.
	// avg does not include those error entries that are infinite.
	std::vector<double> _BucketError_AbsAvg, _BucketError_AbsMin, _BucketError_AbsMax ;
  std::vector<double> _BucketError_Rel;

  // Used for "static" subproblem ordering heuristic
  std::vector<double> _SubtreeError;

  std::vector<int> _Pseudowidth;

	// for each var, distance to the closest descendant with more than 1 MB/BucketErrorQuality>0. INT_MAX means infinite (no descendants).
	// this does not include the node itself; i.e. legal values are >0 & <= INT_MAX (meaning infinity = no such descendants).
	// note that _distToClosestDescendantWithMBs[i] <= _distToClosestDescendantWithLE[i], since the only way a bucket error can occur 
	// is if nMBs>1 (the other way is not neccessarily correct - it could be that nMBs>1, but LE=0).
	std::vector<int> _distToClosestDescendantWithMBs ;
	std::vector<int> _distToClosestDescendantWithLE ;

//	std::vector<MiniBucketElimLHobsoleteError> _LookaheadResiduals ;
//	std::vector<MiniBucketElimLHobsoleteheuristic> _Lookahead ;
	std::vector<MBLHSubtree> _Lookahead ;

public:

	inline std::vector<MBLHSubtree> & LH(void) { return _Lookahead ; }
	inline int64 nLHcalls(int v) { return _nLHcalls[v] ; }
	inline std::set<int> & BucketScope(int v) { return _BucketScopes[v] ; }
	inline std::vector<signed char> & BucketErrorQuality(void) { return _BucketErrorQuality ; }
	inline signed char BucketErrorQuality(int v) { return _BucketErrorQuality[v] ; }
	inline int distToClosestDescendantWithLE(int v) { return _distToClosestDescendantWithLE[v] ; }
	inline double BucketError_AbsAvg(int v) { return _BucketError_AbsAvg[v] ; }
	inline double BucketError_AbsMax(int v) { return _BucketError_AbsMax[v] ; }

	// builds the heuristic, limited to the relevant subproblem, if applicable.
	// if computeTables=false, only returns size estimate (no tables computed)
	virtual size_t build(const std::vector<val_t>* assignment = NULL, bool computeTables = true) ;

	// compute the _distToClosestDescendantWithLE[] array, based on _BucketErrorQuality[] array.
	int Compute_distToClosestDescendantWithLE(void) ;

	virtual void noteOrNodeExpansionBeginning(int Var, std::vector<val_t> & Context, SearchNode *search_node) ;

	// reset the data structures
	virtual void reset() ;

	// computes the heuristic for variable var given a (partial) assignment
	virtual double getHeur(int var, std::vector<val_t>& assignment, SearchNode* n) ;
	virtual double getHeurPerIndSubproblem(int var, std::vector<val_t> & assignment, SearchNode* node, double label, std::vector<double> & subprobH);
	// computes heuristic values for all instantiations of var, given context assignment
	virtual void getHeurAll(int var, std::vector<val_t>& assignment, SearchNode* n, std::vector<double>& out);

	virtual double getOrderingHeur(int var, std::vector<val_t>& assignment,
                                 SearchNode* n) override;

	// computes the lookahead heuristic for variable var given a (partial) assignment.
	// note : this fn returns BucketError[var]. all variables in the output fn of bucket[var] should be instantiated by 'assignment'.
	// note : the return value is in the same representation space (log or normal) as the problem.
	// note : the return values should be >= 0.0, regardless of whether normal/log representation.
	virtual double getLocalError(int var, vector<val_t> & assignment) ;
	// get local error of the 'Child', for each value of the 'Parent'.
	// we assume that 'out' is initialized; furthermore, we will subtract 
	// localerror from current values of 'out'.
	// e.g. typically 'out' is passed in as h/f of 'parent'.
	virtual void getLocalError(int parent, int child, std::vector<val_t>& assignment, std::vector<double>& out) ;

	// compute local MiniBucket error, defined as the different between 
	//		1) the combination of all the mini-bucket output functions,
	//		2) what the normal/full bucket elimination would compute. 
	// note : error_fn/avg_error/avg_exact are in the same representation space (normal/log) as the problem itself.
	// note : localerror should always be >=0, regardless of whether normal/log representation, since item_1 is at least as large as item_2 (above).
	// note : the scope of the table (as fn) is the same as the bucket output fn.
	// note : this fn sets _BucketErrorQuality[var] to 0/1+ when it can be determined; sometimes when bucket error table is all 0, we won't compute the table and will set _BucketErrorQuality[var]=0.
	int computeLocalErrorTable(int var, bool build_table, bool sample_table_if_not_computed, 
		double TableMemoryLimitAsNumElementsLog, 
		double & TableSizeLog, // OUT : bucket error table size, regardless of whether it is actually computed; this is in log scale, i.e. sum_log10(var_domain_size).
		double & avgError, // OUT : avg bucket error; computed when output table is not too large.
		double & avgExact, // OUT : avg value over entire bucket output table; computed when error is computed.
		Function * & error_fn, // OUT : bucket error fn; computed when output table is not too large.
		int64 & nEntriesGenerated) ; // OUT : number of actual table entires computed
	int computeLocalErrorTables(
		bool build_tables, 
		double TotalMemoryLimitAsNumElementsLog, 
		double TableMemoryLimitAsNumElementsLog) ; 

	// compute local MiniBucket error_H, defined as the different between 
	// 1) max-product of m_augmented FNs in the bucket
	// 2) output FNs of all minibuckets of the bucket
	int localError_H(int var, bool computeTable, double & avgError, Function * & errorFn) ;
	int localError_H(double & avgError, bool do_printf) ;

	int deleteLocalErrorFNs(void)
	{
		for (vector<Function*>::iterator it = _BucketErrorFunctions.begin(); it!=_BucketErrorFunctions.end(); ++it)
			if (NULL != *it) delete *it ;
		_BucketErrorFunctions.clear() ;
		return 0 ;
	}

  void InitializeSubtreeErrors(const std::vector<double>& bucket_error);

public :

	MiniBucketElimLH(Problem *p, Pseudotree *pt, ProgramOptions *po, int ib) ;
  void printExtraStats() const;
	virtual ~MiniBucketElimLH(void) ;
} ;

inline MiniBucketElimLH::MiniBucketElimLH(Problem* p, Pseudotree* pt, ProgramOptions* po, int ib) :
    MiniBucketElim(p, pt, po, ib), 
	_BuckerErrorFnTableSizes_Total(-1), 
	_BuckerErrorFnTableSizes_Precomputed(-1), 
	_BuckerErrorFnTableSizes_Ignored(-1), 
	_nBucketsWithNonZeroBuckerError(-1), 
	_nBucketsWithMoreThan1MB(-1) 
{
}

inline void MiniBucketElimLH::printExtraStats() const {
  cout << "Heuristic call counts: " 
       << endl;
  uint64 total_calls_or = 0;
  uint64 total_calls_and = 0;
  for (int i = 0; i < m_problem->getN(); ++i) {
    total_calls_and += var_heur_calls_[i];
    uint64 adjusted_count =
        var_heur_calls_[i] / m_problem->getDomainSize(i);
    cout << " " << adjusted_count;
    total_calls_or += adjusted_count;
  }
  cout << endl;

  cout << "Lookahead node counts: " << endl;

  uint64 total_lookahead_or = 0;
  uint64 total_lookahead_and = 0;
  uint64 count_var_lookahead = 0;
  for (int i = 0; i < m_problem->getN(); ++i) {
    total_lookahead_and += _Stats._NumNodesLookahead[i];
    uint64 adjusted_count =
      _Stats._NumNodesLookahead[i] / m_problem->getDomainSize(i);
    cout << " " << adjusted_count;
    total_lookahead_or += adjusted_count;
    if (adjusted_count > 0) count_var_lookahead += 1;
  }
  cout << endl;
  cout << "Heuristic calls (OR): " << total_calls_or << endl;
  cout << "Heuristic calls (AND): " << total_calls_and << endl;
  cout << "Lookahead calls (OR): " << total_lookahead_or << endl;
  cout << "Lookahead calls (AND): " << total_lookahead_and << endl;
  cout << "Lookahead ratio (OR): " 
       << double(total_lookahead_or) / total_calls_or << endl;
  cout << "Lookahead ratio (AND): " 
       << double(total_lookahead_and) / total_calls_and << endl;
  cout << "Variables w/ lookahead: " << count_var_lookahead << endl;
  cout << "Lookahead total time: " << _Stats._LookaheadTotalTime << endl;
}

inline MiniBucketElimLH::~MiniBucketElimLH(void)
{
  this->reset();
}

/*
int computeMBOutFnArgsVectorPtrMap(
	// IN
	int elim_var, 
	vector<Function*> & Functions, 
	// OUT
	set<int> & scope, 
	int & n, 
	val_t * & tuple, // size is n+1; elim_var is last.
	vector<vector<val_t *>> & idxMap // this array can be used in a call to Function::getValuePtr().
	) ;
*/

class MiniBucketElimLHErrorNode
{
public :
	MiniBucketElimLH *_H ;
public :
	int _v ; // variable of this node.
	int _idxWRTparent ; // idx of this node as child of its parent; this is wrt the childlen list of the parent.
	int _k ; // domain size of this variable.
	int _depth2go ; // depth still to go from this node; as number of edges.
	std::vector<Function *> _RelevantMiniBucketFunctions ; // IF/AF that came from within of the LH-subteee rooted at this node
	std::vector<Function *> _RelevantBucketFunctions ; // OF + AF that came from outside of the LH-subteee rooted at this node
	std::vector<MiniBucketElimLHErrorNode *> _Children ; // children of this node in the bucket tree, that have non-0 bucket error within LH-subteee rooted at this node.
public :
	int Delete(void) 
	{
		for (vector<MiniBucketElimLHErrorNode *>::iterator it = _Children.begin(); it!=_Children.end(); ++it) 
			delete *it ;
		_Children.clear() ;
		_RelevantBucketFunctions.clear() ;
		return 0 ;
	}
	inline double ExactBucketValue(vector<val_t> & assignment) const
	{
/*		int i ;
		vector<Function*>::const_iterator itF_end = _RelevantBucketFunctions.end() ;
		// sum up bucket functions
		vector<double> funVals ;
		vector<double> sumVals ; sumVals.resize((vector<double>::size_type) _k, ELEM_ONE) ;
		for (vector<Function*>::const_iterator itF = _RelevantBucketFunctions.begin() ; itF != itF_end ; ++itF) {
			(*itF)->getValues(assignment, _v, funVals) ;
			for (i = 0 ; i < _k ; ++i) 
				sumVals[i] OP_TIMESEQ funVals[i] ;
			}
		if (_Children.size() > 0) {
			vector<MiniBucketElimLHErrorNode *>::const_iterator itC_end = _Children.end() ;
			for (i = 0 ; i < _k ; ++i) {
				assignment[_v] = i ;
				for (vector<MiniBucketElimLHErrorNode *>::const_iterator itC = _Children.begin() ; itC != itC_end ; ++itC) 
					sumVals[i] OP_TIMESEQ (*itC)->ExactBucketValue(assignment) ;
				}}
		double tableentryB = sumVals[0] ;
		for (i = 1 ; i < _k ; ++i) 
			tableentryB = max(tableentryB, sumVals[i]) ;
		return tableentryB ;*/
		vector<Function*>::const_iterator itF_end = _RelevantBucketFunctions.end() ;
		vector<MiniBucketElimLHErrorNode *>::const_iterator itC_end = _Children.end() ;
		double valueMax = ELEM_ZERO ;
#ifdef DEBUG
    bool assign_is_none = false;
#endif
		for (int i = 0 ; i < _k ; ++i) {
#ifdef DEBUG
      if (assignment[_v] == NONE) {
        assign_is_none = true;
      }
#endif
			assignment[_v] = i ;
			double value = ELEM_ONE ;
			for (vector<Function*>::const_iterator itF = _RelevantBucketFunctions.begin() ; itF != itF_end ; ++itF) 
				value OP_TIMESEQ (*itF)->getValue(assignment) ;
			for (vector<MiniBucketElimLHErrorNode *>::const_iterator itC = _Children.begin() ; itC != itC_end ; ++itC) 
				value OP_TIMESEQ (*itC)->ExactBucketValue(assignment) ;
			valueMax = max(valueMax, value) ;
			}
#ifdef DEBUG
    if (assign_is_none) {
      assignment[_v] = NONE;
    }
#endif
		return valueMax ;
	}
	inline double Error(vector<val_t> & assignment) const
	{
		if (1 == _depth2go) {
			const Function *eFN = (_H->_BucketErrorFunctions)[_v] ;
			if (NULL != eFN) 
				return eFN->getValue(assignment) ;
			}
		double tableentryMB = ELEM_ONE ;
		vector<Function*>::const_iterator itF_end = _RelevantMiniBucketFunctions.end() ;
		for (vector<Function*>::const_iterator itF = _RelevantMiniBucketFunctions.begin() ; itF != itF_end ; ++itF) 
			tableentryMB OP_TIMESEQ (*itF)->getValue(assignment) ;
		if (OUR_OWN_nInfinity == tableentryMB) 
			return 0.0 ; // if must be that tableentryB is also -Infinity and error is 0.
		double tableentryB = ExactBucketValue(assignment) ;
		if (tableentryMB <= tableentryB) return 0.0 ;
		return tableentryMB - tableentryB ;
	}
	// get heuristic error for each value of the variable '_v'.
	// we assume that 'out' is initialized; furthermore, we will subtract heuristic error(s) we compute from current values of 'out'.
	// e.g. typically 'out' is passed in as h/f of 'Parent'.
	inline void Error(int root_var, vector<val_t> & assignment, vector<double> & out) const
	{
		int k, domainsize_root_var = _H->m_problem->getDomainSize(root_var) ;
		vector<double> funVals ;
		if (1 == _depth2go) {
			const Function *eFN = (_H->_BucketErrorFunctions)[_v] ;
			if (NULL != eFN) {
				eFN->getValues(assignment, root_var, funVals) ;
				for (k = 0 ; k < domainsize_root_var ; ++k) 
					out[k] OP_DIVIDEEQ funVals[k] ;
				return ;
				}
			}
		vector<double> MBVals ; MBVals.resize((vector<double>::size_type) domainsize_root_var, ELEM_ONE) ;
		vector<Function*>::const_iterator itF_end = _RelevantMiniBucketFunctions.end() ;
		for (vector<Function*>::const_iterator itF = _RelevantMiniBucketFunctions.begin() ; itF != itF_end ; ++itF) {
			(*itF)->getValues(assignment, root_var, funVals) ;
			for (k = 0 ; k < domainsize_root_var ; ++k) 
				MBVals[k] OP_TIMESEQ funVals[k] ;
			}
		// compute bucket output for the given assignment; this is the exact bucket value, given its functions and assignment.
		// note here we have two unassigned variables : _v and Child (var). we have to enumerate over one of them.
		for (k = 0 ; k < domainsize_root_var ; ++k) {
			if (OUR_OWN_nInfinity == MBVals[k]) 
				continue ; // if must be that tableentryB is also -Infinity and error is 0.
			assignment[root_var] = k ;
			double tableentryB = ExactBucketValue(assignment) ;
			if (MBVals[k] <= tableentryB) 
				continue ; // '<' is bad, '=' is ok. either way, error is 0.
			out[k] OP_DIVIDEEQ (MBVals[k] - tableentryB) ;
			}
	}
	MiniBucketElimLHErrorNode(void) : _H(NULL), _v(-1), _idxWRTparent(-1), _k(-1), _depth2go(-1)
	{
	}
	~MiniBucketElimLHErrorNode(void)
	{
		Delete() ;
	}
} ;

class MiniBucketElimLHError
{
public :
	MiniBucketElimLH *_H ;
public :
	int _v ; // variable of this node.
	int _depth ;
	MiniBucketElimLHErrorNode _RootNode ;
	std::vector<MiniBucketElimLHErrorNode *> _SubtreeNodes ;
public :
	bool IsSubtreeVariable(int v)
	{
		for (vector<MiniBucketElimLHErrorNode *>::iterator itRB = _SubtreeNodes.begin(); itRB!=_SubtreeNodes.end(); ++itRB) {
			if (v == (*itRB)->_v) 
				return true ;
			}
		return false ;
	}
	bool IsSubtreeMBEfunction(Function & f)
	{
		for (vector<MiniBucketElimLHErrorNode *>::iterator itRB = _SubtreeNodes.begin(); itRB!=_SubtreeNodes.end(); ++itRB) {
			MiniBucketElimLHErrorNode *n = *itRB ;
			int u = n->_v ;
			vector<MiniBucket> & minibuckets = _H->_MiniBuckets[u] ;
			for (vector<MiniBucket>::iterator itMB = minibuckets.begin() ; itMB != minibuckets.end() ; ++itMB) {
				Function *fnMB = itMB->output_fn() ;
				if (&f == fnMB) return true ;
				}
			}
		return false ;
	}
	inline double Error(vector<val_t> & assignment) const
	{
		if (0 == _SubtreeNodes.size()) 
			return 0.0 ; // most likely, all the children(buckets) up to depth d have exactly 1 MB.
		// need to sum errors for each child separately, because MB value for some child could be -infinity and hence if we summed all MB, it would also be -infinity, 
		// whereas in reality, some children have 0-error, some non-0-error.
		double e = ELEM_ONE ;
		vector<MiniBucketElimLHErrorNode *>::const_iterator itC_end = _RootNode._Children.end() ;
		for (vector<MiniBucketElimLHErrorNode *>::const_iterator itC = _RootNode._Children.begin() ; itC != itC_end ; ++itC) 
			e OP_TIMESEQ (*itC)->Error(assignment) ;
		return e ;
	}
	inline double ErrorPerIndSubproblem(vector<val_t> & assignment, double current_pruning_gap, std::vector<double> & subprobH) const
	{
		if (0 == _SubtreeNodes.size()) 
			return 0.0 ; // most likely, all the children(buckets) up to depth d have exactly 1 MB.
		// need to sum errors for each child separately, because MB value for some child could be -infinity and hence if we summed all MB, it would also be -infinity, 
		// whereas in reality, some children have 0-error, some non-0-error.
		double e = ELEM_ONE ;
		vector<MiniBucketElimLHErrorNode *>::const_iterator itC_end = _RootNode._Children.end() ;
		for (vector<MiniBucketElimLHErrorNode *>::const_iterator itC = _RootNode._Children.begin() ; itC != itC_end ; ++itC) {
			double eChild = (*itC)->Error(assignment) ; // eChild should be >=0.
			subprobH[(*itC)->_idxWRTparent] OP_DIVIDEEQ eChild ;
			e OP_TIMESEQ eChild ;
			current_pruning_gap OP_DIVIDEEQ eChild ;
			if (current_pruning_gap <= 0.0) 
				break ; // error so far is big enough to cause pruning; no need to improve the error further.
			}
		return e ;
	}
	// get heuristic error for each value of the variable '_v'.
	// we assume that 'out' is initialized; furthermore, we will subtract heuristic error(s) we compute from current values of 'out'.
	// e.g. typically 'out' is passed in as h/f of 'Parent'.
	inline void Error(vector<val_t> & assignment, vector<double> & out) const
	{
		if (0 == _SubtreeNodes.size()) 
			return ; // most likely, all the children(buckets) up to depth d have exactly 1 MB.
#if defined DEBUG || _DEBUG
		vector<double> out_(out) ;
#endif
		// need to sum errors for each child separately, because MB value for some child could be -infinity and hence if we summed all MB, it would also be -infinity, 
		// whereas in reality, some children have 0-error, some non-0-error.
		vector<MiniBucketElimLHErrorNode *>::const_iterator itC_end = _RootNode._Children.end() ;
		for (vector<MiniBucketElimLHErrorNode *>::const_iterator itC = _RootNode._Children.begin() ; itC != itC_end ; ++itC) 
			(*itC)->Error(_v, assignment, out) ;
#if defined DEBUG || _DEBUG
		for (int k = 0 ; k < _H->m_problem->getDomainSize(_v) ; k++) {
			assignment[_v] = k ;
			double e = Error(assignment) ;
			double e_ = out_[k] - e ;
			if (fabs(e_ - out[k]) > 1.0e-32) {
				int bug = 1 ;
				}
			}
#endif
	}
public :
	int Initialize(MiniBucketElimLH & H, int v, int depth)
	{
		Delete() ;
		_H = &H ; _v = v ; _depth = depth ;
		_RootNode._v = v ; _RootNode._k = (H.m_problem)->getDomainSize(v) ; _RootNode._depth2go = depth ;
		// always check if there is any point in computing error
		if (H._distToClosestDescendantWithLE[v] > depth) 
			return 0 ;
		// build tree of nodes
		int memory = depth*5 ; if (memory > H.m_problem->getN()) memory = H.m_problem->getN() ; // assume avg branching factor of 5
		_SubtreeNodes.reserve(memory) ; // try to allocate some memory to make it faster
//		std::deque<MiniBucketElimLHErrorNode *> nodes2process ;
		std::stack<MiniBucketElimLHErrorNode *> nodes2process ; // using stack will turn it into a DFS processing of the tree
		try { nodes2process.push(&_RootNode) ; } catch (...) { return 1 ; }
		int rootchilddepth = depth-1 ;
		MiniBucketElimLHErrorNode *currentrootchild = NULL ;
		while (nodes2process.size() > 0) {
			MiniBucketElimLHErrorNode *N = nodes2process.top() ; nodes2process.pop() ;
			if (rootchilddepth == N->_depth2go) currentrootchild = N ;
			const PseudotreeNode *N_ = H.m_pseudotree->getNode(N->_v) ;
			const vector<PseudotreeNode *> & children = N_->getChildren() ;
			N->_Children.reserve(children.size()) ; // reserver space for all children (some will not be added as children) so that reallocation is not needed
			int depth2go = N->_depth2go - 1 ; int idxChild = children.size() - 1 ;
			for (vector<PseudotreeNode*>::const_reverse_iterator itC = children.rbegin() ; itC != children.rend(); ++itC, --idxChild) {
				int child = (*itC)->getVar() ;
#ifndef USE_FULL_LOOKAHEAD_SUBTREE
				if (H._BucketErrorQuality[child] <= 1 ? H._distToClosestDescendantWithLE[child] > depth2go : false) continue ; // _BucketErrorQuality <= 1 means it is not proven that there is substantial bucket error
#endif // USE_FULL_LOOKAHEAD_SUBTREE
				MiniBucketElimLHErrorNode *n = NULL ; try { n = new MiniBucketElimLHErrorNode ; } catch (...) { return 1 ; } if (NULL == n) return 1 ; n->_H = &H ; n->_v = child ; n->_k = (H.m_problem)->getDomainSize(child) ; n->_depth2go = depth2go ;
				try { _SubtreeNodes.push_back(n) ; } catch (...) { delete n ; return 1 ; }
				try { (N->_Children).push_back(n) ; } catch (...) { delete n ; return 1 ; }
				if (depth2go > 0) { try { nodes2process.push(n) ; } catch (...) { delete n ; return 1 ; }}
				n->_idxWRTparent = N->_Children.size() - 1 ;
				// add all MB output functions of [child] bucket that are in [v] bucket to currentrootchild
				MiniBucketElimLHErrorNode *currentrootchild_ = (NULL == currentrootchild) ? n : currentrootchild ; 
				vector<MiniBucket> & minibuckets = H._MiniBuckets[child] ;
				for (vector<MiniBucket>::iterator itMB = minibuckets.begin() ; itMB != minibuckets.end() ; ++itMB) {
					Function *fnMB = itMB->output_fn() ;
					if (NULL == fnMB) continue ;
					vector<Function *> & vAFlist = H.m_augmented[v] ;
					vector<Function *> & vIFlist = H.m_intermediate[v] ;
					if (vAFlist.end() == std::find(vAFlist.begin(), vAFlist.end(), fnMB) && vIFlist.end() == std::find(vIFlist.begin(), vIFlist.end(), fnMB)) continue ;
					try { (currentrootchild_->_RelevantMiniBucketFunctions).push_back(fnMB) ; } catch (...) { return 1 ; }
					}
				}
			}
		// fill in _RelevantBucketFunctions for each relevant bucket
		for (vector<MiniBucketElimLHErrorNode *>::iterator itRB = _SubtreeNodes.begin(); itRB!=_SubtreeNodes.end(); ++itRB) {
			MiniBucketElimLHErrorNode *n = *itRB ;
			int u = n->_v ;
//			int d = n->_depth2go ;
			std::vector<Function *> & funs = n->_RelevantBucketFunctions ;
			// add OF_n
			const vector<Function *> & OFlist = H.m_pseudotree->getFunctions(u) ;
			funs.insert(funs.end(), OFlist.begin(), OFlist.end()) ;
			// add AF_n^{>d}; i.e. those functions in AF_n that did not come from a subtree bucket
			vector<Function *> & AFlist = H.m_augmented[u] ;
			for (vector<Function*>::iterator itF = AFlist.begin() ; itF != AFlist.end() ; ++itF) {
				Function *fn = *itF ;
//				int vOriginating = abs(fn->getId()) ;
//				if (! IsSubtreeVariable(vOriginating)) 
				if (! IsSubtreeMBEfunction(*fn)) 
					funs.push_back(fn) ;
				}
			}
/*		// fill in a list of relevant functions of the node[v] = all IF/AF in [v] that came from relevant descendant buckets
this code has bug; for all MB output FNs we should check if it is in [v] bucket first
		for (vector<MiniBucketElimLHErrorNode *>::iterator itRB = _SubtreeNodes.begin(); itRB!=_SubtreeNodes.end(); ++itRB) {
			MiniBucketElimLHErrorNode *n = *itRB ;
			int u = n->_v ;
			vector<MiniBucket> & minibuckets = H._MiniBuckets[u] ;
			for (vector<MiniBucket>::iterator itMB = minibuckets.begin() ; itMB != minibuckets.end() ; ++itMB) {
				Function *fnMB = itMB->outputFn() ;
				if (NULL == fnMB) continue ;
				_RootNode._RelevantBucketFunctions.push_back(fnMB) ;
				}
			}*/
		return 0 ;
	}
	int Delete(void) 
	{
		_RootNode.Delete() ;
		_SubtreeNodes.clear() ;
		return 0 ;
	}
	MiniBucketElimLHError(void) : _H(NULL), _v(-1), _depth(-1)
	{
	}
	~MiniBucketElimLHError(void)
	{
	}
} ;

class MiniBucketElimLHheuristicNode
{
public :
	MiniBucketElimLH *_H ;
public :
	int _v ; // variable of this node.
	int _idxWRTparent ; // idx of this node as child of its parent; this is wrt the childlen list of the parent.
	int _k ; // domain size of this variable.
	int _depth2go ; // depth still to go from this node; as number of edges.
	std::vector<Function *> _RelevantBucketFunctions ; // OF + AF (in this bucket) that came from outside of the LH-subteee
	std::vector<MiniBucketElimLHheuristicNode *> _Children ; // children of this node in the bucket tree, that have non-0 bucket error within LH-subteee rooted at this node.
public :
	int Delete(void) 
	{
		for (vector<MiniBucketElimLHheuristicNode *>::iterator it = _Children.begin(); it!=_Children.end(); ++it) 
			delete *it ;
		_Children.clear() ;
		_RelevantBucketFunctions.clear() ;
		return 0 ;
	}
	inline double ExactBucketValue(vector<val_t> & assignment) const
	{
		vector<Function*>::const_iterator itF_end = _RelevantBucketFunctions.end() ;
		vector<MiniBucketElimLHheuristicNode *>::const_iterator itC_end = _Children.end() ;
		double valueMax = ELEM_ZERO ;
#ifdef DEBUG
    bool assign_is_none = false;
#endif
		for (int i = 0 ; i < _k ; ++i) {
#ifdef DEBUG
      if (assignment[_v] == NONE) {
        assign_is_none = true;
      }
#endif
			assignment[_v] = i ;
			double value = ELEM_ONE ;
			for (vector<Function*>::const_iterator itF = _RelevantBucketFunctions.begin() ; itF != itF_end ; ++itF) 
				value OP_TIMESEQ (*itF)->getValue(assignment) ;
			for (vector<MiniBucketElimLHheuristicNode *>::const_iterator itC = _Children.begin() ; itC != itC_end ; ++itC) {
				value OP_TIMESEQ (*itC)->ExactBucketValue(assignment) ;
				if (OUR_OWN_nInfinity == value) 
					goto done_this_iteration ; // no point in continuing with the computation
				}
			valueMax = max(valueMax, value) ;
done_this_iteration : 
			continue ;
			}
#ifdef DEBUG
    if (assign_is_none) {
      assignment[_v] = NONE;
    }
#endif
		return valueMax ;
	}
	MiniBucketElimLHheuristicNode(void) : _H(NULL), _v(-1), _idxWRTparent(-1), _k(-1), _depth2go(-1)
	{
	}
	~MiniBucketElimLHheuristicNode(void)
	{
		Delete() ;
	}
} ;

class MiniBucketElimLHheuristic
{
public :
	MiniBucketElimLH *_H ; // MBE heuristic
public :
	int _v ; // variable of this node.
	int _depth ; // depth of the subtree
	std::vector<MiniBucketElimLHErrorNode *> _SubtreeNodes ; // variables/buckets in the (minimal) subtree
	MiniBucketElimLHErrorNode _RootNode ; // node corresponding to _v
	std::vector<Function *> _IntermediateSubtreeFunctions ; // IF/AF (in bucket of _v) that came from outside of the LH-subteee
public :
	bool IsSubtreeVariable(int v)
	{
		for (vector<MiniBucketElimLHErrorNode *>::iterator itRB = _SubtreeNodes.begin(); itRB!=_SubtreeNodes.end(); ++itRB) {
			if (v == (*itRB)->_v) 
				return true ;
			}
		return false ;
	}
	bool IsSubtreeMBEfunction(Function & f)
	{
		for (vector<MiniBucketElimLHErrorNode *>::iterator itRB = _SubtreeNodes.begin(); itRB!=_SubtreeNodes.end(); ++itRB) {
			MiniBucketElimLHErrorNode *n = *itRB ;
			int u = n->_v ;
			vector<MiniBucket> & minibuckets = _H->_MiniBuckets[u] ;
			for (vector<MiniBucket>::iterator itMB = minibuckets.begin() ; itMB != minibuckets.end() ; ++itMB) {
				Function *fnMB = itMB->output_fn() ;
				if (&f == fnMB) return true ;
				}
			}
		return false ;
	}
	inline double GetHeuristic(vector<val_t> & assignment) const
	{
		if (0 == _SubtreeNodes.size()) 
			return _H->MiniBucketElim::getHeur(_v, assignment, NULL) ; // most likely, all the children(buckets) up to depth d have exactly 1 MB.

		// add up all lambda functions of _v that came from outside the LH subtree; these functions are fully instantiated
		double h = ELEM_ONE ;
		vector<Function*>::const_iterator itF_end = _IntermediateSubtreeFunctions.end() ;
		for (vector<Function*>::const_iterator itF = _IntermediateSubtreeFunctions.begin() ; itF != itF_end ; ++itF) {
			h OP_TIMESEQ (*itF)->getValue(assignment) ;
			if (OUR_OWN_nInfinity == h) 
				return OUR_OWN_nInfinity ; // no point in continuing with the computation
			}

		// add up all children LH values
		vector<MiniBucketElimLHErrorNode *>::const_iterator itC_end = _RootNode._Children.end() ;
		for (vector<MiniBucketElimLHErrorNode *>::const_iterator itC = _RootNode._Children.begin() ; itC != itC_end ; ++itC) {
			h OP_TIMESEQ (*itC)->ExactBucketValue(assignment) ;
			if (OUR_OWN_nInfinity == h) 
				return OUR_OWN_nInfinity ; // no point in continuing with the computation
			}

		return h ;
	}
public :
	int Initialize(MiniBucketElimLH & H, int v, int depth)
	{
		Delete() ;
		_H = &H ; _v = v ; _depth = depth ;
		_RootNode._v = v ; _RootNode._k = (H.m_problem)->getDomainSize(v) ; _RootNode._depth2go = depth ;
		// always check if there is any point in computing error
		if (H._distToClosestDescendantWithLE[v] > depth) 
			return 0 ;
		// build tree of nodes
		int memory = depth*5 ; if (memory > H.m_problem->getN()) memory = H.m_problem->getN() ; // assume avg branching factor of 5
		_SubtreeNodes.reserve(memory) ; // try to allocate some memory to make it faster
//		std::deque<MiniBucketElimLHErrorNode *> nodes2process ;
		std::stack<MiniBucketElimLHErrorNode *> nodes2process ; // using stack will turn it into a DFS processing of the tree
		try { nodes2process.push(&_RootNode) ; } catch (...) { return 1 ; }
		while (nodes2process.size() > 0) {
			MiniBucketElimLHErrorNode *N = nodes2process.top() ; nodes2process.pop() ;
			const PseudotreeNode *N_ = H.m_pseudotree->getNode(N->_v) ;
			const vector<PseudotreeNode *> & children = N_->getChildren() ;
			N->_Children.reserve(children.size()) ; // reserver space for all children (some will not be added as children) so that reallocation is not needed
			int depth2go = N->_depth2go - 1 ; int idxChild = children.size() - 1 ;
			for (vector<PseudotreeNode*>::const_reverse_iterator itC = children.rbegin() ; itC != children.rend(); ++itC, --idxChild) {
				int child = (*itC)->getVar() ;
#ifndef USE_FULL_LOOKAHEAD_SUBTREE
				if (H._BucketErrorQuality[child] <= 1 ? H._distToClosestDescendantWithLE[child] > depth2go : false) continue ; // _BucketErrorQuality <= 1 means it is not proven that there is substantial bucket error
#endif // USE_FULL_LOOKAHEAD_SUBTREE
				MiniBucketElimLHErrorNode *n = NULL ; try { n = new MiniBucketElimLHErrorNode ; } catch (...) { return 1 ; } if (NULL == n) return 1 ; n->_H = &H ; n->_v = child ; n->_k = (H.m_problem)->getDomainSize(child) ; n->_depth2go = depth2go ;
				try { _SubtreeNodes.push_back(n) ; } catch (...) { delete n ; return 1 ; }
				try { (N->_Children).push_back(n) ; } catch (...) { delete n ; return 1 ; }
				if (depth2go > 0) { try { nodes2process.push(n) ; } catch (...) { delete n ; return 1 ; }}
				n->_idxWRTparent = N->_Children.size() - 1 ;
				}
			}
		// fill in _RelevantBucketFunctions for each relevant bucket
		for (vector<MiniBucketElimLHErrorNode *>::iterator itRB = _SubtreeNodes.begin(); itRB!=_SubtreeNodes.end(); ++itRB) {
			MiniBucketElimLHErrorNode *n = *itRB ;
			int u = n->_v ;
//			int d = n->_depth2go ;
			std::vector<Function *> & funs = n->_RelevantBucketFunctions ;
			// add OF_n
			const vector<Function *> & OFlist = H.m_pseudotree->getFunctions(u) ;
			funs.insert(funs.end(), OFlist.begin(), OFlist.end()) ;
			// add AF_n^{>d}; i.e. those functions in AF_n that did not come from a subtree bucket
			vector<Function *> & AFlist = H.m_augmented[u] ;
			for (vector<Function*>::iterator itF = AFlist.begin() ; itF != AFlist.end() ; ++itF) {
				Function *fn = *itF ;
//				int vOriginating = abs(fn->getId()) ;
//				if (! IsSubtreeVariable(vOriginating)) 
				if (! IsSubtreeMBEfunction(*fn)) 
					funs.push_back(fn) ;
				}
			}
		// fill in _IntermediateSubtreeFunctions - IF/AF (in bucket of _v) that came from outside of the LH-subteee
		vector<Function *> & vAFlist = H.m_augmented[v] ;
		for (vector<Function*>::iterator itF = vAFlist.begin() ; itF != vAFlist.end() ; ++itF) {
			Function *fn = *itF ;
			if (! IsSubtreeMBEfunction(*fn)) 
				_IntermediateSubtreeFunctions.push_back(fn) ;
			}
		vector<Function *> & vIFlist = H.m_intermediate[v] ;
		for (vector<Function*>::iterator itF = vIFlist.begin() ; itF != vIFlist.end() ; ++itF) {
			Function *fn = *itF ;
			if (! IsSubtreeMBEfunction(*fn)) 
				_IntermediateSubtreeFunctions.push_back(fn) ;
			}

/*		// fill in a list of relevant functions of the node[v] = all IF/AF in [v] that came from relevant descendant buckets
this code has bug; for all MB output FNs we should check if it is in [v] bucket first
		for (vector<MiniBucketElimLHErrorNode *>::iterator itRB = _SubtreeNodes.begin(); itRB!=_SubtreeNodes.end(); ++itRB) {
			MiniBucketElimLHErrorNode *n = *itRB ;
			int u = n->_v ;
			vector<MiniBucket> & minibuckets = H._MiniBuckets[u] ;
			for (vector<MiniBucket>::iterator itMB = minibuckets.begin() ; itMB != minibuckets.end() ; ++itMB) {
				Function *fnMB = itMB->outputFn() ;
				if (NULL == fnMB) continue ;
				_RootNode._RelevantBucketFunctions.push_back(fnMB) ;
				}
			}*/
		return 0 ;
	}
	int Delete(void) 
	{
		_RootNode.Delete() ;
		_SubtreeNodes.clear() ;
		_IntermediateSubtreeFunctions.clear() ;
		return 0 ;
	}
	MiniBucketElimLHheuristic(void) : _H(NULL), _v(-1), _depth(-1)
	{
	}
	~MiniBucketElimLHheuristic(void)
	{
	}
} ;

}  // namespace daoopt

#endif /* MiniBucketElimLH_H_ */
