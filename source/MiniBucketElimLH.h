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

namespace daoopt {

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
protected :

	MiniBucketElimLHStatistics _Stats ;

	// joint scope of each bucket, including the bucket's var
	std::vector<std::set<int>> _BucketScopes ;
	std::vector<std::vector<Function *>> _BucketFunctions ;
//	std::vector<std::set<int>> _BucketHFnScopes ;

	// a set of minibuckets, one for each var
	std::vector<std::vector<MiniBucket>> _MiniBuckets ;

	std::vector<Function *> _BucketErrorFunctions ;
	std::vector<signed char> _BucketErrorQuality ; // for each var, whether the bucketerror[var]=0; -1 means unknown. 0=bucket error is 0, 1=bucket error is >0.
	double _BuckerErrorFnTableSizes_Total ; // log
	double _BuckerErrorFnTableSizes_Precomputed ; // log
	double _BuckerErrorFnTableSizes_Ignored ; // log
	int64 _nBucketsWithNonZeroBuckerError ;
	int64 _nBucketsWithMoreThan1MB ;

	// for each var, distance to the closest descendant with more than 1 MB/BucketErrorQuality>0. INT_MAX means infinite (no descendants).
	// this does not include the node itself; i.e. legal values are >0 & <= INT_MAX (meaning infinity = no such descendants).
	std::vector<int> _distToClosestDescendantWithMBs ;
	std::vector<int> _distToClosestDescendantWithLE ;

public:

	// builds the heuristic, limited to the relevant subproblem, if applicable.
	// if computeTables=false, only returns size estimate (no tables computed)
	virtual size_t build(const std::vector<val_t>* assignment = NULL, bool computeTables = true) ;

	// reset the data structures
	virtual void reset() ;

	// computes the heuristic for variable var given a (partial) assignment
	virtual double getHeur(int var, std::vector<val_t>& assignment,
                         SearchNode* n);
	// computes heuristic values for all instantiations of var, given context assignment
	virtual void getHeurAll(int var, std::vector<val_t>& assignment, 
                          SearchNode* n, std::vector<double>& out);

	// get heuristic function h error by looking ahead. this allows to improve the value of h.
	// lookaheadDepth==0 means just the bucket of 'var'; each increment of lookaheadDepth means going down one level to children.
	// we assume in the input assignment, the variables in the output fn of bucket[var] are all instantiated.
	virtual double getHeuristicError(int var, vector<val_t> & assignment, int lookaheadDepth) ;

	// computes the lookahead heuristic for variable var given a (partial) assignment.
	// note : this fn returns BucketError[var]. all variables in the output fn of bucket[var] should be instantiated by 'assignment'.
	// note : the return value is in the same representation space (log or normal) as the problem.
	// note : the return values should be >= 0.0, regardless of whether normal/log representation.
	virtual double getLocalError(int var, vector<val_t> & assignment) ;

	// compute local MiniBucket error, defined as the different between 
	//		1) the combination of all the mini-bucket output functions,
	//		2) what the normal/full bucket elimination would compute. 
	// note : error_fn/avg_error/avg_exact are in the same representation space (normal/log) as the problem itself.
	// note : localerror should always be >=0, regardless of whether normal/log representation, since item_1 is at least as large as item_2 (above).
	// note : the scope of the table (as fn) is the same as the bucket output fn.
	// note : this fn sets _BucketErrorQuality[var] to 0/1 when it can be determined; sometimes when bucket error table is all 0, we won't compute the table and will set _BucketErrorQuality[var]=0.
	int computeLocalErrorTable(int var, bool build_table, 
		double TableMemoryLimitAsNumElementsLog, 
		double & TableSizeLog, // OUT : bucket error table size, regardless of whether it is actually computed; this is in log scale, i.e. sum_log10(var_domain_size).
		double & avgError, // OUT : avg bucket error; computed when output table is not too large.
		double & avgExact, // OUT : avg value over entire bucket output table; computed when error is computed.
		Function * & error_fn) ; // OUT : bucket error fn; computed when output table is not too large.
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

public :

	MiniBucketElimLH(Problem *p, Pseudotree *pt, ProgramOptions *po, int ib) ;
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

inline MiniBucketElimLH::~MiniBucketElimLH(void)
{
	deleteLocalErrorFNs() ;
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

}  // namespace daoopt

#endif /* MiniBucketElimLH_H_ */
