/*
 * BoundPropagator.h
 *
 *  Created on: Nov 7, 2008
 *      Author: Lars Otten <lotten@ics.uci.edu>
 */

#ifndef BOUNDPROPAGATOR_H_
#define BOUNDPROPAGATOR_H_

#include "SearchSpace.h"
#include "SearchNode.h"
#include "Pseudotree.h"
#ifdef PARALLEL_STATIC
#include "Statistics.h"
#endif

class BoundPropagator {

protected:

  bool m_doCaching;
  Problem*     m_problem;
  SearchSpace* m_space;

#ifdef PARALLEL_STATIC
  count_t m_subCountCache;
  SubproblemStats m_subStatsCache;
#endif

#ifdef PARALLEL_DYNAMIC
  /* caches the root variable of the last deleted subproblem */
  int m_subRootvarCache;
  /* caches the size of the last deleted subproblem */
  count_t m_subCountCache;
  /* caches the number of leaf nodes in the deleted subproblem */
  count_t m_subLeavesCache;
  /* caches the cumulative leaf depth in the deleted subproblem */
  count_t m_subLeafDCache;
  /* caches the lower/upper bound of the last highest deleted node */
  pair<double,double> m_subBoundsCache;
#endif

public:

  /*
   * propagates the value of the specified search node and removes unneeded nodes
   * returns a pointer to the parent of the highest deleted node
   * @n: the search node to be propagated
   * @reportSolution: should root updates be reported to problem instance?
   */
  SearchNode* propagate(SearchNode* n, bool reportSolution = false, SearchNode* upperLimit = NULL);

#ifdef PARALLEL_STATIC
  const SubproblemStats& getSubproblemStatsCache() const { return m_subStatsCache; }
  count_t getSubCountCache() const { return m_subCountCache; }
  void resetSubCount() { m_subCountCache = 0; }
#endif

#ifdef PARALLEL_DYNAMIC
  int getSubRootvarCache() const { return m_subRootvarCache; }
  count_t getSubCountCache() const { return m_subCountCache; }
  count_t getSubLeavesCache() const { return m_subLeavesCache; }
  count_t getSubLeafDCache() const { return m_subLeafDCache; }
  const pair<double,double>& getBoundsCache() const { return m_subBoundsCache; }
#endif

#ifndef NO_ASSIGNMENT
private:
  void propagateTuple(SearchNode* start, SearchNode* end);
#endif

protected:
  virtual bool isMaster() const { return false; }

public:
  BoundPropagator(Problem* p, SearchSpace* s, bool doCaching = true)
    : m_doCaching(doCaching), m_problem(p), m_space(s)
#if defined PARALLEL_DYNAMIC || defined PARALLEL_STATIC
  , m_subCountCache(0)
#endif
  { /* empty */ }
  virtual ~BoundPropagator() {}
};


#endif /* BOUNDSPROPAGATOR_H_ */
