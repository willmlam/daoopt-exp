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
 *  Created on: May 20, 2010
 *      Author: Lars Otten <lotten@ics.uci.edu>
 */

#ifndef MINIBUCKET_H_
#define MINIBUCKET_H_

#include "Function.h"
#include "Problem.h"


/* A single minibucket, i.e. a collection of functions */
class MiniBucket {

protected:
  int m_bucketVar;               // the bucket variable
  int m_ibound;                  // the ibound
  Problem* m_problem;      // pointer to the bucket elimination structure
  vector<Function*> m_functions; // the functions in the MB
  set<int> m_jointScope;         // keeps track of the joint scope if the functions

public:
  // checks whether the MB has space for a function
  bool allowsFunction(Function*);
  // adds a function to the minibucket and returns its index
  int addFunction(Function*);

  const set<int> &getJointScope() const;

  Function *&getFunctionRef(int i);

  int getVar() const;

  // Joins the MB functions, eliminate the bucket variable, and returns the resulting function
  // set buildTable==false to get only size estimate (table will not be computed)
  Function* eliminate(bool buildTable=true);

  // eliminates the specified set of variables instead of the bucket variables
  Function* eliminate(bool buildTable, const set<int> &elimVars);

  // eliminates the specified set of variables with conditioning
  Function* conditionEliminate(bool buildTable, const map<int,val_t> &cond, const set<int> &elimVars);

  // eliminates the specified set of variables with conditioning
  void conditionEliminateInPlace(bool buildTable, const map<int,val_t> &cond, const set<int> &elimVars, double *table);

  // eliminates the variable while applying max-marginal matching
  Function* eliminateMM(bool buildTable, Function *maxMarginal, Function *avgMaxMarginal);

  // eliminates the variable while applying max-marginal matching with conditioning
  Function* conditionEliminateMM(bool buildTable, const map<int,val_t> &cond, Function *maxMarginal, Function *avgMaxMarginal);

  // Joins the MB functions without elimination
  Function* join(bool buildTable=true);

public:
  MiniBucket(int var, int bound, Problem* p);

};


/* Inline definitions */

inline int MiniBucket::addFunction(Function* f) {
  assert(f);
  // insert function
  m_functions.push_back(f);
  // update joint scope
  m_jointScope.insert(f->getScopeVec().begin(), f->getScopeVec().end() );
  return m_functions.size() - 1;
}

inline Function *&MiniBucket::getFunctionRef(int i) {
    assert(i < m_functions.size());
    return m_functions[i];
}

inline MiniBucket::MiniBucket(int v, int b, Problem* p) :
  m_bucketVar(v), m_ibound(b), m_problem(p) {}

inline const set<int> &MiniBucket::getJointScope() const {
    return m_jointScope;
}

inline int MiniBucket::getVar() const {
    return m_bucketVar;
}


#endif /* MINIBUCKET_H_ */
