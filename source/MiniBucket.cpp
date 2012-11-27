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
 *  Created on: May 20, 2010
 *      Author: Lars Otten <lotten@ics.uci.edu>
 */

#include "MiniBucket.h"

#undef DEBUG


/* checks whether a function 'fits' in this MB */
bool MiniBucket::allowsFunction(Function* f) {

  const set<int>& a=m_jointScope, b=f->getScopeSet();

  set<int>::const_iterator ita=a.begin(),
    itb = b.begin() ;

  int s = 0, d=0;

  while (ita != a.end() && itb != b.end()) {
    d = *ita - *itb;
    if (d>0) {
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

  // accept if no scope increase or new scope not greater than ibound
  return (s == (int) m_jointScope.size()) || (s <= m_ibound+1);
  // new scope would be greater then ibound?
  //return s <= m_ibound+1;

}


/* joins the functions in the MB while marginalizing out the bucket var.,
 * resulting function is returned */
Function* MiniBucket::eliminate(bool buildTable) {

#ifdef DEBUG
  cout << "  Marginalizing together:";
  for (vector<Function*>::iterator it=m_functions.begin();it!=m_functions.end();++it)
    cout << ' ' << (**it);
  cout << endl;
#endif

  set<int> scope;
  int i=0; size_t j=0; // loop variables

  vector<Function*>::const_iterator fit;
  set<int>::const_iterator sit;

  for (fit = m_functions.begin(); fit != m_functions.end(); ++fit) {
    scope.insert( (*fit)->getScopeVec().begin() , (*fit)->getScopeVec().end() );
  }

#ifdef DEBUG
  cout << "   Joint scope: " << scope << endl;
#endif

  // remove elimVar from the new scope
  scope.erase(m_bucketVar);
  int n = scope.size(); // new function arity

#ifdef DEBUG
  cout << "   Target scope: " << scope << endl;
#endif

  // compute new table size and collect domain sizes
  vector<val_t> domains;
  domains.reserve(n);
  size_t tablesize = 1;
  for (sit = scope.begin(); sit!=scope.end(); ++sit) {
//    if (*it != elimVar)
    tablesize *= m_problem->getDomainSize(*sit);
    domains.push_back(m_problem->getDomainSize(*sit));
  }

  double* newTable = NULL;
  if (buildTable) {
    newTable = new double[tablesize];
    for (j=0; j<tablesize; ++j) newTable[j] = ELEM_ZERO;

    // this keeps track of the tuple assignment
    val_t* tuple = new val_t[n+1];
    for (i=0; i<n; ++i) tuple[i] = 0; // i trough n index target variables
    val_t* elimVal = &tuple[n]; // n+1 is elimVar

    // maps each function scope assignment to the full tuple
    vector<vector<val_t*> > idxMap(m_functions.size());

    // holds iterators .begin() and .end() for all function scopes
    vector< pair< set<int>::const_iterator , set<int>::const_iterator > > iterators;
    iterators.reserve(m_functions.size());
    for (fit = m_functions.begin(); fit != m_functions.end(); ++fit) {
      // store begin() and end() for each function scope
      iterators.push_back( make_pair( (*fit)->getScopeSet().begin(), (*fit)->getScopeSet().end() ) );
    }

    // collect pointers to tuple values
    bool bucketVarPassed = false;
    for (i=0, sit=scope.begin(); i<n; ++i, ++sit) {
      if (!bucketVarPassed && *sit > m_bucketVar) { // just went past bucketVar
        for (j=0; j<m_functions.size(); ++j) {
          idxMap[j].push_back(elimVal);
          ++(iterators[j].first); // skip bucketVar in the original function scope
        }
        bucketVarPassed = true;
      }
      for (j=0; j<m_functions.size(); ++j) {
        //      cout << "  f" << funs[j]->getId() << ' ' << *sit << " == " << *(iterators[j].first) << endl;
        if (iterators[j].first != iterators[j].second // scope iterator != end()
            && *sit == *(iterators[j].first)) { // value found
          idxMap[j].push_back(&tuple[i]);
          ++(iterators[j].first);
        }
      }
    }

    if (!bucketVarPassed) { // bucketVar has highest index
      for (j=0; j<m_functions.size(); ++j) {
          if (!m_functions[j]->isConstant())
              idxMap[j].push_back(elimVal);
      }
    }

    // actual computation
    size_t idx; double z;
//    val_t valMax;
    // iterate over all values of elimVar
    for (*elimVal = 0; *elimVal < m_problem->getDomainSize(m_bucketVar); ++(*elimVal) ) {
      idx=0; // go over the full new table
      do {
        z = ELEM_ONE;
        for (j=0; j<m_functions.size(); ++j)
          z OP_TIMESEQ m_functions[j]->getValuePtr(idxMap[j]);
//        if (z>newTable[idx]) valMax = *elimVal;
        newTable[idx] = max(newTable[idx],z);
      } while ( increaseTuple(idx,tuple,domains) );

    }
    // DEBUG
//    cout << "Maximum value: " << int(valMax) << endl;
    // DEBUG
    // clean up
    delete[] tuple;
  }

  return new FunctionBayes(-m_bucketVar,m_problem,scope,newTable,tablesize);
}

/* joins the functions in the MB
 * resulting function is returned */
Function* MiniBucket::join(bool buildTable) {

#ifdef DEBUG
  cout << "  Joining together:";
  for (vector<Function*>::iterator it=m_functions.begin();it!=m_functions.end();++it)
    cout << ' ' << (**it);
  cout << endl;
  for (vector<Function*>::iterator it=m_functions.begin();it!=m_functions.end();++it) {
      Function * fT = *it;
      for (int i=0; i<fT->getTableSize(); ++i) {
          cout << ' ' << fT->getTable()[i] << endl;
      }
      cout << endl;
  }
  cout << endl;
#endif

  set<int> scope;
  int i=0; size_t j=0; // loop variables

  vector<Function*>::const_iterator fit;
  set<int>::const_iterator sit;

  for (fit = m_functions.begin(); fit != m_functions.end(); ++fit) {
    scope.insert( (*fit)->getScopeVec().begin() , (*fit)->getScopeVec().end() );
  }

#ifdef DEBUG
  cout << "   Joint scope: " << scope << endl;
#endif

//  scope.erase(m_bucketVar);
  int n = scope.size(); // new function arity

#ifdef DEBUG
  cout << "   Target scope: " << scope << endl;
#endif

  // compute new table size and collect domain sizes
  vector<val_t> domains;
  domains.reserve(n);
  size_t tablesize = 1;
  for (sit = scope.begin(); sit!=scope.end(); ++sit) {
//    if (*it != elimVar)
    tablesize *= m_problem->getDomainSize(*sit);
    domains.push_back(m_problem->getDomainSize(*sit));
  }

  double* newTable = NULL;
  if (buildTable) {
    newTable = new double[tablesize];
    for (j=0; j<tablesize; ++j) newTable[j] = ELEM_ZERO;

    // this keeps track of the tuple assignment
    val_t* tuple = new val_t[n];
    for (i=0; i<n; ++i) tuple[i] = 0; // i trough n index target variables
//    val_t* elimVal = &tuple[n]; // n+1 is elimVar

    // maps each function scope assignment to the full tuple
    vector<vector<val_t*> > idxMap(m_functions.size());

    // holds iterators .begin() and .end() for all function scopes
    vector< pair< set<int>::const_iterator , set<int>::const_iterator > > iterators;
    iterators.reserve(m_functions.size());
    for (fit = m_functions.begin(); fit != m_functions.end(); ++fit) {
      // store begin() and end() for each function scope
      iterators.push_back( make_pair( (*fit)->getScopeSet().begin(), (*fit)->getScopeSet().end() ) );
    }

    // collect pointers to tuple values
//    bool bucketVarPassed = false;
    for (i=0, sit=scope.begin(); i<n; ++i, ++sit) {
    	/*
      if (!bucketVarPassed && *sit > m_bucketVar) { // just went past bucketVar
        for (j=0; j<m_functions.size(); ++j) {
          idxMap[j].push_back(elimVal);
          ++(iterators[j].first); // skip bucketVar in the original function scope
        }
        bucketVarPassed = true;
      }
      */
      for (j=0; j<m_functions.size(); ++j) {
        //      cout << "  f" << funs[j]->getId() << ' ' << *sit << " == " << *(iterators[j].first) << endl;
        if (iterators[j].first != iterators[j].second // scope iterator != end()
            && *sit == *(iterators[j].first)) { // value found
          idxMap[j].push_back(&tuple[i]);
          ++(iterators[j].first);
        }
      }
    }

    /*
    if (!bucketVarPassed) { // bucketVar has highest index
      for (j=0; j<m_functions.size(); ++j) {
        idxMap[j].push_back(elimVal);
      }
    }
    */

    // actual computation
    size_t idx; double z;
    // iterate over all values of elimVar
    idx=0; // go over the full new table
    do {
    	z = ELEM_ONE;
    	for (j=0; j<m_functions.size(); ++j)
    		z OP_TIMESEQ m_functions[j]->getValuePtr(idxMap[j]);
    	newTable[idx] = z;
    } while ( increaseTuple(idx,tuple,domains) );

    // clean up
    delete[] tuple;
  }

  return new FunctionBayes(-m_bucketVar,m_problem,scope,newTable,tablesize);
}


/* joins the functions in the MB while eliminating the specified set of variables,
 * resulting function is returned */
Function* MiniBucket::eliminate(bool buildTable, const set<int> &elimVars) {

#ifdef DEBUG
  cout << "  Marginalizing together:";
  for (vector<Function*>::iterator it=m_functions.begin();it!=m_functions.end();++it)
    cout << ' ' << (**it);
  cout << endl;
#endif

  set<int> scope;
  int i=0; size_t j=0; // loop variables

  vector<Function*>::const_iterator fit;
  set<int>::const_iterator sit;
  set<int>::const_iterator eit;

  for (fit = m_functions.begin(); fit != m_functions.end(); ++fit) {
    scope.insert( (*fit)->getScopeVec().begin() , (*fit)->getScopeVec().end() );
  }

#ifdef DEBUG
  cout << "   Joint scope: " << scope << endl;
#endif
  vector<val_t> elimDomains;
  elimDomains.reserve(elimVars.size());

  for (eit = elimVars.begin(); eit != elimVars.end(); ++eit) {
      // remove elimVars from the new scope
      scope.erase(*eit);
      elimDomains.push_back(m_problem->getDomainSize(*eit));
  }
  int n = scope.size(); // new function arity

#ifdef DEBUG
  cout << "   Target scope: " << scope << endl;
#endif

  // compute new table size and collect domain sizes
  vector<val_t> domains;
  domains.reserve(n);
  size_t tablesize = 1;
  for (sit = scope.begin(); sit!=scope.end(); ++sit) {
//    if (*it != elimVar)
    tablesize *= m_problem->getDomainSize(*sit);
    domains.push_back(m_problem->getDomainSize(*sit));
  }

  double* newTable = NULL;
  if (buildTable) {
    newTable = new double[tablesize];
    for (j=0; j<tablesize; ++j) newTable[j] = ELEM_ZERO;

    // this keeps track of the tuple assignment
    val_t* tuple = new val_t[n+elimVars.size()];
    val_t* elimTuple = tuple + n;
    for (i=0; i<int(n+elimVars.size()); ++i) tuple[i] = 0; // i trough n index target variables
//    val_t* elimVal = &tuple[n]; // n+1 is elimVar

    // maps each function scope assignment to the full tuple
    vector<vector<val_t*> > idxMap(m_functions.size());

    // holds iterators .begin() and .end() for all function scopes
    vector< pair< set<int>::const_iterator , set<int>::const_iterator > > iterators;
    iterators.reserve(m_functions.size());
    for (fit = m_functions.begin(); fit != m_functions.end(); ++fit) {
      // store begin() and end() for each function scope
      iterators.push_back( make_pair( (*fit)->getScopeSet().begin(), (*fit)->getScopeSet().end() ) );
    }

    // collect pointers to tuple values
//    bool bucketVarPassed = false;
    eit = elimVars.begin();
    int eVarCount = 0;
    for (i=0, sit=scope.begin(); i<n; ++i, ++sit) {
//      if (!bucketVarPassed && *sit > m_bucketVar) { // just went past bucketVar
      for (; eit != elimVars.end() && *sit > *eit; ++eit, ++eVarCount) { 
        for (j=0; j<m_functions.size(); ++j) {
            if (m_functions[j]->getScopeSet().find(*eit) != 
                    m_functions[j]->getScopeSet().end()) {
                idxMap[j].push_back(&tuple[n+eVarCount]);
                ++(iterators[j].first); // skip elimVar in the original function scope
            }
        }
      }
      for (j=0; j<m_functions.size(); ++j) {
        //      cout << "  f" << funs[j]->getId() << ' ' << *sit << " == " << *(iterators[j].first) << endl;
        if (iterators[j].first != iterators[j].second // scope iterator != end()
            && *sit == *(iterators[j].first)) { // value found
          idxMap[j].push_back(&tuple[i]);
          ++(iterators[j].first);
        }
      }
    }

    // one or more of the elimVars have higher indices
    for (; eit != elimVars.end(); ++eit, ++eVarCount)
        for (j=0; j<m_functions.size(); ++j)
            if (m_functions[j]->getScopeSet().find(*eit) != 
                m_functions[j]->getScopeSet().end())
                idxMap[j].push_back(&tuple[n+eVarCount]);

    // actual computation
    size_t idx; double z;
    // iterate over all values of elimVar
    do { 
      idx=0; // go over the full new table
      do {
        z = ELEM_ONE;
        for (j=0; j<m_functions.size(); ++j) {
            /*
            cout << "On function " << *m_functions[j] << endl;
            cout << idxMap[j].size() << endl;
            cout << m_functions[j]->getScopeVec().size() << endl;
            */
          z OP_TIMESEQ m_functions[j]->getValuePtr(idxMap[j]);
        }
#ifdef DEBUG 
        // TUPLE OUTPUT DEBUG
        cout << "Tuple: ";
        for (j=0; j<n+elimVars.size(); ++j) {
            cout << ' ' << int(tuple[j]);
        }
        cout << "  val: " << ELEM_DECODE(z);
        cout << endl;
        // DEBUG
#endif
        newTable[idx] = max(newTable[idx],z);
      } while ( increaseTuple(idx,tuple,domains) );
    } while (increaseTuple(idx,elimTuple,elimDomains));

    // clean up
    delete[] tuple;
  }

  return new FunctionBayes(-m_bucketVar,m_problem,scope,newTable,tablesize);
}

/* joins the functions in the MB while conditoining and eliminating the 
 * specified set of variables, resulting function is returned */
Function* MiniBucket::conditionEliminate(bool buildTable, const map<int,val_t> &cond, const set<int> &elimVarsIn) {

#ifdef DEBUG
  cout << "  Marginalizing together:";
  for (vector<Function*>::iterator it=m_functions.begin();it!=m_functions.end();++it)
    cout << ' ' << (**it);
  cout << endl;
#endif

  set<int> scope;
  int i=0; size_t j=0; // loop variables

  vector<Function*>::const_iterator fit;
  set<int>::const_iterator sit;
  set<int>::const_iterator eit;
  map<int,val_t>::const_iterator cit;

  for (fit = m_functions.begin(); fit != m_functions.end(); ++fit) {
    scope.insert( (*fit)->getScopeVec().begin() , (*fit)->getScopeVec().end() );
  }

#ifdef DEBUG
  cout << "   Joint scope: " << scope << endl;
#endif

  // create a copy and remove conditioning variables
  set<int> elimVars(elimVarsIn);
  for (cit = cond.begin(); cit != cond.end(); ++cit) {
      elimVars.erase(cit->first);
  }

  vector<val_t> elimDomains;
  elimDomains.reserve(elimVars.size());

  // remove elimVars from the new scope
  for (eit = elimVars.begin(); eit != elimVars.end(); ++eit) {
      scope.erase(*eit);
      elimDomains.push_back(m_problem->getDomainSize(*eit));
  }
  // remove conditioned variables from the new scope
  for (cit = cond.begin(); cit != cond.end(); ++cit) {
      scope.erase(cit->first);
  }
  int n = scope.size(); // new function arity

#ifdef DEBUG
  cout << "   Target scope: " << scope << endl;
#endif

  // compute new table size and collect domain sizes
  vector<val_t> domains;
  domains.reserve(n);
  size_t tablesize = 1;
  for (sit = scope.begin(); sit!=scope.end(); ++sit) {
//    if (*it != elimVar)
    tablesize *= m_problem->getDomainSize(*sit);
    domains.push_back(m_problem->getDomainSize(*sit));
  }

  double* newTable = NULL;
  if (buildTable) {
    newTable = new double[tablesize];
    for (j=0; j<tablesize; ++j) newTable[j] = ELEM_ZERO;

    // this keeps track of the tuple assignment
    val_t* tuple = new val_t[n+elimVars.size()+cond.size()];
    val_t* elimTuple = tuple + n;

    for (i=0; i<int(n+elimVars.size()); ++i) tuple[i] = 0; // i trough n index target variables
    for (cit = cond.begin(); cit != cond.end(); ++cit, ++i) tuple[i] = cit->second; 
//    val_t* elimVal = &tuple[n]; // n+1 is elimVar

    // maps each function scope assignment to the full tuple
    vector<vector<val_t*> > idxMap(m_functions.size());

    // holds iterators .begin() and .end() for all function scopes
    vector< pair< set<int>::const_iterator , set<int>::const_iterator > > iterators;
    iterators.reserve(m_functions.size());
    for (fit = m_functions.begin(); fit != m_functions.end(); ++fit) {
      // store begin() and end() for each function scope
      iterators.push_back( make_pair( (*fit)->getScopeSet().begin(), (*fit)->getScopeSet().end() ) );
    }

    // collect pointers to tuple values
//    bool bucketVarPassed = false;
    eit = elimVars.begin();
    cit = cond.begin();
    int eVarCount = 0;
    int cVarCount = 0;
    for (i=0, sit=scope.begin(); i<n; ++i, ++sit) {
//      if (!bucketVarPassed && *sit > m_bucketVar) { // just went past bucketVar

      // ensure we process elim and cond vars first
      while ((eit != elimVars.end() && *sit > *eit) || (cit != cond.end() && *sit > cit->first)) {
          // process only if scope and cond comes after
          for (; eit != elimVars.end() && *sit > *eit && (cit == cond.end() || cit->first > *eit); 
                  ++eit, ++eVarCount) { 
              for (j=0; j<m_functions.size(); ++j) {
                  if (m_functions[j]->getScopeSet().find(*eit) != 
                          m_functions[j]->getScopeSet().end()) {
                      idxMap[j].push_back(&tuple[n+eVarCount]);
                      ++(iterators[j].first); // skip elimVar in the original function scope
                  }
              }
          }
          // process only if scope and elim comes after
          for (; cit != cond.end() && *sit > cit->first && (eit == elimVars.end() || *eit > cit->first); 
                  ++cit, ++cVarCount) {
              for (j=0; j<m_functions.size(); ++j) {
                  if (m_functions[j]->getScopeSet().find(cit->first) !=
                          m_functions[j]->getScopeSet().end()) {
                      idxMap[j].push_back(&tuple[n+elimVars.size()+cVarCount]);
                      ++(iterators[j].first);
                  }
              }
          }
      }
      for (j=0; j<m_functions.size(); ++j) {
        //      cout << "  f" << funs[j]->getId() << ' ' << *sit << " == " << *(iterators[j].first) << endl;
        if (iterators[j].first != iterators[j].second // scope iterator != end()
            && *sit == *(iterators[j].first)) { // value found
          idxMap[j].push_back(&tuple[i]);
          ++(iterators[j].first);
        }
      }
    }

    // elimVars and/or cond variables have higher indices
    while (eit != elimVars.end() || cit != cond.end()) {
        for (; eit != elimVars.end() && (cit == cond.end() || *eit < cit->first); ++eit, ++eVarCount)
            for (j=0; j<m_functions.size(); ++j)
                if (m_functions[j]->getScopeSet().find(*eit) != 
                        m_functions[j]->getScopeSet().end())
                    idxMap[j].push_back(&tuple[n+eVarCount]);
        for (; cit != cond.end() && (eit == elimVars.end() || cit->first < *eit); ++cit, ++cVarCount)
            for (j=0; j<m_functions.size(); ++j)
                if (m_functions[j]->getScopeSet().find(cit->first) != 
                        m_functions[j]->getScopeSet().end())
                    idxMap[j].push_back(&tuple[n+elimVars.size()+cVarCount]);
    }

    // actual computation
    size_t idx; double z;
    // iterate over all values of elimVar
    do { 
      idx=0; // go over the full new table
      do {
        z = ELEM_ONE;
        for (j=0; j<m_functions.size(); ++j) {
            /*
            cout << "On function " << *m_functions[j] << endl;
            cout << idxMap[j].size() << endl;
            cout << m_functions[j]->getScopeVec().size() << endl;
            */
          z OP_TIMESEQ m_functions[j]->getValuePtr(idxMap[j]);
        }
#if 0
        // TUPLE OUTPUT DEBUG
        cout << "Tuple: ";
        for (j=0; j<n+elimVars.size(); ++j) {
            cout << ' ' << int(tuple[j]);
        }
        cout << "  val: " << ELEM_DECODE(z);
        cout << endl;
        // DEBUG
#endif
        newTable[idx] = max(newTable[idx],z);
      } while ( increaseTuple(idx,tuple,domains) );
    } while (increaseTuple(idx,elimTuple,elimDomains));

    // clean up
    delete[] tuple;
  }

  return new FunctionBayes(-m_bucketVar,m_problem,scope,newTable,tablesize);
}

/* joins the functions in the MB while marginalizing out the bucket var. while 
 * also applying max-marginal matching given the max-marginals
 * resulting function is returned */
Function* MiniBucket::eliminateMM(bool buildTable, 
        Function *maxMarginal, Function *avgMaxMarginal) {

#ifdef DEBUG
  cout << "  Marginalizing together:";
  for (vector<Function*>::iterator it=m_functions.begin();it!=m_functions.end();++it)
    cout << ' ' << (**it);
  cout << endl;
#endif

  set<int> scope;
  int i=0; size_t j=0; // loop variables

  vector<Function*>::const_iterator fit;
  set<int>::const_iterator sit;

  for (fit = m_functions.begin(); fit != m_functions.end(); ++fit) {
    scope.insert( (*fit)->getScopeVec().begin() , (*fit)->getScopeVec().end() );
  }

#ifdef DEBUG
  cout << "   Joint scope: " << scope << endl;
#endif

  // remove elimVar from the new scope
  scope.erase(m_bucketVar);
  int n = scope.size(); // new function arity

#ifdef DEBUG
  cout << "   Target scope: " << scope << endl;
#endif

  // compute new table size and collect domain sizes
  vector<val_t> domains;
  domains.reserve(n);
  size_t tablesize = 1;
  for (sit = scope.begin(); sit!=scope.end(); ++sit) {
//    if (*it != elimVar)
    tablesize *= m_problem->getDomainSize(*sit);
    domains.push_back(m_problem->getDomainSize(*sit));
  }

  double* newTable = NULL;
  if (buildTable) {
    newTable = new double[tablesize];
    for (j=0; j<tablesize; ++j) newTable[j] = ELEM_ZERO;

    // this keeps track of the tuple assignment
    val_t* tuple = new val_t[n+1];
    for (i=0; i<n; ++i) tuple[i] = 0; // i trough n index target variables
    val_t* elimVal = &tuple[n]; // n+1 is elimVar

    // maps each function scope assignment to the full tuple
    vector<vector<val_t*> > idxMap(m_functions.size() + 1);

/*
    // maps the maxMarginal scope to the full tuple
    vector<val_t*> mmMap;
    */

    // holds iterators .begin() and .end() for all function scopes
    vector< pair< set<int>::const_iterator , set<int>::const_iterator > > iterators;
    iterators.reserve(m_functions.size() + 1);
    for (fit = m_functions.begin(); fit != m_functions.end(); ++fit) {
      // store begin() and end() for each function scope
      iterators.push_back( make_pair( (*fit)->getScopeSet().begin(), (*fit)->getScopeSet().end() ) );
    }
    iterators.push_back( make_pair( maxMarginal->getScopeSet().begin(),
        maxMarginal->getScopeSet().end() ) );

    // collect pointers to tuple values
    bool bucketVarPassed = false;
    for (i=0, sit=scope.begin(); i<n; ++i, ++sit) {
      if (!bucketVarPassed && *sit > m_bucketVar) { // just went past bucketVar
        for (j=0; j<m_functions.size() + 1; ++j) {
          idxMap[j].push_back(elimVal);
          ++(iterators[j].first); // skip bucketVar in the original function scope
        }
        bucketVarPassed = true;
      }
      for (j=0; j<m_functions.size() + 1; ++j) {
        //      cout << "  f" << funs[j]->getId() << ' ' << *sit << " == " << *(iterators[j].first) << endl;
        if (iterators[j].first != iterators[j].second // scope iterator != end()
            && *sit == *(iterators[j].first)) { // value found
          idxMap[j].push_back(&tuple[i]);
          ++(iterators[j].first);
        }
      }
    }

    if (!bucketVarPassed) { // bucketVar has highest index
      for (j=0; j<m_functions.size() + 1; ++j) {
        idxMap[j].push_back(elimVal);
      }
    }

    // idxmap for max-marginal functions
    vector<val_t*> &mmMap = idxMap[m_functions.size()];

    // actual computation
    size_t idx; double z; 
//    val_t valMax;;
    // iterate over all values of elimVar
    for (*elimVal = 0; *elimVal < m_problem->getDomainSize(m_bucketVar); ++(*elimVal) ) {
      idx=0; // go over the full new table
      do {
        z = ELEM_ONE;
        for (j=0; j<m_functions.size(); ++j) {
    		z OP_TIMESEQ m_functions[j]->getValuePtr(idxMap[j]);
        }
//        z = OP_MINUS(z, maxMarginal->getValuePtr(mmMap));
#ifdef DEBUG
        cout << "avg max-marginal: " << ELEM_DECODE(avgMaxMarginal->getValuePtr(mmMap)) << endl;
#endif
        z OP_TIMESEQ avgMaxMarginal->getValuePtr(mmMap) - maxMarginal->getValuePtr(mmMap);
        /*
        z = ELEM_ENCODE(ELEM_DECODE(z) - 
                ELEM_DECODE(maxMarginal->getValuePtr(mmMap)) +
                ELEM_DECODE(avgMaxMarginal->getValuePtr(mmMap)));
        */
#ifdef DEBUG
        cout << "Final adjusted value: " << ELEM_DECODE(z) << endl;
#endif
        newTable[idx] = max(newTable[idx],z);
//        if (z > newTable[idx]) valMax = *elimVal;
      } while ( increaseTuple(idx,tuple,domains) );

    }
    // DEBUG
//    cout << "Maximum value: " << int(valMax) << endl;
    // DEBUG
    // clean up
    delete[] tuple;
  }

  return new FunctionBayes(-m_bucketVar,m_problem,scope,newTable,tablesize);
}

/* joins the functions in the MB while conditoning a subset of variables, 
 * marginalizing out the bucket var. while also applying max-marginal 
 * matching given the max-marginals
 * resulting function is returned */
Function* MiniBucket::conditionEliminateMM(bool buildTable, 
        const map<int,val_t> &assign, Function *maxMarginal, Function *avgMaxMarginal) {

#ifdef DEBUG
  cout << "  Marginalizing together:";
  for (vector<Function*>::iterator it=m_functions.begin();it!=m_functions.end();++it)
    cout << ' ' << (**it);
  cout << endl;
#endif

  set<int> scope;
  map<int,val_t> cond(assign);
  int i=0; size_t j=0; // loop variables

  vector<Function*>::const_iterator fit;
  set<int>::const_iterator sit;
  map<int,val_t>::iterator cit;

  for (fit = m_functions.begin(); fit != m_functions.end(); ++fit) {
    scope.insert( (*fit)->getScopeVec().begin() , (*fit)->getScopeVec().end() );
  }

#ifdef DEBUG
  cout << "   Joint scope: " << scope << endl;
#endif

  // remove elimVar from the new scope
  scope.erase(m_bucketVar);

  // make the variables in cond and scope be mutually exclusive
  for (cit = cond.begin(); cit != cond.end();) {
      map<int,val_t>::iterator toRemove = cit;
      if (scope.find(cit->first) == scope.end()) {
          cit++;
          cond.erase(toRemove);
      }
      else {
          scope.erase(cit->first);
          cit++;
      }
  }
  int n = scope.size(); // new function arity

#ifdef DEBUG
  cout << "   Conditioning: " << cond << endl;
  cout << "   Target scope: " << scope << endl;
#endif

  // compute new table size and collect domain sizes
  vector<val_t> domains;
  domains.reserve(n);
  size_t tablesize = 1;
  for (sit = scope.begin(); sit!=scope.end(); ++sit) {
//    if (*it != elimVar)
    tablesize *= m_problem->getDomainSize(*sit);
    domains.push_back(m_problem->getDomainSize(*sit));
  }

  double* newTable = NULL;
  if (buildTable) {
    newTable = new double[tablesize];
    for (j=0; j<tablesize; ++j) newTable[j] = ELEM_ZERO;

    // this keeps track of the tuple assignment
    val_t* tuple = new val_t[n+1+cond.size()];
    for (i=0; i<n; ++i) tuple[i] = 0; // i trough n index target variables
    val_t* elimVal = &tuple[n]; // n+1 is elimVar
    for (cit = cond.begin(), i=n+1; cit != cond.end(); ++cit, ++i) tuple[i] = cit->second; 

    // maps each function scope assignment to the full tuple
    vector<vector<val_t*> > idxMap(m_functions.size() + 1);

/*
    // maps the maxMarginal scope to the full tuple
    vector<val_t*> mmMap;
    */

    // holds iterators .begin() and .end() for all function scopes
    vector< pair< set<int>::const_iterator , set<int>::const_iterator > > iterators;
    iterators.reserve(m_functions.size() + 1);
    for (fit = m_functions.begin(); fit != m_functions.end(); ++fit) {
      // store begin() and end() for each function scope
      iterators.push_back( make_pair( (*fit)->getScopeSet().begin(), (*fit)->getScopeSet().end() ) );
    }
    iterators.push_back( make_pair( maxMarginal->getScopeSet().begin(),
        maxMarginal->getScopeSet().end() ) );

#ifdef DEBUG
    cout << "Max-marginal scope: " << maxMarginal->getScopeSet() << endl;
    cout << "Avg Max-marginal scope: " << avgMaxMarginal->getScopeSet() << endl;
#endif

    // collect pointers to tuple values
    bool bucketVarPassed = false;
    cit = cond.begin();
    int cVarCount = 0;
    for (i=0, sit=scope.begin(); i<n; ++i, ++sit) {
      // process conditioning variables as long as they come before the scope var or elim var
      for (; cit != cond.end() && *sit > cit->first && (bucketVarPassed || m_bucketVar > cit->first); 
              ++cit, ++cVarCount) {
            for (j=0; j<m_functions.size(); ++j) {
                if (m_functions[j]->getScopeSet().find(cit->first) != iterators[j].second) {
#ifdef DEBUG
                    cout << j << " : inserted (cond1): " << cit->first << " [" << m_problem->getDomainSize(cit->first) << "]" <<endl;
#endif
                    idxMap[j].push_back(&tuple[n+1+cVarCount]);
                    ++(iterators[j].first);
                }
            }
            if (maxMarginal->getScopeSet().find(cit->first) != iterators[j].second) {
#ifdef DEBUG
                cout << "inserted (cond)" << cit->first << endl;
#endif
                idxMap[j].push_back(&tuple[n+1+cVarCount]);
                ++(iterators[j].first);
            }
      }
      if (!bucketVarPassed && *sit > m_bucketVar) { // just went past bucketVar
        for (j=0; j<m_functions.size() + 1; ++j) {
#ifdef DEBUG
            cout << j << " : inserted : " << m_bucketVar << endl;
#endif
          idxMap[j].push_back(elimVal);
          ++(iterators[j].first); // skip bucketVar in the original function scope
        }
        bucketVarPassed = true;
      }
      for (; cit != cond.end() && *sit > cit->first && (bucketVarPassed || m_bucketVar > cit->first); 
              ++cit, ++cVarCount) {
            for (j=0; j<m_functions.size(); ++j) {
                if (m_functions[j]->getScopeSet().find(cit->first) != iterators[j].second) {
#ifdef DEBUG
                    cout << j << " : inserted (cond2): " << cit->first << " [" << m_problem->getDomainSize(cit->first) << "]" <<endl;
#endif
                    idxMap[j].push_back(&tuple[n+1+cVarCount]);
                    ++(iterators[j].first);
                }
            }
            if (maxMarginal->getScopeSet().find(cit->first) != iterators[j].second) {
#ifdef DEBUG
                cout << "inserted (cond)" << cit->first << endl;
#endif
                idxMap[j].push_back(&tuple[n+1+cVarCount]);
                ++(iterators[j].first);
            }
      }
      for (j=0; j<m_functions.size() + 1; ++j) {
          /*
          if(j<m_functions.size())
              cout << "  f" << m_functions[j]->getId() << ' ' << *sit << " == " << *(iterators[j].first) << endl;
              */
        if (iterators[j].first != iterators[j].second // scope iterator != end()
            && *sit == *(iterators[j].first)) { // value found
#ifdef DEBUG
                cout << j << " : inserted : " << *sit << endl;
#endif
                idxMap[j].push_back(&tuple[i]);
                ++(iterators[j].first);
        }
      }
    }

    // while bucket var or conditioning remains, continue processing
    while (!bucketVarPassed || cit != cond.end()) { 
        for (; cit != cond.end() && (bucketVarPassed || m_bucketVar > cit->first); ++cit, ++cVarCount) {
            for (j=0; j<m_functions.size(); ++j) {
                if (m_functions[j]->getScopeSet().find(cit->first) != iterators[j].second) {
#ifdef DEBUG
                    cout << j << " : inserted (cond3): " << cit->first << " [" << m_problem->getDomainSize(cit->first) << "]" <<endl;
#endif
                    idxMap[j].push_back(&tuple[n+1+cVarCount]);
                }
            }
            if (maxMarginal->getScopeSet().find(cit->first) != iterators[j].second) {
                idxMap[j].push_back(&tuple[n+1+cVarCount]);
            }

        }
        if(!bucketVarPassed) {
            for (j=0; j<m_functions.size() + 1; ++j) {
#ifdef DEBUG
                cout << j << " : inserted : " << m_bucketVar << endl;
#endif
                idxMap[j].push_back(elimVal);
                bucketVarPassed = true;
            }
        }
    }

    // idxmap for max-marginal functions
    vector<val_t*> &mmMap = idxMap[m_functions.size()];

    
    size_t idx; double z; 
//    val_t valMax;;
    // iterate over all values of elimVar
    for (*elimVal = 0; *elimVal < m_problem->getDomainSize(m_bucketVar); ++(*elimVal) ) {
      idx=0; // go over the full new table
      do {
        z = ELEM_ONE;
        for (j=0; j<m_functions.size(); ++j) {
    		z OP_TIMESEQ m_functions[j]->getValuePtr(idxMap[j]);
        }
//        z = OP_MINUS(z, maxMarginal->getValuePtr(mmMap));
#ifdef DEBUG
        cout << "avg max-marginal: " << ELEM_DECODE(avgMaxMarginal->getValuePtr(mmMap)) << endl;
#endif
        z OP_TIMESEQ avgMaxMarginal->getValuePtr(mmMap) - maxMarginal->getValuePtr(mmMap);
        /*
        z = ELEM_ENCODE(ELEM_DECODE(z) - 
                ELEM_DECODE(maxMarginal->getValuePtr(mmMap)) +
                ELEM_DECODE(avgMaxMarginal->getValuePtr(mmMap)));
        */
#ifdef DEBUG
        cout << "Final adjusted value: " << ELEM_DECODE(z) << endl;
#endif
        newTable[idx] = max(newTable[idx],z);
//        if (z > newTable[idx]) valMax = *elimVal;
      } while ( increaseTuple(idx,tuple,domains) );

    }
    // DEBUG
//    cout << "Maximum value: " << int(valMax) << endl;
    // DEBUG
    // clean up
    delete[] tuple;
  }

  return new FunctionBayes(-m_bucketVar,m_problem,scope,newTable,tablesize);
}


