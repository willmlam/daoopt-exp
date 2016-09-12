/*
 * Problem.cpp
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
 *  Created on: Oct 10, 2008
 *      Author: Lars Otten <lotten@ics.uci.edu>
 */

#include "Problem.h"

#include <iomanip>
#include <iostream>
#include <sstream>
#include <fstream>

#include "UAI2012.h"

#undef DEBUG

namespace daoopt {

struct _membuf : std::streambuf {
  _membuf(char* begin, char* end) {
    this->setg(begin, begin, end);
  }
};

//extern time_t time_start;
//
extern string out_bound_file;

void Problem::condition(const map<int,val_t> &cond) {
  assert(m_n!=UNKNOWN);

  vector<Function*> new_funs;

  m_globalConstant = ELEM_ONE;

  /*
  cout << "Conditioning: " << cond << endl;
  cout << "Before: " << endl;
  for (unsigned i = 0; i < m_functions.size(); ++i) {
    cout << *(m_functions[i]) << endl;
  }
  */
  vector<Function*>::iterator fi = m_functions.begin();
  for(; fi != m_functions.end(); ++fi) {
      Function *fn = (*fi);
      // Substitute evidence variables
      Function *new_fn = fn->substitute(cond);
      if (new_fn->isConstant()) {
        m_globalConstant OP_TIMESEQ new_fn->getTable()[0];
        if (fn->getArity() > 0)
        delete new_fn;
      } else {
        new_funs.push_back(new_fn);
      }
      delete fn;
  }

  double *table1 = new double[1];
  table1[0] = m_globalConstant;
  Function *constFun = new FunctionBayes(new_funs.back()->getId()+1, this, set<int>(), table1, 1);
  new_funs.push_back(constFun);
  m_functions = new_funs;

  /*
  cout << "After: " << endl;
  for (unsigned i = 0; i < m_functions.size(); ++i) {
    cout << i << "\t" << *(m_functions[i]) << endl;
  }
  */

  m_c = m_functions.size();

}

void Problem::removeEvidence(bool clearEvid) {

  assert(m_n!=UNKNOWN);

  // record original no. of variables
  m_nOrg = m_n;

  // Declare aux. variables.
  int idx, i, new_r, new_n;
  val_t k, new_k;
  //map<unsigned int, unsigned int> evidence(m_evidence);
  vector<val_t> new_domains;
  vector<Function*> new_funs;

  // eliminateVar[i]==TRUE iff var. is to be eliminated
  vector<bool> eliminateVar(m_n,false);
  for (map<int,val_t>::iterator it = m_evidence.begin(); it != m_evidence.end(); ++it) {
    eliminateVar[it->first] = true;
  }

  // Identify and tag unary domain variables
  for (i = 0; i < m_n; ++i) {
    if (m_domains.at(i) == 1) { // regard unary domains as evidence
      m_evidence.insert(make_pair(i, 0));
      ++m_e;
      eliminateVar.at(i) = true;
    }
  }

  // Identify variables not covered by any function
  vector<bool> covered(m_n, false);
  //cout << m_functions.size() << endl;
  BOOST_FOREACH(Function * f, m_functions) {
    BOOST_FOREACH(int i, f->getScopeVec()) {
      covered.at(i) = true;
    }
  }
  for (i=0; i<m_n; ++i) {
    if (!covered.at(i)) eliminateVar.at(i) = true;
  }

  // Project functions to account for evidence
  m_globalConstant = ELEM_ONE;
  new_r = 0; // max. arity
  vector<Function*>::iterator fi = m_functions.begin();
  for (; fi != m_functions.end(); ++fi) {
    Function *fn = (*fi);
    // Substitute evidence variables
    Function* new_fn = fn->substitute(m_evidence);
    if (new_fn->isConstant()) { // now empty scope
      m_globalConstant OP_TIMESEQ new_fn->getTable()[0];
      delete new_fn;
    } else {
      new_funs.push_back(new_fn); // record new function
      new_r = max(new_r, (int) new_fn->getScopeVec().size());
    }
    if (!m_is_copy) delete fn; // delete old function (only if this is not a copy)

  }
  m_functions = new_funs;

  // Create dummy function for global constant. Technically the field value needs
  // to be reset to ELEM_ONE, but we keep it around for informational purposes --
  // it should never get used in actual computations.
  double* table1 = new double[1];
  table1[0] = m_globalConstant;
  cout << "Applied evidence -- global constant: " << m_globalConstant << endl;
  Function* constFun = new FunctionBayes(m_nOrg, this, set<int>(), table1, 1);
  m_functions.push_back(constFun);

  /*
  // === shrink more by eliminating certain vars ===

  // remember which vars are in which function scope
  map<int, set<Function*> > mapF; // TODO


  // first, build primal graph
  Graph G(m_n);
  for (fi=m_functions.begin(); fi!=m_functions.end(); ++fi)
    G.addClique((*fi)->getScope());

  bool repeatLoop = true;
  while (repeatLoop) {

    for (int i=0; i<m_n; ++i) {
      if (eliminateVar.at(i) || !G.hasNode(i)) continue; // skip this node
    }

    break; // TODO

  }
  */

  // eliminate tagged variables and reorder remaining ones
  idx = new_n = new_k = 0;
  for (int i=0; i<m_n; ++i) {
    if (!eliminateVar.at(i)) {

      m_old2new.insert(make_pair(i,idx));

      k = m_domains.at(i);
      new_domains.push_back(k);
      new_k = max(new_k,k);

      ++idx;
      ++new_n;
    }
  }

  // update variable information
  m_domains = new_domains;
  m_n = new_n;
  m_k = new_k;
#ifndef NO_ASSIGNMENT
//  m_curSolution.resize(m_n,UNKNOWN);
#endif


  // translate scopes of the new functions
  for (fi=m_functions.begin(); fi!=m_functions.end(); ++fi)
    (*fi)->translateScope(m_old2new);

  /*
  cout << "Remapped variables:";
  typedef std::map<int,int>::value_type mtype;
  BOOST_FOREACH( mtype t, m_old2new )
    cout << ' ' << t.first << "->" << t.second;
  cout << endl;
  */

  // update function information
  m_c = m_functions.size();

  if (clearEvid) m_evidence.clear();
}

void Problem::collapseFunctions() {
    // Create a map from scopes to a list of function indexs
    map<set<int>,vector<int>> mapping;
    size_t oldC = m_c;

    for (unsigned int i=0; i<m_functions.size(); ++i) {
        const auto &sc = m_functions[i]->getScopeSet();
        auto itM = mapping.find(sc);
        if(itM == mapping.end()) {
            mapping.insert(make_pair(sc,vector<int>(1,i)));
        }
        else {
            itM->second.push_back(i);
        }
    }

    vector<Function*> newFunctions;
    newFunctions.resize(mapping.size());
    int newid = 0;
    for (auto mapElement : mapping) {
        const auto &fids = mapElement.second;
        Function *f = m_functions[fids[0]];
        f->setId(newid);
        for (unsigned int i=1;i<fids.size();++i) {
            Function *g = m_functions[fids[i]];
            assert(f->getTableSize() == g->getTableSize());
            for (unsigned int k=0;k<f->getTableSize();++k) {
                f->getTable()[k] OP_TIMESEQ g->getTable()[k];
            }
            delete g;
        }
        newFunctions[newid++] = f;
    }
    m_functions = newFunctions;
    m_c = m_functions.size();
    cout << "Collapsed from " << oldC << " to " << m_c << " functions." << endl;
}

void Problem::perturbDeterminism(double epsilon) {
    for (Function *f : m_functions) {
        double *t = f->getTable();
        for (unsigned int i=0; i<f->getTableSize(); ++i) {
            if (t[i] == ELEM_ZERO)
                t[i] = ELEM_ENCODE(epsilon);
        }
    }
}


bool Problem::parseOrdering(const vector<int>& input, vector<int>& elim) const {

  bool fullOrdering = false; // default
  if (input.size() == (size_t) m_nOrg) {
    fullOrdering = true;
  }

  int n=0, x=UNKNOWN;
  vector<bool> check(m_n, false);

  for (vector<int>::const_iterator it=input.begin(); it!=input.end(); ++it) {
    if (!fullOrdering) {

      if (*it < 0 || *it >= m_n) {
        cerr << "Problem reading ordering, variable index " << *it << " out of range" << endl;
        exit(1);
      }

      if (check[*it]) {
        cerr << "Problem reading ordering, variable " << *it << " appears more than once." << endl;
        exit(1);
      } else check[*it] = true;

      elim.push_back(*it); ++n;

    } else { // full order, needs filtering

      if (*it < 0 || *it >= m_nOrg) {
        cerr << "Problem reading ordering, variable index " << *it << " out of range" << endl;
        exit(1);
      }

      map<int,int>::const_iterator it2 = m_old2new.find(*it);
      if (it2 != m_old2new.end()) {
        x = it2->second;
        if (check[x]) {
          cerr << "Problem reading ordering, variable " << *it << " appears more than once." << endl;
          exit(1);
        } else check[x] = true;

        elim.push_back(x); ++n;
      } else { /* evidence */ }

    }

  }

  if (n!=m_n) {
    cerr << "Problem reading ordering, number of variables doesn't match." << endl;
    exit(1);
  }

  return true;
}


bool Problem::parseOrdering(const string& file, vector<int>& elim) const {

  assert(m_n!=UNKNOWN);

  ifstream inTemp(file.c_str());
  inTemp.close();

  if (inTemp.fail()) { // file not existent yet
    return false;
  }

  igzstream in(file.c_str());

  // ignore first line if there's a pound '#' sign (comment)
  if (in.peek() == '#') {
    in.ignore(8192,'\n');
  }

  int nIn;
  in >> nIn; // length of ordering
  if (nIn != m_n && nIn != m_nOrg) {
    cerr << "Problem reading ordering, number of variables doesn't match" << endl;
    in.close(); exit(1);
  }

  // read into buffer first
  list<int> buffer;
  int x=UNKNOWN;
  while(nIn-- && !in.eof()) {
    in >> x;
    buffer.push_back(x);
  }

  bool fullOrdering = false; // default
  if (buffer.size() == (size_t) m_nOrg) {
    fullOrdering = true;
  }

  int n=0;
  vector<bool> check(m_n, false);

  for (list<int>::iterator it=buffer.begin(); it!=buffer.end(); ++it) {
    if (!fullOrdering) {

      if (*it < 0 || *it >= m_n) {
        cerr << "Problem reading ordering, variable index " << *it << " out of range" << endl;
        in.close(); exit(1);
      }

      if (check[*it]) {
        cerr << "Problem reading ordering, variable " << *it << " appears more than once." << endl;
        in.close(); exit(1);
      } else check[*it] = true;

      elim.push_back(*it); ++n;

    } else { // full order, needs filtering

      if (*it < 0 || *it >= m_nOrg) {
        cerr << "Problem reading ordering, variable index " << *it << " out of range" << endl;
        in.close(); exit(1);
      }

      map<int,int>::const_iterator it2 = m_old2new.find(*it);
      if (it2 != m_old2new.end()) {
        x = it2->second;
        if (check[x]) {
          cerr << "Problem reading ordering, variable " << *it << " appears more than once." << endl;
          in.close(); exit(1);
        } else check[x] = true;

        elim.push_back(x); ++n;
      } else { /* evidence */ }

    }

  }

  if (n!=m_n) {
    cerr << "Problem reading ordering, number of variables doesn't match." << endl;
    in.close();
    exit(1);
  }

  in.close();
  return true;
}


void Problem::saveOrdering(const string& file, const vector<int>& elim) const {
  assert( (int) elim.size() == m_n);

  if (file.substr(file.size()-3,3) == ".gz") { // write gzipped

    ogzstream out(file.c_str());
    if ( ! out ) {
      cerr << "Error writing ordering to file " << file << endl;
      exit(1);
    }
    out << "# daoopt ordering for " << m_name << endl << elim.size();
    for (vector<int>::const_iterator it=elim.begin(); it!=elim.end(); ++it)
      out << ' ' << *it;
    out << endl;
    out.close();

  } else { // write straight text

    ofstream out(file.c_str());
    if ( ! out ) {
      cerr << "Error writing ordering to file " << file << endl;
      exit(1);
    }
    out << "# daoopt ordering for " << m_name << endl << elim.size();
    for (vector<int>::const_iterator it=elim.begin(); it!=elim.end(); ++it)
      out << ' ' << *it;
    out << endl;
    out.close();

  }

}

// supports comments and "sparseuai" format.
bool Problem::parseUAI16(char* prob, size_t probN, char* evid, size_t evidN,
                         bool collapse) {
  assert(prob && probN > 0);
  assert(!evid || evidN > 0);

  _membuf probBuf(prob, prob+probN);

  istream in(&probBuf);

  cout << "Reading problem with parseUAI16..." << endl;

  bool read_sparse = false;
  vector<int> arity;
  vector<vector<int>> scopes;
  string s;
  int x;
  int y;
  val_t xs;
  unsigned int z;

  string line;
  while(m_task == UNKNOWN && std::getline(in, s)) {
    // Skip comment lines and empty lines
    if (s[0] == 'c' || s.length() == 0) {
      continue;
    }
    if (s == "BAYES") {
      m_task = TASK_MAX;
      m_prob = PROB_MULT;
    } else if (s == "MARKOV") {
      m_task = TASK_MAX;
      m_prob = PROB_MULT;
    } else if (s == "SPARSEBAYES") {
      m_task = TASK_MAX;
      m_prob = PROB_MULT;
      read_sparse = true;
    } else if (s == "SPARSEMARKOV") {
      m_task = TASK_MAX;
      m_prob = PROB_MULT;
      read_sparse = true;
    } else {
      cerr << "Unsupported problem type \"" << s << "\", aborting." << endl;
      return false;
    }
  }

  in >> x; // No. of variables
  m_n = x;
  m_domains.resize(m_n,UNKNOWN);
#ifndef NO_ASSIGNMENT
//  m_curSolution.resize(m_n,UNKNOWN);
#endif
  m_k = -1;
  for (int i=0; i<m_n; ++i) { // Domain sizes
    in >> x; // read into int first
    if (x > numeric_limits<val_t>::max()) {
      cerr << "Domain size " << x << " out of range for internal representation.\n"
           << "(Recompile with different type for variable values.)" << endl;
      return false;
    }
    xs = (val_t)x;
    m_domains[i] = xs;
    m_k = max(m_k,xs);
  }

  in >> x; // No. of functions
  m_c = x;
  scopes.reserve(m_c);

  // Scope information for functions
  m_r = -1;
  for (int i = 0; i < m_c; ++i)
  {
    vector<int> scope;
    in >> x; // arity

    m_r = max(m_r, x);
    for (int j=0; j<x; ++j) {
      in >> y; // the actual variables in the scope
      if(y>=m_n) {
        cerr << "Variable index " << y << " out of range." << endl;
        return false;
      }
      scope.push_back(y); // preserve order from file
    }
    scopes.push_back(scope);
  }

  // Read functions
  for (int i = 0; i < m_c; ++i) {
    size_t tab_size = 1;

    for (vector<int>::iterator it=scopes[i].begin(); it!=scopes[i].end(); ++it) {
      tab_size *= m_domains[*it];
    }
    // create set version of the scope (ordered)
    set<int> scopeSet(scopes[i].begin(), scopes[i].end());
    unsigned int scope_size = scopeSet.size();

    // compute reindexing map from specified scope to ordered, internal one
    map<int,int> mapping;
    int k = 0;
    for (vector<int>::const_iterator it=scopes[i].begin(); it!=scopes[i].end(); ++it) {
      mapping[*it] = k++;
    }
    vector<int> reidx(scope_size);
    vector<int>::iterator itr = reidx.begin();
    for (set<int>::iterator it=scopeSet.begin(); itr!=reidx.end(); ++it, ++itr) {
      *itr = mapping[*it];
    }

    if (!read_sparse) {
      in >> z; // No. of entries
      assert(tab_size==z); // product of domain sizes matches no. of entries

      // read the full table into an temp. array (to allow reordering)
      vector<double> temp(tab_size);
      for (size_t j=0; j<tab_size;) {
        in >> s ;
        // check sparse-factor encoding "(x:n)" where x is real number and n is an int. meaning of this is, next entries in the table equal x.
        if ('(' == s[0] && ')' == s[s.length()-1]) {
          s.erase(0,1) ; s.erase(s.length()-1, 1) ;
          std::string::size_type pos = s.find(':') ;
          if (std::string::npos == pos)
          return false ;
          string sX = s.substr(0, pos) ;
          if (0 == sX.length()) return false ;
          string sN = s.substr(pos+1, s.length() - pos - 1) ;
          if (0 == sN.length()) return false ;
          double x = atof(sX.c_str()) ;
          size_t nEntriesLeft = tab_size - j ;
          size_t nEntries = atoi(sN.c_str()) ;
          if (nEntries <= 0 || nEntries > nEntriesLeft)
          return false ;
          for (int ni = 0 ; ni < nEntries ; ni++)
          temp[j++] = x ;
        }
        // else token is a table entry
        else {
          temp[j++] = atof(s.c_str()) ;
        }
        //      in >> temp[j];
      }

      // get the variable domain sizes
      vector<val_t> limit; limit.reserve(scope_size);
      for (vector<int>::const_iterator it=scopes[i].begin(); it!=scopes[i].end(); ++it)
      limit.push_back(m_domains[*it]);
      vector<val_t> tuple(scope_size, 0);

      // create the new table (with reordering)
      double* table = new double[tab_size];
      for (size_t j=0; j<tab_size; ) {
        size_t pos=0, offset=1;
        // j is the index in the temp. table
        for (k = scope_size - 1; k >= 0; --k) { // k goes backwards through the ordered scope
          pos += tuple[reidx[k]] * offset;
          offset *= m_domains[scopes[i][reidx[k]]];
        }
        table[pos] = ELEM_ENCODE( temp[j] );
        increaseTuple(j,tuple,limit);
      }

      Function* f = new FunctionBayes(i,this,scopeSet,table,tab_size);
      m_functions.push_back(f);
    } else {
      // Read sparse factor
      double value;
      int num_non_default_entries;
      in >> value;
      in >> num_non_default_entries;

      double* table = new double[tab_size];

      // Fill table with default value first
      for (size_t k = 0; k < tab_size; ++k) {
        table[k] = ELEM_ENCODE(value);
      }
      vector<val_t> tuple(scope_size, 0);
      for (int j = 0; j < num_non_default_entries; ++j) {
        for (int k = 0; k < scope_size; ++k) {
          in >> xs;
          tuple[reidx[k]] = xs;
        }
        size_t pos = 0;
        size_t offset = 1;
        for (int k = scope_size - 1; k >= 0; --k) {
          pos += tuple[k] * offset;
          offset *= m_domains[scopes[i][reidx[k]]];
        }
        in >> value;
        table[pos] = ELEM_ENCODE(value);
      }
      Function* f = new FunctionBayes(i, this, scopeSet, table, tab_size);
      m_functions.push_back(f);
    }
  } // All function tables read

#ifdef DEBUG
  for (Function* f : m_functions) {
    cout << *f << endl;
  }
  cout << endl;
#endif


  if (collapse) collapseFunctions();

  // Read evidence?
  if (!evid) {
    m_e = 0;
    cout << "Problem size (MB): " << (getSize()*sizeof(double) / (1024*1024.0)) << endl;
    return true; // No evidence, return
  }

  cout << "Reading evidence..." << endl;

  _membuf evidBuf(evid, evid+evidN);

  istream in2(&evidBuf);

  /*
  in >> x; // Number of evidence samples
  if (x > 1) {
      myerror("Warning: Ignoring all but one evidence sample.\n");
  }
  */

    in2 >> x;
    m_e = x; // Number of evidence variables

    for (int i=0; i<m_e; ++i) {
        in2 >> x; // Variable index
        in2 >> y; // Variable value
        xs = (val_t) y;
        if (xs >= m_domains[x]) {
        cout << "Variable " << x << " has domain size " << (int) m_domains[x]
            << ", evidence value " << y << " out of range." << endl;
        return false;
        }
        m_evidence.insert(make_pair(x,xs));
    }

  // Compute determinism ratio
  num_tuples_ = 0;
  num_zero_tuples_ = 0;
  for (const Function *f : m_functions) {
    num_tuples_ += f->getTableSize();
    for (int k = 0; k < f->getTableSize(); ++k) {
      if (f->getTable()[k] == ELEM_ZERO) {
        ++num_zero_tuples_;
      }
    }
  }
  determinism_ratio_ = (double) num_zero_tuples_ / num_tuples_;

  cout << "Problem size (MB): " << (getSize()*sizeof(double) / (1024*1024.0)) << endl;
  cout << "Number of tuples:  " << num_tuples_ << endl;
  cout << "Number of zeros:   " << num_zero_tuples_ << endl;
  cout << "Determinism ratio: " << determinism_ratio_ << endl;
  return true;
}



//bool Problem::parseUAI(const string& prob, const string& evid, bool collapse) {
bool Problem::parseUAI(char* prob, size_t probN, char* evid, size_t evidN,
                     bool collapse) {
assert(prob && probN > 0);
assert(!evid || evidN > 0);

  _membuf probBuf(prob, prob+probN);

  istream in(&probBuf);

  cout << "Reading problem..." << endl;

  vector<int> arity;
  vector<vector<int> > scopes;
  string s;
  int x,y;
  val_t xs;
  unsigned int z;

  in >> s; // Problem type
  if (s == "BAYES") {
    m_task = TASK_MAX;
    m_prob = PROB_MULT;
  } else if (s == "MARKOV") {
    m_task = TASK_MAX;
    m_prob = PROB_MULT;
  } else {
    cerr << "Unsupported problem type \"" << s << "\", aborting." << endl;
    return false;
  }

  in >> x; // No. of variables
  m_n = x;
  m_domains.resize(m_n,UNKNOWN);
#ifndef NO_ASSIGNMENT
//  m_curSolution.resize(m_n,UNKNOWN);
#endif
  m_k = -1;
  for (int i=0; i<m_n; ++i) { // Domain sizes
    in >> x; // read into int first
    if (x > numeric_limits<val_t>::max()) {
      cerr << "Domain size " << x << " out of range for internal representation.\n"
           << "(Recompile with different type for variable values.)" << endl;
      return false;
    }
    xs = (val_t)x;
    m_domains[i] = xs;
    m_k = max(m_k,xs);
  }

  in >> x; // No. of functions
  m_c = x;
  scopes.reserve(m_c);

  // Scope information for functions
  m_r = -1;
  for (int i = 0; i < m_c; ++i)
  {
    vector<int> scope;
    in >> x; // arity

    m_r = max(m_r, x);
    for (int j=0; j<x; ++j) {
      in >> y; // the actual variables in the scope
      if(y>=m_n) {
        cerr << "Variable index " << y << " out of range." << endl;
        return false;
      }
      scope.push_back(y); // preserve order from file
    }
    scopes.push_back(scope);
  }

  // Read functions
  for (int i = 0; i < m_c; ++i)
  {
    in >> z; // No. of entries
    size_t tab_size = 1;

    for (vector<int>::iterator it=scopes[i].begin(); it!=scopes[i].end(); ++it) {
      tab_size *= m_domains[*it];
    }

    assert(tab_size==z); // product of domain sizes matches no. of entries

    // create set version of the scope (ordered)
    set<int> scopeSet(scopes[i].begin(), scopes[i].end());
    z = scopeSet.size();

    // compute reindexing map from specified scope to ordered, internal one
    map<int,int> mapping; int k=0;
    for (vector<int>::const_iterator it=scopes[i].begin(); it!=scopes[i].end(); ++it)
      mapping[*it] = k++;
    vector<int> reidx(z);
    vector<int>::iterator itr = reidx.begin();
    for (set<int>::iterator it=scopeSet.begin(); itr!=reidx.end(); ++it, ++itr) {
      *itr = mapping[*it];
    }

    // read the full table into an temp. array (to allow reordering)
    vector<double> temp(tab_size);
    for (size_t j=0; j<tab_size;) {
		in >> s ;
		// check sparse-factor encoding "(x:n)" where x is real number and n is an int. meaning of this is, next entries in the table equal x.
		if ('(' == s[0] && ')' == s[s.length()-1]) {
			s.erase(0,1) ; s.erase(s.length()-1, 1) ;
			std::string::size_type pos = s.find(':') ;
			if (std::string::npos == pos)
				return false ;
			string sX = s.substr(0, pos) ;
			if (0 == sX.length()) return false ;
			string sN = s.substr(pos+1, s.length() - pos - 1) ;
			if (0 == sN.length()) return false ;
			double x = atof(sX.c_str()) ;
			size_t nEntriesLeft = tab_size - j ;
			size_t nEntries = atoi(sN.c_str()) ;
			if (nEntries <= 0 || nEntries > nEntriesLeft)
				return false ;
			for (int ni = 0 ; ni < nEntries ; ni++)
				temp[j++] = x ;
			}
		// else token is a table entry
		else {
			temp[j++] = atof(s.c_str()) ;
			}
//      in >> temp[j];
    }

    // get the variable domain sizes
    vector<val_t> limit; limit.reserve(z);
    for (vector<int>::const_iterator it=scopes[i].begin(); it!=scopes[i].end(); ++it)
      limit.push_back(m_domains[*it]);
    vector<val_t> tuple(z, 0);

    // create the new table (with reordering)
    double* table = new double[tab_size];
    for (size_t j=0; j<tab_size; ) {
      size_t pos=0, offset=1;
      // j is the index in the temp. table
      for (k=z-1; k>=0; --k) { // k goes backwards through the ordered scope
        pos += tuple[reidx[k]] * offset;
        offset *= m_domains[scopes[i][reidx[k]]];
      }
      table[pos] = ELEM_ENCODE( temp[j] );
      increaseTuple(j,tuple,limit);
    }

    Function* f = new FunctionBayes(i,this,scopeSet,table,tab_size);
    m_functions.push_back(f);

  } // All function tables read


  if (collapse) collapseFunctions();

  // Read evidence?
  if (!evid) {
    m_e = 0;
    cout << "Problem size (MB): " << (getSize()*sizeof(double) / (1024*1024.0)) << endl;
    return true; // No evidence, return
  }

  cout << "Reading evidence..." << endl;

  _membuf evidBuf(evid, evid+evidN);

  istream in2(&evidBuf);

  /*
  in >> x; // Number of evidence samples
  if (x > 1) {
      myerror("Warning: Ignoring all but one evidence sample.\n");
  }
  */

    in2 >> x;
    m_e = x; // Number of evidence variables

    for (int i=0; i<m_e; ++i) {
        in2 >> x; // Variable index
        in2 >> y; // Variable value
        xs = (val_t) y;
        if (xs >= m_domains[x]) {
        cout << "Variable " << x << " has domain size " << (int) m_domains[x]
            << ", evidence value " << y << " out of range." << endl;
        return false;
        }
        m_evidence.insert(make_pair(x,xs));
    }


  cout << "Problem size (MB): " << (getSize()*sizeof(double) / (1024*1024.0)) << endl;
  return true;
}


void Problem::outputAndSaveSolution(const string& file, const SearchStats* nodestats,
    const vector<count_t>& nodeProf, const vector<count_t>& leafProf, bool toScreen) const {

  bool writeFile = false;
  if (! file.empty())
    writeFile = true;

  ogzstream out;
  if (writeFile) {
    out.open(file.c_str(), ios::out | ios::trunc | ios::binary);
    if (!out) {
      cerr << "Error writing optimal solution to file " << file << endl;
      writeFile = false;
      out.close();
    }
  }

  oss screen;
  screen << std::setprecision(20);
  screen << "s " << SCALE_LOG(m_curCost);
#ifndef NO_ASSIGNMENT
  int32_t assigSize = UNKNOWN;
  if (m_subprobOnly)
    assigSize = (int32_t) m_curSolution.size();  // no dummy variable included
  else
    assigSize = (int32_t) m_nOrg;
  screen << ' ' << assigSize;
#endif

  if (writeFile) {
    BINWRITE(out, m_curCost); // mpe solution cost
    count_t countOR = 0, countAND = 0;
    if (nodestats) {
      countOR = nodestats->numExpOR;
      countAND = nodestats->numExpAND;
    }
    BINWRITE(out, countOR);
    BINWRITE(out, countAND);
  }

#ifndef NO_ASSIGNMENT
  if (writeFile) {
    BINWRITE(out,assigSize); // no. of variables in opt. assignment
  }

  // generate full assignment (incl. evidence) if necessary and output
  vector<val_t> outputAssig;
  assignmentForOutput(outputAssig);
  BOOST_FOREACH( int32_t v, outputAssig ) {
    screen << ' ' << v;
    if (writeFile) BINWRITE(out, v);
  }
#endif

  // output node profiles in case of subproblem processing
  if (m_subprobOnly) {
    int32_t size = (int32_t) leafProf.size();
    BINWRITE(out, size);
    // leaf nodes first
    for (vector<count_t>::const_iterator it=leafProf.begin(); it!=leafProf.end(); ++it) {
      BINWRITE(out,*it);
    }
    // now full node profile (has same array size)
    for (vector<count_t>::const_iterator it=nodeProf.begin(); it!=nodeProf.end(); ++it) {
      BINWRITE(out,*it);
    }
  }

  screen << endl;
  if (toScreen)
    cout << screen.str();
  if (writeFile)
    out.close();
}


#ifndef NO_ASSIGNMENT
void Problem::assignmentForOutput(vector<val_t>& assg) const {
  assignmentForOutput(m_curSolution, assg);
}

void Problem::assignmentForOutput(const vector<val_t>& inAssg, vector<val_t>& outAssg) const {
  if (m_subprobOnly || inAssg.empty()) {
    outAssg = inAssg;
    // update: no need to remove dummy anymore
  } else {
    outAssg.resize(m_nOrg, UNKNOWN);
    for (int i=0; i<m_nOrg; ++i) {
      map<int,int>::const_iterator itRen = m_old2new.find(i);
      if (itRen != m_old2new.end()) {  // var part of solution
        outAssg.at(i) = inAssg.at(itRen->second);
      } else {
        map<int,val_t>::const_iterator itEvid = m_evidence.find(i);
        if (itEvid != m_evidence.end())  // var part of evidence
          outAssg.at(i) = itEvid->second;
        else  // var had unary domain
          outAssg.at(i) = 0;
      }
    }
  }
}
#endif


void Problem::updateSolution(double cost,
#ifndef NO_ASSIGNMENT
    const vector<val_t>& sol,
#endif
    const SearchStats* nodestats,
    bool output) {

  if (ISNAN(cost))
    return;

  double costCheck = ELEM_ZERO;
#ifndef NO_ASSIGNMENT
  // check for complete assignment first
  for (size_t i = 0; i < sol.size(); ++i) {
    if (sol[i] == NONE) {
      oss ss;
      ss << std::setprecision(20);
      ss << "Warning: skipping incomplete solution, reported " << cost;
      DIAG(ss << " " << sol.size() << " " << sol;)
      ss << endl; myprint(ss.str());
      return;
    }
    if (!m_subprobOnly && sol[i] >= m_domains[i]) {
      oss ss; ss << "Warning: value " << (int)sol[i] << " outside of variable " << i
                 << " domain " << (int)m_domains[i];
      ss << endl; myprint(ss.str());
      return;
    }
  }

  // use Kahan summation to compute exact solution cost
  // TODO (might not work in non-log scale?)
  if (cost != ELEM_ZERO && !m_subprobOnly) {
    costCheck = ELEM_ONE; double comp = ELEM_ONE;  // used across loop iterations
    double y, z;  // reset for each loop iteration
    for (Function* f : m_functions) {
      z = f->getValue(sol);
      y = z OP_DIVIDE comp;
      z = costCheck OP_TIMES y;
      comp = (z OP_DIVIDE costCheck) OP_DIVIDE y;
      costCheck = z;
    }
    if (1e-3 < abs(cost-costCheck)) {
      oss ss;
      ss << std::setprecision(20);
      ss << "Warning: solution cost " << costCheck << " differs significantly"
        << ", reported " << cost << endl;
      myprint(ss.str());
    }
  } else
#endif
  costCheck = cost;

//  if (ISNAN(costCheck) || (!ISNAN(m_curCost) && costCheck <= m_curCost)) { // TODO costCheck =?= ELEM_ZERO )
  if (ISNAN(costCheck) ||
      (!ISNAN(m_curCost) && m_curCost - costCheck >= 1.11e-16)) {
    oss ss;
    ss << std::setprecision(20);
    ss << "Warning: Discarding solution with cost " << costCheck
      << ", reported: " << cost;
#ifndef NO_ASSIGNMENT
    DIAG(ss << " " << sol.size() << " " << sol;)
    vector<val_t> outputAssg;
    assignmentForOutput(sol, outputAssg);
    ss << ' ' << outputAssg.size();
    BOOST_FOREACH( int v, outputAssg ) {
      ss << ' ' << v;
    }
#endif
    ss << endl; myprint(ss.str());
    return;
  }
  m_curCost = costCheck;
  if (costCheck == ELEM_ZERO) output = false;
  ostringstream ss;
  ss << std::setprecision(20);
  if (output) {
    ss << "u ";
    if (nodestats)
      ss << nodestats->numExpOR << ' ' <<  nodestats->numExpAND << ' ';
    else
      ss << "0 0 ";
    ss << SCALE_LOG(costCheck) ;
  }

#ifndef NO_ASSIGNMENT
  // save only the reduced solution
  // NOTE: sol.size() < m_nOrg in conditioned subproblem case
  m_curSolution = sol;
  // output the complete assignment (incl. evidence)
  if (output) {
    vector<val_t> outputAssg;
    assignmentForOutput(outputAssg);
    ss << ' ' << outputAssg.size();
    BOOST_FOREACH( int v, outputAssg ) {
      ss << ' ' << v;
    }
    /*
    UAI2012::outputSolutionValT(outputAssg);
    */
    if (!out_bound_file.empty()) {
      vector<count_t> dummy_node_profile;
      vector<count_t> dummy_leaf_profile;
      outputAndSaveSolution(out_bound_file, nodestats,
                            dummy_node_profile, dummy_leaf_profile, false);
    }

  }
#endif

  if (output) {
    ss << endl;
    myprint(ss.str());
  }
}


void Problem::resetSolution() {
  m_curCost = ELEM_NAN;
#ifndef NO_ASSIGNMENT
  m_curSolution.clear();
#endif
}

void Problem::updateUpperBound(double bound, const SearchStats* nodestats,
    bool output) {
  if (std::isnan(bound)) {
    return;
  }
  if (bound < m_curUpperBound || std::isnan(m_curUpperBound)) {
    m_curUpperBound = bound;
    if (output) {
      oss ss;
      ss << std::setprecision(20);
      ss << "h ";
      if (nodestats)
        ss << nodestats->numExpOR << ' ' <<  nodestats->numExpAND << ' ';
      else
        ss << "0 0 ";
      ss << SCALE_LOG(m_curUpperBound);
      ss << endl;
      myprint(ss.str());
    }
  }
}


void Problem::writeUAI(const string& prob) const {
  assert (prob.size());

  ogzstream out;
  out.open(prob.c_str(), ios::out | ios::trunc);

  if (!out) {
    cerr << "Error writing reduced network to file " << prob << endl;
    exit(1);
  }

  out << "MARKOV" << endl; // TODO hard-coding is not optimal

  // variable info
  out << (m_n) << endl;
  for (vector<val_t>::const_iterator it=m_domains.begin(); it!=m_domains.end(); ++it)
    out << ' ' << ((int) *it) ;

  // function information
  out << endl << m_functions.size() << endl;
  for (vector<Function*>::const_iterator it=m_functions.begin(); it!=m_functions.end(); ++it) {
    const vector<int>& scope = (*it)->getScopeVec();
    out << scope.size() << '\t'; // scope size
    for (vector<int>::const_iterator itS=scope.begin(); itS!=scope.end(); ++itS)
      out << *itS << ' '; // variables in scope
    out << endl;
  }
  out << endl;

  // write the function tables
  for (vector<Function*>::const_iterator it=m_functions.begin(); it!=m_functions.end(); ++it) {
    double * T = (*it)->getTable();
    out << (*it)->getTableSize() << endl; // table size
    for (size_t i=0; i<(*it)->getTableSize(); ++i)
      out << ' ' << SCALE_NORM( T[i] ); // table entries
    out << endl;
  }

  // done
  out << endl;
  out.close();

}


void Problem::addDummy() {
  m_n += 1;
  m_hasDummy = true;
  m_domains.push_back(1); // unary domain
}


void Problem::replaceFunctions(const vector<Function*>& newFunctions, bool asCopy) {
  // delete current functions
  for (vector<Function*>::iterator it = m_functions.begin();
          !m_is_copy && it!= m_functions.end(); ++it) {
    if (*it) delete (*it);
  }
  // store new functions
  if (!asCopy) {
      m_functions = newFunctions;
  }
  else {
      m_functions.clear();
      for (Function *f : newFunctions) {
          m_functions.push_back(f->clone());
      }
  }
  m_c = m_functions.size();
  // update function scopes???
}

/*
void Problem::addEvidence(const map<int,val_t> &evid) {
  map<int,val_t>::const_iterator eit = evid.begin();
  for (; eit != evid.end(); ++eit) {
    int x = eit->first;
    val_t v = eit->second;
    cout << int(v) << " " << int(m_domains[x]) << endl;
    cout << (v >= m_domains[x]) << endl;
    cout << (m_evidence.find(x) != m_evidence.end()) << endl << endl;
    if (v >= m_domains[x] || m_evidence.find(x) != m_evidence.end()) {
        cerr << "WARNING: bad evidence in input, skipping" << endl;
        continue;
    }
    m_e++;
    m_evidence.insert(make_pair(x,v));
  }
}
*/


#ifndef NO_ASSIGNMENT
bool Problem::isEliminated(int i) const {
  map<int,int>::const_iterator itRen = m_old2new.find(i);
  return itRen == m_old2new.end();
}
#endif


size_t Problem::getSize() const {
    size_t S = 0;
    for (const auto &f : m_functions) {
        S += f->getTableSize();
    }
    return S;
}

}  // namespace daoopt
