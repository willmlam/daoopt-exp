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

#include "MiniBucketElim.h"

/* disables DEBUG output */
#undef DEBUG

#ifdef DEBUG

namespace daoopt {
/* ostream operator for debugging */
ostream& operator <<(ostream& os, const vector<Function*>& l) {
  vector<Function*>::const_iterator it = l.begin();
  os << '[';
  while (it!=l.end()) {
    os << (**it);
    if (++it != l.end()) os << ',';
  }
  os << ']';
  return os;
}

}  // namespace daoopt
#endif

namespace daoopt {

/* computes the augmented part of the heuristic estimate */
double MiniBucketElim::getHeur(int var, vector<val_t>& assignment, SearchNode* n)
{
	assert( var >= 0 && var < m_problem->getN());
  ++var_heur_calls_[var];
	double h = ELEM_ONE;

	// go over augmented and intermediate lists and combine all values
	vector<Function*>::const_iterator itF = m_augmented[var].begin();
	for (; itF!=m_augmented[var].end(); ++itF) {
		h OP_TIMESEQ (*itF)->getValue(assignment);
  }

	itF = m_intermediate[var].begin();
	for (; itF!=m_intermediate[var].end(); ++itF) {
		h OP_TIMESEQ (*itF)->getValue(assignment);
  }

	return h;
}

double MiniBucketElim::getHeurPerIndSubproblem(int var, std::vector<val_t> & assignment, SearchNode* node, double label, std::vector<double> & subprobH)
{
	/* note : 
		union of augmented/intermediate functions in this bucket = 
			union of intermediate functions of each child bucket + union of output functions of minibuckets of each child bucket.
		this allows breacking the h into components, one from each child bucket.
	*/
	double h = ELEM_ONE ;
	const PseudotreeNode *n = m_pseudotree->getNode(var) ;
	const vector<PseudotreeNode *> & children = n->getChildren() ;
	subprobH.resize(children.size(), ELEM_ONE) ;
	vector<PseudotreeNode*>::const_iterator itCend = children.end() ;
	int i = 0 ;
	for (vector<PseudotreeNode*>::const_iterator itC = children.begin() ; itC != itCend ; ++itC, ++i) {
		int child = (*itC)->getVar() ;
		double & hSubproblem = subprobH[i] ; hSubproblem = ELEM_ONE ;
		// iterate over subproblem intermediate functions; add to subproblem h
		vector<Function*>::const_iterator itF = m_intermediate[child].begin(), itFend = m_intermediate[child].end() ;
		for (; itF != itFend ; ++itF) 
			hSubproblem OP_TIMESEQ (*itF)->getValue(assignment) ;
		// iterate over subproblem MB output functions; add to subproblem h
		vector<MiniBucket> & minibuckets = _MiniBuckets[child] ;
		for (const MiniBucket& mb : minibuckets) {
			Function *fMB = mb.output_fn() ;
			if (NULL == fMB) continue ;
			hSubproblem OP_TIMESEQ fMB->getValue(assignment) ;
			}
		h OP_TIMESEQ hSubproblem ;
		}
	return h ;
}


void MiniBucketElim::getHeurAll(int var, vector<val_t>& assignment, SearchNode* n, vector<double>& out) {
  out.clear();
  out.resize(m_problem->getDomainSize(var), ELEM_ONE);
  vector<double> funVals;
  vector<Function*>::const_iterator itF;
  for (itF = m_augmented[var].begin(); itF!=m_augmented[var].end(); ++itF) {
    (*itF)->getValues(assignment, var, funVals);
    for (size_t i=0; i<out.size(); ++i)
      out[i] OP_TIMESEQ funVals[i];
  }
  for (itF = m_intermediate[var].begin(); itF!=m_intermediate[var].end(); ++itF) {
    (*itF)->getValues(assignment, var, funVals);
    for (size_t i=0; i<out.size(); ++i)
      out[i] OP_TIMESEQ funVals[i];
  }
}

// In this case, just the functions indexed by the pseudotree
double MiniBucketElim::getLabel(int var, const vector<val_t> &assignment, SearchNode *node) {
    double d = ELEM_ONE;
    for (Function *f : m_pseudotree->getFunctions(var)) {
        d OP_TIMESEQ f->getValue(assignment);
    }
    return d;
}

void MiniBucketElim::getLabelAll(int var, const vector<val_t> &assignment, SearchNode *node, vector<double> &out) {
    vector<double> costTmp(m_problem->getDomainSize(var), ELEM_ONE);
    for (Function *f : m_pseudotree->getFunctions(var)) {
        f->getValues(assignment, var, costTmp);
        for (int i=0; i<m_problem->getDomainSize(var); ++i) {
            out[i] OP_TIMESEQ costTmp[i];
        }
    }
}

// For the baseline minibucket heuristic, the ordering heuristic is identical
// to the original heuristic. This is always used after getHeur().
double MiniBucketElim::getOrderingHeur(int var, std::vector<val_t>& assignment,
    SearchNode* node) {
  return node->getHeur();
}


void MiniBucketElim::reset() {

  for (vector<MiniBucket> & mb : _MiniBuckets) 
	  mb.clear() ;
  _MiniBuckets.clear() ;

//  vector<vector<Function*> > empty;
//  m_augmented.swap(empty);
//  vector<vector<Function*> > empty2;
//  m_intermediate.swap(empty2);
	for (vector<vector<Function*> >::iterator itA = m_augmented.begin(); itA!=m_augmented.end(); ++itA)
		for (vector<Function*>::iterator itB = itA->begin(); itB!=itA->end(); ++itB)
			delete *itB;
	m_augmented.clear();
	m_intermediate.clear();
}

size_t MiniBucketElim::build(int task, const vector<val_t> * assignment, bool computeTables){
	return build(assignment, computeTables);
}

size_t MiniBucketElim::build(const vector<val_t> * assignment, bool computeTables) {

#ifdef DEBUG
  cout << "$ Building MBE(" << m_ibound << ")" << endl;
#endif

  this->reset();
  if (computeTables) {
    LPReparameterization();
  }

  vector<int> elimOrder; // will hold dfs order
  findDfsOrder(elimOrder); // computes dfs ordering of relevant subtree

  m_augmented.resize(m_problem->getN());
  m_intermediate.resize(m_problem->getN());
  _MiniBuckets.resize(m_problem->getN());

  // keep track of total memory consumption
  size_t memSize = 0;

  // ITERATES OVER BUCKETS, FROM LEAVES TO ROOT
  for (vector<int>::reverse_iterator itV=elimOrder.rbegin(); itV!=elimOrder.rend(); ++itV) {
	  int v = *itV;  // this is the variable being eliminated
	  // partition functions into minibuckets
	  vector<MiniBucket>& minibuckets = _MiniBuckets[v];
	  minibuckets.clear();

#ifdef DEBUG
    cout << "$ Bucket for variable " << *itV << endl;
#endif

    // collect relevant functions in funs
    vector<Function*> funs;
    const vector<Function*>& fnlist = m_pseudotree->getFunctions(*itV);
    funs.insert(funs.end(), fnlist.begin(), fnlist.end());
    funs.insert(funs.end(), m_augmented[*itV].begin(), m_augmented[*itV].end());
#ifdef DEBUG
    for (vector<Function*>::iterator itF=funs.begin(); itF!=funs.end(); ++itF)
      cout << ' ' << (**itF);
    cout << endl;
#endif

    // compute global upper bound for root (dummy) bucket
    if (*itV == elimOrder[0]) {// variable is dummy root variable
      if (computeTables && assignment) { // compute upper bound if assignment is given
        m_globalUB = ELEM_ONE;
        for (vector<Function*>::iterator itF=funs.begin(); itF!=funs.end(); ++itF)
          m_globalUB OP_TIMESEQ (*itF)->getValue(*assignment);
        cout << "    MBE-ALL  = " << SCALE_LOG(m_globalUB) << " (" << SCALE_NORM(m_globalUB) << ")" << endl;
        m_globalUB OP_DIVIDEEQ m_problem->globalConstInfo();  // for backwards compatibility of output
        cout << "    MBE-ROOT = " << SCALE_LOG(m_globalUB) << " (" << SCALE_NORM(m_globalUB) << ")" << endl;
      }
      continue; // skip the dummy variable's bucket
    }

    // sort functions by decreasing scope size
    sort(funs.begin(), funs.end(), scopeIsLarger);

    // partition functions into minibuckets
//    vector<Function*>::iterator itF; bool placed;
    for (vector<Function*>::iterator itF = funs.begin(); itF!=funs.end();
         ++itF) {
      bool placed = false;
      for (vector<MiniBucket>::iterator itB=minibuckets.begin();
            !placed && itB!=minibuckets.end(); ++itB)
      {
        if (itB->allowsFunction(*itF)) { // checks if function fits into bucket
          itB->addFunction(*itF);
          placed = true;
        }
      }
      if (!placed) { // no fit, need to create new bucket
        MiniBucket mb(*itV,m_ibound,m_problem);
        mb.addFunction(*itF);
        minibuckets.push_back(mb);
      }
    }


    // Moment-matching step, performed only if we have partitioning.
    vector<Function*> max_marginals;
    std::unique_ptr<Function> average_mm_function;
    if (computeTables && m_options->match && minibuckets.size() > 1) {

      set<int> scope_intersection;
      bool first_mini_bucket = true;
      for (const MiniBucket& mini_bucket : minibuckets) {
        if (first_mini_bucket) {
          scope_intersection = mini_bucket.getJointScope();
          first_mini_bucket = false;
        } else {
          scope_intersection =
              intersection(scope_intersection, mini_bucket.getJointScope());
        }
      }
      for (MiniBucket& mini_bucket : minibuckets) {
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

      double* average_mm_table = new double[table_size];
      for (size_t i = 0; i < table_size; ++i) {
        average_mm_table[i] = ELEM_ONE;
      }
      for (const Function* max_marginal : max_marginals) {
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

    // minibuckets for current bucket are now ready, process each
    // and place resulting function
    int bucket_idx = 0;
    for (MiniBucket& mini_bucket : minibuckets) {
      Function* new_function;
      if (!computeTables || !m_options->match || minibuckets.size() <= 1) {
        new_function = mini_bucket.eliminate(computeTables); 
      } else {
        new_function = mini_bucket.eliminateMM(computeTables,
                                               max_marginals[bucket_idx++],
                                               average_mm_function.get());
      }

      const set<int>& new_scope = new_function->getScopeSet();
      memSize += new_function->getTableSize();
      // go up in tree to find target bucket
      PseudotreeNode* n = m_pseudotree->getNode(*itV)->getParent();
      while (new_scope.find(n->getVar()) == new_scope.end() &&
             n != m_pseudotree->getRoot()) {
        m_intermediate[n->getVar()].push_back(new_function);
        n = n->getParent();
      }
      // matching bucket found OR root of pseudo tree reached
      m_augmented[n->getVar()].push_back(new_function);
    }
    // all minibuckets processed and resulting functions placed
  }

#ifdef DEBUG
  // output augmented and intermediate buckets
  if (computeTables)
    for (int i=0; i<m_problem->getN(); ++i) {
      cout << "$ AUG" << i << ": " << m_augmented[i] << " + " << m_intermediate[i] << endl;
    }
#endif

  // clean up for estimation mode
  if (!computeTables) {
	  this->reset() ;
/*    for (vector<vector<Function*> >::iterator itA = m_augmented.begin(); itA!=m_augmented.end(); ++itA)
      for (vector<Function*>::iterator itB = itA->begin(); itB!=itA->end(); ++itB)
        delete *itB;
    m_augmented.clear();
    m_intermediate.clear();*/
  }

  return memSize;
}

// Re-parameterize problem using FGLP on original factors
bool MiniBucketElim::DoFGLP() {
  bool changed_functions = false;

  if (m_options && (m_options->mplp > 0 || m_options->mplps > 0)) {
    if (m_options->usePriority) {
      m_fglpRoot = new PriorityFGLP(m_problem, m_options->useNullaryShift);
    } else {
      m_fglpRoot = new FGLP(m_problem, m_options->useNullaryShift);
    }
    m_fglpRoot->Run(m_options->mplp < 0 ? 5 : m_options->mplp,
                    m_options->mplps, m_options->mplpt);
    m_fglpRoot->set_owns_factors(false);
    m_problem->replaceFunctions(m_fglpRoot->factors(), true);
    changed_functions = true;
  }
  return changed_functions;
}

bool MiniBucketElim::DoJGLP() {
  assert(m_pseudotree);
  bool changed_functions = false;

  if (m_options && (m_options->jglp > 0 || m_options->jglps >0))  {
    mex::mbe _jglp(CopyFactors());
    mex::VarOrder var_order(m_pseudotree->getElimOrder().begin(),
                            --m_pseudotree->getElimOrder().end());
    _jglp.setOrder(var_order);

    mex::VarOrder parents(m_problem->getN() - 1); // copy pseudotree information
    for (int i = 0; i < m_problem->getN() - 1; ++i) {
      int parent_var = m_pseudotree->getNode(i)->getParent()->getVar();
      parents[i] = parent_var == m_pseudotree->getRoot()->getVar() 
                   ? -1 : parent_var;
    }
    _jglp.setPseudotree(parents);
    _jglp.setIBound(m_options->jglpi);
    _jglp.setProperties("DoMatch=1,DoFill=1,DoJG=1,DoMplp=0,DoHeur=0");
    _jglp.init();

    _jglp.tighten(m_options->jglp > 0 ? m_options->jglp : 100,
                  m_options->jglps);
    RewriteFactors(_jglp.factors());
    changed_functions = true;
  }
  return changed_functions;
}

// Copy daoopt Function class into mex::Factor class structures
mex::vector<mex::Factor> MiniBucketElim::CopyFactors() {
  mex::vector<mex::Factor> functions(m_problem->getC());
  for (int i = 0; i < m_problem->getC(); ++i) {
    functions[i] = m_problem->getFunctions()[i]->asFactor().exp();
  }
  return functions;
}

void MiniBucketElim::RewriteFactors(const vector<mex::Factor>& factors) {
  vector<Function*> new_functions;
  for (size_t function_idx = 0; function_idx < factors.size(); ++function_idx) {
    const mex::Factor& factor = factors[function_idx];
    double* table_ptr = new double[factor.nrStates()];
    std::set<int> scope;
    for (mex::VarSet::const_iterator var = factor.vars().begin();
         var != factor.vars().end(); ++var) {
      scope.insert(var->label());
    }
    if (scope.size() > 0 && factor.nrStates() == 1) {
      continue;
    }
    new_functions.push_back(
        new FunctionBayes(function_idx, m_problem, scope, table_ptr,
                          factor.nrStates()));
    // Current Mex code leaves factors in log form
    new_functions.back()->fromFactor(factor);
  }
  // Replace the problem definition with the new functions.
  m_problem->replaceFunctions(new_functions);

  cout << "Rewrote factors, problem size (MB) now: " <<
                 m_problem->getSize() * sizeof(double) /
                     (1024 * 1024.0) << endl;
}

void MiniBucketElim::LPReparameterization() {
  if (m_options->mplp > 0 || m_options->mplps > 0) {
    cout << "Running FGLP" << endl;
    DoFGLP();
    m_pseudotree->resetFunctionInfo(m_problem->getFunctions());
  }

  if (m_options->jglp > 0 || m_options->jglps > 0) {
    if (m_options->jglpi > 0 && m_options->memlimit != NONE) {
      LimitJGLPIBound(m_options->memlimit, NULL);
      m_options->jglpi /= 2;
      cout << "Adjusted JGLP i-bound to maximum i-bound / 2" << endl;
    } else {
      m_options->jglpi = m_ibound / 2;
      cout << "Setting JGLP i-bound to the current i-bound / 2" << endl;
    }
    cout << "Running JGLP with i-bound " << m_options->jglpi << endl;
    DoJGLP();
    m_pseudotree->resetFunctionInfo(m_problem->getFunctions());

    // May need to readjust i-bound based on new parameterization
    cout << "Readjusting minibucket i-bound" << endl;
    limitSize(m_options->memlimit, NULL);
  }
}

/* finds a dfs order of the pseudotree (or the locally restricted subtree)
 * and writes it into the argument vector */
void MiniBucketElim::findDfsOrder(vector<int>& order) const {
  order.clear();
  stack<PseudotreeNode*> dfs;
  dfs.push(m_pseudotree->getRoot());
  PseudotreeNode* n = NULL;
  while (!dfs.empty()) {
    n = dfs.top();
    dfs.pop();
    order.push_back(n->getVar());
    for (vector<PseudotreeNode*>::const_iterator it=n->getChildren().begin();
          it!=n->getChildren().end(); ++it) {
      dfs.push(*it);
    }
  }
}


size_t MiniBucketElim::limitSize(size_t memlimit, const vector<val_t> * assignment) {

  // convert to bits
  memlimit *= 1024 *1024 / sizeof(double);

  int decreaseIbd = 0;
  int ibound = m_options->ibound;
  if (ibound <= 0)  { // if ibound were provided, skip this step (finding max. i bound)
	  decreaseIbd = m_options->ibound;
	  ibound = m_options->jglpi; // set ibound = induced width
	  m_options->jglpi = 0; // reset jglpi to be 0 (jglpi used to pass induced width)
  }
  
  cout << "Adjusting mini bucket i-bound..." << endl;
  this->setIbound(ibound);
  size_t mem = this->build(assignment, false);
  cout << " i=" << ibound << " -> " << ((mem / (1024*1024.0)) * sizeof(double) )
       << " MBytes" << endl;

  while (mem > memlimit && ibound > 1) {
    this->setIbound(--ibound);
    mem = this->build(assignment, false);
    cout << " i=" << ibound << " -> " << ((mem / (1024*1024.0)) * sizeof(double) )
         << " MBytes" << endl;
  }

  if (decreaseIbd < 0) {
	  ibound = ibound + decreaseIbd;
	  cout << " decreased i bound by " << decreaseIbd << ", i= " << ibound << endl;
  }
  m_options->ibound = ibound;
  return mem;
}

size_t MiniBucketElim::LimitJGLPIBound(size_t memlimit,
                                       const vector<val_t>* assignment) {
  memlimit *= 1024 * 1024 / sizeof(double);

  int ibound = m_options->jglpi;
  int original_ibound = this->getIbound();
  this->setIbound(ibound);

  size_t mem = this->build(assignment, false);
  cout << " i=" << ibound << " -> " << (mem / 1024 * 1024.0) * sizeof(double)
       << " MBytes" << endl;

  while (mem > memlimit && ibound > 1 ) {
    this->setIbound(--ibound);
    mem = this->build(assignment, false);
    cout << " i=" << ibound << " -> " << (mem / 1024 * 1024.0) * sizeof(double)
         << " MBytes" << endl;
  }
  this->setIbound(original_ibound);
  m_options->jglpi = ibound;

  return mem;
}


size_t MiniBucketElim::getSize() const {
  size_t S = 0;
  for (vector<vector<Function*> >::const_iterator it=m_augmented.begin(); it!= m_augmented.end(); ++it) {
    for (vector<Function*>::const_iterator itF=it->begin(); itF!=it->end(); ++itF)
      S += (*itF)->getTableSize();
  }
  return S;
}


/*
 * mini bucket file format (all data in binary):
 * - size_t: no. of variables
 * - int: i-bound
 * - double: global upper bound
 * for every variable:
 *   - size_t: number of functions in bucket structure
 *   for every such function:
 *     - int: function ID
 *     - size_t: scope size
 *     for every scope variable:
 *       - int: variable index
 *     - size_t: table size
 *     for every table entry:
 *       - double: CPT entry
 * for every variable:
 *   - size_t: number of intermediate function pointers
 *   for every function pointer:
 *     - size_t: function index (implicit order from above)
 */

bool MiniBucketElim::writeToFile(string fn) const {

  ogzstream out(fn.c_str());
  if ( ! out ) {
    cerr << "Error writing mini buckets to file " << fn << endl;
    return false;
  }

  // used later
  int x = NONE;
  size_t y = NONE;

  // number of variables
  size_t sz = m_augmented.size();
  out.write((char*)&( sz ), sizeof( sz ));

  // i-bound
  out.write((char*)&( m_ibound ), sizeof( m_ibound ));

  // global UB
  out.write((char*)&( m_globalUB ), sizeof( m_globalUB ));

  map<const Function*,size_t> funcMap;

  // over m_augmented
  for (size_t i=0; i<sz; ++i) {
    size_t sz2 = m_augmented[i].size();
    out.write((char*)&( sz2 ), sizeof( sz2 ));

    vector<Function*>::const_iterator itF = m_augmented[i].begin();
    for (size_t j=0; j<sz2; ++j, ++itF) {
      const Function* f = *itF;
      funcMap.insert(make_pair(f,funcMap.size()));

      // function ID
      int id = f->getId();
      out.write((char*)&( id ), sizeof( id ));


      // scope
      size_t sz3 = f->getScopeVec().size();
      out.write((char*)&( sz3 ), sizeof( sz3 ));
 // scope size
      for (vector<int>::const_iterator it=f->getScopeVec().begin(); it!=f->getScopeVec().end(); ++it) {
        x = *it;
        out.write((char*)&( x ), sizeof( x ));
 // vars from scope
      }

      // table size
      sz3 = f->getTableSize();
      out.write((char*)&( sz3 ), sizeof( sz3 ));

      // table
      out.write((char*) ( f->getTable() ), sizeof( double ) * sz3);

    }


  }

  // over m_intermediate
  for (size_t i=0; i<sz; ++i) {
    size_t sz2 = m_intermediate[i].size();
    out.write((char*)&( sz2 ), sizeof( sz2 ));

    vector<Function*>::const_iterator itF = m_intermediate[i].begin();
    for (size_t j=0; j<sz2; ++j, ++itF) {
      y = funcMap.find(*itF)->second;
      out.write((char*) &( y ), sizeof(y));
    }

  }

  out.close();

  return true;

}


bool MiniBucketElim::readFromFile(string fn) {

  ifstream inTemp(fn.c_str());
  inTemp.close();
  if (inTemp.fail()) { // file not existent yet
    return false;
  }

  igzstream in(fn.c_str());

  this->reset();

  // used later
  int x = NONE;
  size_t y = NONE;
  vector<Function*> allFuncs;

  // no. of variables
  size_t sz;
  in.read((char*) &( sz ), sizeof( sz ));

  if (sz != (size_t) m_problem->getN()) {
    cerr << "Number of variables in mini bucket file doesn't match" << endl;
    return false;
  }

  m_augmented.resize(sz);
  m_intermediate.resize(sz);

  // i-bound
  int ibound;
  in.read((char*) &(ibound), sizeof(ibound));
  m_ibound = ibound;

  // global UB
  double ub;
  in.read((char*) &(ub), sizeof( ub ));
  m_globalUB = ub;

  // over variables for m_augmented
  for (size_t i=0; i<sz; ++i) {
    size_t sz2;
    in.read((char*) &(sz2), sizeof(sz2));

    // over functions
    for (size_t j=0; j<sz2; ++j) {
      int id;
      in.read((char*) &( id ), sizeof(id));

      // scope
      size_t sz3;
      in.read((char*) &(sz3), sizeof(sz3));
      set<int> scope;
      for (size_t k=0; k<sz3; ++k) {
        in.read((char*) &( x ), sizeof(x));
        scope.insert(x);
      }

      // table size and table
      in.read((char*) &( sz3 ), sizeof(sz3));
      double* T = new double[sz3];
      in.read((char*) ( T ), sizeof(double)*sz3);

      // create function and store it
      Function* f = new FunctionBayes(id,m_problem,scope,T,sz3);
      m_augmented[i].push_back(f);
      allFuncs.push_back(f);
    }
  }

  for (size_t i=0; i<sz; ++i) {
    // no. of function pointers
    size_t sz2;
    in.read((char*) &(sz2), sizeof(sz2));

    for (size_t j=0; j<sz2; ++j) {
      // function index
      in.read((char*) &(y), sizeof(y));
      m_intermediate[i].push_back(allFuncs.at(y));
    }
  }

  in.close();
  cout << "Read mini bucket with i-bound " << ibound << " from file " << fn << endl;
  return true;
}

}  // namespace daoopt
