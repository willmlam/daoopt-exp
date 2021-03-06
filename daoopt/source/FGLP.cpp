#include "FGLP.h"
#include <chrono>

namespace daoopt {
using namespace std;
using namespace std::chrono;

FGLP::FGLP(Problem *p, bool use_nullary_shift)
  : problem_(p),
    owns_factors_(true),
    factors_by_variable_(p->getN(), vector<Function*>()),
    global_const_factor_(nullptr),
    bound_contribs_(problem_->getN(), ELEM_ONE),
    max_marginals_(p->getN(), vector<double*>()),
    conditioned_cost_(ELEM_ONE),
    ub_(ELEM_ONE),
    use_nullary_shift_(use_nullary_shift),
    use_cost_shift_reversal_(false),
    verbose_(true) {

      for (Function *f : problem_->getFunctions()) {
        if (f->getArity() == 0) {
          global_const_factor_ = f->clone();
          continue;
        }
        factors_.push_back(f->clone());
        for (int vs : factors_.back()->getScopeVec()) {
          factors_by_variable_[vs].push_back(factors_.back());
        }
      }

      // Allocate total shift
      total_shift_.resize(factors_.back()->getId()+1,
                          map<int,vector<double>>());
      for (Function *f : factors_) {
        for (int v : f->getScopeVec()) {
          total_shift_[f->getId()].insert(
              make_pair(v, vector<double>(problem_->getDomainSize(v),
                                          ELEM_ONE)));
        }
      }
      if (!global_const_factor_) {
        double *table_const = new double[1];
        table_const[0] = ELEM_ONE;
        global_const_factor_ = new FunctionBayes(factors_.back()->getId()+1, 
            problem_, set<int>(), table_const, 1);
      }
      factors_.push_back(global_const_factor_);
      conditioned_cost_ = global_const_factor_->getTable()[0];

      // allocate max marginal storage
      for (int i = 0; i < problem_->getN(); ++i) {
        max_marginals_[i].resize(factors_by_variable_[i].size(), nullptr);
      }

    }

FGLP::FGLP(FGLP *parent_fglp, const map<int,val_t> &assignment, 
    const set<int> &subVars, int condition_var)
: 
  problem_(parent_fglp->problem_),
  owns_factors_(true),
  factors_by_variable_(problem_->getN(), vector<Function*>()),
  global_const_factor_(nullptr),
  bound_contribs_(parent_fglp->bound_contribs_),
  total_shift_(parent_fglp->total_shift_),
  max_marginals_(problem_->getN(), vector<double*>()),
  conditioned_cost_(parent_fglp->conditioned_cost_),
  ub_(ELEM_ONE),
  use_nullary_shift_(parent_fglp->use_nullary_shift_),
  use_cost_shift_reversal_(parent_fglp->use_cost_shift_reversal_),
  verbose_(false) {

    Condition(parent_fglp->factors(), assignment, subVars, condition_var);

    for (Function *f : factors_) {
      for (int vs : f->getScopeVec()) {
        factors_by_variable_[vs].push_back(f);
      }
    }

    // allocate max marginal storage
    for (size_t i = 0; i < factors_by_variable_.size(); ++i) {
      if (factors_by_variable_[i].size() == 0) {
        continue;
      }
      max_marginals_[i].resize(factors_by_variable_[i].size(), NULL);
    }
  }

void FGLP::Run(int max_iter, double max_time, double tolerance) {
  double diff = numeric_limits<double>::max();
  UpdateUB();
  if (verbose_) {
    cout << "Initial UB: " << ub_ << endl;
  }
  steady_clock::time_point time_start = steady_clock::now();
  steady_clock::time_point time_end;

  int iter;
  for (iter = 0; iter < max_iter || (max_iter == -1 && max_time > 0); ++iter) {
    time_end = steady_clock::now();
    if (max_time > 0 && 
        duration_cast<seconds>(time_end - time_start).count() >= max_time) 
      break;

    if (fabs(diff) < tolerance || std::isnan(diff)) break;

    // Basic round robin update
    for (int v = 0; v < problem_->getN(); ++v) {
      vars_updated_.insert(v);
      time_end = steady_clock::now();
      if (max_time > 0 && 
          duration_cast<seconds>(time_end - time_start).count() >= max_time) 
        break;

      if (v >= int(factors_by_variable_.size())) continue;

      // Skip this variable, no factors have this variable
      if (factors_by_variable_[v].size() == 0) continue;

      // To store max marginals of each of the functions

      vector<double*> &var_max_marginals = max_marginals_[v];

      vector<double> mm_max(var_max_marginals.size(),
          -numeric_limits<double>::max());

      size_t table_size = problem_->getDomainSize(v);
      double total_nullary_shift = ELEM_ONE;
      // Compute max marginals and shifted versions
      for (size_t i = 0; i < factors_by_variable_[v].size(); ++i) {
        Function *fv = factors_by_variable_[v][i];
        var_max_marginals[i] = MaxMarginal(fv, v);
      }

      // Compute average max marginal
      double *avg_mm = new double[table_size];

      for (size_t i = 0; i < table_size; ++i) {
        avg_mm[i] = ELEM_ONE;
      }

      for (double *mm : var_max_marginals) {
        for (size_t i = 0; i < table_size; ++i) {
          avg_mm[i] OP_TIMESEQ mm[i];
        }
      }
      for (size_t i = 0; i < table_size; ++i) {
        avg_mm[i] = OP_ROOT(avg_mm[i],var_max_marginals.size());
      }

      // Compute nullary shift by shifting out minimum of avg_mm
      if (use_nullary_shift_) {
        double max_avg_mm = -numeric_limits<double>::max();
        for (size_t i = 0; i< table_size; ++i) {
          max_avg_mm = max(max_avg_mm, avg_mm[i]);
        }
        for (size_t i = 0; i< table_size; ++i) {
          avg_mm[i] OP_DIVIDEEQ max_avg_mm;
        }
        total_nullary_shift = var_max_marginals.size() * max_avg_mm;
        // Shift max into the global constant
        global_const_factor_->getTable()[0] OP_TIMESEQ total_nullary_shift;
        bound_contribs_[v] OP_TIMESEQ total_nullary_shift;
      }

      // Reparameterize
      for (size_t i = 0; i < factors_by_variable_[v].size(); ++i) {
        Function *fv = factors_by_variable_[v][i];
        Reparameterize(fv, var_max_marginals[i], avg_mm, v);

        if (use_cost_shift_reversal_) {
          // Record reparameterization
          for (size_t j = 0; j < table_size; ++j) {
            total_shift_[fv->getId()][v][j] OP_TIMESEQ 
              avg_mm[j] OP_DIVIDE var_max_marginals[i][j];
          }
        }
        delete var_max_marginals[i];
        var_max_marginals[i] = NULL;
      }
      delete avg_mm;
    }
    diff = UpdateUB();
    if (verbose_) {
      cout << "UB: " << ub_ << " (d=" << diff << ")" << endl;
      /*
         cout << "constant: " << m_globalConstFactor->getTable()[0] << endl;
         cout << "non-const UB: " << m_UBNonConstant << endl;
         */
    }
  }
  time_end = steady_clock::now();
  runtime_ = double(duration_cast<milliseconds>(
        time_end - time_start).count()) / 1000;
  runiters_ = iter;

  if (verbose_) {
    cout << "FGLP (" << iter << " iter, " << runtime_ << " sec): " << ub_ 
         << endl;
  }

}

void FGLP::GetLabelAll(int var, vector<double> &out) {
  out.clear();
  out.resize(problem_->getDomainSize(var), ELEM_ONE);

  // For each function, check if it's unary and contains the variable
  for (Function *f : factors_) {
    if (f->getArity() == 1 && f->hasInScope(var)) {
      for (size_t i = 0; i < out.size(); ++i) {
        out[i] OP_TIMESEQ f->getTable()[i];
      }
    }
  }
}

size_t FGLP::GetSize() const {
  size_t s = 0;
  for (Function *f : factors_)
    s += f->getTableSize();
  return s;
}

void FGLP::Condition(const vector<Function*>& fns, const map<int,val_t>& assn, 
    const set<int>& sub_vars, int condition_var) {
  double *table_const = new double[1];

  // Move the bound contributed by the variable about to be conditioned
  // into the conditioned cost (g), but only if it will not be reversed.
  if (condition_var >= 0 && !use_cost_shift_reversal_) {
    conditioned_cost_ OP_TIMESEQ bound_contribs_[condition_var];
  }
  // Calculate global constant
  table_const[0] = conditioned_cost_;

  // Add in contributions from each variable
  for (int v : sub_vars) {
    table_const[0] OP_TIMESEQ bound_contribs_[v];
  }

  // The constant is now representing g+h
  for (Function *f : fns) {
    if (f->getArity() == 0) {
      continue;
    }

    Function* new_f = nullptr;
    // Figure out whether we need to drop this function
    // We do so if it's not part of the subproblem and that 
    // it's not a unary function on the conditioning variable
    if (intersectionEmpty(f->getScopeSet(), sub_vars) && 
        !(ContainsKey(f->getScopeSet(), condition_var) &&
          f->getArity() == 1)) {
      continue;
    } else if (ContainsKey(f->getScopeSet(), condition_var)) {
    // At this point we know that the function needs to be processed
    // (conditioned or copied)
      new_f = f->substitute(assn);
      assert(f->getArity() - new_f->getArity() == 1);

      if (use_cost_shift_reversal_) {
        // Compute function reversal shift needed
        // We compute the total cost shifted into the function over 
        // the variable/value that was just conditioned.
        double shift = ELEM_ONE;
        val_t value = assn.find(condition_var)->second;
        shift OP_TIMESEQ total_shift_[f->getId()][condition_var][value];

        for (size_t i = 0; i < new_f->getTableSize(); ++i) {
          new_f->getTable()[i] OP_DIVIDEEQ shift;
        }
      }
    } else {
      // Clone it if it is not affected by the conditioning
      new_f = f->clone();
    } 

    if (new_f->isConstant()) {
      table_const[0] OP_TIMESEQ new_f->getTable()[0];
      conditioned_cost_ OP_TIMESEQ new_f->getTable()[0];
      delete new_f;
    } else {
      factors_.push_back(new_f);
    }
  }

  // The resulting constant is g+h+additional g from conditioning.
  global_const_factor_ = new FunctionBayes(factors_.back()->getId()+1, 
      problem_, set<int>(), table_const, 1);
  factors_.push_back(global_const_factor_);
}

double FGLP::UpdateUB() {
  double old_ub = ub_;
  if (!use_nullary_shift_) {
    ub_ = ELEM_ONE;
    for (Function *f : factors_) {
      double z = ELEM_ZERO;
      for (size_t i = 0; i < f->getTableSize(); ++i)
        z = max(z,f->getTable()[i]);
      ub_ OP_TIMESEQ z;
    }
  }
  // If using nullary shift, the bound is collected in the nullary function
  else {
    ub_ = global_const_factor_->getTable()[0];
  } 
  return ub_ OP_DIVIDE old_ub;
}

void FGLP::GetVarUB(int var, vector<double> &out) {
  // Fill out the output with the current g+h. Then figure out what we need 
  // to subtract.
  for (int i = 0; i < problem_->getDomainSize(var); ++i)
    out[i] OP_TIMESEQ global_const_factor_->getTable()[0];
  /*
     cout << "var: " << var << endl;
     cout << "out: " << out << endl;
     */
  for (Function *f : factors_) {
    if (f->getArity() == 0) {
      continue;
    }
    // If the factor does not have the variable, take the maximum
    if (!f->hasInScope(var)) {
      double z = ELEM_ZERO;
      for (size_t i = 0; i < f->getTableSize(); ++i)
        z = max(z, f->getTable()[i]);
      for (int i = 0; i < problem_->getDomainSize(var); ++i)
        out[i] OP_TIMESEQ z;
      continue;
    }

    // Read off table directly if the factor is unary over the variable
    if (f->getArity() == 1 && f->hasInScope(var)) {
      for (int i = 0; i < problem_->getDomainSize(var); ++i) {
        out[i] OP_TIMESEQ f->getTable()[i];
      }
      continue;
    }

    // Otherwise, take the maximum of the the subtable
    val_t *tuple = new val_t[f->getArity()];
    val_t *unfixed_tuple = tuple + 1;

    for (int i = 0; i < f->getArity(); ++i) tuple[i] = 0;

    vector<val_t> unfixed_domains;
    for (int v : f->getScopeSet())
      if (v != var) unfixed_domains.push_back(problem_->getDomainSize(v));

    vector<val_t> var_domain;
    var_domain.push_back(problem_->getDomainSize(var));

    vector<val_t*> idx_map;
    int i = 0;
    for (int v : f->getScopeSet()) {
      if (v == var)
        idx_map.push_back(&tuple[0]);
      else
        idx_map.push_back(&unfixed_tuple[i++]);
    }

    size_t idx = 0;
    do {
      double z = ELEM_ZERO;
      do {
        z = max(z, f->getValuePtr(idx_map));
      } while (increaseTuple(idx, unfixed_tuple, unfixed_domains));
      out[tuple[0]] OP_TIMESEQ z;
    } while (increaseTuple(idx, tuple, var_domain));
    delete [] tuple;
  }
}


double *FGLP::MaxMarginal(Function *f, int v) {
  assert(f->hasInScope(v));
  size_t table_size = problem_->getDomainSize(v);
  double *new_table = new double[table_size];
  for (size_t i = 0; i < table_size; ++i) new_table[i] = ELEM_ZERO;

  val_t *tuple = new val_t[f->getArity()];
  val_t *elim_tuple = tuple + 1;

  vector<val_t> elim_domains;
  for (int vs : f->getScopeSet()) {
    if (vs != v) elim_domains.push_back(problem_->getDomainSize(vs));
  }

  vector<val_t> var_domain;
  var_domain.push_back(problem_->getDomainSize(v));

  for (int i = 0; i < f->getArity(); ++i) tuple[i] = 0;

  vector<val_t*> idx_map;

  int i = 0;
  for (int vs : f->getScopeSet()) {
    if (vs == v)
      idx_map.push_back(&tuple[0]);
    else
      idx_map.push_back(&elim_tuple[i++]);
  }

  size_t idx;
  // iterate over all values non-v variables
  do {
    idx = 0;
    do {
      new_table[idx] = max(new_table[idx], f->getValuePtr(idx_map));
    } while (increaseTuple(idx, tuple, var_domain));
  } while (increaseTuple(idx, elim_tuple, elim_domains));
  delete [] tuple;

  return new_table;
}

void FGLP::Reparameterize(Function *f, double *mm, double *avg_mm, int var) {
  assert(f->hasInScope(var));
  vector<val_t> domains;

  val_t *tuple = new val_t[f->getArity()];
  for (int i = 0; i < f->getArity(); ++i) tuple[i] = 0;

  val_t *mm_val = NULL;

  int i = 0;
  for (int vs : f->getScopeSet()) {
    domains.push_back(problem_->getDomainSize(vs));
    if (vs == var) mm_val = &tuple[i];
    i++;
  }

  assert(int(domains.size()) == f->getArity());

  size_t idx = 0;
  do {
    if (std::isinf(avg_mm[*mm_val])) {
      f->getTable()[idx] = ELEM_ZERO;
    }
    else
      f->getTable()[idx] OP_TIMESEQ 
        (avg_mm[*mm_val] OP_DIVIDE mm[*mm_val]);
  } while(increaseTuple(idx, tuple, domains));
  delete [] tuple;
}

void FGLP::PrintAllFactors() const { 
  for (Function *f : factors_) {
    cout << *f << endl;
    for (size_t k = 0; k < f->getTableSize(); ++k) {
      cout << " " << f->getTable()[k] << endl;
    }
    cout << endl;
  }
}

}  // namespace daoopt
