#include "PriorityFGLP.h"
#include <chrono>

namespace daoopt {

using namespace std;
using namespace std::chrono;

PriorityFGLP::PriorityFGLP(Problem* p, bool use_nullary_shift)
    : FGLP(p, use_nullary_shift) {}

PriorityFGLP::PriorityFGLP(PriorityFGLP* parent_fglp,
                           const map<int, val_t>& assignment,
                           const set<int>& sub_vars,
                           int condition_var) {
  problem_ = parent_fglp->problem_;
  owns_factors_ = true;
  factors_by_variable_.resize(problem_->getN(), vector<Function *>());
  global_const_factor_ = nullptr;
  bound_contribs_ = parent_fglp->bound_contribs_;
  total_shift_ = parent_fglp->total_shift_;
  max_marginals_.resize(problem_->getN(), vector<double*>());
  conditioned_cost_ = parent_fglp->conditioned_cost_;
  ub_ = ELEM_ONE;
  use_nullary_shift_ = parent_fglp->use_nullary_shift_;
  use_cost_shift_reversal_ = parent_fglp->use_cost_shift_reversal_;
  verbose_ = false;

  // Copy over only relevant priorities.
  for (const int var : sub_vars) {
    var_priority_.insert(parent_fglp->var_priority_.getP(var), var);
  }
  Condition(parent_fglp->factors(), assignment, sub_vars, condition_var);

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
    max_marginals_[i].resize(factors_by_variable_[i].size(), nullptr);
  }
}

void PriorityFGLP::Run(int max_updates, double max_time, double tolerance) {
  double diff = numeric_limits<double>::max();
  UpdateUB();
  if (verbose_) {
    cout << "Initial UB: " << ub_ << endl;
  }

  // Initialize all priorities to max if nothing is in the queue
  if (var_priority_.size() == 0) {
    for (int i = 0; i < problem_->getN(); ++i)
      var_priority_.insert(numeric_limits<double>::max(), i);
  }

  steady_clock::time_point time_start = steady_clock::now();
  steady_clock::time_point time_end;

  int iter;
  for (iter = 0; iter < max_updates || (max_updates == -1 && max_time > 0) ||
       var_priority_.top().first == numeric_limits<double>::max();
       ++iter) {
    time_end = steady_clock::now();
    if (max_time > 0 &&
        duration_cast<seconds>(time_end - time_start).count() >= max_time)
      break;

    if (var_priority_.size() == 0) {
      break;
    }

    // Prioritized update
    double p = var_priority_.top().first;
    if (fabs(p) < tolerance || std::isnan(p)) break;
    int v = var_priority_.top().second;
    var_priority_.pop();

    if (v >= int(factors_by_variable_.size())) continue;

    // Skip this variable, no factors have this variable
    if (factors_by_variable_[v].size() == 0) continue;

    vars_updated_.insert(v);

    // To store max marginals of each of the functions

    vector<double *> &var_max_marginals = max_marginals_[v];

    vector<double> mm_max(var_max_marginals.size(),
                          -numeric_limits<double>::max());

    size_t table_size = problem_->getDomainSize(v);
    double total_nullary_shift = ELEM_ONE;
    // Compute max marginals and shifted versions
    for (size_t i = 0; i < factors_by_variable_[v].size(); ++i) {
      Function *f = factors_by_variable_[v][i];
      var_max_marginals[i] = MaxMarginal(f, v);

      if (use_nullary_shift_) {
        for (size_t j = 0; j < table_size; ++j) {
          mm_max[i] = max(mm_max[i], var_max_marginals[i][j]);
        }
        total_nullary_shift OP_TIMESEQ mm_max[i];

        // Shift max into the global constant
        global_const_factor_->getTable()[0] OP_TIMESEQ mm_max[i];
        bound_contribs_[v] OP_TIMESEQ mm_max[i];
      }
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
      if (use_nullary_shift_) {
        avg_mm[i] OP_DIVIDEEQ total_nullary_shift;
      }
      avg_mm[i] = OP_ROOT(avg_mm[i], var_max_marginals.size());
    }

    // Compute new priority values for each function
    for (size_t i = 0; i < factors_by_variable_[v].size(); ++i) {
      double p_new =
          MessageDist(var_max_marginals[i], avg_mm,
                      total_nullary_shift / var_max_marginals.size(), v);
      Function *f = factors_by_variable_[v][i];
      for (int vs : f->getScopeVec()) {
        if (vs == v) continue;
        double p_max = max(p_new, var_priority_.getP(vs));
        var_priority_.insert(p_max, vs);
      }
    }

    // Reparameterize
    for (size_t i = 0; i < factors_by_variable_[v].size(); ++i) {
      Function* fv = factors_by_variable_[v][i];
      Reparameterize(fv, var_max_marginals[i], avg_mm, v);

      if (use_cost_shift_reversal_) {
        // Record reparameterization
        for (size_t j = 0; j < table_size; ++j) {
          total_shift_[fv->getId()][v][j] OP_TIMESEQ
            avg_mm[j] OP_DIVIDE var_max_marginals[i][j];
        }
      }
      delete var_max_marginals[i];
      var_max_marginals[i] = nullptr;
    }
    delete avg_mm;

    if (iter > 0 && iter % problem_->getN() == 0) {
      diff = UpdateUB();
      if (verbose_) {
        cout << "UB: " << ub_ << " (d=" << diff << ")" << endl;
        /*
           cout << "constant: " << m_globalConstFactor->getTable()[0] << endl;
           cout << "non-const UB: " << m_UBNonConstant << endl;
           */
      }
    }
  }
  UpdateUB();
  time_end = steady_clock::now();
  runtime_ =
      double(duration_cast<milliseconds>(time_end - time_start).count()) / 1000;
  runiters_ = iter;

  if (verbose_) {
    cout << "FGLP (" << iter << " iter, " << runtime_ << " sec): " << ub_
         << endl;
  }
}

void PriorityFGLP::Condition(const vector<Function*>& fns,
                             const map<int, val_t>& assn,
                             const set<int>& sub_vars,
                             int condition_var) {

  // Should be identical to FGLP::Condition except for a step that updates
  // the priorities.
  double *table_const = new double[1];

  if (condition_var >= 0 && !use_cost_shift_reversal_) {
    conditioned_cost_ OP_TIMESEQ bound_contribs_[condition_var];
  }

  table_const[0] = conditioned_cost_;

  for (const int v : sub_vars) {
    table_const[0] OP_TIMESEQ bound_contribs_[v];
  }

  for (Function *f : fns) {
    if (f->getArity() == 0) {
      continue;
    }

    Function* new_f = nullptr;
    if (intersectionEmpty(f->getScopeSet(), sub_vars) &&
        !(ContainsKey(f->getScopeSet(), condition_var) &&
          f->getArity() == 1)) {
      continue;
    } else if (ContainsKey(f->getScopeSet(), condition_var)) {
      new_f = f->substitute(assn);
      assert(f->getArity() - new_f->getArity() == 1);

      if (use_cost_shift_reversal_) {
        double shift = ELEM_ONE;
        val_t value = assn.find(condition_var)->second;
        shift OP_TIMESEQ total_shift_[f->getId()][condition_var][value];

        for (size_t i = 0; i < new_f->getTableSize(); ++i) {
          new_f->getTable()[i] OP_DIVIDEEQ shift;
        }
      }
    } else {
      new_f = f->clone();
    }

    if (new_f->isConstant()) {
      table_const[0] OP_TIMESEQ new_f->getTable()[0];
      conditioned_cost_ OP_TIMESEQ new_f->getTable()[0];
      delete new_f;
    } else {
      // Otherwise, increase the priority of all its variables
      // to max
      for (int vs : new_f->getScopeVec()) {
        var_priority_.insert(numeric_limits<double>::max(), vs);
      }
      factors_.push_back(new_f);
    }
  }
  global_const_factor_ = new FunctionBayes(
      factors_.back()->getId() + 1, problem_, set<int>(), table_const, 1);
  factors_.push_back(global_const_factor_);
}

double PriorityFGLP::MessageDist(double *m1, double *m2, double ns, int var) {
  double dist = 0;
  for (int i = 0; i < problem_->getDomainSize(var); ++i) {
    //        dist = max(dist, fabs(m1[i] - m2[i]));
    double m2_orig = m2[i] + ns;
    double diff;
    if (std::isinf(m1[i]) && std::isinf(m2_orig)) {
      diff = 0;
    } else {
      diff = fabs(m1[i] - m2_orig);
    }
    dist += diff;
  }
  return dist;
}

}  // namespace daoopt
