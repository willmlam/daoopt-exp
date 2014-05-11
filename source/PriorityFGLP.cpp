#include "PriorityFGLP.h"
#include <chrono>
using namespace std;
using namespace std::chrono;

PriorityFGLP::PriorityFGLP(Problem *p, bool use_nullary_shift) 
    : FGLP(p, use_nullary_shift) {

}

PriorityFGLP::PriorityFGLP(PriorityFGLP *parent_fglp, 
        const map<int,val_t> &assignment)
: var_priority_(parent_fglp->var_priority_) {
    problem_ = parent_fglp->problem_;
    owns_factors_ = true;
    factors_by_variable_.resize(problem_->getN(), vector<Function*>());
    global_const_factor_ = nullptr;
    max_marginals_.resize(problem_->getN(), vector<double*>());
    ub_ = ELEM_ONE;
    use_nullary_shift_ = parent_fglp->use_nullary_shift_;
    verbose_ = false;
    
    Condition(parent_fglp->factors(), assignment);

    for (Function *f : factors_) {
        for (int vs : f->getScopeVec()) {
            factors_by_variable_[vs].push_back(f);
        }
    }

    // allocate max marginal storage
    for (size_t i = 0; i < factors_by_variable_.size(); ++i) {
        if (factors_by_variable_[i].size() == 0) continue;
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
    for (iter = 0; iter < max_updates || (max_updates == -1 && max_time > 0); ++iter) {
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

        vector<double*> &var_max_marginals = max_marginals_[v];

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
            }
        }

        // Compute average max marginal
        double *avg_mm = new double[table_size];

        for (size_t i = 0; i < table_size; ++i) avg_mm[i] = ELEM_ONE;

        for (double *mm : var_max_marginals) {
            for (size_t i = 0; i < table_size; ++i)
                avg_mm[i] OP_TIMESEQ mm[i];
        }
        for (size_t i = 0; i < table_size; ++i) {
            if (use_nullary_shift_)
                avg_mm[i] OP_DIVIDEEQ total_nullary_shift;
            avg_mm[i] = OP_ROOT(avg_mm[i],var_max_marginals.size());
        }

        // Compute new priority values for each function
        for (size_t i = 0; i < factors_by_variable_[v].size(); ++i) {
            double p_new = MessageDist(var_max_marginals[i], avg_mm, total_nullary_shift / var_max_marginals.size(), v);
            Function *f = factors_by_variable_[v][i];
            for (int vs : f->getScopeVec()) {
                if (vs == v) continue;
                double p_max = max(p_new, var_priority_.getP(vs));
                var_priority_.insert(p_max, vs);
            }
        }

        // Reparameterize
        for (size_t i = 0; i < factors_by_variable_[v].size(); ++i) {
            Reparameterize(factors_by_variable_[v][i],
                    var_max_marginals[i],avg_mm,v);
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
    runtime_ = double(duration_cast<milliseconds>(
                time_end - time_start).count()) / 1000;
    runiters_ = iter;

    if (verbose_) {
        cout << "FGLP (" << iter << " iter, " << runtime_ << " sec): " << ub_ << endl;
    }

}

void PriorityFGLP::Condition(const vector<Function*> &fns,
        const map<int,val_t> &assn) {

    double *table_const = new double[1];
    table_const[0] = ELEM_ONE;

    for (Function *f : fns) {
        if (f->getArity() == 0) {
            table_const[0] OP_TIMESEQ f->getTable()[0];
            continue;
        }

        Function *new_f = f->substitute(assn);

        // If function changes, we will need to recompute the updates
        // for all of the variables in this function

        if (new_f->isConstant()) {
            table_const[0] OP_TIMESEQ new_f->getTable()[0];
            delete new_f;
        }
        else {
            // If function is changed, increase the priority of all its variables
            // to max
            if (new_f->getArity() != f->getArity()) {
                for (int vs : new_f->getScopeVec()) {
                    var_priority_.insert(numeric_limits<double>::max(), vs);
                }
            }
            factors_.push_back(new_f);
        }
    }
    global_const_factor_ = new FunctionBayes(factors_.back()->getId()+1, 
                problem_, set<int>(), table_const, 1);
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
        }
        else {
            diff = fabs(m1[i] - m2_orig);
        }
        dist += diff;
    }
    return dist;
}
