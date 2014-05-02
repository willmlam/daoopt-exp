#include "ResidualFGLP.h"
#include <chrono>

using namespace std;
using namespace std::chrono;

ResidualFGLP::ResidualFGLP(Problem *p,
        bool use_nullary_shift) : 
    problem_(p),
    owns_factors_(true),
    factors_by_variable_(p->getN(), vector<Function*>()),
    global_const_factor_(nullptr),
    ub_(ELEM_ONE),
    use_nullary_shift_(use_nullary_shift),
    verbose_(true),
    message_updates_(vector<FGLPVariableUpdate*>(p->getN(),nullptr)) {

    for (Function *f : p->getFunctions()) {
        if (f->getArity() == 0) {
            global_const_factor_ = f->clone(); 
            continue;
        }
        factors_.push_back(f->clone());
        for (int vs : factors_.back()->getScopeVec()) {
            factors_by_variable_[vs].push_back(factors_.back());
        }
    }

    if (global_const_factor_) {
        factors_.push_back(global_const_factor_);
    }

    for (int i = 0; i < p->getN(); ++i)
        ComputeVariableMessage(i);
//    message_queue_.debug();

}

ResidualFGLP::ResidualFGLP(ResidualFGLP *parent_fglp, 
        const map<int,val_t> &assignment) :
    problem_(parent_fglp->problem_),
    owns_factors_(true),
    factors_by_variable_(problem_->getN(), vector<Function*>()),
    global_const_factor_(nullptr),
    ub_(ELEM_ONE),
    use_nullary_shift_(parent_fglp->use_nullary_shift_),
    verbose_(false),
    message_updates_(parent_fglp->message_updates_),
    message_queue_(parent_fglp->message_queue_) {

    Condition(parent_fglp->factors(), assignment);

    for (Function *f : factors_) {
        for (int vs : f->getScopeVec()) {
            factors_by_variable_[vs].push_back(f);
        }
    }
    for (int v : vars_to_update_)
        ComputeVariableMessage(v, true);
}

void ResidualFGLP::run(int max_updates, double max_time, double tolerance) {
    double diff = numeric_limits<double>::max();
    UpdateUB();
    if (verbose_) {
        cout << "Initial UB: " << ub_ << endl;
    }

    steady_clock::time_point time_start = steady_clock::now();
    steady_clock::time_point time_end;

    int iter;
    for (iter = 0; iter < max_updates || (max_updates = -1 && max_time > 0); 
            ++iter) {
        time_end = steady_clock::now();
        if (max_time > 0 && 
                duration_cast<seconds>(time_end - time_start).count() >= max_time)
            break;

        double dist = message_queue_.top().first;

        if (fabs(dist) < tolerance || std::isnan(dist)) break;

        int var_to_update = message_queue_.top().second;

        message_queue_.pop();

        set<int> to_recompute;
        to_recompute.insert(var_to_update);

        // Apply update from stored messages
        FGLPVariableUpdate *message_update = 
                message_updates_[var_to_update];
        for (Function *f : factors_by_variable_[var_to_update]) {
            for (int v : f->getScopeVec()) 
                to_recompute.insert(v);
            Reparameterize(f, message_update->factor_to_var_[f->getId()],
                message_update->var_to_f_,var_to_update);
        }
        global_const_factor_->getTable()[0] OP_TIMESEQ 
            message_update->nullary_shift();

        // Update messages affected by applied update
        // Replace the current message without modifying
        // the previous if the variable was not actually updated
        for (int v : to_recompute) ComputeVariableMessage(v, v != var_to_update);
        if (iter % problem_->getN() == 0) {
            diff = UpdateUB();
            if (verbose_) 
                cout << "UB: " << ub_ << " (d=" << diff << ")" << endl;
        }
    }
    time_end = steady_clock::now();
    runtime_ = double(
            duration_cast<milliseconds>(time_end - time_start).count()) / 1000;
    runiters_ = iter;

    if (verbose_) {
        cout << "RFGLP (" << iter << " updates, " << runtime_ << " sec): " 
            << ub_ << endl;
    }
}

void ResidualFGLP::GetLabelAll(int var, vector<double> &out) {
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

size_t ResidualFGLP::GetSize() const {
    size_t s = 0;
    for (Function *f : factors_)
        s += f->getTableSize();
    return s;
}

void ResidualFGLP::ComputeVariableMessage(int v, bool replace_mode) {

    int domain_size = problem_->getDomainSize(v);

    // Release old update and create new message update owned by this class
    /*
    if (v < int(message_queue_.size()))
        message_queue_.erase(v);
        */

    if (message_updates_[v] && message_updates_[v]->owner_ != this) {
        message_updates_[v] = message_updates_[v]->Clone();
    }
    else if (!message_updates_[v]) {
        message_updates_[v] = new FGLPVariableUpdate(v, domain_size);
    }
    message_updates_[v]->set_owner(this);

    // If using nullary shift, find the max of each of those messages
    double nullary_shift = ELEM_ONE;

    double *avg_mm = new double[domain_size];
    for (int i = 0; i < domain_size; ++i) avg_mm[i] = ELEM_ONE;

    // Find all max-marginals / accumulate total max-marginal
    for (Function *f : factors_by_variable_[v]) {
        double *mm = MaxMarginal(f, v);
        /*
        for (int i = 0; i < domain_size; ++i) {
            cout << mm[i] << endl;
        }
        cout << endl;
        */
        message_updates_[v]->UpdateFactorMessage(f->getId(), mm);

        for (int i = 0; i < domain_size; ++i) avg_mm[i] OP_TIMESEQ mm[i];

        if (use_nullary_shift_) {
            double m_max = -numeric_limits<double>::max();
            for (int i = 0; i < domain_size; ++i) m_max = max(m_max,mm[i]);
            nullary_shift OP_TIMESEQ m_max;
        }
    }

    // Shift out maximum (if needed) and average
    for (int i = 0; i < domain_size; ++i) {
        if (use_nullary_shift_) avg_mm[i] OP_DIVIDEEQ nullary_shift;
        avg_mm[i] = OP_ROOT(avg_mm[i], factors_by_variable_[v].size());
    }

    if (replace_mode) {
        message_updates_[v]->ReplaceVariableMessage(avg_mm);        
        if (use_nullary_shift_) message_updates_[v]->ReplaceNullaryShift(nullary_shift);
    }
    else {
        message_updates_[v]->UpdateVariableMessage(avg_mm);        
        if (use_nullary_shift_) message_updates_[v]->UpdateNullaryShift(nullary_shift);
    }

    message_queue_.insert(message_updates_[v]->UpdateResidual(), v);
}

void ResidualFGLP::Condition(const vector<Function*> &fns, 
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
            if (new_f->getArity() != f->getArity()) {
                for (int v : new_f->getScopeVec())
                    vars_to_update_.insert(v); 
            }
            factors_.push_back(new_f);
        }
    }
    global_const_factor_ = new FunctionBayes(factors_.back()->getId()+1, 
                problem_, set<int>(), table_const, 1);
    factors_.push_back(global_const_factor_);

}

double ResidualFGLP::UpdateUB() {
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
    else {
        /*
        ub_ = ELEM_ONE;
        for (Function *f : factors_) {
            if (f == global_const_factor_) continue;
            double z = ELEM_ZERO;
            for (size_t i = 0; i < f->getTableSize(); ++i)
                z = max(z,f->getTable()[i]);
            ub_ OP_TIMESEQ z;
        }
        cout << "rest of system = " << ub_ << endl;
        */
        ub_ = global_const_factor_->getTable()[0];
    }
    return ub_ OP_DIVIDE old_ub;
}

void ResidualFGLP::GetVarUB(int var, vector<double> &out) {
    for (Function *f : factors_) {
        if (f->getArity() == 0) continue;
        // If the factor does not have the variable, take the maximum
        if (!f->hasInScope(var)) {
            double z = ELEM_ZERO;
            for (size_t i = 0; i < f->getTableSize(); ++i)
                z = max(z, f->getTable()[i]);
            for (int i = 0; i < problem_->getDomainSize(var); ++i)
                out[i] OP_TIMESEQ z;
            continue;
        }

        if (f->getArity() == 1 && f->hasInScope(var)) continue;

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

double *ResidualFGLP::MaxMarginal(Function *f, int v) {
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

void ResidualFGLP::Reparameterize(Function *f, double *mm, double *avg_mm, 
        int var) {
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
        if (std::isinf(avg_mm[*mm_val]))
            f->getTable()[idx] = ELEM_ZERO;
        else
            f->getTable()[idx] OP_TIMESEQ 
                (avg_mm[*mm_val] OP_DIVIDE mm[*mm_val]);
    } while(increaseTuple(idx, tuple, domains));
    delete [] tuple;
}


