#ifndef _FGLP_H_
#define _FGLP_H_

#include "Problem.h"
#include "Function.h"


// FGLP class, holds a reference to the source problem to be reparmaeterized and 
// generates/stores a reparameterized version of the problem

class FGLP {
public:
    static constexpr double DEFAULT_TOLERANCE = 1e-7;

    // Default constructor: does nothing
    FGLP() { }

    // Constructor for reparameterizing the original problem
    FGLP(Problem *p, bool use_nullary_shift = false);

    // Constructor for reparameterizing dynamically during search
    // Additionally takes in the set of variables
    FGLP(FGLP *parent_fglp, const map<int,val_t> &assignment,
         const set<int> &subVars, int root_var);

    virtual void Run(int max_iter, double max_time, double tolerance=DEFAULT_TOLERANCE);

    // Get the upper bounds with respect to all assignments to a variable
    void GetVarUB(int var, vector<double> &out);

    // Get the upper bound produced directly from reparameterization
    inline double ub() const { return ub_; }

    // Get the constant value (nullary function value)
    inline double global_constant() const { 
        return global_const_factor_->getTable()[0]; 
    }

    void GetLabelAll(int var, vector<double> &out);

    inline void set_verbose(bool v) { verbose_ = v; }
    inline void set_owns_factors(bool o) { owns_factors_ = o; }

    inline const vector<Function*> &factors() const { return factors_; }

    inline const set<int> vars_updated() const { return vars_updated_; }

    inline double runtime() const { return runtime_; }
    inline double runiters() const { return runiters_; }

    inline void set_use_cost_shift_reversal(bool r) { 
      use_cost_shift_reversal_ = r;
    }
    inline const vector<double> &bound_contribs() const {
      return bound_contribs_;
    }

    inline const vector<map<int,vector<double>>> &total_shift() const {
      return total_shift_;
    }

    size_t GetSize() const;

    virtual ~FGLP() {
        if (owns_factors_) {
            for (Function *f : factors_) delete f;
        }
    }
protected:
    // condition the functions in fns according to the assignment and fill them 
    // into m_factors. Also removes factors not in the subproblem
    virtual void Condition(const vector<Function*> &fns, 
            const map<int,val_t> &assignment, const set<int> &subVars,
            int condition_var);

    // Updates UB and returns the amount it changed
    // Not used when using nullary shift version, which updates the UB on the fly
    // Simply sets the UB value instead.
    double UpdateUB();

    // Utility function to compute the max marginal of a function
    double *MaxMarginal(Function *f, int v);

    // Utility function to reparameterize a function given its max marginal and average max marginal
    void Reparameterize(Function *f, double *mm, double *avg_mm, int v);

    // (For debugging) Prints out all of the factors
    virtual void PrintAllFactors() const;

    // Problem to reparameterize
    Problem *problem_;

    // Does this class own its factors?
    bool owns_factors_;

    // Local copy of all of the functions of the problem
    vector<Function*> factors_;

    // Mapping of variables to factors
    vector<vector<Function*> > factors_by_variable_;

    // Pointer to global constant factor
    Function *global_const_factor_;

    // vector of individual bound contributions from each variable
    vector<double> bound_contribs_;

    // storage of total cost shift across each function
    vector<map<int,vector<double>>> total_shift_;

    // Storage of max marginals during tightening
    vector<vector<double*> > max_marginals_;

    set<int> vars_updated_;

    // Stores the portion of the nullary function from conditioning
    double conditioned_cost_;

    // Stores the current upper bound of the problem
    // (computed by taking the product of the max values of each factor)
    // (or if using nullary shift, taking the nullary function value)
    double ub_;

    // Use version of update that shifts the max-marginals into the nullary function?
    bool use_nullary_shift_;

    // Use cost-shift reversal (required for AOBB with caching)
    bool use_cost_shift_reversal_;

    // Use verbose output (show bound progression)
    bool verbose_;

    // store the time it took to run
    double runtime_;

    // store the number of iterations it took
    int runiters_;





};

#endif
