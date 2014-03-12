#ifndef _FGLP_H_
#define _FGLP_H_

#include "Function.h"

class SearchNode;


// FGLP class, holds a reference to the source problem to be reparmaeterized and 
// generates/stores a reparameterized version of the problem

class FGLP {
private:
    static constexpr double DEFAULT_TOLERANCE = 1e-6;
    /*
    // Problem to reparameterize
    Problem *m_problem;

    // Local copy of the problem that will have its factors overwritten
    Problem *m_problemLocal;
    */

    // Domain sizes of variables
    const vector<val_t> &m_domains;

    // Ordering of variables to use
    const vector<int> &m_ordering;


    // Local copy of all of the functions of the problem
    vector<Function*> m_factors;

    // Mapping of variables to factors
    vector<vector<Function*> > m_factorsByVariable;

    // Map to unary factors
    vector<Function*> m_unaryFactors;

    // Pointer to global constant factor
    Function *m_globalConstFactor;

    // Storage of max marginals during tightening
    vector<vector<double*> >m_maxMarginals;

    // Stores the current upper bound of the problem
    // (computed by taking the product of the max values of each factor)
    double m_UB;

    // Stores the upper bound for the non constant portion of the problem
    double m_UBNonConstant;

    // Accumulates the cost for functions already conditioned
    double m_ancestorCost;

    // Use verbose output (show bound progression)
    bool m_verbose;


    // condition the functions in fns according to the assignment and fill them 
    // into m_factors. Also removes factors not in the subproblem
    void condition(const vector<Function*> &fns, const map<int,val_t> &assignment);


    // store the time it took to run
    double m_runtime;

    // store the number of iterations it took
    int m_runiters;

    // Updates UB and returns the amount it changed
    double updateUB();

    // Utility function to compute the max marginal of a function
    double *maxMarginal(Function *f, int v);

    // Utility function to reparameterize a function given its max marginal and average max marginal
    void reparameterize(Function *f, double *maxMarginal, double *averageMaxMarginal, int v);

public:
    // Constructor for reparameterizing the original problem
    FGLP(int nVars, const vector<val_t> &domains, const vector<Function*> &fns, const vector<int> &ordering);

    // Constructor for reparameterizing dynamically during search
    FGLP(int nVars, const vector<val_t> &domains, const vector<Function*> &fns, const vector<int> &ordering, const map<int,val_t> &assignment);

    void run(int maxIter, double maxTime, double tolerance=DEFAULT_TOLERANCE);

    // Get the upper bounds with respect to all assignments to a variable
    void getVarUB(int var, vector<double> &out);

    // Get the upper bound produced directly from reparameterization
    inline double getUB() const { return m_UB; }

    // Get the upper bound produced directly from reparameterization for only 
    // the non constant portion of the problem
    inline double getUBNonConstant() const { return m_UBNonConstant; }

    // Get the constant value (nullary function value)
    inline double getConstant() const { return m_globalConstFactor->getTable()[0]; }

    void getLabelAll(int var, vector<double> &out);

    inline void setVerbose(bool v) { m_verbose = v; }

    inline const vector<Function*> &getFactors() const { return m_factors; }

    inline double getRuntime() const { return m_runtime; }
    inline double getRunIters() const { return m_runiters; }

    size_t getSize() const;

    ~FGLP() {
        for (size_t i = 0; i < m_factors.size(); ++i)
            delete m_factors[i];
    }

};

#endif
