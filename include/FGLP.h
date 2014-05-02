#ifndef _FGLP_H_
#define _FGLP_H_

#include "Function.h"

class SearchNode;


// FGLP class, holds a reference to the source problem to be reparmaeterized and 
// generates/stores a reparameterized version of the problem

class FGLP {
protected:
    static constexpr double DEFAULT_TOLERANCE = 1e-7;
//    static constexpr int DEFAULT_UPDATE_DISTANCE = std::numeric_limits<int>::max();
    /*
    // Problem to reparameterize
    Problem *m_problem;

    // Local copy of the problem that will have its factors overwritten
    Problem *m_problemLocal;
    */

    // Use version of update that shifts the max-marginals into the nullary function?
    bool m_useNullaryShift;

    // Number of variables
    int m_nVars;

    // Domain sizes of variables
    const vector<val_t> &m_domains;

    // Local copy of all of the functions of the problem
    vector<Function*> m_factors;

    // Does this class own its factors?
    bool m_ownsFactors;

    // Mapping of variables to factors
    vector<vector<Function*> > m_factorsByVariable;

    // Map to unary factors
    vector<Function*> m_unaryFactors;

    // Pointer to global constant factor
    Function *m_globalConstFactor;


    // Store the variables that must be updated
    // (Conditioning definitely affects these variables in general)
    // (only needed for first iteration?) (rest can do zero check?)
    vector<bool> m_forceUpdateVars;

    // Stores the current upper bound of the problem
    // (computed by taking the product of the max values of each factor)
    double m_UB;

    // Stores the upper bound for the non constant portion of the problem
    double m_UBNonConstant;

    // Use verbose output (show bound progression)
    bool m_verbose;

    // store the time it took to run
    double m_runtime;

    // store the number of iterations it took
    int m_runiters;

    // condition the functions in fns according to the assignment and fill them 
    // into m_factors. Also removes factors not in the subproblem
    void condition(const vector<Function*> &fns, const map<int,val_t> &assignment);

    // Updates UB and returns the amount it changed
    // Not used when using nullary shift version, which updates the UB on the fly
    // Simply sets the UB value instead.
    double updateUB();

    // Utility function to compute the max marginal of a function
    double *maxMarginal(Function *f, int v);

    // Utility function to reparameterize a function given its max marginal and average max marginal
    void reparameterize(Function *f, double *maxMarginal, double *averageMaxMarginal, int v);


private:
    // Update ordering for each iteration
    const vector<int> &m_updateOrdering;

    // Storage of max marginals during tightening
    vector<vector<double*> >m_maxMarginals;

    // Skip variable updates if zero?
    bool m_skipVars;

    // Flags for each variable to determine if all of their functions are zero
    vector<bool> m_varFactorMaxIsZero;

    // Goes through each function to check if the function maximum is 0
    // if not, mark it in the vector
    void updateVarFactorMaxIsZero();

public:
    // Constructor for reparameterizing the original problem
    FGLP(int nVars, const vector<val_t> &domains, const vector<Function*> &fns, const vector<int> &ordering, bool useNullaryShift = false);

    // Constructor for reparameterizing dynamically during search
    FGLP(int nVars, const vector<val_t> &domains, const vector<Function*> &fns, const vector<int> &ordering, const map<int,val_t> &assignment, const vector<bool> &zeroFlags, bool useNullaryShift = false);

    virtual void run(int maxIter, double maxTime, double tolerance=DEFAULT_TOLERANCE);

    // Get the upper bounds with respect to all assignments to a variable
    void getVarUB(int var, vector<double> &out);

    // Get the upper bound produced directly from reparameterization
    inline double getUB() const { return m_UB; }

    // Get the upper bound produced directly from reparameterization for only 
    // the non constant portion of the problem
    inline double getUBNonConstant() const { return m_UBNonConstant; }

    // Get the constant value (nullary function value)
    inline double getConstant() const { return m_globalConstFactor->getTable()[0]; }

    // Get the flags stating whether all of the functions containing a variable
    // have maximums at 0.
    inline const vector<bool> &getVarFactorMaxIsZero() const {
        return m_varFactorMaxIsZero;
    }
    void getLabelAll(int var, vector<double> &out);

    /*
    inline bool addToUpdateOrdering(int v) {
        if (!m_inUpdateOrdering[v]) {
            m_updateOrdering.push_back(v);
            m_inUpdateOrdering[v] = true;
            return true;
        }
        return false;
    }
    */
//    void populateOrdering(int dist=DEFAULT_UPDATE_DISTANCE);

    inline void setVerbose(bool v) { m_verbose = v; }
    inline void setOwnsFactors(bool o) { m_ownsFactors = o; }

    inline const vector<Function*> &getFactors() const { return m_factors; }

    inline double getRuntime() const { return m_runtime; }
    inline double getRunIters() const { return m_runiters; }

    size_t getSize() const;

    ~FGLP() {
        if (m_ownsFactors) {
            for (size_t i = 0; i < m_factors.size(); ++i)
                delete m_factors[i];
        }
    }

};

#endif
