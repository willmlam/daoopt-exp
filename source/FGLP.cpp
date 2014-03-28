#include "FGLP.h"
#include <chrono>

using namespace std;
using namespace std::chrono;

FGLP::FGLP(int nVars, const vector<val_t> &domains, const vector<Function*> &fns, const vector<int> &ordering) 
    :  
    m_domains(domains),
    m_updateOrdering(ordering),
    m_factorsByVariable(nVars, vector<Function*>()), 
    m_unaryFactors(nVars, NULL),
    m_globalConstFactor(NULL),
    m_maxMarginals(nVars, vector<double*>()),
    m_UB(ELEM_ONE),
    m_UBNonConstant(ELEM_ONE),
    m_ancestorCost(ELEM_ONE),
    m_verbose(true) {


    // copy factors and create a mapping of variables to factors
    //for (int i=0; i< fns.size(); ++i) cout << *fns[i] << endl;


    for (size_t i = 0; i < fns.size(); ++i) {

        // Check if factor exists in this subproblem
        bool factorExists = false;
        vector<int>::const_iterator it = m_updateOrdering.begin();
        
        for (; it != m_updateOrdering.end(); ++it) {
            if (fns[i]->hasInScope(*it)) factorExists = true;
        }
        if (fns[i]->getArity() == 0) m_globalConstFactor = fns[i]->clone();
        
        if (!factorExists) {
            continue;
        }

        m_factors.push_back(fns[i]->clone());
        const vector<int> &scope = m_factors.back()->getScopeVec();
        vector<int>::const_iterator itS = scope.begin();
        for (; itS != scope.end(); ++itS) {
            m_factorsByVariable[*itS].push_back(m_factors.back());    
            if (scope.size() == 1) m_unaryFactors[*itS] = m_factors.back();
        }
    }
    if (m_globalConstFactor != NULL)
        m_factors.push_back(m_globalConstFactor);

    // allocate max marginal storage
    for (size_t v = 0; v < m_factorsByVariable.size(); ++v) {
        m_maxMarginals[v].resize(m_factorsByVariable[v].size(), NULL);
    }

//    addToUpdateOrdering(m_ordering[1]);
//    populateOrdering();
}

FGLP::FGLP(int nVars, const vector<val_t> &domains, const vector<Function*> &fns, const vector<int> &ordering, const map<int,val_t> &assignment) 
    : 
    m_domains(domains),
    m_updateOrdering(ordering), 
    m_factorsByVariable(nVars, vector<Function*>()), 
    m_unaryFactors(nVars, NULL),
    m_globalConstFactor(NULL),
    m_maxMarginals(nVars, vector<double*>()),
    m_UB(ELEM_ONE),
    m_ancestorCost(ELEM_ONE),
    m_verbose(false) {

    condition(fns,assignment);


    for (size_t i = 0; i < m_factors.size(); ++i) {

//        m_factors[i] = fns[i]->clone();
        const vector<int> &scope = m_factors[i]->getScopeVec();
        vector<int>::const_iterator itS = scope.begin();
        for (; itS != scope.end(); ++itS) {
            m_factorsByVariable[*itS].push_back(m_factors[i]);    
            if (scope.size() == 1) m_unaryFactors[*itS] = m_factors[i];
        }
    }

    // allocate max marginal storage
    for (size_t v = 0; v < m_factorsByVariable.size(); ++v) {
        m_maxMarginals[v].resize(m_factorsByVariable[v].size(), NULL);
    }
//    populateOrdering();

}

double FGLP::updateUB() {
    double oldUB = m_UB;
    m_UB = ELEM_ONE;
    m_UBNonConstant = ELEM_ONE;
    vector<Function*>::const_iterator itF = m_factors.begin();
    for (; itF != m_factors.end(); ++itF) {
       // Find the maximum value of the factor
       double z = ELEM_ZERO;
       for (size_t i=0; i<(*itF)->getTableSize(); ++i) {
           z = max(z,(*itF)->getTable()[i]);
       }
       m_UB OP_TIMESEQ z;
       if ((*itF)->getArity() > 0) m_UBNonConstant OP_TIMESEQ z;
    }
    return m_UB OP_DIVIDE oldUB;
}

void FGLP::getVarUB(int var, vector<double> &out) {
    vector<Function*>::const_iterator itF = m_factors.begin();
    set<int>::const_iterator itS;
    int i;

//    cout << endl;
    for (; itF != m_factors.end(); ++itF) {
        if ((*itF)->getArity() == 0) continue;
        // If var is not in the factor, take the maximum value
        if (!(*itF)->hasInScope(var)) {
            // Find the maximum value of the factor
            double z = ELEM_ZERO;
            for (size_t i=0; i<(*itF)->getTableSize(); ++i) {
                z = max(z,(*itF)->getTable()[i]);
            }
            for (int i=0; i<m_domains[var]; ++i) {
                out[i] OP_TIMESEQ z;
            }
            continue;
        }

        // Skip if unary factor with the current variable scope (would be pushed to the global constant in actual conditioning)
        if ((*itF)->getArity() == 1 && (*itF)->hasInScope(var)) {
            /*
            cout << " unary [ ";
            for (size_t k=0; k<(*itF)->getTableSize(); ++k) {
                cout << (*itF)->getTable()[k] << " ";
            }
            cout << "]" << endl;
            */

            continue;
        }

//        cout << **itF << endl;
        

        // Generate a tuple to iterate over the factor
        val_t *tuple = new val_t[(*itF)->getArity()];
        val_t *unfixedTuple = tuple + 1;

        for (i=0;i<(*itF)->getArity();++i) tuple[i] = 0;

        vector<val_t> unfixedDomains;
        for (itS = (*itF)->getScopeSet().begin(); itS != (*itF)->getScopeSet().end(); ++itS) {
            if (*itS != var) unfixedDomains.push_back(m_domains[*itS]);    
        }

        vector<val_t> varDomain;
        varDomain.push_back(m_domains[var]);
        
        vector<val_t*> idxMap;
        for (i=0, itS=(*itF)->getScopeSet().begin(); itS != (*itF)->getScopeSet().end(); ++itS) {
            if (*itS == var) 
                idxMap.push_back(&tuple[0]);
            else 
                idxMap.push_back(&unfixedTuple[i++]);
        }
        
        size_t idx=0;
        do {
            // Find the best configuration for current varTuple
            double z = ELEM_ZERO;
            do {
                z = max(z,(*itF)->getValuePtr(idxMap));
            } while (increaseTuple(idx,unfixedTuple,unfixedDomains));
//            cout << z << endl;
            out[tuple[0]] OP_TIMESEQ z;
//            cout << "Out tuple " << int(tuple[0]) << " " << out[tuple[0]] << endl;
        } while (increaseTuple(idx,tuple,varDomain));
        delete [] tuple;
    }

    if(m_globalConstFactor) {
//        cout << "GConst: " << m_globalConstFactor->getTable()[0] << endl;
    }
}

void FGLP::run(int maxIter, double maxTime, double tolerance) {

    double diff = updateUB();
    if (m_verbose) {
        cout << "Initial UB: " << m_UB << endl;
        /*
        for(Function *f : m_factors) {
            cout << *f << endl;
            cout << f->getTableSize() << endl;
            cout << " ";
            for (size_t i=0; i<f->getTableSize(); ++i) {
                cout << f->getTable()[i] << " ";
            }
            cout << endl;
        }
        */
    }
    steady_clock::time_point timeStart = steady_clock::now();
    steady_clock::time_point timeEnd;

    int iter;
    for (iter = 0; iter < maxIter || (maxIter == -1 && maxTime > 0); ++iter) {
        timeEnd = steady_clock::now();
        if (maxTime > 0 && duration_cast<seconds>(timeEnd-timeStart).count() >= maxTime) break;
        if (fabs(diff) < tolerance || std::isnan(diff)) break;
//        vector<int>::const_reverse_iterator rit = m_ordering.rbegin();
//        m_updateOrdering = m_ordering;
//        cout << m_updateOrdering << endl;
        vector<int>::const_iterator it = m_updateOrdering.begin();
        for (; it != m_updateOrdering.end(); ++it) {
//        for (; rit != m_ordering.rend(); ++rit) {
            timeEnd = steady_clock::now();
            if (maxTime > 0 && duration_cast<seconds>(timeEnd-timeStart).count() >= maxTime) break;
            // update variable v this iteration
            int v = *it;

            if (v >= int(m_factorsByVariable.size())) continue;

            // Skip this variable, no factors have this variable
            if (m_factorsByVariable[v].size() == 0) continue;

            // *May not be needed*
            /*
            // Check if a unary factor exists for it and create it if it does not
            if (!m_unaryFactors[v]) {
                int fid = m_factors.size();
                set<int> scope;
                scope.insert(v);
                double *newTable = new double[m_domains[v]];
                m_factors.push_back(
                        new FunctionBayes(m_factors.size(),
                            m_problem,scope,newTable,m_problem->getDomainSize(v)));
                m_unaryFactors[v] = m_factors.back();
            }
            */

            // To store max marginals of each of the functions

            vector<double*> &varMaxMarginals = m_maxMarginals[v];

            // Compute max marginals
            for (size_t i = 0; i < m_factorsByVariable[v].size(); ++i) {
                varMaxMarginals[i] = maxMarginal(m_factorsByVariable[v][i], v);
            }

            // Compute average max marginal
            size_t tableSize = m_domains[v];
            double *avgMaxMarginalTable = new double[tableSize];

            for (unsigned int i=0; i<tableSize; ++i) avgMaxMarginalTable[i] = ELEM_ONE;
            for (vector<double*>::iterator itMM=varMaxMarginals.begin();
                    itMM!=varMaxMarginals.end(); ++itMM) {
                for (unsigned int i=0; i<tableSize; ++i)
                    avgMaxMarginalTable[i] OP_TIMESEQ (*itMM)[i];
            }
            for (unsigned int i=0; i<tableSize; ++i) {
                avgMaxMarginalTable[i] = OP_ROOT(avgMaxMarginalTable[i],varMaxMarginals.size());
            }
//            for (int i=0; i<int(tableSize);++i) cout << avgMaxMarginalTable[i] << endl;

            // Reparameterize
            for (size_t i=0; i<m_factorsByVariable[v].size(); ++i) {
                reparameterize(m_factorsByVariable[v][i],varMaxMarginals[i],avgMaxMarginalTable,v);
                delete varMaxMarginals[i];
                varMaxMarginals[i] = NULL;
            }
            delete avgMaxMarginalTable;
        }
        diff = updateUB();
        if (m_verbose) {
            cout << "UB: " << m_UB << " (d=" << diff << ")" << endl;
            /*
            cout << "constant: " << m_globalConstFactor->getTable()[0] << endl;
            cout << "non-const UB: " << m_UBNonConstant << endl;
            */
        }
    }
    timeEnd = steady_clock::now();
    m_runtime = double(duration_cast<milliseconds>(timeEnd-timeStart).count()) / 1000;
    m_runiters = iter;
    if (m_verbose) {
        cout << "FGLP (" << iter << " iter, " << m_runtime << " sec): " << m_UB << endl;
        /*
        for(Function *f : m_factors) {
            cout << *f << endl;
            cout << f->getTableSize() << endl;
            cout << " ";
            for (size_t i=0; i<f->getTableSize(); ++i) {
                cout << f->getTable()[i] << " ";
            }
            cout << endl;
        }
        */
    }

}
    
double *FGLP::maxMarginal(Function *f, int v) {
    assert(f->hasInScope(v));
    size_t tableSize = m_domains[v];
    double *newTable = new double[tableSize];
    for (size_t i=0; i<tableSize; ++i) newTable[i] = ELEM_ZERO;

    // Tuple corresponds to the same vars, but v is moved to the front
    val_t *tuple = new val_t[f->getScopeSet().size()];
    val_t *elimTuple = tuple + 1;

    set<int>::const_iterator itS;

    vector<val_t> elimDomains;
    for (itS = f->getScopeSet().begin(); itS != f->getScopeSet().end(); ++itS) {
        if (*itS != v) elimDomains.push_back(m_domains[*itS]);    
    }

    vector<val_t> varDomain;
    varDomain.push_back(m_domains[v]);

    for (int i=0; i<f->getArity(); ++i) tuple[i] = 0;

    vector<val_t*> idxMap;

    int i;
    for (i=0, itS = f->getScopeSet().begin(); itS != f->getScopeSet().end(); ++itS) {
        if (*itS == v) 
            idxMap.push_back(&tuple[0]);
        else
            idxMap.push_back(&elimTuple[i++]);
    }

    size_t idx;
    // iterate over all values non-v variables
    do {
        idx=0;
        do {
            //cout << idx << endl;
            newTable[idx] = max(newTable[idx],f->getValuePtr(idxMap));
        } while (increaseTuple(idx,tuple,varDomain));
    } while (increaseTuple(idx,elimTuple,elimDomains));

    delete [] tuple;

//    for (int i=0; i<int(tableSize);++i) cout << newTable[i] << endl;
    return newTable;
}

void FGLP::reparameterize(Function *f, double *maxMarginal, double *averageMaxMarginal, int v) {
    assert(f->hasInScope(v));
    vector<val_t> domains;
    domains.reserve(f->getArity());
    set<int>::const_iterator itS;
    int i;

    val_t *tuple = new val_t[f->getArity()];
    for (i=0;i<f->getArity();++i) tuple[i] = 0;

    val_t *mmVal = NULL;

    for (i=0,itS=f->getScopeSet().begin(); itS!=f->getScopeSet().end(); ++i,++itS) {
        domains.push_back(m_domains[*itS]);
        if (*itS == v) mmVal = &tuple[i];
    }

    size_t idx = 0;
    do {
        // if the max-marginal value passed is -inf (0), then this is a hard constraint
        if (std::isinf(averageMaxMarginal[*mmVal]))
            f->getTable()[idx] = ELEM_ZERO;
        else
            f->getTable()[idx] OP_TIMESEQ (averageMaxMarginal[*mmVal] OP_DIVIDE maxMarginal[*mmVal]);
    } while(increaseTuple(idx,tuple,domains));
    delete [] tuple;

}

void FGLP::condition(const vector<Function*> &fns, const map<int,val_t> &assignment) {

    double globalConstant = ELEM_ONE;
    int arityZeroFnIdx = -1;
//    cout << assignment << endl;
    for (size_t i = 0; i < fns.size(); ++i) {
        if (fns[i]->getArity() == 0) {
            globalConstant OP_TIMESEQ fns[i]->getTable()[0];
            arityZeroFnIdx = i;
            continue;
        }

        /*
        // Check if factor exists in this subproblem
        bool factorExists = false;
        vector<int>::const_iterator it = m_ordering.begin();
        for (; it!=m_ordering.end(); ++it) {
            if (fns[i]->hasInScope(*it)) factorExists = true;
        }
        
        if (!factorExists) {
            continue;
        }
        */

        /*
        for (int j=0; j<fns[i]->getTableSize(); ++j) {
            cout << " " << fns[i]->getTable()[j] << endl;
        }
        cout << endl;
        */
        Function *new_fn = fns[i]->substitute(assignment);

        // Function changed, place other variables in update ordering
        /*
        if (new_fn->getArity() != fns[i]->getArity()) {
//            cout << "new arity: " << new_fn->getArity() << " , old arity:" << fns[i]->getArity() << endl;
            for (int v : new_fn->getScopeVec()) {
//                cout << "Adding to initial update ordering: " << v << endl;
                addToUpdateOrdering(v);
            }
        }
        */

        if (new_fn->isConstant()) {
//            cout << "Changed function is constant (nothing to add)" << endl;
            globalConstant OP_TIMESEQ new_fn->getTable()[0];
            delete new_fn;
        } else {
            m_factors.push_back(new_fn);
        }
    }

    Function *constFun = fns[arityZeroFnIdx]->clone();
    constFun->getTable()[0] = globalConstant;
    m_globalConstFactor = constFun;
    m_factors.push_back(m_globalConstFactor);
}

void FGLP::getLabelAll(int var, vector<double> &out) {
    out.clear();
    out.resize(m_domains[var],ELEM_ONE);
    // For each function, check if its unary and contains the variable
    for (const auto &f : m_factors) {
        if (f->getArity() == 1 && f->hasInScope(var)) {
            for (size_t i=0; i < out.size(); ++i) {
                out[i] OP_TIMESEQ f->getTable()[i];
            }
        }
    }
}

/*
void FGLP::populateOrdering(int dist) {
    queue<int> newVars;
    for (int k : m_updateOrdering) {
        newVars.push(k);
    }
    int firstAdded;
    for (int k=0; k<dist; ++k) {
//        cout << k << endl;
        if (newVars.empty()) {
            break;
        }
        firstAdded = NONE;
        while(!newVars.empty() && newVars.front() != firstAdded) {
            int currentVar = newVars.front();
            newVars.pop();
            for (Function *f : m_factorsByVariable[currentVar]) {
                for (int k : f->getScopeVec()) {
                    if (addToUpdateOrdering(k)) {
                        newVars.push(k);
                        if (firstAdded == NONE) firstAdded = k;
                    }
                }
            }
        }
    }
}
*/

size_t FGLP::getSize() const {
    size_t S = 0;
    for (const auto &f : m_factors) {
        S += f->getTableSize();
    }
    return S;
}
