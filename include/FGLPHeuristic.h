#ifndef FGLPHEURISTIC_H_
#define FGLPHEURISTIC_H_

#include "Heuristic.h"
#include "Function.h"
#include "Problem.h"
#include "ProgramOptions.h"
#include "Pseudotree.h"
#include "utils.h"
#include "FGLP.h"
#include "ExtraNodeInfo.h"

// Utility class to contain all the FGLP processed problems
// corresponding to each variable assignment
// Also stores the cost to the node based on the original functions
class FGLPNodeInfo : public ExtraNodeInfo {
    std::vector<FGLP*> fglpStore;
    std::vector<double> origCostToNode;
    public:
    const std::vector<FGLP*> &getFGLPStore() const {
        return fglpStore;
    }
    void addToStore(FGLP *fglp) {
        fglpStore.push_back(fglp);
    }

    const std::vector<double> &getOrigCostToNode() const {
        return origCostToNode;
    }
    void addToCosts(double v) {
        origCostToNode.push_back(v);
    }
    ~FGLPNodeInfo() {
        for (FGLP *f : fglpStore) {
            delete f;
        }
    }
};

class FGLPHeuristic : public Heuristic {
protected:
    double m_globalUB;  
    FGLP *rootFGLP;
    vector<vector<int>> m_ordering;

    map<int,val_t> tempAssn;

    vector<set<int>> m_subproblemFunIds;

    // Cannot do this with AO search
    // Whenever getHeur is called, this is popped until the top of the stack
    // is the parent of the current search node or it is empty
//    std::stack<pair<SearchNode*,vector<FGLP*>>> fglpStore;
//
    void computeSubproblemFunIds();

    
public:
    FGLPHeuristic(Problem *p, Pseudotree *pt, ProgramOptions *po);

    // No memory adjustments for FGLP
    size_t limitSize(size_t limit, const std::vector<val_t> *assignment) { 
        return 0;
    }

    // No memory size for FGLP
    size_t getSize() const { return 0; }

    // Get the initial bound for the problem via FGLP
    size_t build(const std::vector<val_t> *assignment = NULL, bool computeTables = true);

    bool readFromFile(std::string filename) { return false; } 
    bool writeToFile(std::string filename) const { return false; } 

    double getGlobalUB() const { return m_globalUB; }

    double getHeur(int var, const std::vector<val_t> &assignment, SearchNode *node);

    void getHeurAll(int var, const std::vector<val_t> &assignment, SearchNode *node, 
            std::vector<double> &out);

    // Readjusts the heuristic value so it is consistent with the original functions
    // that are already assigned
    void getHeurAllAdjusted(int var, const std::vector<val_t> &assignment, SearchNode *node, std::vector<double> &out);

    double getLabel(int var, const std::vector<val_t> &assignment, SearchNode *node);
    void getLabelAll(int var, const std::vector<val_t> &assignment, SearchNode *node, std::vector<double> &out);

    void printExtraStats() const {}

    virtual ~FGLPHeuristic() { 
        if (rootFGLP)
            delete rootFGLP;
    }

protected:
    void findDfsOrder(vector<int>& ordering, int var) const;


};

#endif
