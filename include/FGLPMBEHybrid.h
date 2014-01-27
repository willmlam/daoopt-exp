#ifndef FGLPMBEHYBRID_H_
#define FGLPMBEHYBRID_H_

#include "Heuristic.h"
#include "FGLPHeuristic.h"
#include "MiniBucketElim.h"
#include "Problem.h"
#include "ProgramOptions.h"
#include "Pseudotree.h"
#include "utils.h"
#include "FGLP.h"
#include "ExtraNodeInfo.h"

// A wrapper to take the tightest heuristics between FGLP and MBE

class FGLPMBEHybrid : public Heuristic {
    FGLPHeuristic *fglpHeur;
    MiniBucketElim *mbeHeur;

    unsigned long timesFGLPUsed;
    unsigned long timesMBEUsed;

public:
    FGLPMBEHybrid(Problem *p, Pseudotree *pt, ProgramOptions *po); 

    size_t limitSize(size_t limit, const std::vector<val_t> *assignment) {
        return mbeHeur->limitSize(limit,assignment);
    }

    size_t getSize() const { return mbeHeur->getSize(); }

    size_t build(const std::vector<val_t> *assignment = NULL, bool computeTables = true);
    bool readFromFile(std::string filename) { return mbeHeur->readFromFile(filename); } 
    bool writeToFile(std::string filename) const { return mbeHeur->writeToFile(filename); } 

    double getGlobalUB() const { return min(fglpHeur->getGlobalUB(), mbeHeur->getGlobalUB()); }

    double getHeur(int var, const std::vector<val_t> &assignment, SearchNode *node);

    void getHeurAll(int var, const std::vector<val_t> &assignment, SearchNode *node, 
            std::vector<double> &out);


    double getLabel(int var, const std::vector<val_t> &assignment, SearchNode *node);
    void getLabelAll(int var, const std::vector<val_t> &assignment, SearchNode *node, std::vector<double> &out);

    void printExtraStats() const {
        cout << "Times FGLP better:" << timesFGLPUsed << endl;
        cout << "Times MBE better:" << timesMBEUsed << endl;
    }

    virtual ~FGLPMBEHybrid() { 
        if (fglpHeur)
            delete fglpHeur;
        if (mbeHeur)
            delete mbeHeur;
    }
};

#endif
