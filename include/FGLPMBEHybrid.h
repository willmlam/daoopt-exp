#ifndef FGLPMBEHYBRID_H_
#define FGLPMBEHYBRID_H_

#include "Heuristic.h"
#include "FGLPHeuristic.h"
#include "MiniBucketElim.h"
#include "Problem.h"
#include "ProgramOptions.h"
#include "Pseudotree.h"
#include "utils.h"
#include "mex/mbe.h"
#include "FGLP.h"
#include "ExtraNodeInfo.h"

// A wrapper to take the tightest heuristics between FGLP and MBE

class FGLPMBEHybrid : public Heuristic {
    FGLPHeuristic *fglpHeur;
    MiniBucketElim *mbeHeur;

    vector<unsigned long> timesFGLPUsed;
    vector<unsigned long> timesMBEUsed;

    vector<unsigned long> timesFGLPPruned;
    vector<unsigned long> timesMBEPruned;
    vector<unsigned long> timesBothPruned;

public:
    FGLPMBEHybrid(Problem *p, Pseudotree *pt, ProgramOptions *po); 

    inline size_t limitSize(size_t limit, const std::vector<val_t> *assignment) {
        if (m_options->fglpMBEHeur)
            return mbeHeur->limitSize(limit,assignment);
        else
            return 0;
    }

    inline size_t getSize() const { 
        if (m_options->fglpMBEHeur) 
            return mbeHeur->getSize(); 
        else
            return 0;
    }

    size_t build(const std::vector<val_t> *assignment = NULL, bool computeTables = true);
    bool readFromFile(std::string filename) { return mbeHeur->readFromFile(filename); } 
    bool writeToFile(std::string filename) const { return mbeHeur->writeToFile(filename); } 

    inline double getGlobalUB() const { return min(fglpHeur->getGlobalUB(), (mbeHeur ? mbeHeur->getGlobalUB() : ELEM_ONE)); }

    double getHeur(int var, const std::vector<val_t> &assignment, SearchNode *node);

    void getHeurAll(int var, const std::vector<val_t> &assignment, SearchNode *node, 
            std::vector<double> &out);


    double getLabel(int var, const std::vector<val_t> &assignment, SearchNode *node);
    void getLabelAll(int var, const std::vector<val_t> &assignment, SearchNode *node, std::vector<double> &out);

    bool calculatePruning(int var, SearchNode *node, double curPSTVal);

    mex::vector<mex::Factor> copyFactors();
    void rewriteFactors( const vector<mex::Factor> &factors);

    inline bool isAccurate() {
        assert(m_pseudotree);
        return mbeHeur && (m_pseudotree->getWidthCond() <= mbeHeur->getIbound());
    }

    size_t computeMBEMemory(int ibound);

    size_t limitJGLPIBound(size_t memlimit);

    inline void printExtraStats() const {
        cout << "depth,FGLPBetter,MBEBetter,OnlyFGLPPruned,OnlyMBEPruned,BothPruned" << endl;
        for (size_t i=0; i<timesFGLPUsed.size(); ++i) {
            cout << i << "," 
                << timesFGLPUsed[i] << "," 
                << timesMBEUsed[i] << ","
                << timesFGLPPruned[i] << ","
                << timesMBEPruned[i] << ","
                << timesBothPruned[i] << endl;
        }
    }

    inline virtual ~FGLPMBEHybrid() { 
        if (fglpHeur)
            delete fglpHeur;
        if (mbeHeur)
            delete mbeHeur;
    }
};

#endif
