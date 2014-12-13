#ifndef FGLPMBEHYBRID_H_
#define FGLPMBEHYBRID_H_

#include "ExtraNodeInfo.h"
#include "FGLP.h"
#include "FGLPHeuristic.h"
#include "Heuristic.h"
#include "MiniBucketElim.h"
#include "Problem.h"
#include "ProgramOptions.h"
#include "Pseudotree.h"
#include "utils.h"
#include "mex/mbe.h"

// A wrapper to take the tightest heuristics between FGLP and MBE

namespace daoopt {

class FGLPMBEHybrid : public Heuristic {
  std::unique_ptr<FGLPHeuristic> fglpHeur;
  std::unique_ptr<MiniBucketElim> mbeHeur;

  std::unique_ptr<FGLP> fglp;

  unordered_map<int32, uint32> timesFGLPUsed;
  unordered_map<int32, uint32> timesMBEUsed;

  unordered_map<int32, uint32> timesFGLPPruned;
  unordered_map<int32, uint32> timesMBEPruned;
  unordered_map<int32, uint32> timesBothPruned;


public:
    FGLPMBEHybrid(Problem *p, Pseudotree *pt, ProgramOptions *po); 

    inline virtual ~FGLPMBEHybrid() { }

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

    double getHeur(int var, std::vector<val_t> &assignment, SearchNode *node);

    void getHeurAll(int var, std::vector<val_t> &assignment, SearchNode *node, 
            std::vector<double> &out);


    double getLabel(int var, const std::vector<val_t> &assignment, SearchNode *node);
    void getLabelAll(int var, const std::vector<val_t> &assignment, SearchNode *node, std::vector<double> &out);

    bool calculatePruning(int var, SearchNode *node, double curPSTVal);

    mex::vector<mex::Factor> copyFactors();
    void rewriteFactors( const vector<mex::Factor> &factors);

    inline bool isAccurate() {
        assert(m_pseudotree);
        return false;
//        return mbeHeur && (m_pseudotree->getWidthCond() <= mbeHeur->getIbound());
    }

    size_t computeMBEMemory(int ibound);

    size_t limitJGLPIBound(size_t memlimit);

    inline void printExtraStats() const {
        cout << "depth,FGLPBetter,MBEBetter,OnlyFGLPPruned,OnlyMBEPruned,BothPruned" << endl;
        for (int32 i = -1; i <= m_pseudotree->getHeight(); ++i) {
            cout << i << "," 
                << timesFGLPUsed.at(i) << "," 
                << timesMBEUsed.at(i) << ","
                << timesFGLPPruned.at(i) << ","
                << timesMBEPruned.at(i) << ","
                << timesBothPruned.at(i) << endl;
        }

        unsigned long long totalIterationsRun = fglpHeur->getTotalIterationsRun();
        unsigned long long totalInitiated = fglpHeur->getTotalInitiated();

        const vector<unsigned long> &countVars = fglpHeur->getCountVars();
        const vector<unsigned long> &varsUpdated = fglpHeur->getVarsUpdated();

        cout << "Total iterations run: " << totalIterationsRun << endl;
        cout << "Total initiated: " << totalInitiated << endl;

        cout << "Average iterations per node: " << totalIterationsRun / double(totalInitiated) << endl;

        cout << "Average variables updated at each variable: " << endl;
        cout << "var,depth,average" << endl;
        for (int i = 0; i < m_problem->getN(); ++i) {
            cout << i << "," << m_pseudotree->getNode(i)->getDepth() << ","
                << varsUpdated[i]/double(countVars[i]) << endl;
        }
    }

};

}  // namespace daoopt

#endif
