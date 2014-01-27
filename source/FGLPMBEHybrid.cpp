#include "FGLPMBEHybrid.h"
using namespace std;

FGLPMBEHybrid::FGLPMBEHybrid(Problem *p, Pseudotree *pt, ProgramOptions *po)
    : Heuristic(p,pt,po) {
    fglpHeur = new FGLPHeuristic(new Problem(p),pt,po);
    mbeHeur = new MiniBucketElim(new Problem(p),pt,po,po->ibound);
    timesFGLPUsed = 0;
    timesMBEUsed = 0;
}

size_t FGLPMBEHybrid::build(const std::vector<val_t> *assignment, bool computeTables) {
    fglpHeur->build(assignment,computeTables);
    return mbeHeur->build(assignment,computeTables);
}

double FGLPMBEHybrid::getHeur(int var, const vector<val_t> &assignment, SearchNode *node) {
    // TO DO
    return ELEM_ZERO;
}

void FGLPMBEHybrid::getHeurAll(int var, const vector<val_t> &assignment, SearchNode *node, 
        vector<double> &out) {
    vector<double> fglpOut(out.size(),ELEM_ONE);
    vector<double> mbeOut(out.size(),ELEM_ONE);
    fglpHeur->getHeurAllAdjusted(var,assignment,node,fglpOut);
    mbeHeur->getHeurAll(var,assignment,node,mbeOut);

    for (unsigned int i=0; i<out.size(); ++i) {
        if (fglpOut[i] < mbeOut[i]) {
            out[i] = fglpOut[i];
            timesFGLPUsed++;
        }
        else {
            out[i] = mbeOut[i];
            timesMBEUsed++;
        }
    }
}

double FGLPMBEHybrid::getLabel(int var, const vector<val_t> &assignment, SearchNode *node) {
    double d = ELEM_ONE;
    for (Function *f : m_pseudotree->getFunctions(var)) {
        d OP_TIMESEQ f->getValue(assignment);
    }
    return d;
}

void FGLPMBEHybrid::getLabelAll(int var, const vector<val_t> &assignment, SearchNode *node,
        vector<double> &out) {
    vector<double> costTmp(m_problem->getDomainSize(var), ELEM_ONE);
    for (Function *f : m_pseudotree->getFunctions(var)) {
        f->getValues(assignment, var, costTmp);
        for (int i=0; i<m_problem->getDomainSize(var); ++i) {
            out[i] OP_TIMESEQ costTmp[i];
        }
    }
}
