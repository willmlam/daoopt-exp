#include "FGLPMBEHybrid.h"
using namespace std;

FGLPMBEHybrid::FGLPMBEHybrid(Problem *p, Pseudotree *pt, ProgramOptions *po)
    : Heuristic(p,pt,po) {
    // do any original preprocessing first
    if (m_options!=NULL && (m_options->mplp > 0 || m_options->mplps > 0)) {
        FGLP *fglp = new FGLP(m_problem->getNOrg(), m_problem->getDomains(), m_problem->getFunctions(),m_pseudotree->getElimOrder());
        fglp->run(m_options->mplp < 0 ? 5 : m_options->mplp, m_options->mplps);
        m_problem->replaceFunctions(fglp->getFactors());
        m_pseudotree->addFunctionInfo(m_problem->getFunctions()); 
        m_options->mplp = 0;
        m_options->mplps = 0;
    }
    fglpHeur = new FGLPHeuristic(p,pt,po);
    mbeHeur = new MiniBucketElim(p,pt,po,po->ibound);
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
    cout << "FGLP:"<< fglpOut << endl;
    cout << "MBE :"<< mbeOut << endl;

    for (unsigned int i=0; i<out.size(); ++i) {
        if (fglpOut[i] < mbeOut[i]) {
            cout << "(var=" << var << ", val=" << i << ")" << endl;
            cout << "Found " << fglpOut[i] << " < " << mbeOut[i] << endl;
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
