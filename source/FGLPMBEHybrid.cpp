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
    if (m_options->fglpMBEHeur)
        mbeHeur = new MiniBucketElim(p,pt,po,po->ibound);
    else
        mbeHeur = nullptr;
    timesFGLPUsed.resize(p->getN(),0);
    timesMBEUsed.resize(p->getN(),0);
    timesFGLPPruned.resize(p->getN(),0);
    timesMBEPruned.resize(p->getN(),0);
    timesBothPruned.resize(p->getN(),0);
}

size_t FGLPMBEHybrid::build(const std::vector<val_t> *assignment, bool computeTables) {
    fglpHeur->build(assignment,computeTables);
    if (m_options->fglpMBEHeur)
        return mbeHeur->build(assignment,computeTables);
    else
        return 0;
}

double FGLPMBEHybrid::getHeur(int var, const vector<val_t> &assignment, SearchNode *node) {
    // TO DO
    return ELEM_ZERO;
}

void FGLPMBEHybrid::getHeurAll(int var, const vector<val_t> &assignment, SearchNode *node, 
        vector<double> &out) {
    vector<double> fglpOut(out.size(),ELEM_ONE);
    vector<double> mbeOut(out.size(),ELEM_ONE);
    if (m_options->useShiftedLabels) {
        fglpHeur->getHeurAll(var,assignment,node,fglpOut);
    }
    else {
        fglpHeur->getHeurAllAdjusted(var,assignment,node,fglpOut);
    }

    if (m_options->fglpMBEHeur)
        mbeHeur->getHeurAll(var,assignment,node,mbeOut);

//    cout << "FGLP:"<< fglpOut << endl;
//    cout << "MBE :"<< mbeOut << endl;

    // Count the number of possible prunings for each heuristics wrt to the 
    // values (only if the other heuristic doesn't prune)
//    if (m_options->comparePruning) {
      if (m_options->fglpMBEHeur) {
        vector<double> labels(out.size(),ELEM_ONE);
        getLabelAll(var,assignment,node,labels);
        for (unsigned int i=0; i<fglpOut.size(); ++i) {
            double fglpH = fglpOut[i] OP_TIMES labels[i];
            double mbeH = mbeOut[i] OP_TIMES labels[i];
            bool fglpPruned = calculatePruning(var,node,fglpH);
            bool mbePruned = calculatePruning(var,node,mbeH);
            if (fglpPruned && !mbePruned) {
                timesFGLPPruned[m_pseudotree->getNode(var)->getDepth()]++;
            }
            else if (!fglpPruned && mbePruned) {
                timesMBEPruned[m_pseudotree->getNode(var)->getDepth()]++;
            }
            else if (fglpPruned && mbePruned) {
                timesBothPruned[m_pseudotree->getNode(var)->getDepth()]++;
            }
        }
        for (unsigned int i=0; i<out.size(); ++i) {
            if (fglpOut[i] < mbeOut[i]) {
                //            cout << "(var=" << var << ", val=" << i << ")" << endl;
                //            cout << "Found " << fglpOut[i] << " < " << mbeOut[i] << endl;
                out[i] = fglpOut[i];
                timesFGLPUsed[m_pseudotree->getNode(var)->getDepth()]++;
            }
            else {
                out[i] = mbeOut[i];
                timesMBEUsed[m_pseudotree->getNode(var)->getDepth()]++;
            }
        }
      }
      else {
          for (unsigned int i=0; i<out.size(); ++i) {
              out[i] = fglpOut[i];
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
    if (m_options->useShiftedLabels) {
        fglpHeur->getLabelAll(var, assignment, node, out);
    }
    else {
        vector<double> costTmp(m_problem->getDomainSize(var), ELEM_ONE);
        for (Function *f : m_pseudotree->getFunctions(var)) {
            f->getValues(assignment, var, costTmp);
            for (int i=0; i<m_problem->getDomainSize(var); ++i) {
                out[i] OP_TIMESEQ costTmp[i];
            }
        }
    }
}

bool FGLPMBEHybrid::calculatePruning(int var, SearchNode *node, double curPSTVal) {
    if (curPSTVal == ELEM_ZERO) return true;

    SearchNode *curOR = (node->getType() == NODE_OR) ? node : node->getParent();

    if (curPSTVal <= curOR->getValue()) return true;

    SearchNode *curAND = NULL;
    
    while (curOR->getParent()) {
        curAND = curOR->getParent();
        curPSTVal OP_TIMESEQ curAND->getLabel();
        curPSTVal OP_TIMESEQ curAND->getSubSolved();

        NodeP* children = curAND->getChildren();
        for(size_t i=0; i<curAND->getChildCountFull(); ++i) {
            if (!children[i] || children[i] == curOR) continue;
            else curPSTVal OP_TIMESEQ children[i]->getHeur();
        }
        curOR = curAND->getParent();

        if (curPSTVal <= curOR->getValue()) {
            return true;
        }
    }
    return false;
}
