#include "FGLPHeuristic.h"
using namespace std;

FGLPHeuristic::FGLPHeuristic(Problem *p, Pseudotree *pt, ProgramOptions *po) 
: Heuristic(p,pt,po), rootFGLP(nullptr) {
    // Precompute lists of variables for each subproblem
    m_ordering.resize(p->getN());
    for (int i = 0; i < p->getN(); ++i) {
        findDfsOrder(m_ordering[i], i);
    }
    computeSubproblemFunIds();

}

size_t FGLPHeuristic::build(const std::vector<val_t> *assignment, bool computeTables) {
    rootFGLP = new FGLP(m_problem->getNOrg(), m_problem->getDomains(), 
            m_problem->getFunctions(), m_ordering.back());
//    rootFGLP->setVerbose(true);
    rootFGLP->run(m_options->mplp, m_options->mplps);
    m_globalUB = rootFGLP->getUB();
    return 0;
}

double FGLPHeuristic::getHeur(int var, const vector<val_t> &assignment, SearchNode *node) {
    // TO DO
    return ELEM_ZERO;
}

void FGLPHeuristic::getHeurAll(int var, const vector<val_t> &assignment, SearchNode *node,
        vector<double> &out) {

    FGLP *parentFGLP;
    double parentCost;
//    double parentCostShifted;
    // Is the root node?
    if (!node->getParent()) {
        parentFGLP = rootFGLP;
        parentCost = ELEM_ONE;
//        parentCostShifted = ELEM_ONE;
    }
    else {
        SearchNode *parentOR = node->getParent()->getParent();
        FGLPNodeInfo *parentInfo = 
            static_cast<FGLPNodeInfo*>(parentOR->getExtraNodeInfo().get());

        /*
        cout << "parent var/val: " 
            << node->getParent()->getVar() << "/" << int(node->getParent()->getVal())
            << endl;
        cout << "current var: " << node->getVar() << endl;
        */

        parentFGLP = parentInfo->getFGLPStore().get();
        parentCost = parentInfo->getOrigCostToNode();
        parentCost OP_TIMESEQ node->getParent()->getLabel();
//        parentCostShifted = parentFGLP->getConstant();


        // verify the computed cost
        SearchNode *cur = node->getParent();
        double vCost = ELEM_ONE;
        while (cur) {
            if (cur->getType() == NODE_AND) {
                vCost OP_TIMESEQ cur->getLabel();
            }
            cur = cur->getParent();
        }
//        cout << "Computed original cost vs cost from space: ";
        if (fabs(parentCost-vCost) >= 1e-8) {
            cout << parentCost << " != " << vCost << " --WARNING!" << endl;
        }
    }

    m_tempAssn.clear();
    const vector<int> &relVars = m_pseudotree->getNode(var)->getFullContextVec();
    for (int v : relVars) {
        m_tempAssn[v] = assignment[v];
    }

    FGLPNodeInfo *info = new FGLPNodeInfo();
    node->setExtraNodeInfo(info);
    vector<Function*> funs;
    for (Function *f : parentFGLP->getFactors()) {
        if (m_subproblemFunIds[var].find(f->getId()) != m_subproblemFunIds[var].end() ||
                f->getArity() == 0) {
            funs.push_back(f);
        }
    }

    FGLP *varFGLP = new FGLP(m_problem->getN(),
            m_problem->getDomains(),
            funs,
            m_ordering[var],
            m_tempAssn);

    varFGLP->run(m_options->ndfglp, m_options->ndfglps);

//    cout << "FGLP size (MB): " << (varFGLP->getSize()*sizeof(double)) / (1024*1024.0)  << endl;
    
    varFGLP->getVarUB(var, out);
    info->setFGLPStore(varFGLP);
    info->setOrigCostToNode(parentCost);

}

void FGLPHeuristic::getHeurAllAdjusted(int var, const vector<val_t> &assignment, SearchNode *node, vector<double> &out) {
    getHeurAll(var, assignment, node, out);
    // Calculate the difference between the path costs and adjust if needed
    FGLPNodeInfo* info = static_cast<FGLPNodeInfo*>(node->getExtraNodeInfo().get());
    double adjustment;
    double originalCost;

    vector<double> costTmp(m_problem->getDomainSize(var), ELEM_ONE);
    m_tempLabels.clear();
    m_tempLabels.resize(m_problem->getDomainSize(var), ELEM_ONE);
    for (Function *f : m_pseudotree->getFunctions(var)) {
        f->getValues(assignment, var, costTmp);
        for (int i=0; i<m_problem->getDomainSize(var); ++i) {
            m_tempLabels[i] OP_TIMESEQ costTmp[i];
        }
    }
    // Need to get original labels into m_tempLabels
    for (size_t i=0; i<out.size(); ++i) {
        originalCost = info->getOrigCostToNode() OP_TIMES m_tempLabels[i];
        adjustment = info->getFGLPStore()->getConstant() OP_DIVIDE originalCost;
        if (adjustment != 0 && !std::isnan(adjustment)) {
            out[i] OP_TIMESEQ adjustment;
        }
    }
}

double FGLPHeuristic::getLabel(int var, const vector<val_t> &assignment, SearchNode *node) {
    assert(assignment[var] != NONE);
    const vector<Function*> &functions = static_cast<FGLPNodeInfo*>(node->getExtraNodeInfo().get())->getFGLPStore()->getFactors();
    double labelValue = ELEM_ONE;
    for (auto f : functions)
        if (f->getArity() == 1 && f->getScopeVec()[0] == var)
            labelValue OP_TIMESEQ f->getTable()[assignment[var]];
    return labelValue;
}

void FGLPHeuristic::getLabelAll(int var, const vector<val_t> &assignment, SearchNode *node, vector<double> &out) {
    const vector<Function*> &functions = static_cast<FGLPNodeInfo*>(node->getExtraNodeInfo().get())->getFGLPStore()->getFactors();
    double labelValue;
    for (int i=0; i<m_problem->getDomainSize(var); ++i) {
        labelValue = ELEM_ONE;
        for (auto f : functions)
            if (f->getArity() == 1 && f->getScopeVec()[0] == var)
                labelValue OP_TIMESEQ f->getTable()[i];
        out[i] = labelValue;
    }
}

void FGLPHeuristic::findDfsOrder(vector<int> &order, int var) const {
    order.clear();
    stack<PseudotreeNode*> dfs;
    if (m_pseudotree->getRoot()->getVar() != var) {
        order.push_back(m_pseudotree->getRoot()->getVar());
    }
    dfs.push(m_pseudotree->getNode(var));
    PseudotreeNode* n = NULL;
    while (!dfs.empty()) {
        n = dfs.top();
        dfs.pop();
        order.push_back(n->getVar());
        for (vector<PseudotreeNode*>::const_iterator it=n->getChildren().begin();
                it!=n->getChildren().end(); ++it) {
            dfs.push(*it);
        }
    }
}

void FGLPHeuristic::computeSubproblemFunIds() {
    m_subproblemFunIds.clear();
    m_subproblemFunIds.resize(m_problem->getN(), set<int>());

    // For each node, add all of its current function ids to its parent
    for (auto itV = m_ordering.back().rbegin(); itV != m_ordering.back().rend(); ++itV) {
        if (*itV == m_pseudotree->getRoot()->getVar()) continue;
        for (Function *f : m_pseudotree->getNode(*itV)->getFunctions()) {
            m_subproblemFunIds[*itV].insert(f->getId());
        }
        for (int id : m_subproblemFunIds[*itV]) {
            m_subproblemFunIds[m_pseudotree->getNode(*itV)->
                getParent()->getVar()].insert(id);
        }
    }
}
