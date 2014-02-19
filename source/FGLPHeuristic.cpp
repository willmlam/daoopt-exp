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
            static_cast<FGLPNodeInfo*>(parentOR->getExtraNodeInfo());

        /*
        cout << "parent var/val: " 
            << node->getParent()->getVar() << "/" << int(node->getParent()->getVal())
            << endl;
        cout << "current var: " << node->getVar() << endl;
        */

        parentFGLP = parentInfo->getFGLPStore()[node->getParent()->getVal()];
        parentCost = parentInfo->getOrigCostToNode()[node->getParent()->getVal()];
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

    tempAssn.clear();
    const vector<int> &relVars = m_pseudotree->getNode(var)->getFullContextVec();
    for (int v : relVars) {
        tempAssn[v] = assignment[v];
    }

    vector<double> costTmp(m_problem->getDomainSize(var), ELEM_ONE);
    vector<double> costs(m_problem->getDomainSize(var), ELEM_ONE);
    for (Function *f : m_pseudotree->getFunctions(var)) {
        f->getValues(assignment, var, costTmp);
        for (int i=0; i<m_problem->getDomainSize(var); ++i) {
            costs[i] OP_TIMESEQ costTmp[i];
        }
    }


    FGLPNodeInfo *info = new FGLPNodeInfo();
    node->setExtraNodeInfo(info);
    // For each value, get factors of parentFGLP and run FGLP on them with conditioning
    for (int i=0; i<m_problem->getDomainSize(var); ++i) {
//        cout << "val: " << i << endl;
        tempAssn[var] = i;
        /*
        for (Function *f : parentFGLP->getFactors()) {
            cout << *f << endl;
        }
        cout << endl;
        */
        // work off shifted version or original problem?
//        const vector<Function*> &funsOrig = m_problem->getFunctions();
        vector<Function*> funs;
        for (Function *f : parentFGLP->getFactors()) {
            if (m_subproblemFunIds[var].find(f->getId()) != m_subproblemFunIds[var].end() ||
                    f->getArity() == 0) {
                funs.push_back(f);
            }
        }

        FGLP *valFGLP = new FGLP(m_problem->getN(),
                m_problem->getDomains(),
                funs,
//                parentFGLP->getFactors(),
                m_ordering[var],
                tempAssn);
//        valFGLP->setVerbose(true);
        valFGLP->run(m_options->ndfglp, m_options->ndfglps);

        // Add in original cost for traversing this edge and add to parent cost
        info->addToCosts(parentCost + costs[i]);
        /*
        cout << "Cost (original):" << parentCost + costs[i] << endl;
        cout << "Cost (shifted) :" << parentCostShifted + valFGLP->getLabel() << endl;
        */

        out[i] = valFGLP->getUBNonConstant();
        info->addToStore(valFGLP);

    }
}

void FGLPHeuristic::getHeurAllAdjusted(int var, const vector<val_t> &assignment, SearchNode *node, vector<double> &out) {
    getHeurAll(var, assignment, node, out);
    for (unsigned int i=0; i<out.size(); ++i) {
        FGLPNodeInfo* info = static_cast<FGLPNodeInfo*>(node->getExtraNodeInfo());
//        cout << info->getFGLPStore()[i]->getConstant() << endl;
//        cout << info->getOrigCostToNode()[i] << endl;
        double adjustment = info->getFGLPStore()[i]->getConstant() OP_DIVIDE info->getOrigCostToNode()[i];
        if (adjustment != 0 && !std::isnan(adjustment)) {
            /*
            cout << "(var=" << var << ", val=" << i << ")" << endl;
            cout << "Adjusted:  " << out[i] << " -> ";
            */
            out[i] OP_TIMESEQ adjustment;
//            cout << out[i] << endl;
        }
    }

}

double FGLPHeuristic::getLabel(int var, const vector<val_t> &assignment, SearchNode *node) {
    return static_cast<FGLPNodeInfo*>(node->getExtraNodeInfo())->
        getFGLPStore()[assignment[var]]->getLabel();
}

void FGLPHeuristic::getLabelAll(int var, const vector<val_t> &assignment, SearchNode *node, vector<double> &out) {
    for (int i=0; i<m_problem->getDomainSize(var); ++i) {
        out[i] = static_cast<FGLPNodeInfo*>(node->getExtraNodeInfo())->
            getFGLPStore()[i]->getLabel();
    }
    /*
    vector<double> costTmp(m_problem->getDomainSize(var), ELEM_ONE);
    for (Function *f : m_pseudotree->getFunctions(var)) {
        f->getValues(assignment, var, costTmp);
        for (int i=0; i<m_problem->getDomainSize(var); ++i) {
            out[i] OP_TIMESEQ costTmp[i];
        }
    }
    */
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
