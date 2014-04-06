#include "FGLPHeuristic.h"
#include "Graph.h"
using namespace std;

FGLPHeuristic::FGLPHeuristic(Problem *p, Pseudotree *pt, ProgramOptions *po) 
: Heuristic(p,pt,po), rootFGLP(nullptr), totalIterationsRun(0), totalInitiated(0) {
    // Precompute lists of variables for each subproblem
    m_ordering.resize(p->getN());
    for (int i = 0; i < p->getN(); ++i) {
        findDfsOrder(m_ordering[i], i);
    }

    Graph g(p->getN());
    for (Function *f : p->getFunctions()) {
        g.addClique(f->getScopeVec());
    }

    m_updateOrdering.resize(p->getN());
    for (int v : m_ordering.back()) {

        // keeps track of which node has been visited
        SETCLASS<int> nodes(++m_ordering[v].begin(), m_ordering[v].end());
#if defined HASH_GOOGLE_SPARSE || defined HASH_GOOGLE_DENSE
        nodes.set_deleted_key(UNKNOWN);
#endif
#ifdef HASH_GOOGLE_DENSE
        nodes.set_empty_key(UNKNOWN-1);
#endif

        queue<int> bfs;
        if (v == pt->getRoot()->getVar()) {
            assert(m_ordering.back().size() > 1);
            bfs.push(m_ordering.back()[1]);
        }
        else {
            bfs.push(v);
        }
        nodes.erase(bfs.front());
        while (!bfs.empty()) {
            int currentVar = bfs.front(); 
            m_updateOrdering[v].push_back(currentVar);
            bfs.pop();
            const set<int> &nb = g.getNeighbors(currentVar);
            for (int vn : nb) {
                if (nodes.find(vn) != nodes.end()) {
                    bfs.push(vn);
                    nodes.erase(vn);
                }
            }
            if (v == pt->getRoot()->getVar() && !nodes.empty()) {
                for (int k : m_ordering.back()) {
                    if (nodes.find(k) != nodes.end()) {
                        bfs.push(k);
                        nodes.erase(k);
                        break;
                    }
                }
            }
        }
        g.removeNode(v);
    }
    computeSubproblemFunIds();

}

size_t FGLPHeuristic::build(const std::vector<val_t> *assignment, bool computeTables) {
    rootFGLP = new FGLP(m_problem->getNOrg(), m_problem->getDomains(), 
            m_problem->getFunctions(), m_updateOrdering.back());
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

//    cout << "var: " << var << endl;
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
    /*
    vector<Function*> funs;
    for (Function *f : parentFGLP->getFactors()) {
        if (m_subproblemFunIds[var].find(f->getId()) != m_subproblemFunIds[var].end() ||
                f->getArity() == 0) {
            funs.push_back(f);
        }
        else {
            cout << "not included" << endl;
            cout << *f;
        }
    }
    */


    FGLP *varFGLP = new FGLP(m_problem->getN(),
            m_problem->getDomains(),
//       funs,
            parentFGLP->getFactors(),
            (m_options->useFglpBfs ? m_updateOrdering[var] : m_ordering[var]),
            m_tempAssn);

//    varFGLP->setVerbose(true);
    varFGLP->run(m_options->ndfglp, m_options->ndfglps, m_options->ndfglpt);

    totalIterationsRun += varFGLP->getRunIters();
    if (varFGLP->getRunIters() > 0) {
        totalInitiated++;
    }

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
    double edgeCostDiff;

    vector<double> costTmp(m_problem->getDomainSize(var), ELEM_ONE);
    m_tempLabelsFGLP.clear();
    m_tempLabelsFGLP.resize(m_problem->getDomainSize(var), ELEM_ONE);
    info->getFGLPStore()->getLabelAll(var,m_tempLabelsFGLP);

    m_tempLabels.clear();
    m_tempLabels.resize(m_problem->getDomainSize(var), ELEM_ONE);

    for (Function *f : m_pseudotree->getFunctions(var)) {
        f->getValues(assignment, var, costTmp);
        for (int i=0; i<m_problem->getDomainSize(var); ++i) {
            m_tempLabels[i] OP_TIMESEQ costTmp[i];
        }
    }

    // Compute the difference in the costs going to the current node
//    cout << "Difference in costs to node: " << endl;
//    cout << info->getFGLPStore()->getConstant() << " , " << info->getOrigCostToNode() << endl;
    double restCostDiff = info->getFGLPStore()->getConstant() OP_DIVIDE
            info->getOrigCostToNode();

    // Compute the difference in costs in selecting the next value
    for (size_t i=0; i<out.size(); ++i) {
//        cout << "Difference in assignment (" << i << "): " << endl;
//        cout << m_tempLabelsFGLP[i] << " , " << m_tempLabels[i] << endl;
        edgeCostDiff = m_tempLabelsFGLP[i] OP_DIVIDE m_tempLabels[i];
        adjustment = restCostDiff OP_TIMES edgeCostDiff;
        if (adjustment != 0 && !std::isnan(adjustment)) {
            out[i] OP_TIMESEQ adjustment;
//            cout << "Adjusted: " << adjustment << endl;
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

void FGLPHeuristic::findBfsOrder(vector<int> &order, int var) const {
    order.clear();
    queue<PseudotreeNode*> bfs;
    if (m_pseudotree->getRoot()->getVar() != var) {
        order.push_back(m_pseudotree->getRoot()->getVar());
    }
    bfs.push(m_pseudotree->getNode(var));
    PseudotreeNode* n = NULL;
    while (!bfs.empty()) {
        n = bfs.front();
        bfs.pop();
        order.push_back(n->getVar());
        for (vector<PseudotreeNode*>::const_iterator it=n->getChildren().begin();
                it!=n->getChildren().end(); ++it) {
            bfs.push(*it);
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
        for (Function *f : m_pseudotree->getNode(*itV)->getParent()->getFunctions()) {
            m_subproblemFunIds[*itV].insert(f->getId());
        }
        for (int id : m_subproblemFunIds[*itV]) {
            m_subproblemFunIds[m_pseudotree->getNode(*itV)->
                getParent()->getVar()].insert(id);
        }
    }
}
