#include "FGLPHeuristic.h"
#include "Graph.h"

namespace daoopt {

using namespace std;

FGLPHeuristic::FGLPHeuristic(Problem *p, Pseudotree *pt, ProgramOptions *po)
    : Heuristic(p, pt, po),
      rootFGLP(nullptr),
      totalIterationsRun(0),
      totalInitiated(0) {

  m_countVars.resize(p->getN(), 0);
  m_varsUpdated.resize(p->getN(), 0);
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
    nodes.set_empty_key(UNKNOWN - 1);
#endif

    queue<int> bfs;
    if (v == pt->getRoot()->getVar()) {
      assert(m_ordering.back().size() > 1);
      bfs.push(m_ordering.back()[1]);
    } else {
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
  computeSubproblemVars();
}

size_t FGLPHeuristic::build(const std::vector<val_t> *assignment,
                            bool computeTables) {
  int root_iterations = m_options->mplp;
  if (m_options->usePriority) {
    rootFGLP = new PriorityFGLP(m_problem, m_options->useNullaryShift);
    root_iterations *= m_problem->getN();
  } else {
    rootFGLP = new FGLP(m_problem, m_options->useNullaryShift);
  }
  //    rootFGLP->setVerbose(true);

  rootFGLP->Run(root_iterations, m_options->mplps, m_options->mplpt);
  // Need to turn on cost shift reversal if we are doing AO with caching
  if (!m_options->orSearch) {
    rootFGLP->set_use_cost_shift_reversal(true);
  }
  m_globalUB = rootFGLP->ub();
  return 0;
}

double FGLPHeuristic::getHeur(int var, vector<val_t> &assignment,
                              SearchNode *node) {
  // TO DO
  return ELEM_ZERO;
}

double FGLPHeuristic::getHeurPerIndSubproblem(int var,
                                              std::vector<val_t> &assignment,
                                              SearchNode *node, double label,
                                              std::vector<double> &subprobH) {
  // TO DO
  return ELEM_ZERO;
}

void FGLPHeuristic::getHeurAll(int var, vector<val_t> &assignment,
                               SearchNode *node, vector<double> &out) {

  //    cout << "var: " << var << endl;
  FGLP *parentFGLP;
  double parentCost;
  //    double parentCostShifted;
  // Is the root node?
  if (!node->getParent()) {
    parentFGLP = rootFGLP;
    parentCost = ELEM_ONE;
    //        parentCostShifted = ELEM_ONE;
  } else {
    SearchNode *parentOR = node->getParent()->getParent();
    FGLPNodeInfo *parentInfo =
        static_cast<FGLPNodeInfo *>(parentOR->getExtraNodeInfo().get());

    /*
    cout << "parent var/val: "
        << node->getParent()->getVar() << "/" <<
    int(node->getParent()->getVal())
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
    if (fabs(parentCost - vCost) >= 1e-8) {
      cout << parentCost << " != " << vCost << " --WARNING!" << endl;
    }
  }

  m_tempAssn.clear();
  //    const vector<int> &relVars =
  // m_pseudotree->getNode(var)->getFullContextVec();
  /*
  for (int v : relVars) {
      m_tempAssn[v] = assignment[v];
  }
  */
  int condition_var = node->getParent()
                          ? m_pseudotree->getNode(var)->getParent()->getVar()
                          : -1;
  m_tempAssn[condition_var] = assignment[condition_var];

  FGLPNodeInfo *info = new FGLPNodeInfo();
  node->setExtraNodeInfo(info);
  /*
  vector<Function*> funs;
  for (Function *f : parentFGLP->getFactors()) {
      if (m_subproblemFunIds[var].find(f->getId()) !=
  m_subproblemFunIds[var].end() ||
              f->getArity() == 0) {
          funs.push_back(f);
      }
      else {
          cout << "not included" << endl;
          cout << *f;
      }
  }
  */

  FGLP *varFGLP;
  if (m_options->usePriority) {
    varFGLP =
        new PriorityFGLP(dynamic_cast<PriorityFGLP *>(parentFGLP), m_tempAssn,
                         m_subproblemVars[var], condition_var);
  } else {
    varFGLP =
        new FGLP(parentFGLP, m_tempAssn, m_subproblemVars[var], condition_var);
  }

  //    varFGLP->set_verbose(true);
  varFGLP->Run(m_options->ndfglp, m_options->ndfglps, m_options->ndfglpt);

  totalIterationsRun += varFGLP->runiters();
  if (varFGLP->runiters() > 0) {
    totalInitiated++;
  }

  if (node->getParent()) {
    int parentVar = node->getParent()->getVar();
    m_countVars[parentVar]++;
    m_varsUpdated[parentVar] += varFGLP->vars_updated().size();
  }

  //    cout << "FGLP size (MB): " << (varFGLP->getSize()*sizeof(double)) /
  // (1024*1024.0)  << endl;

  varFGLP->GetVarUB(var, out);
  //    cout << "var UB: " << out << endl;
  info->setFGLPStore(varFGLP);
  info->setOrigCostToNode(parentCost);

  if (!m_options->useShiftedLabels) {
    AdjustHeurAll(var, assignment, node, out);
  }
}

void FGLPHeuristic::AdjustHeurAll(int var, const vector<val_t> &assignment,
                                  SearchNode *node, vector<double> &out) {
  // Calculate the difference between the path costs and adjust if needed
  FGLPNodeInfo *info =
      static_cast<FGLPNodeInfo *>(node->getExtraNodeInfo().get());
  double adjustment;
  //  double edgeCostDiff;

  vector<double> costTmp(m_problem->getDomainSize(var), ELEM_ONE);
  /*
  m_tempLabelsFGLP.clear();
  m_tempLabelsFGLP.resize(m_problem->getDomainSize(var), ELEM_ONE);
  info->getFGLPStore()->GetLabelAll(var,m_tempLabelsFGLP);
  */

  m_tempLabels.clear();
  m_tempLabels.resize(m_problem->getDomainSize(var), ELEM_ONE);

  for (Function *f : m_pseudotree->getFunctions(var)) {
    f->getValues(assignment, var, costTmp);
    for (int i = 0; i < m_problem->getDomainSize(var); ++i) {
      m_tempLabels[i] OP_TIMESEQ costTmp[i];
    }
  }

  // Compute the difference in the costs going to the current node
  //    cout << "Difference in costs to node: " << endl;
  //    cout << info->getFGLPStore()->global_constant() << " , " <<
  // info->getOrigCostToNode() << endl;
  //    double restCostDiff = info->getFGLPStore()->global_constant() OP_DIVIDE
  // info->getOrigCostToNode();

  // Compute the difference in costs in selecting the next value
  for (size_t i = 0; i < out.size(); ++i) {
    //        cout << "Difference in assignment (" << i << "): " << endl;
    //        cout << m_tempLabelsFGLP[i] << " , " << m_tempLabels[i] << endl;
    //    edgeCostDiff = m_tempLabelsFGLP[i] OP_DIVIDE m_tempLabels[i];
    // Combine the original cost up until before conditioning and the
    // original cost for conditioning.
    adjustment = info->getOrigCostToNode();
    adjustment OP_TIMESEQ m_tempLabels[i];
    if (adjustment != 0 && !std::isnan(adjustment)) {
      if (!std::isinf(out[i])) {
        out[i] OP_DIVIDEEQ adjustment;
      }
    }
  }
}

double FGLPHeuristic::getLabel(int var, const vector<val_t> &assignment,
                               SearchNode *node) {
  assert(assignment[var] != NONE);
  const vector<Function *> &functions =
      static_cast<FGLPNodeInfo *>(node->getExtraNodeInfo().get())
          ->getFGLPStore()
          ->factors();
  double labelValue = ELEM_ONE;
  for (auto f : functions)
    if (f->getArity() == 1 && f->getScopeVec()[0] == var)
      labelValue OP_TIMESEQ f->getTable()[assignment[var]];
  return labelValue;
}

void FGLPHeuristic::getLabelAll(int var, const vector<val_t> &assignment,
                                SearchNode *node, vector<double> &out) {
  if (m_options->useShiftedLabels) {
    const vector<Function *> &functions =
        static_cast<FGLPNodeInfo *>(node->getExtraNodeInfo().get())
            ->getFGLPStore()
            ->factors();
    double labelValue;
    for (int i = 0; i < m_problem->getDomainSize(var); ++i) {
      labelValue = ELEM_ONE;
      for (auto f : functions) {
        if (f->getArity() == 1 && f->getScopeVec()[0] == var) {
          labelValue OP_TIMESEQ f->getTable()[i];
        }
      }
      out[i] = labelValue;
    }
  } else {
    vector<double> cost_tmp(m_problem->getDomainSize(var), ELEM_ONE);
    for (Function *f : m_pseudotree->getFunctions(var)) {
      f->getValues(assignment, var, cost_tmp);
      for (int i = 0; i < m_problem->getDomainSize(var); ++i) {
        out[i] OP_TIMESEQ cost_tmp[i];
      }
    }
  }
}

void FGLPHeuristic::findDfsOrder(vector<int> &order, int var) const {
  order.clear();
  stack<PseudotreeNode *> dfs;
  if (m_pseudotree->getRoot()->getVar() != var) {
    order.push_back(m_pseudotree->getRoot()->getVar());
  }
  dfs.push(m_pseudotree->getNode(var));
  PseudotreeNode *n = NULL;
  while (!dfs.empty()) {
    n = dfs.top();
    dfs.pop();
    order.push_back(n->getVar());
    for (vector<PseudotreeNode *>::const_iterator it = n->getChildren().begin();
         it != n->getChildren().end(); ++it) {
      dfs.push(*it);
    }
  }
}

void FGLPHeuristic::findBfsOrder(vector<int> &order, int var) const {
  order.clear();
  queue<PseudotreeNode *> bfs;
  if (m_pseudotree->getRoot()->getVar() != var) {
    order.push_back(m_pseudotree->getRoot()->getVar());
  }
  bfs.push(m_pseudotree->getNode(var));
  PseudotreeNode *n = NULL;
  while (!bfs.empty()) {
    n = bfs.front();
    bfs.pop();
    order.push_back(n->getVar());
    for (vector<PseudotreeNode *>::const_iterator it = n->getChildren().begin();
         it != n->getChildren().end(); ++it) {
      bfs.push(*it);
    }
  }
}

void FGLPHeuristic::computeSubproblemFunIds() {
  m_subproblemFunIds.clear();
  m_subproblemFunIds.resize(m_problem->getN(), set<int>());

  // For each node, add all of its current function ids to its parent
  for (auto itV = m_ordering.back().rbegin(); itV != m_ordering.back().rend();
       ++itV) {
    if (*itV == m_pseudotree->getRoot()->getVar()) continue;
    for (Function *f : m_pseudotree->getNode(*itV)->getFunctions()) {
      m_subproblemFunIds[*itV].insert(f->getId());
    }
    for (int id : m_subproblemFunIds[*itV]) {
      m_subproblemFunIds[m_pseudotree->getNode(*itV)->getParent()->getVar()]
          .insert(id);
    }
  }
}

void FGLPHeuristic::computeSubproblemVars() {
  m_subproblemVars.clear();
  m_subproblemVars.resize(m_problem->getN(), set<int>());

  // For each node, add all of its current variables to its parent
  for (auto itV = m_ordering.back().rbegin(); itV != m_ordering.back().rend();
       ++itV) {
    if (*itV == m_pseudotree->getRoot()->getVar()) continue;
    m_subproblemVars[*itV].insert(*itV);
    for (int v : m_subproblemVars[*itV]) {
      m_subproblemVars[m_pseudotree->getNode(*itV)->getParent()->getVar()]
          .insert(v);
    }
  }
}

}  // namespace daoopt
