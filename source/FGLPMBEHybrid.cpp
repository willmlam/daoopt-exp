#include "FGLPMBEHybrid.h"
using namespace std;

FGLPMBEHybrid::FGLPMBEHybrid(Problem *p, Pseudotree *pt, ProgramOptions *po)
    : Heuristic(p, pt, po), fglp(nullptr) {
  // do any original preprocessing first
  if (m_options != NULL && (m_options->mplp > 0 || m_options->mplps > 0)) {
    if (m_options->usePriority)
      fglp = new PriorityFGLP(m_problem, po->useNullaryShift);
    else
      fglp = new FGLP(m_problem, po->useNullaryShift);
    fglp->Run(m_options->mplp < 0 ? 5 : m_options->mplp, m_options->mplps,
              m_options->mplpt);
    fglp->set_owns_factors(false);
    m_problem->replaceFunctions(fglp->factors());
    m_pseudotree->addFunctionInfo(m_problem->getFunctions());
    m_options->mplp = 0;
    m_options->mplps = 0;

    cout << "Problem size (MB) after FGLP preprocessing: "
         << (m_problem->getSize() * sizeof(double) / (1024 * 1024.0)) << endl;
  }

  // Construct the fglp heuristic here on the smaller functions
  fglpHeur = new FGLPHeuristic(p, pt, po);

  if (m_options->fglpMBEHeur) {
    if (m_options != NULL && (m_options->jglp > 0 || m_options->jglps > 0)) {
      mex::mbe _jglp(copyFactors());
      mex::VarOrder ord(m_pseudotree->getElimOrder().begin(),
                        --m_pseudotree->getElimOrder().end());
      // cout << m_pseudotree->getElimOrder().size() << endl;
      _jglp.setOrder(ord);

      // cout << m_problem->getN() << endl;
      mex::VarOrder parents(m_problem->getN() -
                            1);  // copy pseudotree information
      for (int i = 0; i < m_problem->getN() - 1; ++i) {
        int par = m_pseudotree->getNode(i)->getParent()->getVar();
        parents[i] = (par == m_pseudotree->getRoot()->getVar()) ? -1 : par;
      }
      _jglp.setPseudotree(parents);

      limitJGLPIBound(m_options->memlimit);

      _jglp.setIBound(m_options->jglpi);

      _jglp.setProperties("DoMatch=1,DoFill=0,DoJG=1,DoMplp=0");
      cout << "Running JGLP(" << m_options->jglpi << ")" << endl;

      _jglp.init();

      int iter;
      if (m_options->jglp > 0)
        iter = m_options->jglp;
      else
        iter = 100;
      _jglp.tighten(iter, m_options->jglps);
      rewriteFactors(_jglp.factors());

      // JGLP code duplicates scopes, recollapse functions
      m_problem->collapseFunctions();
      m_pseudotree->addFunctionInfo(m_problem->getFunctions());

      // Make sure preprocessing doesn't run again.
      m_options->jglp = 0;
      m_options->jglps = 0;

      cout << "Problem size (MB) after JGLP preprocessing: "
           << (m_problem->getSize() * sizeof(double) / (1024 * 1024.0)) << endl;
    }
    mbeHeur = new MiniBucketElim(p, pt, po, po->ibound);
  } else {
    mbeHeur = nullptr;
  }

  timesFGLPUsed.resize(p->getN(), 0);
  timesMBEUsed.resize(p->getN(), 0);
  timesFGLPPruned.resize(p->getN(), 0);
  timesMBEPruned.resize(p->getN(), 0);
  timesBothPruned.resize(p->getN(), 0);
}

size_t FGLPMBEHybrid::build(const std::vector<val_t> *assignment,
                            bool computeTables) {
  fglpHeur->build(assignment, computeTables);
  if (m_options->usePriority) {
    PriorityFGLP *pfglp = dynamic_cast<PriorityFGLP *>(fglpHeur->getRootFGLP());
    pfglp->set_var_priority(dynamic_cast<PriorityFGLP *>(fglp)->var_priority());
  }
  if (m_options->fglpMBEHeur)
    return mbeHeur->build(assignment, computeTables);
  else
    return 0;
}

double FGLPMBEHybrid::getHeur(int var, const vector<val_t> &assignment,
                              SearchNode *node) {
  return ELEM_ZERO;
}

void FGLPMBEHybrid::getHeurAll(int var, const vector<val_t> &assignment,
                               SearchNode *node, vector<double> &out) {
  vector<double> fglpOut(out.size(), ELEM_ONE);
  vector<double> mbeOut(out.size(), ELEM_ONE);

  // Do not use FGLP at all if MBE is accurate
  if (!isAccurate()) {
    // If cost reversal is on, we may not need to adjust? (TODO)
    if (m_options->useShiftedLabels) {
      fglpHeur->getHeurAll(var, assignment, node, fglpOut);
    } else {
      fglpHeur->getHeurAllAdjusted(var, assignment, node, fglpOut);
    }
  }

  if (m_options->fglpMBEHeur)
    mbeHeur->getHeurAll(var, assignment, node, mbeOut);

  /*
     cout << "FGLP:"<< fglpOut << endl;
     cout << "MBE :"<< mbeOut << endl;
     */

  // Count the number of possible prunings for each heuristics wrt to the
  // values (only if the other heuristic doesn't prune)
  //    if (m_options->comparePruning) {
  if (m_options->fglpMBEHeur) {
    vector<double> labels(out.size(), ELEM_ONE);
    getLabelAll(var, assignment, node, labels);
    for (unsigned int i = 0; i < fglpOut.size(); ++i) {
      double fglpH = fglpOut[i] OP_TIMES labels[i];
      double mbeH = mbeOut[i] OP_TIMES labels[i];
      bool fglpPruned = calculatePruning(var, node, fglpH);
      bool mbePruned = calculatePruning(var, node, mbeH);
      if (fglpPruned && !mbePruned) {
        timesFGLPPruned[m_pseudotree->getNode(var)->getDepth()]++;
      } else if (!fglpPruned && mbePruned) {
        timesMBEPruned[m_pseudotree->getNode(var)->getDepth()]++;
      } else if (fglpPruned && mbePruned) {
        timesBothPruned[m_pseudotree->getNode(var)->getDepth()]++;
      }
    }
    for (unsigned int i = 0; i < out.size(); ++i) {
      if (fglpOut[i] < mbeOut[i]) {
        //            cout << "(var=" << var << ", val=" << i << ")" << endl;
        //            cout << "Found " << fglpOut[i] << " < " << mbeOut[i] <<
        // endl;
        out[i] = fglpOut[i];
        timesFGLPUsed[m_pseudotree->getNode(var)->getDepth()]++;
      } else {
        out[i] = mbeOut[i];
        timesMBEUsed[m_pseudotree->getNode(var)->getDepth()]++;
      }
    }
  } else {
    for (unsigned int i = 0; i < out.size(); ++i) {
      out[i] = fglpOut[i];
    }
  }
}

double FGLPMBEHybrid::getLabel(int var, const vector<val_t> &assignment,
                               SearchNode *node) {
  double d = ELEM_ONE;
  for (Function *f : m_pseudotree->getFunctions(var)) {
    d OP_TIMESEQ f->getValue(assignment);
  }
  return d;
}

void FGLPMBEHybrid::getLabelAll(int var, const vector<val_t> &assignment,
                                SearchNode *node, vector<double> &out) {
  if (m_options->useShiftedLabels) {
    fglpHeur->getLabelAll(var, assignment, node, out);
  } else {
    vector<double> costTmp(m_problem->getDomainSize(var), ELEM_ONE);
    for (Function *f : m_pseudotree->getFunctions(var)) {
      f->getValues(assignment, var, costTmp);
      for (int i = 0; i < m_problem->getDomainSize(var); ++i) {
        out[i] OP_TIMESEQ costTmp[i];
      }
    }
  }
}

bool FGLPMBEHybrid::calculatePruning(int var, SearchNode *node,
                                     double curPSTVal) {
  if (curPSTVal == ELEM_ZERO) return true;

  SearchNode *curOR = (node->getType() == NODE_OR) ? node : node->getParent();

  if (curPSTVal <= curOR->getValue()) return true;

  SearchNode *curAND = NULL;

  while (curOR->getParent()) {
    curAND = curOR->getParent();
    curPSTVal OP_TIMESEQ curAND->getLabel();
    curPSTVal OP_TIMESEQ curAND->getSubSolved();

    NodeP *children = curAND->getChildren();
    for (size_t i = 0; i < curAND->getChildCountFull(); ++i) {
      if (!children[i] || children[i] == curOR)
        continue;
      else
        curPSTVal OP_TIMESEQ children[i]->getHeur();
    }
    curOR = curAND->getParent();

    if (curPSTVal <= curOR->getValue()) {
      return true;
    }
  }
  return false;
}

// Copy DaoOpt Function class into mex::Factor class structures
mex::vector<mex::Factor> FGLPMBEHybrid::copyFactors(void) {
  mex::vector<mex::Factor> fs(m_problem->getC());
  for (int i = 0; i < m_problem->getC(); ++i)
    fs[i] = m_problem->getFunctions()[i]->asFactor().exp();
  return fs;
}

// Mini-bucket may have re-parameterized the original functions; if so, replace
// them
void FGLPMBEHybrid::rewriteFactors(const vector<mex::Factor> &factors) {
  //  vector<Function*> newFunctions(factors.size()); // to hold replacement,
  // reparameterized functions
  vector<Function *> newFunctions;
  for (size_t f = 0; f < factors.size();
       ++f) {  // allocate memory, copy variables into std::set
    double *tablePtr = new double[factors[f].nrStates()];
    std::set<int> scope;
    for (mex::VarSet::const_iterator v = factors[f].vars().begin();
         v != factors[f].vars().end(); ++v)
      scope.insert(v->label());
    if (scope.size() > 0 && factors[f].nrStates() == 1) continue;
    newFunctions.push_back(new FunctionBayes(f, m_problem, scope, tablePtr,
                                             factors[f].nrStates()));
    newFunctions[f]
        ->fromFactor(log(factors[f]));  // write in log factor functions
  }
  double *table1 = new double[1];
  table1[0] = m_problem->globalConstInfo();
  Function *constFun = new FunctionBayes(newFunctions.back()->getId() + 1,
                                         m_problem, set<int>(), table1, 1);
  newFunctions.push_back(constFun);

  m_problem->replaceFunctions(
      newFunctions);  // replace them in the problem definition
}

size_t FGLPMBEHybrid::computeMBEMemory(int ibound) {
  const vector<int> &elimOrder = m_pseudotree->getElimOrder();
  size_t memSize = 0;

  vector<vector<Function *> > augmented(m_problem->getN(),
                                        vector<Function *>());

  for (vector<int>::const_iterator itV = elimOrder.begin();
       itV != elimOrder.end(); ++itV) {
    vector<Function *> funs;
    const vector<Function *> &fnlist = m_pseudotree->getFunctions(*itV);
    funs.insert(funs.end(), fnlist.begin(), fnlist.end());
    funs.insert(funs.end(), augmented[*itV].begin(), augmented[*itV].end());

    if (*itV == elimOrder.back()) {
      continue;
    }

    sort(funs.begin(), funs.end(), scopeIsLarger);

    vector<MiniBucket> minibuckets;
    for (vector<Function *>::iterator itF = funs.begin(); itF != funs.end();
         ++itF) {
      bool placed = false;
      for (vector<MiniBucket>::iterator itB = minibuckets.begin();
           !placed && itB != minibuckets.end(); ++itB) {
        if (itB->allowsFunction(*itF)) {
          itB->addFunction(*itF);
          placed = true;
        }
      }
      if (!placed) {
        MiniBucket mb(*itV, ibound, m_problem);
        mb.addFunction(*itF);
        minibuckets.push_back(mb);
      }
    }

    for (vector<MiniBucket>::iterator itB = minibuckets.begin();
         itB != minibuckets.end(); ++itB) {
      Function *newf = itB->eliminate(false);
      const set<int> &newscope = newf->getScopeSet();
      memSize += newf->getTableSize();

      PseudotreeNode *n = m_pseudotree->getNode(*itV)->getParent();
      while (newscope.find(n->getVar()) == newscope.end() &&
             n != m_pseudotree->getRoot()) {
        n = n->getParent();
      }
      augmented[n->getVar()].push_back(newf);
    }
  }

  for (vector<vector<Function *> >::iterator itA = augmented.begin();
       itA != augmented.end(); ++itA)
    for (vector<Function *>::iterator itB = itA->begin(); itB != itA->end();
         ++itB)
      delete *itB;
  augmented.clear();

  return memSize;
}

size_t FGLPMBEHybrid::limitJGLPIBound(size_t memlimit) {
  memlimit *= 1024 * 1024 / sizeof(double);

  double divideFactor = m_options->jglpi == -1 ? 2 : 1;

  if (m_options->jglpi == -1) m_options->jglpi = m_options->ibound;
  int ibound = m_options->jglpi;

  size_t mem = computeMBEMemory(ibound);
  cout << " i=" << ibound << " -> "
       << ((mem / (1024 * 1024.0)) * sizeof(double)) << " MBytes" << endl;

  while (mem > memlimit && ibound > 1) {
    ibound--;
    mem = computeMBEMemory(ibound);
    cout << " i=" << ibound << " -> "
         << ((mem / (1024 * 1024.0)) * sizeof(double)) << " MBytes" << endl;
  }

  m_options->jglpi = ibound / divideFactor;

  return mem;
}
