#include "FGLPMBEHybrid.h"

namespace daoopt {

using namespace std;

FGLPMBEHybrid::FGLPMBEHybrid(Problem *p, Pseudotree *pt, ProgramOptions *po)
    : Heuristic(p, pt, po) {
  // do any original preprocessing first
  int root_iterations = m_options->mplp < 0 ? 5 : m_options->mplp;
  if (m_options != NULL && (m_options->mplp > 0 || m_options->mplps > 0)) {
    if (m_options->usePriority) {
      fglp.reset(new PriorityFGLP(m_problem, po->useNullaryShift));
      root_iterations *= m_problem->getN();
    } else {
      fglp.reset(new FGLP(m_problem, po->useNullaryShift));
    }
    fglp->Run(root_iterations, m_options->mplps, m_options->mplpt);
    fglp->set_owns_factors(false);
    m_problem->replaceFunctions(fglp->factors());
    m_pseudotree->resetFunctionInfo(m_problem->getFunctions());
    m_options->mplp = 0;
    m_options->mplps = 0;

    cout << "Problem size (MB) after FGLP preprocessing: "
         << (m_problem->getSize() * sizeof(double) / (1024 * 1024.0)) << endl;
  }

  // Construct the fglp heuristic here on the smaller functions
  fglpHeur.reset(new FGLPHeuristic(p, pt, po));

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
      m_pseudotree->resetFunctionInfo(m_problem->getFunctions());

      // Make sure preprocessing doesn't run again.
      m_options->jglp = 0;
      m_options->jglps = 0;

      cout << "Problem size (MB) after JGLP preprocessing: "
           << (m_problem->getSize() * sizeof(double) / (1024 * 1024.0)) << endl;
    }
    mbeHeur.reset(new MiniBucketElim(p, pt, po, po->ibound));
  } 

  for (int32 i = -1; i <= m_pseudotree->getHeight(); ++i) {
    timesFGLPUsed[i] = 0;
    timesMBEUsed[i] = 0;
    timesFGLPPruned[i] = 0;
    timesMBEPruned[i] = 0;
    timesBothPruned[i] = 0;
  }
}

size_t FGLPMBEHybrid::build(const std::vector<val_t> *assignment,
                            bool computeTables) {
  fglpHeur->build(assignment, computeTables);
  if (m_options->fglpMBEHeur && mbeHeur)
    return mbeHeur->build(assignment, computeTables);
  else
    return 0;
}

double FGLPMBEHybrid::getHeur(int var, vector<val_t> &assignment,
                              SearchNode *node) {
  return ELEM_ZERO;
}

double FGLPMBEHybrid::getHeurPerIndSubproblem(int var, std::vector<val_t> & assignment, SearchNode* node, double label, std::vector<double> & subprobH) {
  return ELEM_ZERO;
}

void FGLPMBEHybrid::getHeurAll(int var, vector<val_t> &assignment,
                               SearchNode *node, vector<double> &out) {
  vector<double> fglpOut(out.size(), ELEM_ONE);
  vector<double> mbeOut(out.size(), ELEM_ONE);

  // Do not use FGLP at all if MBE is accurate
  if (!isAccurate()) {
    fglpHeur->getHeurAll(var, assignment, node, fglpOut);
  }

  if (m_options->fglpMBEHeur && mbeHeur) {
    mbeHeur->getHeurAll(var, assignment, node, mbeOut);
  }

  /*
  cout << "FGLP:" << fglpOut << endl;
  cout << "MBE: " << mbeOut << endl;
  cout << endl;
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
//        cin.get();
      } else if (!fglpPruned && mbePruned) {
        timesMBEPruned[m_pseudotree->getNode(var)->getDepth()]++;
      } else if (fglpPruned && mbePruned) {
        timesBothPruned[m_pseudotree->getNode(var)->getDepth()]++;
      }
    }
    for (unsigned int i = 0; i < out.size(); ++i) {
      if (fglpOut[i] < mbeOut[i]) {
        out[i] = fglpOut[i];
        timesFGLPUsed[m_pseudotree->getNode(var)->getDepth()]++;
        /*
        cout << StrCat("var: ", var) << endl;
        cout << StrCat("val: ", i) << endl;
        cout << StrCat(fglpOut[i], "<", mbeOut[i]) << endl; 
        cout << assignment << endl;
        cin.get();
        */
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

}  // namespace daoopt
