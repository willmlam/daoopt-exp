/*
 * Main.cpp
 *
 *  Copyright (C) 2008-2012 Lars Otten
 *  This file is part of DAOOPT.
 *
 *  DAOOPT is free software: you can redistribute it and/or modify it
 *  under the terms of the GNU General Public License as published by
 *  the Free Software Foundation, either version 3 of the License, or
 *  (at your option) any later version.
 *
 *  DAOOPT is distributed in the hope that it will be useful, but
 *  WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU General Public License for more details.
 *
 *  You should have received a copy of the GNU General Public License
 *  along with DAOOPT.  If not, see <http://www.gnu.org/licenses/>.
 *
 *  Created on: Oct 18, 2011
 *      Author: Lars Otten <lotten@ics.uci.edu>
 */

#include "Main.h"
#include "DaooptInterface.h"
#include "MiniBucketElim.h"
#include "MiniBucketElimLH.h"
#include "MiniBucketElimInc.h"

#include "UAI2012.h"

#include "ARP/ARPall.hxx"

#include <chrono>
using namespace std::chrono;

#define VERSIONINFO "1.1.2"

namespace daoopt {

string UAI2012::filename = "";
string out_bound_file = "";

high_resolution_clock::time_point _time_start, _time_pre;

bool Main::parseOptions(int argc, char** argv) {
  // Reprint command line
  for (int i = 0; i < argc; ++i) cout << argv[i] << ' ';
  cout << endl;
  // parse command line
  ProgramOptions* opt = parseCommandLine(argc, argv);
  if (!opt) {
    err_txt("Error parsing command line.");
    return false;
  }

  if (opt->seed == NONE) opt->seed = time(0);
  rand::seed(opt->seed);

  m_options.reset(opt);

  size_t idx = m_options->in_problemFile.find_last_of("/");
  UAI2012::filename = m_options->in_problemFile.substr(idx + 1) + ".MPE";
  out_bound_file = m_options->out_boundFile;

  return true;
}

bool Main::setOptions(const ProgramOptions& options) {
  ProgramOptions* opt = new ProgramOptions(options);
  if (opt->seed == NONE) {
    opt->seed = time(0);
  }
  rand::seed(opt->seed);

  m_options.reset(opt);
  out_bound_file = m_options->out_boundFile;

  return true;
}

bool Main::setSLSOptions(int slsIter, int slsTimePerIter) {
  m_options.get()->slsIter = slsIter;
  m_options.get()->slsTime = slsTimePerIter;
  return true;
}

bool Main::setFGLPOptions(int mplp, int mplps) {
  m_options.get()->mplp = mplp;
  m_options.get()->mplps = mplps;
  return true;
}

bool Main::setJGLPOptions(int jglp, int jglps) {
  m_options.get()->jglp = jglp;
  m_options.get()->jglps = jglps;
  return true;
}

bool Main::setIboundOptions(int ibound) {
  m_options.get()->ibound = ibound;
  return true;
}

bool Main::loadProblem() {
  m_problem.reset(new Problem);

  // load problem file
  assert(m_options->in_problemFile != "" || m_options->problemSpec);
  string evid_string;
  if (!m_options->evidSpec) {
    evid_string = "0\n";
    m_options->evidSpec = &evid_string[0];
    m_options->evidSpec_len = evid_string.size();
  }
  string problem_string;
  if (!m_options->problemSpec) {
    problem_string = getFileContents(m_options->in_problemFile.c_str());
    m_options->problemSpec = &problem_string[0];
    m_options->problemSpec_len = problem_string.size();
  }
  //  if (!m_problem->parseUAI(m_options->in_problemFile,
  // m_options->in_evidenceFile, m_options->collapse))
  if (!m_problem->parseUAI16(m_options->problemSpec, m_options->problemSpec_len,
                           m_options->evidSpec, m_options->evidSpec_len,
                           m_options->collapse))
    return false;

  if (m_options->perturb > 0) {
    m_problem->perturbDeterminism(m_options->perturb);
  }
  cout << "Created problem with " << m_problem->getN() << " variables and "
       << m_problem->getC() << " functions." << endl;

  // Remove evidence variables
  m_problem->removeEvidence();
  cout << "Removed evidence, now " << m_problem->getN() << " variables and "
       << m_problem->getC() << " functions." << endl;

#if defined PARALLEL_DYNAMIC || defined PARALLEL_STATIC
  if (!m_options->par_solveLocal) {
    // Re-output problem file for parallel processing
    m_options->out_reducedFile = string("temp_prob.") + m_options->problemName +
                                 string(".") + m_options->runTag +
                                 string(".gz");
    m_problem->writeUAI(m_options->out_reducedFile);
    cout << "Saved problem to file " << m_options->out_reducedFile << endl;
  }
#else
  // Output reduced network?
  if (!m_options->out_reducedFile.empty()) {
    cout << "Writing reduced network to file " << m_options->out_reducedFile
         << endl;
    m_problem->writeUAI(m_options->out_reducedFile);
  }
#endif

  // Some statistics
  cout << "Global constant:\t" << SCALE_LOG(m_problem->globalConstInfo())
       << endl;
  cout << "Max. domain size:\t" << (int)m_problem->getK() << endl;
  cout << "Max. function arity:\t" << m_problem->getR() << endl;

  // Compute average tightness ratio
  double sum_tightness = 0.0;
  for (const Function* f : m_problem->getFunctions()) {
    sum_tightness += f->getTightness() / f->getTableSize();
  }
  cout << "Average function tightness:\t" << sum_tightness / m_problem->getC()
       << endl;

  return true;
}

bool Main::findOrLoadOrdering() {
  // Create primal graph of *reduced* problem
  Graph g(m_problem->getN());
  const vector<Function*>& fns = m_problem->getFunctions();
  for (vector<Function*>::const_iterator it = fns.begin(); it != fns.end();
       ++it) {
    g.addClique((*it)->getScopeVec());
  }
  cout << "Graph with " << g.getStatNodes() << " nodes and " << g.getStatEdges()
       << " edges created." << endl;

#ifdef PARALLEL_STATIC
  // for static parallelization post-processing mode, look for
  // ordering from preprocessing step
  if (m_options->par_postOnly) {
    m_options->in_orderingFile = string("temp_elim.") + m_options->problemName +
                                 string(".") + m_options->runTag +
                                 string(".gz");
    m_options->order_iterations = 0;
  }
#endif

  // Find variable ordering
  vector<int> elim;
  int w = numeric_limits<int>::max();
  bool orderFromFile = false;
  if (!m_options->in_orderingFile.empty()) {
    orderFromFile = m_problem->parseOrdering(m_options->in_orderingFile, elim);
  }

  // Init. pseudo tree
  m_pseudotree.reset(new Pseudotree(m_problem.get(), m_options->subprobOrder));

  if (NULL != m_options->varOrder) {  // Order passed in through program options
    m_problem->parseOrdering(*m_options->varOrder, elim);
    m_pseudotree->build(g, elim, m_options->cbound);
    w = m_pseudotree->getWidth();
    cout << "Using provided elimination ordering (" << w << '/'
         << m_pseudotree->getHeight() << ")." << endl;
  } else if (orderFromFile) {  // Reading from file succeeded (i.e. file exists)
    m_pseudotree->build(g, elim, m_options->cbound);
    w = m_pseudotree->getWidth();
    cout << "Read elimination ordering from file " << m_options->in_orderingFile
         << " (" << w << '/' << m_pseudotree->getHeight() << ")." << endl;
  } else {
    if (m_options->order_timelimit == NONE)
      // compute at least one
      m_options->order_iterations = max(1, m_options->order_iterations);
  }

  high_resolution_clock::time_point time_order_start, time_order_cur;
  double timediff = 0.0;
  time_order_start = high_resolution_clock::now();

  unique_ptr<ARE::Graph> cvo_graph;
  unique_ptr<ARE::Graph> cvo_master_graph;
  int temp_adj_var_space_exists = 0;
  ARE::AdjVar** temp_adj_var_space_size_extra_array = nullptr;

  if (m_options->order_cvo) {
    vector< const vector<int>*> fn_signatures;
    for (Function* f : m_problem->getFunctions()) {
      fn_signatures.push_back(&f->getScopeVec());
    }

    cvo_master_graph.reset(new ARE::Graph());
    cvo_master_graph->Create(m_problem->getN(), fn_signatures);
    cvo_master_graph->RNG().seed(m_options->seed);
    if (!cvo_master_graph->_IsValid) {
      return false;
    }
    temp_adj_var_space_size_extra_array = new ARE::AdjVar*[1000];

    cvo_master_graph->ComputeVariableEliminationOrder_Simple(0, INT_MAX, false,
        DBL_MAX, false, true, 1, 1, 0.0, temp_adj_var_space_exists,
        temp_adj_var_space_size_extra_array);
    cvo_master_graph->ReAllocateEdges();
    
  }
  // Search for variable elimination ordering, looking for min. induced
  // width, breaking ties via pseudo tree height
  cout << "Searching for elimination ordering,";
  if (m_options->order_cvo) {
    cout << " CVO,";
  }

  if (m_options->order_iterations != NONE)
    cout << " " << m_options->order_iterations << " iterations";
  if (m_options->order_timelimit != NONE)
    cout << " " << m_options->order_timelimit << " seconds";
  cout << ":" << flush;

  int iterCount = 0, sinceLast = 0;
  int remaining = m_options->order_iterations;

  while (true) {

    if (m_options->order_iterations != NONE && remaining == 0) break;

    vector<int> elimCand;   // new ordering candidate
    bool improved = false;  // improved in this iteration?
    int new_w;
    if (m_options->order_cvo) {
      cvo_graph.reset(new ARE::Graph());
      *cvo_graph = *cvo_master_graph;
      // seed based on daoopt's RNG; this lets us have deterministic iterations
      // subject to daoopt's seed.
      cvo_graph->RNG().seed(rand::next());
      int width_limit = INT_MAX;
      double space_limit = DBL_MAX;
      new_w = cvo_graph->ComputeVariableEliminationOrder_Simple(
          0, width_limit, false, space_limit, false, false, 10, 1, 0.0,
          temp_adj_var_space_exists, temp_adj_var_space_size_extra_array);
      if (new_w != 0) {
        cout << "WARNING! (errorcode: " << new_w << ")" << endl;
        new_w = INT_MAX;
      } else {
        new_w = cvo_graph->_VarElimOrderWidth;
        elimCand.assign(cvo_graph->_VarElimOrder,
                        cvo_graph->_VarElimOrder + cvo_graph->_nNodes);
      }
    } else {
      new_w =
          m_pseudotree->eliminate(g, elimCand, w, m_options->order_tolerance);
    }
    if (new_w < w) {
      elim = elimCand;
      w = new_w;
      improved = true;
      m_pseudotree->build(g, elimCand, m_options->cbound);
      cout << " " << iterCount << ':' << w << '/' << m_pseudotree->getHeight()
           << flush;
    } else if (new_w == w) {
      Pseudotree ptCand(m_problem.get(), m_options->subprobOrder);
      ptCand.build(g, elimCand, m_options->cbound);
      if (ptCand.getHeight() < m_pseudotree->getHeight()) {
        elim = elimCand;
        improved = true;
        m_pseudotree->build(g, elim, m_options->cbound);
        cout << " " << iterCount << ':' << w << '/' << m_pseudotree->getHeight()
             << flush;
      }
    }
    ++iterCount, ++sinceLast, --remaining;

    // Adaptive ordering scheme
    if (improved && m_options->autoIter && remaining > 0) {
      remaining = max(remaining, sinceLast + 1);
      sinceLast = 0;
    }

    // check termination conditions
    time_order_cur = high_resolution_clock::now();
    timediff = duration_cast<duration<double>>(time_order_cur -
                                               time_order_start).count();
    if (m_options->order_timelimit != NONE &&
        timediff > m_options->order_timelimit)
      break;
  }
  if (temp_adj_var_space_size_extra_array) {
    for (int i = temp_adj_var_space_exists - 1; i >= 0; --i) {
      delete[] temp_adj_var_space_size_extra_array[i];
    }
  }
  time_order_cur = high_resolution_clock::now();
  timediff = duration_cast<duration<double>>(time_order_cur - time_order_start)
                 .count();
  cout << endl << "Ran " << iterCount << " iterations (" << int(timediff)
       << " seconds), lowest width/height found: " << w << '/'
       << m_pseudotree->getHeight() << '\n';

  // Save order to file?
  if (!m_options->in_orderingFile.empty() && !orderFromFile) {
    m_problem->saveOrdering(m_options->in_orderingFile, elim);
    cout << "Saved ordering to file " << m_options->in_orderingFile << endl;
  }
#if defined PARALLEL_DYNAMIC || defined PARALLEL_STATIC
#if defined PARALLEL_STATIC
  if (!m_options->par_solveLocal &&
      !m_options->par_postOnly)  // no need to write ordering
#endif
  {
    m_options->in_orderingFile = string("temp_elim.") + m_options->problemName +
                                 string(".") + m_options->runTag +
                                 string(".gz");
    m_problem->saveOrdering(m_options->in_orderingFile, elim);
    cout << "Saved ordering to file " << m_options->in_orderingFile << endl;
  }
#endif

  // OR search?
  if (m_options->orSearch) {
    cout << "Rebuilding pseudo tree as chain." << endl;
    m_pseudotree->buildChain(g, elim, m_options->cbound);
  }

  // Pseudo tree has dummy node after build(), add to problem
  m_problem->addDummy();  // add dummy variable to problem, to be in sync with
                          // pseudo tree
  m_pseudotree->resetFunctionInfo(m_problem->getFunctions());
  m_pseudotree->addDomainInfo(m_problem->getDomains());

#if defined PARALLEL_STATIC || TRUE
  m_pseudotree->computeSubprobStats();
#endif
#if defined PARALLEL_DYNAMIC  //|| defined PARALLEL_STATIC
  int cutoff = m_pseudotree->computeComplexities(m_options->threads);
  cout << "Suggested cutoff:\t" << cutoff << " (ignored)" << endl;
//  if (opt.autoCutoff) {
//    cout << "Auto cutoff:\t\t" << cutoff << endl;
//    opt.cutoff_depth = cutoff;
//  }
#endif

  // Output pseudo tree to file for plotting?
  if (!m_options->out_pstFile.empty()) {
    m_pseudotree->outputToFile(m_options->out_pstFile);
    cout << "Saved pseudotree for plotting to file " << m_options->out_pstFile
         << "." << endl;
  }

  if (m_options->maxWidthAbort != NONE &&
      m_options->maxWidthAbort < m_pseudotree->getWidth()) {
    oss msg;
    msg << "Problem instance with width w=" << m_pseudotree->getWidth()
        << " is too complex, aborting (limit is w=" << m_options->maxWidthAbort
        << ")";
    err_txt(msg.str());
    return false;
  }
  return true;
}

bool Main::runSLS() {
  // allows to stop sls4mpe computation.
  sls4mpe::global_abort = false;
#ifdef PARALLEL_STATIC
  if (m_options->par_postOnly) return true;  // skip SLS for static post mode
#endif
#ifdef ENABLE_SLS
  if (!m_options->in_subproblemFile.empty())
    return true;  // no SLS in case of subproblem processing
  if (m_options->slsIter <= 0) return true;
  oss ss;
  ss << "Running SLS " << m_options->slsIter << " times for "
     << m_options->slsTime << " seconds" << endl;
  myprint(ss.str());

  switch (m_options->slsAlgo) {
    case 1:
      m_slsWrapper.reset(new SLSWrapper());
      break;
    case 2:
      m_slsWrapper.reset(new SLSWrapperHybrid());
      break;
    default:
      m_slsWrapper.reset(new SLSWrapperInc(m_options->slsConvergeRate));
      break;
  }
  // m_slsWrapper.reset(new SLSWrapper());

  m_slsWrapper->init(m_problem.get(), m_options->slsIter, m_options->slsTime);
  m_slsWrapper->run();
  myprint("SLS finished.\n");
#endif
  return true;
}

bool Main::stopSLS() {
  sls4mpe::global_abort = true;
  return true;
}

Heuristic* Main::newHeuristic(Problem* p, Pseudotree* pt, ProgramOptions* po) {
#ifdef NO_HEURISTIC
  return new Unheuristic;
#else
  if (po->fglpHeur || po->fglpMBEHeur) {
    bool useFGLP = true;

    /*
    // Decide if we shouldn't be using dynamic FGLP.
    if (po->fglpMBEHeurChoice) {
      // Try to see the ibound possible first
      // Temporarily set mplp to 0
      int store_mplp = po->mplp;
      int store_mplps = po->mplps;
      int store_ibound = po->ibound;
      po->mplp = 0;
      po->mplps = 0;
      unique_ptr<MiniBucketElim> test_heuristic(
          new MiniBucketElim(p, pt, po, po->ibound));
      test_heuristic->limitSize(po->memlimit, nullptr);
      int ib = test_heuristic->getIbound();

      // limitSize changes the ibound in po: restore it.
      po->ibound = store_ibound;

      po->mplp = store_mplp;
      po->mplps = store_mplps;
      if (ib < pt->getWidth() / 2) {
        cout << "ibound < w/2, using dynamic FGLP" << endl;
      } else {
        cout << "ibound >= w/2, using regular MBE-MM" << endl;
        useFGLP = false;
      }
    }
    */

    if (useFGLP) {
      if (po->fglpHeur) {
        return new FGLPHeuristic(p, pt, po);
      } else if (po->fglpMBEHeur) {
        return new FGLPMBEHybrid(p, pt, po);
      }
    }
  }

  if (po->lookaheadDepth > 0 || po->lookaheadSubtreeSizeLimit > 0 ||
      po->aobf_subordering == "static_be" ||
      po->aobf_subordering == "sampled_be" ||
      po->aobf_subordering == "sampled_st_be") {
    return new MiniBucketElimLH(p, pt, po, po->ibound);
  }
  if (po->incrementalJG)  // default value of the incrementalJG is false
    return new MiniBucketElimInc(p, pt, po, po->ibound);
  else
    return new MiniBucketElim(p, pt, po, po->ibound);
#endif
}

bool Main::initDataStructs() {

// The main search space
#ifdef PARALLEL_DYNAMIC
  m_space.reset(new SearchSpaceMaster(m_pseudotree.get(), m_options.get()));
#else
  if (m_options->algorithm == "aobb") {
    m_space.reset(new SearchSpace(m_pseudotree.get(), m_options.get()));
  } else if (m_options->algorithm == "aobf" ||
      m_options->algorithm == "aaobf") {
    m_space.reset(new BFSearchSpace(m_pseudotree.get(), m_options.get(),
                                    m_problem->getN()));
  }
  m_space->stats.numORVar.resize(m_pseudotree->getN(), 0);
  m_space->stats.numANDVar.resize(m_pseudotree->getN(), 0);
  m_space->stats.numProcORVar.resize(m_pseudotree->getN(), 0);
  m_space->stats.numProcANDVar.resize(m_pseudotree->getN(), 0);
#endif

  // Heuristic is initialized here, built later in compileHeuristic()
  m_heuristic.reset(
      newHeuristic(m_problem.get(), m_pseudotree.get(), m_options.get()));

  m_prop.reset(new BoundPropagator(m_problem.get(), m_space.get(),
                                   !m_options->nocaching));

// Main search engine
#if defined PARALLEL_DYNAMIC
  m_search.reset(new BranchAndBoundMaster(m_problem.get(), m_pseudotree.get(),
                                          m_space.get(),
                                          m_heuristic.get()));  // TODO
#elif defined PARALLEL_STATIC
  m_search.reset(new ParallelManager(m_problem.get(), m_pseudotree.get(),
                                     m_space.get(), m_heuristic.get()));
#else
  if (m_options->algorithm == "aobb") {
    if (m_options->rotate) {
      m_search.reset(new BranchAndBoundRotate(
          m_problem.get(), m_pseudotree.get(), m_space.get(), m_heuristic.get(),
          m_prop.get(), m_options.get()));
    } else {
      m_search.reset(new BranchAndBound(m_problem.get(), m_pseudotree.get(),
                                        m_space.get(), m_heuristic.get(),
                                        m_prop.get(), m_options.get()));
    }
  } else if (m_options->algorithm == "aobf") {
    m_search.reset(new AOStar(m_problem.get(), m_pseudotree.get(),
                              m_space.get(), m_heuristic.get(), m_prop.get(),
                              m_options.get()));
  } else if (m_options->algorithm == "aaobf") {
    m_search.reset(new AnytimeAOStar(m_problem.get(), m_pseudotree.get(),
          m_space.get(), m_heuristic.get(), m_prop.get(),
          m_options.get()));
  } else {
    cout << "Invalid algorithm option." << endl;
    return false;
  }
#endif

  // Subproblem specified? If yes, restrict.
  if (!m_options->in_subproblemFile.empty()) {
    if (m_options->in_orderingFile.empty()) {
      err_txt("Subproblem specified but no ordering given.");
      return false;
    } else {
      m_problem->setSubprobOnly();
      m_options->order_iterations = 0;
      cout << "Reading subproblem from file " << m_options->in_subproblemFile
           << '.' << endl;
      if (!m_search->restrictSubproblem(m_options->in_subproblemFile)) {
        err_txt("Subproblem restriction failed.");
        return false;
      }
    }
  }

  cout << "Induced width:\t\t" << m_pseudotree->getWidthCond() << " / "
       << m_pseudotree->getWidth() << endl;
  if (m_options->orSearch) {
    cout << "Pathwidth:\t\t" << m_pseudotree->getPathwidth() << endl;
  }
  cout << "Pseudotree depth:\t" << m_pseudotree->getHeightCond() << " / "
       << m_pseudotree->getHeight() << endl;
  cout << "Problem variables:\t" << m_pseudotree->getSizeCond() << " / "
       << m_pseudotree->getSize() << endl;
#ifdef PARALLEL_STATIC
  cout << "State space bound:\t" << m_pseudotree->getStateSpaceCond() << endl;
#endif
  cout << "Disconn. components:\t" << m_pseudotree->getComponentsCond() << " / "
       << m_pseudotree->getComponents() << endl;

#ifdef MEMDEBUG
  malloc_stats();
#endif

  return true;
}

bool Main::preprocessHeuristic() {
  Pseudotree* cur_pt = m_pseudotree.get();  // could be null
  const vector<val_t>* curAsg =
      (m_search.get()) ? &m_search->getAssignment() : NULL;

  if (cur_pt) {
    m_options->ibound = min(m_options->ibound, m_pseudotree->getWidthCond());
  }
  size_t sz = 0;

  // pseudo tree can be NULL right now!
  m_heuristic.reset(newHeuristic(m_problem.get(), cur_pt, m_options.get()));

  if (m_options->memlimit != NONE && cur_pt) {
    sz = m_heuristic->limitSize(m_options->memlimit, curAsg);
    sz *= sizeof(double) / (1024 * 1024.0);
    cout << "Enforcing memory limit resulted in i-bound " << m_options->ibound
         << " with " << sz << " MByte." << endl;
  }
  if (m_options->nosearch && !m_options->force_compute_tables) {
    cout << "Skipping heuristic preprocessing..." << endl;
    return false;
  }

  if (m_heuristic->preprocess(curAsg)) {
    if (cur_pt) {
      m_pseudotree->resetFunctionInfo(m_problem->getFunctions());
    }
  }
  return true;
}

bool Main::compileHeuristicBeforeBuild() {
  // m_options->ibound = min(m_options->ibound, m_pseudotree->getWidthCond());
  // m_options->jglpi = min(m_options->jglpi, m_pseudotree->getWidthCond());
  m_options->jglpi =
      m_pseudotree->getWidth();  // min(m_options->jglpi,
                                 // m_pseudotree->getWidthCond());
  size_t sz = 0;

  // do not limit ibound in before running classifier. it should be remained as
  // 0.
  if (m_options->memlimit != NONE) {
    sz =
        m_heuristic->limitSize(m_options->memlimit, &m_search->getAssignment());
    sz *= sizeof(double) / (1024 * 1024.0);
    cout << "Enforcing memory limit resulted in i-bound " << m_options->ibound
         << " with " << sz << " MByte." << endl;
  }

  if (m_options->nosearch && !m_options->force_compute_tables) {
    cout << "Simulating mini bucket heuristic..." << endl;
    sz = m_heuristic->build(&m_search->getAssignment(),
                            false);  // false = just compute memory estimate
  }

  return true;
}

bool Main::compileHeuristicAfterBuild() {

  // heuristic might have changed problem functions, pseudotree needs remapping
  m_pseudotree->resetFunctionInfo(m_problem->getFunctions());

  // set initial lower bound if provided (but only if no subproblem was
  // specified)
  if (m_options->in_subproblemFile.empty()) {
    if (m_options->in_boundFile.size()) {
      cout << "Loading initial lower bound from file "
           << m_options->in_boundFile << '.' << endl;
      if (!m_search->loadInitialBound(m_options->in_boundFile)) {
        err_txt("Loading initial bound failed");
        return false;
      }
    } else if (!ISNAN(m_options->initialBound)) {
#ifdef NO_ASSIGNMENT
      cout << "Setting external lower bound " << m_options->initialBound
           << endl;
      m_search->updateSolution(m_options->initialBound);
#else
      err_txt("Compiled with tuple support, value-based bound not possible.");
      return false;
#endif
    }
  }

#ifndef NO_HEURISTIC
  if (m_search->getCurOptValue() >= m_heuristic->getGlobalUB()) {
    m_solved = true;
    cout << endl << "--------- Solved during preprocessing ---------" << endl;
  } else if (m_heuristic->isAccurate()) {
    cout << endl << "Heuristic is accurate!" << endl;
    m_options->lds = 0;  // set LDS to 0 (sufficient given perfect heuristic)
    m_solved = true;
  }
#endif

  return true;
}

// added method for separating FGLP/JGLP/MBE compilation step
// assuming that we are using MiniBucketLimInc object,
// MiniBucketElimInc::build(int taks) method handles
// FGLP, JGLP, MBEMM separately,
// task = 0 run FGLP, task = 1 run JGLP, task =2 MBEMM, otherwise return false
// Each of them will be called by DaooptInterface methods as
// DaooptInterface::compileFGLP(), DaooptInterface::compileJGLP(),
// DaooptInterface::compileMbeAndFinishPreprocessing()
// ans such methods calls, _daoopt.compileFGLP() etc will be triggerd inside
// daoopt trhead function
bool Main::compileHeuristic(int task) {
  size_t sz = 0;
  bool mbFromFile = false;
  mbFromFile = m_heuristic->readFromFile(m_options->in_minibucketFile);

  if (m_options->nosearch) {
    return true;
  } else {

    if (task == 2) {
      _time_pre = high_resolution_clock::now();
      if (!m_options->in_minibucketFile.empty()) sz = m_heuristic->getSize();

      if (!mbFromFile) {
        cout << "Computing mini bucket heuristic..." << endl;
        cout << "(Moment matching: " << (m_options->match ? "yes" : "no") << ")"
             << endl;
      }
    }

    if (!mbFromFile)
      sz = m_heuristic->build(task, &m_search->getAssignment(),
                              true);  // true =  actually compute heuristic

    if (task == 2) {
      high_resolution_clock::time_point cur_time = high_resolution_clock::now();
      double time_passed =
          duration_cast<duration<double>>(cur_time - _time_pre).count();
      cout << "\tMini bucket finished in " << time_passed << " seconds" << endl;
    }

    if (task == 2 && !mbFromFile && !m_options->in_minibucketFile.empty()) {
      cout << "\tWriting mini bucket to file " << m_options->in_minibucketFile
           << " ..." << flush;
      m_heuristic->writeToFile(m_options->in_minibucketFile);
      cout << " done" << endl;
    }

    if (task == 2) {
      cout << '\t' << (sz / (1024 * 1024.0)) * sizeof(double) << " MBytes"
           << endl;
    }
    return true;
  }
}

// orignial method for compiling Heuristic that calls Build method for the
// Heuristic obejct
bool Main::compileHeuristic() {
  m_options->ibound = min(m_options->ibound, m_pseudotree->getWidthCond());
  m_options->jglpi = min(m_options->jglpi, m_pseudotree->getWidthCond());
  size_t sz = 0;
  if (m_options->memlimit != NONE) {
    sz =
        m_heuristic->limitSize(m_options->memlimit, &m_search->getAssignment());
    sz *= sizeof(double) / (1024 * 1024.0);
    cout << "Enforcing memory limit resulted in i-bound " << m_options->ibound
         << " with " << sz << " MByte." << endl;
  }

  if (m_options->nosearch && !m_options->force_compute_tables) {
    cout << "Simulating mini bucket heuristic..." << endl;
    sz = m_heuristic->build(&m_search->getAssignment(),
                            false);  // false = just compute memory estimate
  } else {
    _time_pre = high_resolution_clock::now();
    bool mbFromFile = false;
    if (!m_options->in_minibucketFile.empty()) {
      mbFromFile = m_heuristic->readFromFile(m_options->in_minibucketFile);
      sz = m_heuristic->getSize();
    }
    if (!mbFromFile) {
      cout << "Computing mini bucket heuristic..." << endl;
      cout << "(Moment matching: " << (m_options->match ? "yes" : "no") << ")"
           << endl;
      sz = m_heuristic->build(&m_search->getAssignment(),
                              true);  // true =  actually compute heuristic
      high_resolution_clock::time_point cur_time = high_resolution_clock::now();
      double time_passed =
          duration_cast<duration<double>>(cur_time - _time_pre).count();
      cout << "\tMini bucket finished in " << time_passed << " seconds" << endl;
    }
    if (!mbFromFile && !m_options->in_minibucketFile.empty()) {
      cout << "\tWriting mini bucket to file " << m_options->in_minibucketFile
           << " ..." << flush;
      m_heuristic->writeToFile(m_options->in_minibucketFile);
      cout << " done" << endl;
    }
  }
  cout << '\t' << (sz / (1024 * 1024.0)) * sizeof(double) << " MBytes" << endl;

  // heuristic might have changed problem functions, pseudotree needs remapping
  m_pseudotree->resetFunctionInfo(m_problem->getFunctions());

  // set initial lower bound if provided (but only if no subproblem was
  // specified)
  if (m_options->in_subproblemFile.empty()) {
    if (m_options->in_boundFile.size()) {
      cout << "Loading initial lower bound from file "
           << m_options->in_boundFile << '.' << endl;
      if (!m_search->loadInitialBound(m_options->in_boundFile)) {
        err_txt("Loading initial bound failed");
        return false;
      }
    } else if (!ISNAN(m_options->initialBound)) {
#ifdef NO_ASSIGNMENT
      cout << "Setting external lower bound " << m_options->initialBound
           << endl;
      m_search->updateSolution(m_options->initialBound);
#else
      err_txt("Compiled with tuple support, value-based bound not possible.");
      return false;
#endif
    }
  }

#ifndef NO_HEURISTIC
  if (m_search->getCurOptValue() >= m_heuristic->getGlobalUB()) {
    m_solved = true;
    cout << endl << "--------- Solved during preprocessing ---------" << endl;
  } else if (m_heuristic->isAccurate()) {
    cout << endl << "Heuristic is accurate!" << endl;
    m_options->lds = 0;  // set LDS to 0 (sufficient given perfect heuristic)
    m_solved = true;
  }
#endif

  return true;
}

bool Main::runLDS() {
#ifdef PARALLEL_STATIC
  if (m_options->par_postOnly) return true;  // skip LDS for static post mode
#endif
  // Run LDS if specified
  if (m_options->lds != NONE) {
    cout << "Running LDS with limit " << m_options->lds << endl;
    scoped_ptr<SearchSpace> spaceLDS(
        new SearchSpace(m_pseudotree.get(), m_options.get()));
    spaceLDS->stats.numORVar.resize(m_pseudotree->getN());
    spaceLDS->stats.numANDVar.resize(m_pseudotree->getN());
    spaceLDS->stats.numProcORVar.resize(m_pseudotree->getN(), 0);
    spaceLDS->stats.numProcANDVar.resize(m_pseudotree->getN(), 0);
    unique_ptr<BoundPropagator> propLDS(new BoundPropagator(
        m_problem.get(), spaceLDS.get(), false));  // doCaching = false
    LimitedDiscrepancy lds(m_problem.get(), m_pseudotree.get(), spaceLDS.get(),
                           m_heuristic.get(), propLDS.get(), m_options.get(),
                           m_options->lds);
    if (!m_options->in_subproblemFile.empty()) {
      if (!lds.restrictSubproblem(m_options->in_subproblemFile)) {
        err_txt("Subproblem restriction for LDS failed.");
        return false;
      }
    }

    // load current best solution into LDS
    if (lds.updateSolution(m_problem->getSolutionCost()
#ifndef NO_ASSIGNMENT
                           ,
                           m_problem->getSolutionAssg()
#endif
                           ))
      cout << "LDS: Initial solution loaded." << endl;

    lds.finalizeHeuristic();
    /*
    SearchNode* n = lds.nextLeaf();
    while (n) {
      propLDS.propagate(n,true); // true = report solution
      n = lds.nextLeaf();
    }
    */
    lds.solve(0);
    cout << "LDS: explored " << spaceLDS->stats.numExpOR << '/'
         << spaceLDS->stats.numExpAND << " OR/AND nodes" << endl;
    cout << "LDS: solution cost " << lds.getCurOptValue() << endl;

#ifndef NO_HEURISTIC
    if (m_search->getCurOptValue() >= m_heuristic->getGlobalUB()) {
      m_solved = true;
      cout << endl << "--------- Solved by LDS ---------" << endl;
    }
#endif
  }

  return true;
}

bool Main::finishPreproc() {

  // load current best solution from preprocessing into search instance
  if (m_search->updateSolution(m_problem->getSolutionCost()
#ifndef NO_ASSIGNMENT
                               ,
                               m_problem->getSolutionAssg()
#endif
                               ))
    cout << "Initial problem lower bound: " << m_search->curLowerBound()
         << endl;

#ifndef NO_HEURISTIC
  if (!m_options->nosearch || m_options->force_compute_tables)
    m_search->finalizeHeuristic();
#endif

#ifdef PARALLEL_STATIC
  if (m_options->par_preOnly) {
    m_search->storeLowerBound();
  } else if (m_options->par_postOnly) {
    m_search->loadLowerBound();
  }
#endif

  // Record time after preprocessing
  _time_pre = high_resolution_clock::now();
  double time_passed =
      duration_cast<duration<double>>(_time_pre - _time_start).count();
  cout << "Preprocessing complete: " << time_passed << " seconds" << endl;

  return true;
}

#ifdef PARALLEL_DYNAMIC
/* dynamic master mode for distributed execution */
bool Main::runSearchDynamic() {

  if (!m_solved) {
    // The propagation engine
    BoundPropagatorMaster prop(m_problem.get(), m_space.get());

    // take care of signal handling
    sigset_t new_signal_mask;  //, old_signal_mask;
    sigfillset(&new_signal_mask);
    pthread_sigmask(SIG_BLOCK, &new_signal_mask, NULL);  // block all signals

    try {
      CondorSubmissionEngine cse(m_space.get());

      boost::thread thread_prop(boost::ref(prop));
      boost::thread thread_bab(boost::ref(*m_search));
      boost::thread thread_cse(boost::ref(cse));

      // start signal handler
      SigHandler sigH(&thread_bab, &thread_prop, &thread_cse);
      boost::thread thread_sh(boost::ref(sigH));

      thread_bab.join();
      thread_prop.join();
      thread_cse.interrupt();
      thread_cse.join();
      cout << endl;
    }
    catch (...) {
      myerror("Caught signal during master execution, aborting.\n");
      return false;
    }

    // unblock signals
    pthread_sigmask(SIG_UNBLOCK, &new_signal_mask, NULL);

  } else {  // !m_solved
    BoundPropagator prop(m_problem.get(), m_space.get());
    SearchNode* n = m_search->nextLeaf();
    while (n) {
      prop.propagate(n, true);
      n = m_search->nextLeaf();
    }
  }
  return true;
}
#endif

#ifdef PARALLEL_STATIC
/* static master mode for distributed execution */
bool Main::runSearchStatic() {
  bool preOnly = m_options->par_preOnly, postOnly = m_options->par_postOnly,
       local = m_options->par_solveLocal;

  bool success = true;
  if (!postOnly) {
    success = success && m_search->doLearning();
  }
  if (true) {  // TODO
    /* find frontier from scratch */
    success = success && m_search->findFrontier();
  }
  if (!postOnly) {
    /* writes CSV with subproblem stats */
    success = success && m_search->writeSubprobStats();
  }
  if (!postOnly && !local) {
    /* generate files for subproblems */
    success = success && m_search->writeJobs();
    if (m_search->getSubproblemCount() == 0) m_solved = true;
  }
  if (local && !preOnly) {
    /* solve external subproblems locally */
    success = success && m_search->extSolveLocal();
  }
  if (!local && !preOnly && !postOnly) {
    /* run Condor and wait for results */
    success = success && m_search->runCondor();
  }
  if (!local && !preOnly) {
    /* read external results */
    success = success && m_search->readExtResults();
  }

  if (!success) {
    myprint("!!! Search failed. !!!\n");
    err_txt("Main search routine failed.");
    return false;
  }

  return true;
}
#endif

/* sequential mode or worker mode for distributed execution */
bool Main::runSearchWorker(size_t nodeLimit) {
  m_solved = m_search->solve(nodeLimit);
  return m_solved;
}

bool Main::outputStats() const {
  if (m_options->nosearch) {
    cout << "Found '-no_search', full search skipped, exiting." << endl;
    return true;
  }

  // Output cache statistics
  if (m_space->cache) {
    m_space->cache->printStats();
  }
  // Output search stats
  if (!m_search->printStats()) {

    cout << endl << "--------- Search done ---------" << endl;
    cout << "Problem name:       " << m_options->problemName << endl;
#if defined PARALLEL_DYNAMIC || defined PARALLEL_STATIC
    cout << "Condor jobs:        " << m_search->getSubproblemCount() << endl;
#endif
    cout << "OR nodes:           " << m_space->stats.numExpOR << endl;
    cout << "AND nodes:          " << m_space->stats.numExpAND << endl;
#if defined PARALLEL_STATIC
    cout << "OR external:        " << m_space->stats.numORext << endl;
    cout << "AND external:       " << m_space->stats.numANDext << endl;
#endif
    cout << "OR processed:       " << m_space->stats.numProcOR << endl;
    cout << "AND processed:      " << m_space->stats.numProcAND << endl;
    cout << "Leaf nodes:         " << m_space->stats.numLeaf << endl;
    cout << "Pruned nodes:       " << m_space->stats.numPruned << endl;
    cout << "Deadend nodes:      " << m_space->stats.numDead << endl;
    cout << "Deadend nodes (CP): " << m_space->stats.numDeadCP << endl;

#ifdef PARALLEL_STATIC
    if (m_options->par_preOnly && m_solved) {
      ofstream slvd("SOLVED");
      slvd << "SOLVED" << endl;
      slvd.close();
    }
#endif
  }

  high_resolution_clock::time_point time_end = high_resolution_clock::now();
  double time_passed =
      duration_cast<duration<double>>(time_end - _time_start).count();
  cout << "Time elapsed:       " << time_passed << " seconds" << endl;
  time_passed =
      duration_cast<duration<double>>(_time_pre - _time_start).count();
  cout << "Preprocessing:      " << time_passed << " seconds" << endl;
  time_passed = duration_cast<duration<double>>(time_end - _time_pre).count();
  cout << "Search:             " << time_passed << " seconds" << endl;
  cout << "-------------------------------" << endl;

  cout << endl;
  m_heuristic->printExtraStats();
  cout << endl;

#ifdef PARALLEL_STATIC
  if (!m_options->par_preOnly ||
      m_solved) {  // parallel static: only output if solved
#endif

    double mpeCost = m_problem->getSolutionCost();
    cout << setprecision(6);
    cout << SCALE_LOG(mpeCost) << " (" << SCALE_NORM(mpeCost) << ')' << endl;

    // Output node and leaf profiles per depth
    const vector<count_t>& prof = m_search->getNodeProfile();
    cout << endl << "p " << prof.size();
    for (vector<count_t>::const_iterator it = prof.begin(); it != prof.end();
         ++it) {
      cout << ' ' << *it;
    }
    const vector<count_t>& leaf = m_search->getLeafProfile();
    cout << endl << "l " << leaf.size();
    for (vector<count_t>::const_iterator it = leaf.begin(); it != leaf.end();
         ++it) {
      cout << ' ' << *it;
    }
    cout << endl;

    //  pair<size_t,size_t> noNodes = make_pair(m_space->nodesOR,
    // m_space->nodesAND);
    m_problem->outputAndSaveSolution(
        m_options->out_solutionFile, &m_space->stats,
        m_search->getNodeProfile(), m_search->getLeafProfile());
#ifdef PARALLEL_STATIC
  }
#endif

  cout << endl;
  return true;
}

int Main::outputStatsToFile() const {
  FILE* fp = m_options->_fpLogFile;
  if (!fp) {
    return 1;
  }

  if (m_options->nosearch) {
    fprintf(fp, "\nNo search option used, full search skipped, exiting.");
    return 0;
  }

  fprintf(fp, "\n--------- Search done ---------");
  fprintf(fp, "\nProblem name:  %s ", m_options->problemName.c_str());
  fprintf(fp, "\nOR nodes:      %lld", (int64)m_space->stats.numExpOR);
  fprintf(fp, "\nAND nodes:     %lld", (int64)m_space->stats.numExpAND);
  fprintf(fp, "\nOR processed:  %lld", (int64)m_space->stats.numProcOR);
  fprintf(fp, "\nAND processed: %lld", (int64)m_space->stats.numProcAND);
  fprintf(fp, "\nLeaf nodes:    %lld", (int64)m_space->stats.numLeaf);
  fprintf(fp, "\nPruned nodes:  %lld", (int64)m_space->stats.numPruned);
  fprintf(fp, "\nDeadend nodes: %lld", (int64)m_space->stats.numDead);

  high_resolution_clock::time_point time_end = high_resolution_clock::now();
  double time_passed =
      duration_cast<duration<double>>(time_end - _time_start).count();
  fprintf(fp, "\nTime elapsed:  %g seconds", time_passed);
  time_passed =
      duration_cast<duration<double>>(_time_pre - _time_start).count();
  fprintf(fp, "\nPreprocessing: %g seconds", time_passed);
  fprintf(fp, "\n-------------------------------");

  double mpeCost = m_problem->getSolutionCost();
  fprintf(fp, "\n%g (%g)\n", (double)(SCALE_LOG(mpeCost)),
          (double)(SCALE_NORM(mpeCost)));

  return 0;
}

bool Main::start() const {

  _time_start = high_resolution_clock::now();

  // compile version string
  string version = "DAOOPT ";
  version += VERSIONINFO;

#ifdef PARALLEL_DYNAMIC
  version += " PARALLEL-DYNAMIC (";
#elif defined PARALLEL_STATIC
  version += " PARALLEL-STATIC (";
#else
  version += " STANDALONE (";
#endif
  version += boost::lexical_cast<std::string>(sizeof(val_t) * 8);

#ifdef NO_ASSIGNMENT
  version += ") w/o assig.";
#else
  version += ") w. assig.";
#endif

#if defined(ANYTIME_DEPTH)
  version += " /dive";
#endif

#if defined LIKELIHOOD
  version += " /likelihood";
#endif

  ostringstream oss;

  oss << "------------------------------------------------------------------"
      << endl << version << endl
      << "  by Lars Otten, UC Irvine <lotten@ics.uci.edu>" << endl
      << "------------------------------------------------------------------"
      << endl << "(FGLP Dynamic Heuristic Components)" << endl
      << "  by William Lam, UC Irvine <willmlam@ics.uci.edu>" << endl
      << "------------------------------------------------------------------"
      << endl;

  cout << oss.str();

  return true;
}

bool Main::outputInfo() const {
  assert(m_options.get());

#if defined PARALLEL_DYNAMIC || defined PARALLEL_STATIC
  if (m_options->runTag == "") {
    m_options->runTag = "notag";
  }
#endif

  ostringstream oss;
  oss << "+ algorithm:\t" << m_options->algorithm << endl;
  oss << "+ i-bound:\t" << m_options->ibound << endl << "+ j-bound:\t"
      << m_options->cbound << endl << "+ Memory limit:\t" << m_options->memlimit
      << endl << "+ Suborder:\t" << m_options->subprobOrder << " ("
      << subprob_order[m_options->subprobOrder] << ")" << endl
      << "+ Random seed:\t" << m_options->seed << endl;
#if defined PARALLEL_DYNAMIC || defined PARALLEL_STATIC
  oss << "+ Cutoff depth:\t" << m_options->cutoff_depth << endl
      << "+ Cutoff size:\t" << m_options->cutoff_size << endl
      << "+ Max. workers:\t" << m_options->threads << endl << "+ Run tag:\t"
      << m_options->runTag << endl;
#else
  if (m_options->rotate)
    oss << "+ rotate:\ton (" << m_options->rotateLimit << ")" << endl;
  else
    oss << "+ rotate:\toff" << endl;
#endif

  cout << oss.str();
  return true;
}

double Main::getSolution() const { return m_problem->getSolutionCost(); }

#ifndef NO_ASSIGNMENT
const vector<val_t>& Main::getSolutionAssg() const {
  return m_problem->getSolutionAssg();
}

void Main::getSolutionAssgOrg(vector<val_t>& org_sol) const {
  m_problem->assignmentForOutput(m_problem->getSolutionAssg(), org_sol);
}
#endif

double computeAvgDepth(const vector<count_t>& before,
                       const vector<count_t>& after, int offset) {
  size_t leafCount = 0;
  offset = max(offset, 0);  // offset at least 0.
  for (size_t d = offset; d < before.size(); ++d)
    leafCount += after[d] - before[d];
  double avg = 0.0;
  for (size_t d = offset; d < before.size(); ++d)
    avg += (after[d] - before[d]) * (d - offset) * 1.0 / leafCount;
  return avg;
}

double Main::evaluate(SearchNode* node) const {
  assert(node && node->getType() == NODE_OR);

  int var = node->getVar();
  PseudotreeNode* ptnode = m_pseudotree->getNode(var);
  const SubprobStats* stats = ptnode->getSubprobStats();
  const SubprobFeatures* feats = node->getSubprobFeatures();

  // "dynamic subproblem features
  double ibnd = m_options->ibound;
  double ub = node->getHeur(), lb = node->getInitialBound();
  double rPruned = feats->ratioPruned, rLeaf = feats->ratioLeaf,
         rDead = feats->ratioLeaf;
  double avgNodeD = feats->avgNodeDepth, avgLeafD = feats->avgLeafDepth,
         avgBraDg = feats->avgBranchDeg;

  if (lb == -INFINITY) {
    myprint("Warning: No initial lower bound, will yield infinite estimate.\n");
  }

  // "static" subproblem properties
  double D = node->getDepth(), Vars = ptnode->getSubprobSize(),
         Leafs = stats->getLeafCount();
  double Wmax = stats->getClusterStats(stats->MAX),
         Wavg = stats->getClusterStats(stats->AVG),
         Wsdv = stats->getClusterStats(stats->SDV),
         Wmed = stats->getClusterStats(stats->MED);
  double WCmax = stats->getClusterCondStats(stats->MAX),
         WCavg = stats->getClusterCondStats(stats->AVG),
         WCsdv = stats->getClusterCondStats(stats->SDV),
         WCmed = stats->getClusterCondStats(stats->MED);
  double Kmax = stats->getDomainStats(stats->MAX),
         Kavg = stats->getDomainStats(stats->AVG),
         Ksdv = stats->getDomainStats(stats->SDV);
  double Hmax = stats->getDepthStats(stats->MAX),
         Havg = stats->getDepthStats(stats->AVG),
         Hsdv = stats->getDepthStats(stats->SDV),
         Hmed = stats->getDepthStats(stats->MED);

  /*
  double z = (lb == ELEM_ZERO) ? 0.0 : 2*(ub-lb);
  z += Hmax + 0.5 * log10(Vars);
  */

  double z =
      // LassoLars(alpha=0.01), degree=1, stats4.csv (10603 samples)
      +(1.64620e-04 * (ub)) + (3.83243e-01 * (ub - lb)) +
      (6.77123e-02 * (avgNodeD)) - (4.05910e-02 * (D)) +
      (3.78208e-03 * (Vars)) - (7.44384e-03 * (Leafs)) +
      (4.63889e-01 * (Wavg)) + (2.47618e-01 * (Wsdv)) -
      (1.68938e-01 * (WCmax)) + (9.91102e-02 * (Havg)) - (2.57769e-01 * (Hsdv));
  // LassoLars(alpha=0.01), degree=2, stats4.csv (10603 samples)
  //- ( 1.06230e-03 * (lb)*(avgNodeD))  + ( 1.69370e-03 * (lb)*(avgLeafD))  + (
  //1.39504e-03 * (lb)*(Vars))  - ( 7.41810e-03 * (lb)*(Leafs))  - ( 2.40857e-04
  //* (lb)*(Hmax))  - ( 3.04937e-03 * (ub)*(ub))  - ( 7.36406e-04 *
  //(ub)*(ub-lb))  + ( 9.39502e-04 * (ub)*(D))  - ( 8.06125e-04 * (ub)*(WCsdv))
  //- ( 1.54377e-02 * (ub-lb)*(ub-lb))  + ( 3.90821e-04 * (ub-lb)*(avgNodeD))  +
  //( 1.20563e-02 * (ub-lb)*(D))  - ( 6.50204e-04 * (ub-lb)*(Vars))  + (
  //2.60533e-04 * (ub-lb)*(Leafs))  - ( 3.72410e-02 * (ub-lb)*(WCmax))  + (
  //1.00336e-02 * (ub-lb)*(Hmax))  + ( 7.41045e-03 * (ub-lb)*(Hsdv))  + (
  //8.60623e-03 * (ub-lb)*(Hmed))  - ( 9.21925e-04 * (rPruned)*(Vars))  + (
  //1.46400e-02 * (rDead)*(Vars))  + ( 2.89068e-03 * (rLeaf)*(Vars))  + (
  //1.43511e-03 * (avgNodeD)*(avgNodeD))  - ( 5.70504e-04 * (avgNodeD)*(Vars))
  //+ ( 1.68691e-03 * (avgNodeD)*(Leafs))  - ( 1.79182e-03 * (avgNodeD)*(Hmed))
  //+ ( 3.72780e-04 * (avgLeafD)*(D))  + ( 4.70759e-04 * (avgLeafD)*(Vars))  - (
  //1.65965e-03 * (avgLeafD)*(Leafs))  + ( 1.51711e-03 * (avgLeafD)*(Wmax))  + (
  //5.29760e-03 * (avgLeafD)*(Wsdv))  - ( 1.72148e-03 * (avgLeafD)*(Hmax))  + (
  //1.15015e-02 * (avgLeafD)*(Havg))  - ( 7.82195e-03 * (avgLeafD)*(Hmed))  - (
  //5.89913e-03 * (avgBraDg)*(Vars))  + ( 4.15688e-03 * (D)*(D))  + (
  //1.90502e-04 * (D)*(Vars))  - ( 4.01007e-05 * (D)*(Leafs))  + ( 1.07401e-02 *
  //(D)*(WCmax))  + ( 9.35496e-05 * (Vars)*(Vars))  - ( 1.28385e-03 *
  //(Vars)*(Wmax))  - ( 1.09807e-04 * (Vars)*(Wmed))  + ( 1.14059e-03 *
  //(Vars)*(WCmax))  + ( 3.61452e-03 * (Vars)*(WCavg))  + ( 2.92006e-04 *
  //(Vars)*(WCsdv))  - ( 1.12599e-02 * (Vars)*(Ksdv))  + ( 1.38700e-04 *
  //(Vars)*(Havg))  - ( 4.01527e-04 * (Vars)*(Hsdv))  - ( 4.28666e-04 *
  //(Vars)*(Hmed))  - ( 1.19686e-03 * (Leafs)*(Leafs))  - ( 8.13292e-03 *
  //(Leafs)*(WCmax))  + ( 3.42008e-03 * (Leafs)*(Havg))  - ( 5.46289e-03 *
  //(Wmax)*(Hmax))  + ( 1.62756e-02 * (WCmax)*(WCmax))  + ( 4.64474e-05 *
  //(Hmed)*(Hmed));

  if (z < 0.0 || z > 100.0) {
    oss ss;
    ss << "evaluate: unreasonable estimate for node " << *node << ": " << z
       << endl;
    myprint(ss.str());
  }

  node->setComplexityEstimate(z);
  DIAG(oss ss; ss << "eval " << *node << " : " << z << endl; myprint(ss.str()));

  return z;
}

double Main::runEstimation(size_t nodeLimit) {
  SearchNode* node = m_space->getTrueRoot();
  int var = node->getVar();
  PseudotreeNode* ptnode = m_pseudotree->getNode(var);
  const SubprobStats* stats =
      ptnode->getSubprobStats();  // TODO: where would this get set?

  // Copy current node profiles and search stats
  vector<count_t> startNodeP = m_search->getNodeProfile();
  vector<count_t> startLeafP = m_search->getLeafProfile();
  SearchStats startStats = m_space->stats;

  node->setInitialBound(m_search->getNodeLowerBound(node));

  SearchStats& currStats = m_space->stats;
  if (runSearchWorker(nodeLimit)) {
    // Search finished, return 0.0
    return 0.0;
  }

  // Compute features
  SubprobFeatures* features = node->getSubprobFeatures();
  features->ratioPruned =
      (currStats.numPruned - startStats.numPruned) * 1.0 / nodeLimit;
  features->ratioDead =
      (currStats.numDead - startStats.numDead) * 1.0 / nodeLimit;
  features->ratioLeaf =
      (currStats.numLeaf - startStats.numLeaf) * 1.0 / nodeLimit;

  features->avgNodeDepth =
      computeAvgDepth(startNodeP, m_search->getNodeProfile(), node->getDepth());
  features->avgLeafDepth =
      computeAvgDepth(startLeafP, m_search->getLeafProfile(), node->getDepth());
  features->avgBranchDeg = pow(nodeLimit, 1.0 / features->avgLeafDepth);

  // Compute estimate
  double estimate = evaluate(node);

  return estimate;
}

bool Main::exportProblemStats(int& fN, int& fF, int& fK, int& fS, int& fW,
                              int& fH) {
  fN = m_problem->getN();
  fF = m_problem->getC();
  fK = m_problem->getK();
  fS = m_problem->getR();
  fW = m_pseudotree->getWidth();
  fH = m_pseudotree->getHeight();
  return true;
}

}  // namespace daoopt
// EoF
