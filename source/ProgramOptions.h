/*
 * ProgramOptions.h
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
 *  Created on: Nov 15, 2008
 *      Author: Lars Otten <lotten@ics.uci.edu>
 */

#ifndef PROGRAMOPTIONS_H_
#define PROGRAMOPTIONS_H_

#include <gflags/gflags.h>
#include "_base.h"


#include <string>
#include <iostream>

namespace daoopt {

struct ProgramOptions {
public:
  bool nosearch; // abort before starting the actual search
  bool nocaching; // disable caching
  bool autoCutoff; // enable automatic cutoff
  bool autoIter; // enable adaptive ordering limit
  bool orSearch; // use OR search (builds pseudo tree as chain)
  bool par_solveLocal; // solve all parallel subproblems locally
  bool par_preOnly; // static parallel: preprocessing only (generate subproblems)
  bool par_postOnly; // static parallel: postprocessing only (read solution files)
  bool rotate; // enables breadth-rotating AOBB
  bool match; // uses moment matching during MBE
  int mplp;  // enables MPLP in Alex Ihler's MBE library (# iters)
  double mplps;  // enables MPLP in Alex Ihler's MBE library (# sec)
  double mplpt;  // convergence tolerance for MPLP
  int jglp;  // enables JGLP tightening in Alex Ihler's MBE library (# iters)
  double jglps;  // enables JGLP tightening in Alex Ihler's MBE library (# sec)
  int jglpi;  // specifies the i-bound used for JGLP
  int ibound; // bucket elim. i-bound
  int cbound; // cache context size bound
  int cbound_worker; // cache bound for worker processes
  int threads; // max. number of parallel subproblems
  int order_iterations; // no. of randomized order finding iterations
  int order_timelimit; // no. of seconds to look for variable ordering
  int order_tolerance; // allowed range of deviation from suggested optimal minfill heuristic
  int cutoff_depth; // fixed cutoff depth for central search
  int cutoff_width; // fixed width for central cutoff
  int nodes_init; // number of nodes for local initialization (times 10^6)
  int memlimit; // memory limit (in MB)
  int cutoff_size; // fixed cutoff subproblem size (times 10^6)
  int local_size; // lower bound for problem size to be solved locally (times 10^6)
  int maxSubprob; // only generate this many subproblems, then abort (for testing)
  int lds;  // run initial LDS with this limit (-1: enabled)
  int seed; // the seed for the random number generator
  int rotateLimit; // how many nodes to expand per subproblem stack before rotating
  int subprobOrder; // subproblem ordering, integers defined in _base.h
  int sampleDepth; // max. depth for randomness in sampler (will follow heuristic otherwise)
  int sampleScheme; // sampling scheme (TBD)
  int sampleRepeat; // how many times to repeat the sample size sequence
  int maxWidthAbort; // upper bound for induced width, abort if above this
  int slsIter; // number of SLS iterations for initial lower bound
  int slsTime; // time per SLS iteration (in seconds)
  int aobbLookahead;  // max. number of nodes for parallel static AOBB subproblem lookahead

  /* CVO */
  bool order_cvo; // Use Kalev Kask's variable ordering code
  int cvo_n_random_pick;
  double cvo_e_random_pick;

  /* DYNAMIC COST SHIFTING HEURISTIC OPTIONS */
  bool fglpHeur; // use pure FGLP heuristic
  bool fglpMBEHeur; // use FGLP/MBE hybrid heuristic
  bool fglpMBEHeurChoice; // use FGLP/MBE choice heuristic

  bool useShiftedLabels; // use shifted labels induced by FGLP
  bool useNullaryShift; // use FGLP update that shifts maximums into a nullary function
  bool usePriority; // use prioritized FGLP update schedule
  int ndfglp; // enables FGLP computation at every node dynamic MBE is used (# iters per node)
  double ndfglps; // enables FGLP computation at every node dynamic MBE is used (# seconds per node)
  double ndfglpt; // convergence tolerance for FGLP

  /* LOOKAHEAD HEURISTIC OPTIONS */
  int lookaheadDepth; // depth of lookahead when computing the h (heuristic) function; 0=no lookahead
  double lookahead_LE_SingleTableLimit; // as number of entries; log10.
  double lookahead_LE_AllTablesTotalLimit; // as number of entries; log10.
  double lookahead_LE_IgnoreThreshold ; // if bucket error is larger than this, it will be ignored for by lookahead. default is DBL_MIN.

  /* MISC OPTIONS */
  bool collapse; // collapse functions with identical scopes onto each other
  double perturb; // sets zeros to this value in the functions

  // Added to allow program to terminate itself.
  int maxTime; // timeout threshold (seconds)

  double initialBound; // initial lower bound

  char* problemSpec;  // problem specification in UAI format
  size_t problemSpec_len;  // length of the problemSpec buffer
  char* evidSpec;  // evidence specification in UAI format
  size_t evidSpec_len;  // length of the evidSpec buffer

  const std::vector<int>* varOrder;  // externally provided variable elimination
                                     // ordering

  std::string executableName; // name of the executable
  std::string problemName; // name of the problem
  std::string runTag; // string tag of this particular run
  std::string sampleSizes; // Sequence of sizes for subproblem samples (for prediction),
                           // comma-separated list of doubles.
  std::string in_problemFile; // problem file path
  std::string in_evidenceFile; // evidence file path
  std::string in_orderingFile; // ordering file path
  std::string in_minibucketFile; // minibucket file path
  std::string in_subproblemFile; // subproblem file path
  std::string in_boundFile; // file with initial lower bound (from SLS, e.g.)
  std::string out_solutionFile; // file path to write solution to
  std::string out_boundFile; // file to output best lower bound to
  std::string out_reducedFile; // file to save reduced network to
  std::string out_pstFile; // file to output pseudo tree description to (for plotting)

  std::string log_file; // logging: currently only used by lookahead
  FILE* _fpLogFile; // logging

public:
  ProgramOptions();
};

ProgramOptions* parseCommandLine(int argc, char** argv);

inline ProgramOptions::ProgramOptions() 
  : nosearch(false), nocaching(false), autoCutoff(false), autoIter(false),
    orSearch(false), par_solveLocal(false), par_preOnly(false),
    par_postOnly(false), rotate(false), ibound(0), cbound(0), cbound_worker(0),
    threads(0), order_iterations(0), order_timelimit(0), order_tolerance(0),
    cutoff_depth(NONE), cutoff_width(NONE),
    nodes_init(NONE), memlimit(NONE),
    cutoff_size(NONE), local_size(NONE), maxSubprob(NONE),
    lds(NONE), seed(NONE), rotateLimit(0), subprobOrder(NONE),
    sampleDepth(NONE), sampleScheme(NONE), sampleRepeat(NONE),
    maxWidthAbort(NONE), slsIter(0), slsTime(5),
    aobbLookahead(0),
    initialBound(ELEM_NAN),
    problemSpec(NULL), problemSpec_len(0),
    evidSpec(NULL), evidSpec_len(0), varOrder(NULL),
    order_cvo(false), cvo_n_random_pick(-1), cvo_e_random_pick(0.0),
    fglpHeur(false), fglpMBEHeur(false),
    mplp(0), mplps(0),
    jglp(0), jglps(0),
    useShiftedLabels(false), useNullaryShift(false), usePriority(false),
    ndfglp(0), ndfglps(0),
    match(false),
    collapse(false),
    perturb(0),
    lookaheadDepth(0), lookahead_LE_SingleTableLimit(-1.0),
    lookahead_LE_AllTablesTotalLimit(-1.0),
	lookahead_LE_IgnoreThreshold(DBL_MIN), 
    maxTime(kint32max),

    _fpLogFile(nullptr) {
}

}  // namespace daoopt

#endif /* PROGRAMOPTIONS_H_ */
