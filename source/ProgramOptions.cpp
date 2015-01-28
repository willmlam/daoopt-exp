/*
 * ProgramOptions.cpp
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

#include "ProgramOptions.h"
#include "google-utils/base/integral_types.h"
#include "google-utils/strings/strcat.h"

DEFINE_string(input_file, "", "path to problem file (required)");
DEFINE_string(evid_file, "", "path to evidence file (optional)");
DEFINE_string(ordering_file, "",
              "read elimination ordering from this file (first to last), or"
              "write elimination ordering to this file if it does not exist");
DEFINE_bool(adaptive, false, "enable adaptive ordering scheme");
DEFINE_int32(max_time, kint32max, "timeout threshold in seconds");
DEFINE_string(minibucket_file, "", "path to read/store minibucket heuristic");
DEFINE_string(subproblem_file, "",
              "path to subproblem specification to limit search to");
DEFINE_int32(suborder, 0, 
             "subproblem order "
             "(0:width-inc 1:width-dec 2:heur-inc 3:heur-dec)");
DEFINE_string(sol_file, "", "path to output optimal solution to");
DEFINE_string(out_bound_file, "", "path to output current best solution to");

DEFINE_int32(order_iterations, 25, "iterations for finding ordering");
DEFINE_int32(order_time, -1, "maximum time for finding ordering");
DEFINE_int32(order_tolerance, 0,
             "allowed deviation from minfill suggested optimal");
DEFINE_int32(max_width, -1, "max. induced width to process, abort otherwise");

DEFINE_bool(cvo, false, "use Kalev Kask's CVO ordering");
DEFINE_int32(cvo_n_random_pick, -1, "CVO: pool size");
DEFINE_double(cvo_e_random_pick, 0.0,
             "CVO: parameter for randomly choosing from the pool "
             "0 - uniform distribution, <0 - prefer lower cost, "
             ">0 - prefer higher cost");

DEFINE_int32(ibound, 10, "i-bound for minibucket heuristics");
DEFINE_int32(cbound, 1000, "context size bound for caching");

DEFINE_int32(fglp_iterations, -1, "do FGLP preprocessing (# iterations)");
DEFINE_double(fglp_time, -1.0, "do FGLP preprocessing (# seconds)");
DEFINE_double(fglp_tolerance, 1e-7, "convergence tolerance for FGLP");
DEFINE_int32(jglp_iterations, -1, "do JGLP preprocessing (# iterations)");
DEFINE_double(jglp_time, -1.0, "do JGLP preprocessing (# seconds)");
DEFINE_int32(jglp_ibound, -1, "i-bound for JGLP");

DEFINE_bool(rotate, false, "use breadth-rotating AOBB");
DEFINE_int32(rotate_limit, 1000,
             "nodes per subproblem stack rotation (0: disabled)");

DEFINE_string(bound_file, "", "file with initial lower bound on solution cost");
DEFINE_double(initial_bound, ELEM_NAN, "initial lower bound on solution cost");


#ifdef ENABLE_SLS
DEFINE_int32(sls_iterations, 0, "Number of initial SLS iterations");
DEFINE_int32(sls_time, 5, "Time per SLS iteration");
#endif
DEFINE_int32(lds_limit, -1,
             "run initial LDS search with given limit (-1: disabled)");

DEFINE_int32(mem_limit, -1, "approximate memory limit for minibuckets (in MB)");
DEFINE_int32(seed, -1, "seed for random number generator, time() otherwise");

DEFINE_bool(or_search, false, "use OR search (build pseudotree as a chain)");
DEFINE_bool(no_caching, false, "disable context-based caching during search");
DEFINE_bool(no_search, false, "perform preprocessing, output stats, and exit");
DEFINE_bool(match, false, "use moment-matching during MBE");

DEFINE_string(reduce_file, "",
              "path to output reduced network to (removes evidence and unary "
              "variables)");
DEFINE_bool(collapse, false,
            "collapse functions with identical scopes onto each other");
DEFINE_double(zero_perturb, 0, "set all zero values to this value");

DEFINE_bool(heuristic_fglp, false, "use pure FGLP dyanmic heuristic");
DEFINE_bool(heuristic_fglp_mbe_hybrid, false,
            "use FGLP dynamic/MBE hybrid heuristic");
DEFINE_bool(heuristic_fglp_mbe_choice, false,
            "use pure FGLP dyanmic heuristic if the i-bound is less "
            "than half the width; otherwise, use MBE");

DEFINE_bool(dfglp_shifted_labels, false, "use shifted labels induced by FGLP");
DEFINE_bool(dfglp_nullary_shift, false,
            "use FGLP update that shifts maximums into a nullary function");
DEFINE_bool(fglp_schedule_priority, false,
            "use FGLP update schedule with priorities");


DEFINE_int32(dfglp_iterations, -1, "# iterations of FGLP at every node");
DEFINE_double(dfglp_time, -1.0, "time for FGLP at every node");
DEFINE_double(dfglp_tolerance, 1e-7, "convergence tolerance for dynamic FGLP");

// Lookahead optinas
DEFINE_int32(lookahead_depth, 0, "depth of lookahead when computing the h"
             "(heuristic) function; 0=no lookahead");
DEFINE_double(lookahead_local_error_single_table_limit, 7.0, 
              "lookahead: limit as number of entries for a single local error" 
              "table (in log10)");
DEFINE_double(lookahead_local_error_all_tables_total_limit, 8.0, 
              "lookahead: limit as number of entries for all local error" 
              "tables (in log10)");
DEFINE_double(lookahead_local_error_ignore_threshold, -DBL_MIN, 
              "lookahead: ignore threshold");

DEFINE_string(pst_file, "", "path to output the pseudotree to, for plotting");

namespace daoopt {

ProgramOptions* parseCommandLine(int ac, char** av) {

  ProgramOptions* opt = new ProgramOptions;

  // executable name
  opt->executableName = av[0];

  try{
    string usage_message(opt->executableName);
    usage_message.append(" -input_file <problem.uai>");
    gflags::SetUsageMessage(usage_message);
    gflags::ParseCommandLineFlags(&ac, &av, true);

    if (FLAGS_input_file.empty()) {
      cout << "No input file given." << endl;
      return NULL;
    }

    opt->in_problemFile = FLAGS_input_file;
    opt->in_evidenceFile = FLAGS_evid_file;
    opt->in_orderingFile = FLAGS_ordering_file;
    opt->autoIter = FLAGS_adaptive;
    opt->maxTime = FLAGS_max_time;
    opt->in_subproblemFile = FLAGS_subproblem_file;
    opt->out_solutionFile = FLAGS_sol_file;
    opt->out_boundFile = FLAGS_out_bound_file;
    opt->in_minibucketFile = FLAGS_minibucket_file;
    opt->subprobOrder = FLAGS_suborder;
    if (opt->subprobOrder < 0 || opt->subprobOrder > 3) {
      cout << "Invalid subproblem order" << endl;
      exit(0);
    }

    opt->ibound = FLAGS_ibound;
    opt->cbound = FLAGS_cbound;
    opt->cbound_worker = FLAGS_cbound;
    opt->mplp = FLAGS_fglp_iterations;
    opt->mplps = FLAGS_fglp_time;
    opt->mplpt = FLAGS_fglp_tolerance;
    opt->jglp = FLAGS_jglp_iterations;
    opt->jglps = FLAGS_jglp_time;
    opt->jglpi = FLAGS_jglp_ibound;

    opt->fglpHeur = FLAGS_heuristic_fglp;
    opt->fglpMBEHeur = FLAGS_heuristic_fglp_mbe_hybrid;
    opt->fglpMBEHeurChoice = FLAGS_heuristic_fglp_mbe_choice;

    opt->useShiftedLabels = FLAGS_dfglp_shifted_labels;
    opt->useNullaryShift = FLAGS_dfglp_nullary_shift;
    opt->usePriority = FLAGS_fglp_schedule_priority;

    opt->ndfglp = FLAGS_dfglp_iterations;
    opt->ndfglps = FLAGS_dfglp_time;
    opt->ndfglpt = FLAGS_dfglp_tolerance;

    opt->lookaheadDepth = FLAGS_lookahead_depth;
    opt->lookahead_LE_SingleTableLimit =
        FLAGS_lookahead_local_error_single_table_limit;
    opt->lookahead_LE_AllTablesTotalLimit =
        FLAGS_lookahead_local_error_all_tables_total_limit;
    opt->lookahead_LE_IgnoreThreshold = 
        FLAGS_lookahead_local_error_ignore_threshold;

    opt->order_iterations = FLAGS_order_iterations;
    opt->order_timelimit = FLAGS_order_time;
    opt->order_tolerance = FLAGS_order_tolerance;
    opt->maxWidthAbort = FLAGS_max_width;
    opt->order_cvo = FLAGS_cvo;
    opt->cvo_n_random_pick = FLAGS_cvo_n_random_pick;
    opt->cvo_e_random_pick = FLAGS_cvo_e_random_pick;

    opt->in_boundFile = FLAGS_bound_file;
    opt->initialBound = FLAGS_initial_bound;

    opt->lds = FLAGS_lds_limit;

#ifdef ENABLE_SLS
    opt->slsIter = FLAGS_sls_iterations;
    opt->slsTime = FLAGS_sls_time;
#endif

    opt->memlimit = FLAGS_mem_limit;

    opt->orSearch = FLAGS_or_search;
    opt->nosearch = FLAGS_no_search;
    opt->nocaching = FLAGS_no_caching;
    opt->match = FLAGS_match;

    opt->rotate = FLAGS_rotate;
    opt->rotateLimit = FLAGS_rotate_limit;

    opt->seed = FLAGS_seed;

    opt->out_reducedFile = FLAGS_reduce_file;

    opt->collapse = FLAGS_collapse;
    
    opt->perturb = FLAGS_zero_perturb;

    opt->out_pstFile = FLAGS_pst_file;


    if (!FLAGS_subproblem_file.empty() && FLAGS_ordering_file.empty()) {
      cerr << "Error: Specifying a subproblem requires reading a fixed ordering from file." << endl;
      exit(1);
    }

    {
      //Extract the problem name
      string& fname = opt->in_problemFile;
      size_t len, start, pos1, pos2;
#if defined(WIN32)
      pos1 = fname.find_last_of("\\");
#elif defined(LINUX)
      pos1 = fname.find_last_of("/");
#endif
      pos2 = fname.find_last_of(".");
      if (pos1 == string::npos) { len = pos2; start = 0; }
      else { len = (pos2-pos1-1); start = pos1+1; }
      opt->problemName = fname.substr(start, len);
    }

  } catch (exception& e) {
    cerr << "error: " << e.what() << endl;
    return NULL;
  } catch (...) {
    cerr << "Exception of unknown type!" << endl;
    exit(1);
  }

  return opt;

}

}  // namespace daoopt

