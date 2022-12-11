#include "AOStar.h"
#include "BFSearchSpace.h"

#include <chrono>
using namespace std::chrono;

//#undef DEBUG

namespace daoopt {

extern high_resolution_clock::time_point _time_start; // From Main.cpp

bool AOStar::solve(size_t nodeLimit) {
  // for when AOStar search runs out of memory
  // allocate 64K memory.
  char* emergency_memory = new char[65536];

  bool solved = false;

  // Initial best-first search phase
  try {
    InitBFSearchSpace();
    solved = DoSearch();
  } catch (std::bad_alloc& ba_exception) {
    delete[] emergency_memory;
    emergency_memory = nullptr;
    best_first_limit_reached_ = true;
    solved = false;
  } catch (...) {
    std::cout << "Non bad alloc error?" << std::endl;
    solved = false;
  }

  // TODO(lamw): If the problem is not already solved, use AOBB the rest of
  // the way (would have returned earlier)

  if (emergency_memory) delete[] emergency_memory;

  if (timed_out_ || best_first_limit_reached_) {
    high_resolution_clock::time_point time_now = high_resolution_clock::now();
    double time_elapsed =
      duration_cast<duration<double>>(time_now - _time_start).count();
    if (timed_out_) {
      cout << "TIMED OUT at ";
    } else if (best_first_limit_reached_) {
      cout << "OUT OF MEMORY at ";
    }
    cout << time_elapsed << " seconds." << endl;
    printStats();
  }

  return solved;
}

bool AOStar::DoSearch() {
  BFSearchNode* root = dynamic_cast<BFSearchNode*>(search_space_->getRoot());
  double h = m_heuristic->getGlobalUB();
  double w = ELEM_ONE;
  double* dv = new double[2];
  dv[0] = h;
  dv[1] = w;

  root->setHeur(h);
  root->setValue(h);
  root->set_terminal(false);
  root->set_fringe(true);
  root->setHeurCache(dv);

  solution_cost_ = ELEM_ZERO;
  heuristic_bound_ = h;
  prev_reported_time_ = -1;

  while(!root->is_solved()) {
    assert(tip_nodes_.size() > 0);
    high_resolution_clock::time_point time_now = high_resolution_clock::now();
    double time_elapsed =
      duration_cast<duration<double>>(time_now - _time_start).count();
    if (time_elapsed > m_options->maxTime) {
      timed_out_ = true;
      return false;
    }
    ArrangeTipNodes();
    BFSearchNode* n = ChooseTipNode();
    ExpandAndRevise(n);
  }

  if (root->is_solved()) {
    solution_cost_ = root->getValue();
    m_problem->updateUpperBound(heuristic_bound_, &(search_space_->stats));
#ifdef NO_ASSIGNMENT
    m_problem->updateSolution(solution_cost_, &(search_space_->stats));
#else
    m_problem->updateSolution(solution_cost_, m_assignment, &(search_space_->stats));
#endif
  }

  return true;
}

void AOStar::ExpandAndRevise(BFSearchNode* node) {
  assert(node->is_fringe() && !node->is_solved());
  Expand(node);

#ifdef DEBUG
  cout << "Expanded " << node->ToString() << endl;
#endif

  // Revise
  std::multiset<BFSearchNode*, CompNodeIndexAsc> revise_set;

  revise_set.insert(node);
  node->set_visited(true);

  while (!revise_set.empty()) {
    BFSearchNode *e = *revise_set.begin();
    revise_set.erase(revise_set.begin());
    e->set_visited(false);

    assert(e->index() == 0);
    bool change = Revise(e);

#ifdef DEBUG
    cout << "Revised " << e->ToString() << (change ? " (changed)" : "") << endl;
#endif

    if (change) {
      // Step 12 Nilssons
      if (e->getType() == NODE_AND) {
        assert(e->parents().size() == 1);
        BFSearchNode* parent = e->parents().front();

        // Count children of e still in the revise_set
        // The index will place this node later the more children there are.
        size_t index = 0;
        for (BFSearchNode* child : parent->children()) {
          if (child->is_visited()) {
            ++index;
          }
        }
        BFSearchNode* best = parent->best_child();
        bool found = (best == e); // Is 'e' the marked AND child of the parent?
        if (parent->is_visited()) {
          // decrease parent index and replace
          auto si = revise_set.begin();
          for (; si != revise_set.end(); ++si) {
            if (parent == *si) {
              if (parent->index() > 0) {
                revise_set.erase(si);
                parent->decrement_index();
                revise_set.insert(parent);
              }
              break;
            }
          }
        } else if (found) {
          parent->set_index(index);
          parent->set_visited(true);
          revise_set.insert(parent);
        }
      } else if (e->getType() == NODE_OR) {
        // Multiple parents in the CMAO graph.
        for (BFSearchNode* parent : e->parents()) {
          size_t index = 0;
          for (BFSearchNode* child : parent->children()) {
            if (child->is_visited()) {
              ++index;
            }
          }

          if (parent->is_visited()) {
            auto si = revise_set.begin();
            for (; si != revise_set.end(); ++si) {
              if (parent == *si) {
                if (parent->index() > 0) {
                  revise_set.erase(si);
                  parent->decrement_index();
                  revise_set.insert(parent);
                }
                break;
              }
            }
          } else {
            parent->set_index(index);
            parent->set_visited(true);
            revise_set.insert(parent);
          }
        }
      }
    } else {
      assert(e->getType() != NODE_AND || e->parents().size() == 1);
      for (BFSearchNode* parent : e->parents()) {
        if (parent->is_visited()) {
          auto si = revise_set.begin();
          for (; si != revise_set.end(); ++si) {
            if (parent == *si) {
              if (parent->index() > 0) {
                revise_set.erase(si);
                parent->decrement_index();
                revise_set.insert(parent);
              }
              break;
            }
          }
        }
      }
    }
  }
  assert(revise_set.empty());

  FindBestPartialTree();

  assert(dynamic_cast<BFSearchNode*>(search_space_->getRoot())->is_solved() ||
         tip_nodes_.size() > 0);
}

bool AOStar::Expand(BFSearchNode* node) {
  bool no_children = true;
  int var = node->getVar();
  int depth = node->getDepth();
  PseudotreeNode* pt_node = m_pseudotree->getNode(var);
  if (node->getType() == NODE_AND) {
    const std::vector<PseudotreeNode*>& children = pt_node->getChildren();
    for (PseudotreeNode* ch : children) {
      int var_child = ch->getVar();
      PseudotreeNode* pt_child = m_pseudotree->getNode(var_child);
      std::string str = Context(NODE_OR, pt_child->getFullContext());
      BFSearchState state(NODE_OR, str);

#ifdef DEBUG
      cout << var_child << " OR context: " << str.c_str() << endl;
      cout << pt_node->getFullContext() << endl;
      cout.flush();
#endif

      BFSearchNodeOR* c;
      if (!m_options->nocaching && search_space_->find_node(var_child, state)) {
        c = (BFSearchNodeOR*) search_space_->node(var_child, state);
      } else {
        c = new BFSearchNodeOR(node, var_child, depth + 1);
        double h = assignCostsOR(c);
        c->setHeur(h);
        c->setValue(h);
        c->setFeasibleValue(ELEM_ZERO);

        search_space_->add_node(state, c);
      }
      c->add_parent(node);
      node->add_child(c);
      no_children = false;
    }

    node->set_expanded(true);
    node->set_fringe(false);
    node->set_terminal(children.empty());
    search_space_->IncNodesExpanded(NODE_AND);
  } else {
    assert(node->getHeurCache());
    for (int val = 0; val < m_problem->getDomainSize(var); ++val) {
      std::string str = Context(NODE_AND, pt_node->getFullContext());
      BFSearchState state(NODE_AND, str);
      BFSearchNodeAND* c = new BFSearchNodeAND(node, var, val, depth + 1);
      node->add_child(c);
      c->add_parent(node);
      double* heur_cache = node->getHeurCache();
      double* ordering_heur_cache = node->getOrderingHeurCache();
      double h = heur_cache[2 * val];
      double w = heur_cache[2 * val + 1];
      double oh = ordering_heur_cache[val];
      double h_child;
      if (h == ELEM_ZERO) {
        h_child = ELEM_ZERO;
      } else {
        h_child = h - w;
      }
      c->setHeur(h_child);
      c->setValue(h_child);
      c->setOrderingHeur(oh);
      search_space_->add_node(state, c);

      no_children = false;
    }
    node->set_expanded(true);
    node->set_fringe(false);
    node->set_terminal(false);
    search_space_->IncNodesExpanded(NODE_OR);
  }
  return no_children;
}

bool AOStar::Revise(BFSearchNode* node) {
  assert(node);
//  cout << "Revise BFS: (" << node << ") "<< node->ToString() << endl;


  bool change = true;
  bool tightening_threshold_reached = false;
  if (node->getType() == NODE_AND) {
    if (node->is_terminal()) {
      node->setValue(ELEM_ONE);
      node->set_solved(true);
      node->set_fringe(false);

      change = true;
    } else {
      double old_value = node->getValue();
      bool solved = true;
      double q_value = ELEM_ONE;
      for (BFSearchNode* child : node->children()) {
        solved = solved && child->is_solved();
        q_value OP_TIMESEQ child->getValue();
      }

      node->setValue(q_value);
      if (solved) {
        node->set_solved(true);
        node->set_fringe(false);
      }

      change = solved || (q_value != old_value);
    }
  } else if (node->getType() == NODE_OR) {
    double old_value = node->getValue();
    double q_value = ELEM_ZERO;
    BFSearchNode* best = nullptr;

    for (BFSearchNode* child : node->children()) {
      int val = child->getVal();
      double w = node->getWeight(val);
      double q = w OP_TIMES child->getValue();

      if (q - 1e-10 > q_value) {
        q_value = q;
        best = child;
      } else if (q == q_value) {
        if (!best || child->is_solved()) {
          best = child;
        }
      }
    }

    assert(best);

    node->setValue(q_value);
    node->set_best_child(best);
    bool solved = best->is_solved();
    if (solved) {
      node->set_solved(true);
      node->set_fringe(false);

      if (node != search_space_->getRoot()) {
        node->getParent()->addSubSolved(node->getValue());
      }
    }
    change = solved || (q_value != old_value);
    if (change && node == search_space_->getRoot()) {
      if (heuristic_bound_ - q_value > 1e-10) {
        heuristic_bound_ = q_value;
        if (m_options->algorithm != "aaobf") {
          high_resolution_clock::time_point now = high_resolution_clock::now();
          double t = duration_cast<duration<double>>(now - _time_start).count();
          if (prev_reported_time_ < 0 || t - prev_reported_time_ > 5) {
            m_problem->updateLowerUpperBound(solution_cost_, heuristic_bound_,
                &(search_space_->stats));
            prev_reported_time_ = t;
          }
        }
      }
    }
  }

  return change;
}

bool AOStar::FindBestPartialTree() {
  tip_nodes_.clear();
  BFSearchNode* root = reinterpret_cast<BFSearchNode*>(search_space_->getRoot());

  std::fill(m_assignment.begin(), m_assignment.end(), UNKNOWN);

  std::stack<BFSearchNode*> dfs_stack;
  dfs_stack.push(root);
  while (!dfs_stack.empty()) {
    BFSearchNode* e = dfs_stack.top();
    dfs_stack.pop();

    if (e->getType() == NODE_AND) {
      int var = e->getVar();
      int val = e->getVal();
      m_assignment[var] = val;

      if (!e->children().empty()) {
        for (BFSearchNode* child : e->children()) {
          child->set_current_parent(e);
          if (!child->is_solved()) {
            dfs_stack.push(child);
          }
        }
      } else {
        e->set_fringe(true);
        tip_nodes_.push_back(e);
      }
    } else {
      if (!e->children().empty()) {
        BFSearchNode* best_child = e->best_child();
        assert(best_child);
        best_child->set_current_parent(e);
        dfs_stack.push(best_child);
      } else {
        e->set_fringe(true);
        tip_nodes_.push_back(e);
      }
    }
  }
  return !tip_nodes_.empty();
}

void AOStar::ArrangeTipNodes() {
  std::sort(tip_nodes_.begin(), tip_nodes_.end(),
            comp_node_ordering_heur_desc_fn_);
}

BFSearchNode* AOStar::ChooseTipNode() {
  if (!tip_nodes_.empty()) {
    return reinterpret_cast<BFSearchNode*>(*tip_nodes_.begin());
  } else {
    return nullptr;
  }
}

std::string AOStar::Context(int node_type, const std::set<int>& ctxt) {
  std::ostringstream oss;
  if (ctxt.empty() || node_type == NODE_AND) {
    oss << "s" << global_search_index_++;
  } else {
    for (int cvar : ctxt) {
      assert(m_assignment[cvar] != UNKNOWN);
      oss << "x" << cvar << "=" << m_assignment[cvar] << ";";
    }
  }
  return oss.str();
}


void AOStar::reset(SearchNode* p) {
  search_space_->clear_all_nodes();
}

SearchNode* AOStar::initSearch() {
  reset(nullptr);
  int var_root = m_pseudotree->getRoot()->getVar();

  BFSearchNode* root = new BFSearchNodeOR(nullptr, var_root, 0);
  root->set_terminal(false);
  root->set_fringe(true);

  search_space_->setRoot(root);
  return root;
}

void AOStar::InitBFSearchSpace() {
  BFSearchNode* root = reinterpret_cast<BFSearchNode*>(
      search_space_->getRoot());
  BFSearchState state(NODE_OR, "s-2");
  search_space_->add_node(state, root);
  tip_nodes_.push_back(root);
}

bool AOStar::printStats() const {
  cout << "Search Stats: " << endl;
  cout << "============= " << endl;
  cout << "OR nodes:           " << m_space->stats.numExpOR << endl;
  cout << "AND nodes:          " << m_space->stats.numExpAND << endl;
  /*
  cout << "Deadend nodes:      " << m_space->stats.numDead << endl;
  cout << "Deadend nodes (CP): " << m_space->stats.numDeadCP << endl;
  */
  return true;
}

AOStar::AOStar(Problem* p, Pseudotree* pt, SearchSpace* space,
                     Heuristic* heur, BoundPropagator* prop,
                     ProgramOptions* po)
: Search(p, pt, space, heur, prop, po), global_search_index_(0),
  best_first_limit_reached_(false),
  timed_out_(false),
  comp_node_ordering_heur_desc_fn_(NodeOrderingHeurDesc()) {
  search_space_ = dynamic_cast<BFSearchSpace*>(space);
  initSearch();
}

}  // namespace daoopt
