#include "BestFirst.h"
#include "BFSearchSpace.h"

#undef DEBUG

namespace daoopt {


bool BestFirst::solve(size_t nodeLimit) {
  // for when BestFirst search runs out of memory
  // allocate 64K memory.
  char* emergency_memory = new char[65536]; 

  bool solved = false;

  // Initial best-first search phase
  try {
    solved = AOStar();
  } catch (std::bad_alloc& ba_exception) {
    delete [] emergency_memory;
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

  return solved;
}

bool BestFirst::AOStar() {
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

  heuristic_bound_ = h;

  while(!root->is_solved()) {
    assert(tip_nodes_.size() > 0);
    ArrangeTipNodes();
    BFSearchNode* n = ChooseTipNode();
    ExpandAndRevise(n);
  }

  if (root->is_solved()) {
    solution_cost_ = root->getValue();
    m_problem->updateSolution(solution_cost_, &(search_space_->stats));
  }

  return true;
}

void BestFirst::ExpandAndRevise(BFSearchNode* node) {
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

    assert(e->get_index() == 0);
    bool change = Revise(e);

#ifdef DEBUG
    cout << "Revised " << e->ToString() << (change ? " (changed)" : "") << endl;
#endif
    
    if (change) {
      // Step 12 Nilssons
      if (e->getType() == NODE_AND) {
        assert(e->get_parents().size() == 1);
        BFSearchNode* parent = e->get_parents().front();
        
        // Count children of e still in the revise_set
        // The index will place this node later the more children there are.
        size_t index = 0;
        for (BFSearchNode* child : parent->get_children()) {
          if (child->is_visited()) {
            ++index;
          }
        }
        BFSearchNode* best = parent->get_best_child();
        bool found = (best == e); // Is 'e' the marked AND child of the parent?
        if (parent->is_visited()) {
          // decrease parent index and replace
          auto si = revise_set.begin();
          for (; si != revise_set.end(); ++si) {
            if (parent == *si) {
              if (parent->get_index() > 0) {
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
        for (BFSearchNode* parent : e->get_parents()) {
          size_t index = 0;
          for (BFSearchNode* child : parent->get_children()) {
            if (child->is_visited()) {
              ++index;
            }
          }

          if (parent->is_visited()) {
            auto si = revise_set.begin();
            for (; si != revise_set.end(); ++si) {
              if (parent == *si) {
                if (parent->get_index() > 0) {
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
      assert(e->getType() != NODE_AND || e->get_parents().size() == 1);
      for (BFSearchNode* parent : e->get_parents()) {
        if (parent->is_visited()) {
          auto si = revise_set.begin();
          for (; si != revise_set.end(); ++si) {
            if (parent == *si) {
              if (parent->get_index() > 0) {
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

  tip_nodes_.clear();
  FindBestPartialTree();

  assert(dynamic_cast<BFSearchNode*>(search_space_->getRoot())->is_solved() ||
         tip_nodes_.size() > 0);
}

void BestFirst::Expand(BFSearchNode* node) {
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
      if (search_space_->find_node(var_child, state)) {
        c = (BFSearchNodeOR*) search_space_->get_node(var_child, state);
      } else {
        c = new BFSearchNodeOR(node, var_child, depth + 1);
        double h = assignCostsOR(c);
        c->setHeur(h);
        c->setValue(h);

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
      double h = heur_cache[2 * val];
      double w = heur_cache[2 * val + 1];
      c->setHeur(h - w);
      c->setValue(h - w);
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

bool BestFirst::Revise(BFSearchNode* node) {
  assert(node);

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
      for (BFSearchNode* child : node->get_children()) {
        solved = solved && child->is_solved();
        q_value OP_TIMESEQ child->getValue();
      }

      node->setValue(q_value);
      if (solved) {
        node->set_solved(true);
        node->set_fringe(false);
      }

      change = solved || q_value != old_value;
    }
  } else if (node->getType() == NODE_OR) {
    double old_value = node->getValue();
    double q_value = -std::numeric_limits<double>::infinity();
    BFSearchNode* best = nullptr;
    
    for (BFSearchNode* child : node->get_children()) {
      int val = child->getVal();
      double w = node->getWeight(val);
      double q = w OP_TIMES child->getValue();

      if (q > q_value) {
        q_value = q;
        best = child;
      } else if (q == q_value) {
        if (!best) {
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
    }
    change = solved || q_value != old_value;
    if (change && node == search_space_->getRoot()) {
      if (heuristic_bound_ - q_value > 1e-5) {
        heuristic_bound_ = q_value;
        m_problem->updateUpperBound(heuristic_bound_, &(search_space_->stats));
      }
    }
  }

  return change;
}

bool BestFirst::FindBestPartialTree() {
  BFSearchNode* root = search_space_->getRoot();

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

      if (!e->get_children().empty()) {
        for (BFSearchNode* child : e->get_children()) {
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
      if (!e->get_children().empty()) {
        BFSearchNode* best_child = e->get_best_child();
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

void BestFirst::ArrangeTipNodes() {
  std::sort(tip_nodes_.begin(), tip_nodes_.end(), CompNodeOrderingHeurDescFn);
}

BFSearchNode* BestFirst::ChooseTipNode() {
  if (!tip_nodes_.empty()) {
    return *tip_nodes_.begin();
  } else {
    return nullptr;
  }
}

std::string BestFirst::Context(int node_type, const std::set<int>& ctxt) {
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


void BestFirst::reset(SearchNode* p) {
  search_space_->clear_all_nodes();
}

SearchNode* BestFirst::initSearch() {
  reset(nullptr);
  int var_root = m_pseudotree->getRoot()->getVar();

  BFSearchNode* root = new BFSearchNodeOR(nullptr, var_root, 0);
  root->set_terminal(false);
  root->set_fringe(true);

  BFSearchState state(NODE_OR, "s-2");
  search_space_->setRoot(root);
  search_space_->add_node(state, root);
  tip_nodes_.push_back(root);
  return root;
}


// Empty implementations for unused functions.
bool BestFirst::doCompleteProcessing(SearchNode* n) {
  assert(false);
  return false;
}

bool BestFirst::doExpand(SearchNode* n) {
  assert(false);
  return false;
}

SearchNode* BestFirst::nextNode() {
  assert(false);
  return nullptr;
}

BestFirst::BestFirst(Problem* p, Pseudotree* pt, SearchSpace* space,
                     Heuristic* heur, BoundPropagator* prop,
                     ProgramOptions* po)
: Search(p, pt, space, heur, prop, po), global_search_index_(0) {
  search_space_ = dynamic_cast<BFSearchSpace*>(space);
  this->initSearch();
  CompNodeOrderingHeurDescFn = CompNodeOrderingHeurDesc(pt);
}

}  // namespace daoopt
