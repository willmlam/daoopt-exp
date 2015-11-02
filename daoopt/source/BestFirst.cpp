#include "BestFirst.h"
#include "BFSearchSpace.h"

namespace daoopt {

bool BestFirst::solve(size_t nodeLimit) {
  // for when BestFirst search runs out of memory
  // allocate 64K memory.
  char* emergency_memory = new char[65536]; 

  // Initial best-first search phase
  try {
    m_solved = AOStar();
  } catch (std::bad_alloc& ba_exception) {
    delete [] emergency_memory;
    best_first_limit_reached_ = true;
    m_solved = false;
  } catch (...) {
    std::cout << "Non bad alloc error?" << std::endl;
    m_solved = false;
  }

  // If the problem is not already solved, use AOBB the rest of the way
  // (would have returned earlier)

  if (emergency_memory) delete[] emergency_memory;

  return m_solved;

}

bool BestFirst::AOStar() {
  int var_root = m_pseudotree->getRoot()->getVar();
  double h = ELEM_ENCODE(m_heuristic->getGlobalUB());
  double w = 0.0;
  double* dv = new double[2];
  dv[0] = h;
  dv[1] = w;

  BFSearchNode* root = new BFSearchNodeOR(var_root, 0);
  root->setHeur(h);
  root->setValue(h);
  root->setHeurCache(dv);
  root->set_terminal(false);
  root->set_fringe(true);

  BFSearchState state(NODE_OR, "s-2");
  search_space_->setRoot(root);
  search_space_->add_node(state, root);
  tip_nodes_.push_back(root);
  while(!root->is_solved()) {
    assert(tip_nodes_.size() > 0);
    BFSearchNode* n = *tip_nodes_.begin();
    ExpandAndRevise(n);
  }

  if (root->is_solved()) {
    m_solutionCost = root->getValue();
  }
  return true;
}

void BestFirst::ExpandAndRevise(BFSearchNode* node) {
  assert(node->is_fringe() && !node->is_solved());
  Expand(node);

  // Revise
  std::multiset<BFSearchNode*> revise_set;

  s.insert(node);
  node->set_visited(true);

  while (!revise_set.empty()) {
    BFSearchNode *e = *revise_set.begin();
    revise_set.erase(revise_set.begin());
    e->set_visited(false);
    bool change = Revise(e);

    if (change) {
      if (e->getType() == NODE_AND) {
        for (BFSearchNode* parent : e->get_parents()) {
          size_t index = 0;
          BFSearchNode* best = parent->get_best_child();
          bool found = (best == e);

          for (BFSearchNode* child : p->get_children) {
            if (child->is_visited()) {
              ++index;
            }
          }

          if (parent->is_visited()) {
            revise_set.erase(parent);
            parent->decrement_index();
          } else if (found) {
            parent->set_index(index);
            parent->set_visited(true);
          }
          revise_set.insert(parent);
        }
      } else {
        for (BFSearchNode* parent : e->get_parents()) {
          size_t index = 0;
          for (BFSearchNode* child : parent->get_children()) {
            if (child->is_visited()) {
              ++index;
            }
          }
          parent->set_index(index);
          parent->set_visited(true);
          revise_set.insert(parent);
        }
      }
    } else {
      for (BFSearchNode* parent : e->get_parents()) {
        if (parent->is_visited()) {
          revise_set.erase(parent);
          parent->decrease_index();
          revise_set.insert(parent);
        }
      }
    }
  }
  assert(revise_set.empty());

  tip_nodes_.clear();
  FindBestPartialTree();

  assert(search_space_->getRoot()->isSolved() || tip_nodes_.size() > 0);
}

void BestFirst::Expand(BFSearchNode* node) {
  bool no_children = true;
  int var = node->getVar();
  int depth = node->getDepth();
  PseudotreeNode* pt_node = m_pseudotree->getNode(var);
  if (node->getType == NODE_AND) {
    const std::vector<PseudotreeNode*>& children = pt_node->getChildren();
    for (PseudotreeNode* ch : children) {
      int var_child = ch->getVar();
      BFSearchNodeOR* c = new BFSearchNodeOR(var_child, depth + 1);
      c->add_parent(node);
      node->add_child(c);

      int old_value = m_assignment[c->getVar()];
      h = assignCostsOR(c);
      m_assignment[c->getVar] = old_value;
      
      c->setHeur(h);
      c->setValue(h);

      std::string str = Context(var_child, NONE, pt_node->getFullContext());
      BFSearchState state(NODE_OR, str);
      search_space_->add_node(state, c);
      no_children = false;
    }

    node->set_expanded(true);
    node->set_fringe(false);
    node->set_terminal(children.empty());
    search_space_->incNodesExpanded(NODE_AND);
  } else {
    for (int val = 0; val < m_problem->getDomainSize(var); ++val) {
      std::string str = Context(var, val, pt_node->getAndContext());
      BFSearchState(NODE_AND, str);
      if (search_space_->find_node(state)) {
        BFSearchNodeAND* c = (BFSearchNodeAND*) search_space_->get(state);
        node->add_child(c);
        c->add_parent(node);
      } else {
        BFSearchNodeAND* c = new BFSearchNodeAND(var, val, depth + 1);
        node->add_child(c);
        c->add_parent(node);
        double* heur_cache = node->getHeurCache();
        double h = heur_cache[2 * val];
        double w = heur_cache[2 * val + 1];
        c->setHeur(h - w);
        c->setValue(h - w);
        search_space_->add_node(state, c);
      }
      no_children = false;
    }
    node->set_expanded(true);
    node->set_fringe(false);
    node->set_terminal(false);
    search_space_->incNodesExpanded(NODE_OR);
  }
  return no_children;
}

bool BestFirst::Revise(BFSearchNode* node) {
  bool change = true;
  if (node->getType() == NODE_AND) {
    if (node->is_terminal()) {
      node->setValue(0.0);
      node->set_solved(true);
      node->set_fringe(false);

      change = true;
    } else {
      double old_value = node->getValue();
      bool solved = true;
      double q_value = 0.0;
      for (BFSearchNode* child : node->get_children()) {
        solved &= child->is_solved();
        q_value += child->getValue();
      }

      node->setValue(q_value);
      if (solved) {
        node->set_solved(true);
        node->set_fringe(false);
      }
      change = solved || q_value != old_value;
    }
  } else {
    double old_value = node->getValue();
    double q_value = std::numeric_limits<double>::max();
    BFSearchNode* best = nullptr;
    
    for (BFSearchNode* child : node->get_children()) {
      int val = child->getVal();
      double w = node->getWeight(val);
      double q = w + child->getValue();

      if (q + 1e-15 < q_value) {
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
    s.pop();

    if (e->getType() == NODE_AND) {
      int var = e->getVar();
      int val = e->getVal();
      m_assignment[var] = val;

      if (!e->get_children().empty()) {
        for (BFSearchNode* child : e->get_children()) {
          child->set_current_parent(e);
          if (!c->is_solved()) {
            dfs_stack.push(c);
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

void BestFirst::Context(int var, int val, const std::set<int>& ctxt) {
  std::ostringstream oss;
  if (ctxt.empty() || val == NONE) {
    oss << "s" << m_globalSearchIndex++;
  } else {
    for (int cvar : ctxt) {
      if (cvar == var) {
        oss << "x" << var << "=" << val << ";";
      } else {
        assert(m_assignment[cvar] != UNKNOWN);
        oss << "x" << cvar << "=" << m_assignment[cvar] << ";";
      }
    }
  }
  return oss.str();
}


void BestFirst::reset(SearchNode* p) {
  
}

}  // namespace daoopt
