#include "AnytimeAOStar.h"
#include "BFSearchSpace.h"

#include <boost/range/adaptor/reversed.hpp>
#include <chrono>
using namespace std::chrono;

namespace daoopt {

extern high_resolution_clock::time_point _time_start; // From Main.cpp


bool AnytimeAOStar::solve(size_t nodeLimit) {
  char* emergency_memory = new char[65536];

  bool solved = false;

  try {
    InitBFSearchSpace();
    solved = DoSearch();
  } catch (std::bad_alloc& ba_exception) {
    delete[] emergency_memory;
    emergency_memory = nullptr;
    best_first_limit_reached_ = true;
    solved = false;
  } catch (...) {
    std::cerr << "Non bad alloc error?" << std::endl;
    solved = false;
  }

  if (emergency_memory) delete[] emergency_memory;

  return solved;
}

bool AnytimeAOStar::DoSearch() {
  BFSearchNode* root = dynamic_cast<BFSearchNode*>(search_space_->getRoot());
  double h = m_heuristic->getGlobalUB();
  double w = ELEM_ONE;
  double* dv = new double[2];
  dv[0] = h;
  dv[1] = w;

  root->setHeur(h);
  root->setValue(h);
  root->setFeasibleValue(ELEM_ZERO);
  root->set_terminal(false);
  root->set_fringe(true);
  root->setHeurCache(dv);

  solution_cost_ = ELEM_ZERO;
  heuristic_bound_ = h;
  m_problem->updateLowerUpperBound(solution_cost_, heuristic_bound_, &search_space_->stats);

  prev_reported_time_ = -1;

  assignment_feasible_ = m_assignment;

  bool repair_needed = false;
  while (!root->is_solved()) {
    BFSearchNode* next = nullptr;
    if (!tip_nodes_feasible_.empty()) {
      m_assignment = assignment_feasible_;
      ArrangeTipNodesFeasible();
      next = ChooseTipNodeFeasible();
      ExpandAndRevise(next, BoundType::kFeasible);
      IncExpNodeCount(BoundType::kFeasible, next->getType());
      FindFeasiblePartialTree();
      repair_needed = true;
    } else {
      if (root->getFeasibleValue() > solution_cost_) {
        solution_cost_ = root->getFeasibleValue();
        m_problem->updateLowerUpperBound(solution_cost_,
                                         heuristic_bound_,
                                         &(search_space_->stats));
      }

      if (repair_needed) {
        Repair();
        repair_needed = false;
        FindBestPartialTree();
      }

      // no tip nodes in BST means the problem is solved
      if (tip_nodes_.empty()) {
        break;
      }

      ArrangeTipNodes();
      next = ChooseTipNode();
      ExpandAndRevise(next, BoundType::kHeuristic);
      IncExpNodeCount(BoundType::kHeuristic, next->getType());
      FindFeasiblePartialTree();
    }
  }

//  Repair();
//  solution_cost_ = root->getFeasibleValue();
  solution_cost_ = heuristic_bound_;
  m_problem->updateLowerUpperBound(solution_cost_,
                                   heuristic_bound_,
                                   &(search_space_->stats));
  cout << "FST exp (OR): " << exp_depth_first_or_ << endl;;
  cout << "FST exp (AND): " << exp_depth_first_and_ << endl;;
  cout << "BST exp (OR): " << exp_best_first_or_ << endl;;
  cout << "BST exp (AND): " << exp_best_first_and_ << endl;;
  return true;
}

void AnytimeAOStar::InitBFSearchSpace() {
  // Setup layers for repariring step
  // number of layers is 2*(h+1): height + dummy and AND/OR layers
  int num_layers = 2 * (m_pseudotree->getHeight() + 2);
  dynamic_cast<BFSearchSpace*>(m_space)->initialize_layers(num_layers);

  BFSearchNode* root = reinterpret_cast<BFSearchNode*>(
      search_space_->getRoot());
  BFSearchState state(NODE_OR, "s-2");
  root->setFeasibleValue(ELEM_ZERO);
  search_space_->add_node(state, root);
  tip_nodes_feasible_.push_back(root);
}

bool AnytimeAOStar::MarkFeasibleChild(BFSearchNode* node) {
  // OR nodes only.
  assert(node->getType() == NODE_OR);

  BFSearchNode* m = nullptr;
  // If there is already a feasible value, then we need to update it to the
  // best one, which is based on the lower-upper bound gap.
  double feasible_value = node->getFeasibleValue();
  if (feasible_value != ELEM_ZERO) {
    double max_gap = ELEM_ZERO;
    for (BFSearchNode* c : node->children()) {
      int val = c->getVal();
      double w = node->getWeight(val);
      double q = c->getValue();

      // pruning check
      double heur_bound = q OP_TIMES w;
      if (heur_bound < feasible_value) {
#ifdef DEBUG
        cout << "heur bound (" << heur_bound<< ") <= feasible value (" <<
            feasible_value << ")" << endl;
        cout << "val: " << val << endl;
        cout << node->ToString() << endl;
        cout << c->ToString() << endl;
#endif
        continue; // the current node already has a better feasible solution
      }

      // At this point, we are checking each child for whichever one has a
      // heuristic estimate that is much better.
      // the gap is the ratio of the estimate over the known feasible value.
      double child_feasible = feasible_value OP_DIVIDE w;
      double gap = q / child_feasible;
      /*
      if (!c->is_terminal()) {
        gap = q OP_DIVIDE gap;
      } 
      */

      if (gap + 1e-10 > max_gap) {
        max_gap = gap;
        m = c;
      }
    }
#ifdef DEBUG
    cout << "Marked: " << node << " -> " << m << endl;
#endif
    node->set_feasible_child(m);
  } else {
    double max_q = ELEM_ZERO;
    if (!node->children().empty()) {
      for (BFSearchNode* c : node->children()) {
        int val = c->getVal();
        double w = node->getWeight(val);
        double q = c->getValue();

        if (w OP_TIMES q + 1e-10 > max_q) {
          max_q = w OP_TIMES q;
          m = c;
        }
      }
#ifdef DEBUG
      cout << "Marked: " << node << " -> " << m << endl;
#endif
      node->set_feasible_child(m);
    }
  }
  return true;
}

void AnytimeAOStar::ExpandAndRevise(BFSearchNode* node, BoundType bound_type) {
  assert(node->is_fringe());
  assert(node->children().empty());

#ifdef DEBUG
  cout << "ExpandAndRevise: " << bound_type << " ";
  cout << "node: " << node << ", " << node->ToString() << endl;
#endif

  node->set_fringe(false);
  Expand(node);
  node->set_expanded(true);
  if (node->getType() == NODE_OR) {
    MarkFeasibleChild(node);
  }

  switch (bound_type) {
    case BoundType::kHeuristic:
      UpdateBound(node, BoundType::kHeuristic);
      UpdateBound(node, BoundType::kFeasible);
      heuristic_bound_ = search_space_->getRoot()->getValue();
      break;
    case BoundType::kFeasible:
      UpdateBound(node, BoundType::kFeasible);
      break;
  }
}

void AnytimeAOStar::UpdateBound(BFSearchNode* node, BoundType bound_type) {
#ifdef DEBUG
  cout << "UpdateBound ";
  switch (bound_type) {
    case BoundType::kHeuristic:
      cout << "(kHeuristic): ";
      break;
    case BoundType::kFeasible:
      cout << "(kFeaible): ";
      break;
  }
  cout << node->ToString() << endl;
#endif
  std::multiset<BFSearchNode*, CompNodeIndexAsc> revise_set;

  revise_set.insert(node);
  node->set_visited(true);

  while (!revise_set.empty()) {
    BFSearchNode* e = *revise_set.begin();
    revise_set.erase(revise_set.begin());
    e->set_visited(false);
    
    assert(e->index() == 0);
    bool change = false;
  
    switch (bound_type) {
      case BoundType::kHeuristic: {
        change = Revise(e);
        break;
      }
      case BoundType::kFeasible: {
        change = ReviseFeasible(e);
        break;
      }
    }

    if (change) {
      if (e->getType() == NODE_AND) {
        assert(e->parents().size() == 1);
        BFSearchNode* parent = e->parents().front();

        size_t index = 0;
        for (BFSearchNode* child : parent->children()) {
          if (child->is_visited()) {
            ++index;
          }
        }
        bool found = false;
        if (bound_type == BoundType::kHeuristic) {
          found = (parent->best_child() == e);
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
        } else if (found) {
          assert(bound_type == BoundType::kHeuristic);
          parent->set_index(index);
          parent->set_visited(true);
          revise_set.insert(parent);
        }
      } else if (e->getType() == NODE_OR) {
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

  if (bound_type == BoundType::kHeuristic) {
    FindBestPartialTree();

    assert(dynamic_cast<BFSearchNode*>(
            search_space_->getRoot())->is_solved() ||
            tip_nodes_.size() > 0);
  }
}

bool AnytimeAOStar::ReviseFeasible(BFSearchNode* node) {
  bool result = true;
//  cout << "Revise DFS: (" << node << ") "<< node->ToString() << endl;
  switch (node->getType()) {
    case NODE_AND: {
      if (node->is_terminal()) {
        node->setFeasibleValue(ELEM_ONE);
        node->set_feasible(true);
        node->set_solved(true);
        node->set_fringe(false);
        result = true;
      } else {
        double old_feasible_value = node->getFeasibleValue();
        double feasible_value = ELEM_ONE;
        bool feasible = true;
        for (BFSearchNode* c : node->children()) {
          assert(c->getType() == NODE_OR);
          feasible = feasible && (c->is_feasible() || c->is_solved());

          if(c->is_feasible() || c->is_solved()) {
            feasible_value OP_TIMESEQ c->getFeasibleValue();
          }
        }

        if (feasible) {
          node->setFeasibleValue(feasible_value);
          node->set_feasible(true);
          node->set_fringe(false);
        }
        result = feasible || (feasible_value != old_feasible_value);
      }
      break;
    }
    case NODE_OR: {
      double old_feasible_value = node->getFeasibleValue();
      double feasible_value = ELEM_ZERO;
      bool feasible = false;
      for (BFSearchNode* c : node->children()) {
        assert(c->getType() == NODE_AND);
        feasible = feasible || c->is_feasible();
        int val = c->getVal();
        double w = node->getWeight(val);
        double f = c->is_feasible() ?
          w OP_TIMES c->getFeasibleValue() :
          ELEM_ZERO;
        feasible_value = std::max(feasible_value, f);
      }

      node->setFeasibleValue(feasible_value);
      if (feasible) {
        node->set_feasible(true);
        node->set_fringe(false);
      }

      result = feasible || (feasible_value != old_feasible_value);
      break;
    }
  }
  return result;
}


bool AnytimeAOStar::Repair() {
  AOLayers& layers = search_space_->layers();
  for (auto layer_list : boost::adaptors::reverse(layers)) {
    for (BFSearchNode* node : layer_list) {
      int var = node->getVar();
      node->set_solved(false);
      if (!node->is_expanded()) {
        double heur = node->getHeur();
        node->setValue(heur);
        if (!node->is_terminal()) {
          node->set_solved(false);
        }
        node->set_fringe(true);
      } else {
        Revise(node);
        ReviseFeasible(node);
      }

      if (node->getType() == NODE_OR && !node->is_solved()) {
        MarkFeasibleChild(node);
      }
    }
  }
  return true;
}

bool AnytimeAOStar::FindFeasiblePartialTree() {
  tip_nodes_feasible_.clear();
  size_t feasible_tree_node_count = 0;

  BFSearchNode* root = dynamic_cast<BFSearchNode*>(search_space_->getRoot());
  if (!root->is_solved()) {
    double g = ELEM_ONE;
    double h = ELEM_ONE;

    std::fill(assignment_feasible_.begin(),
        assignment_feasible_.end(), UNKNOWN);

    std::stack<BFSearchNode*> dfs;
    dfs.push(root);

    while(!dfs.empty()) {
      BFSearchNode* node = dfs.top();
      dfs.pop();
      ++feasible_tree_node_count;
      switch(node->getType()) {
        case NODE_AND: {
          int var = node->getVar();
          int val = node->getVal();
          assignment_feasible_[var] = val;

          if (!node->children().empty()) {
            for (BFSearchNode* c : node->children()) {
              dfs.push(c);
            }
          } else {
            if (!node->is_solved()) {
              node->set_fringe(true);
              tip_nodes_feasible_.push_back(node);
              h OP_TIMESEQ node->getHeur();
            } else {
              g OP_TIMESEQ node->getFeasibleValue();
            }
          }
          break;
        }
        case NODE_OR: {
          if (!node->children().empty()) {
            BFSearchNode* m = node->feasible_child();
            if (!m) {
              tip_nodes_feasible_.clear();
//              cout << node->ToString() << endl;
              return false;
            } else {
              dfs.push(m);
              g OP_TIMESEQ node->getWeight(m->getVal());
            }
          } else {
            if (!node->is_solved()) {
              node->set_fringe(true);
              tip_nodes_feasible_.push_back(node);
              h OP_TIMESEQ node->getHeur();
            } else {
              g OP_TIMESEQ node->getFeasibleValue();
            }
          }
          break;
        }
      }
    }

    if (g OP_TIMES h <= solution_cost_) {
      tip_nodes_feasible_.clear();
      return false;
    }
#ifdef DEBUG
    cout << "Feasible tree size: " << feasible_tree_node_count << endl;
    cout << "assignment: " << assignment_feasible_ << endl;
#endif
    return true;
  }
  return false; // already solved
}

void AnytimeAOStar::ArrangeTipNodesFeasible() {
  std::sort(tip_nodes_feasible_.begin(),
      tip_nodes_feasible_.end(), comp_node_ordering_gap_fn_);
}

BFSearchNode* AnytimeAOStar::ChooseTipNodeFeasible() {
  if (!tip_nodes_feasible_.empty()) {
    return reinterpret_cast<BFSearchNode*>(*tip_nodes_feasible_.begin());
  } else {
    return nullptr;
  }
}

AnytimeAOStar::AnytimeAOStar(Problem* p, Pseudotree* pt, SearchSpace* space,
    Heuristic* heur, BoundPropagator* prop, ProgramOptions* po)
: Search(p, pt, space, heur, prop, po),
  AOStar(p, pt, space, heur, prop, po),
  exp_depth_first_or_(0), exp_depth_first_and_(0),
  exp_best_first_or_(0), exp_best_first_and_(0),
  comp_node_ordering_gap_fn_(NodeHeurDesc()) {
}

}  // namespace daoopt
