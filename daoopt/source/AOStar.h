#ifndef AOSTAR_H_
#define AOSTAR_H_

#include "BFSearchSpace.h"
#include "Search.h"

namespace daoopt {
class NodeComp {
 public:
  bool operator() (const SearchNode* a, const SearchNode* b) const {
    return a->getHeur() < b->getHeur();
  }
};

class NodeOrderingGap {
 public:
  bool operator()(const SearchNode* x, const SearchNode* y) const {
    double x_gap = x->getValue() OP_DIVIDE x->getFeasibleValue();
    double y_gap = y->getValue() OP_DIVIDE y->getFeasibleValue();
    if (std::isinf(x_gap) && std::isinf(y_gap)) {
      return x->getHeur() > y->getHeur();
    }
    return x_gap > y_gap;
  }
};

class NodeOrderingHeurDesc {
 public:
  bool operator()(const SearchNode* x, const SearchNode* y) const {
    if (x->getOrderingHeur() == y->getOrderingHeur()) {
      return x->getHeur() > y->getHeur();
    }
    return x->getOrderingHeur() > y->getOrderingHeur();
  }
};

class AOStar : virtual public Search {
 public:
  AOStar(Problem* p, Pseudotree* pt, SearchSpace* space,
      Heuristic* heur, BoundPropagator* prop, ProgramOptions* po);
  bool solve(size_t nodeLimit);
  bool printStats() const override;

 protected:
  virtual SearchNode* initSearch();
  virtual void InitBFSearchSpace();
  virtual bool DoSearch();
  void ExpandAndRevise(BFSearchNode* node);
  bool Expand(BFSearchNode* node);
  bool Revise(BFSearchNode* node);
  bool FindBestPartialTree();
  void ArrangeTipNodes();

  BFSearchNode* ChooseTipNode();


  std::string Context(int node_type, const std::set<int>& ctxt);

  struct CompNodeHeurAsc : public std::binary_function<BFSearchNode*,
    BFSearchNode*, bool> {
    bool operator()(const BFSearchNode* x, const BFSearchNode* y) const {
      return x->getHeur() < y->getHeur();
    }
  };

  struct CompNodeHeurDesc : public std::binary_function<BFSearchNode*,
    BFSearchNode*, bool> {
    bool operator()(const BFSearchNode* x, const BFSearchNode* y) const {
      return x->getHeur() > y->getHeur();
    }
  };

  struct CompNodeIndexAsc : public std::binary_function<BFSearchNode*,
    BFSearchNode*, bool> {
    bool operator()(const BFSearchNode* x, const BFSearchNode* y) const {
      return x->index() < y->index();
    }
  };

  // Ordering heuristic, tie breaks in favor of smaller depth.
  struct CompNodeOrderingHeurDesc : public std::binary_function<BFSearchNode*,
    BFSearchNode*, bool> {
    bool operator()(const BFSearchNode* x, const BFSearchNode* y) const {
      if (x->getOrderingHeur() == y->getOrderingHeur()) {
        return x->getHeur() > y->getHeur();
      }
      return x->getOrderingHeur() > y->getOrderingHeur();
    }
  };

  void reset(SearchNode* p);


  bool isDone() const {
    return dynamic_cast<BFSearchNode*>(search_space_->getRoot())->is_solved();
  }
  bool isMaster() const { return false; }

  // Empty implementation for unused functions.
  bool doCompleteProcessing(SearchNode* n) {
    assert(false);
    return false;
  }
  bool doExpand(SearchNode* n) {
    assert(false);
    return false;
  }
  SearchNode* nextNode() {
    assert(false);
    return false;
  }

  std::vector<SearchNode*> tip_nodes_;

  bool best_first_limit_reached_;
  bool timed_out_;

  BFSearchSpace* search_space_;
  size_t global_search_index_;
  double solution_cost_;
  double heuristic_bound_;
  double prev_reported_time_;

  // function objects
  std::function<bool(const SearchNode*, const SearchNode*)>
      comp_node_ordering_heur_desc_fn_;
};




}  // namespace daoopt

#endif
