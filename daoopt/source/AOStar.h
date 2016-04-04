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
class AOStar : virtual public Search {
 public:
  AOStar(Problem* p, Pseudotree* pt, SearchSpace* space,
      Heuristic* heur, BoundPropagator* prop, ProgramOptions* po);
  bool solve(size_t nodeLimit);

 protected:
  virtual SearchNode* initSearch();
  bool DoSearch();
  void ExpandAndRevise(BFSearchNode* node);
  void Expand(BFSearchNode* node);
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
      return x->get_index() < y->get_index();
    }
  };

  // Ordering heuristic, tie breaks in favor of smaller depth.
  struct CompNodeOrderingHeurDesc : public std::binary_function<BFSearchNode*,
    BFSearchNode*, bool> {
    bool operator()(const BFSearchNode* x, const BFSearchNode* y) const {
      if (x->getOrderingHeur() == y->getOrderingHeur()) {
        return x->getHeur() < y->getHeur();
      } 
      return x->getOrderingHeur() > y->getOrderingHeur();
    }
  };
  
  void reset(SearchNode* p);
  SearchNode* nextNode();

  bool doCompleteProcessing(SearchNode* n);
  bool doExpand(SearchNode* n);
  bool isDone() const {
    return dynamic_cast<BFSearchNode*>(search_space_->getRoot())->is_solved();
  }
  bool isMaster() const { return false; }

  std::vector<SearchNode*> tip_nodes_;

  bool best_first_limit_reached_;

  BFSearchSpace* search_space_;
  size_t global_search_index_;
  double solution_cost_;
  double heuristic_bound_;
  double prev_reported_time_;
};




}  // namespace daoopt

#endif
