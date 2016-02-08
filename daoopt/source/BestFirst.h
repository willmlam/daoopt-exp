#ifndef BESTFIRST_H_
#define BESTFIRST_H_

#include "BFSearchSpace.h"
#include "Search.h"

#include <queue>

namespace daoopt {
class NodeComp {
 public:
  bool operator() (const SearchNode* a, const SearchNode* b) const {
    return a->getHeur() < b->getHeur();
  }
};
class BestFirst : virtual public Search {
 public:
  BestFirst(Problem* p, Pseudotree* pt, SearchSpace* space,
      Heuristic* heur, BoundPropagator* prop, ProgramOptions* po);
  bool solve(size_t nodeLimit);

 protected:
  virtual SearchNode* initSearch();
  bool AOStar();
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
  
  struct CompNodeOrderingHeurDesc {
    CompNodeOrderingHeurDesc(Pseudotree* pt) : pt_(pt) { }
    bool operator()(const BFSearchNode* x, const BFSearchNode* y) const {
      double x_ordering_heur = pt_->getNode(x->getVar())
        ->getOrderingHeuristic();
      double y_ordering_heur = pt_->getNode(y->getVar())
        ->getOrderingHeuristic();
      return x_ordering_heur > y_ordering_heur;
    }
    Pseudotree* pt_;
  };

  std::function<bool(const BFSearchNode*,
                     const BFSearchNode*)> CompNodeOrderingHeurDescFn;

  struct CompNodeOrderingHeurDesc2 {
    CompNodeOrderingHeurDesc2(Pseudotree* pt) : pt_(pt) { }
    bool operator()(const BFSearchNode* x, const BFSearchNode* y) const {
      double x_ordering_heur = pt_->getNode(x->getVar())
        ->getOrderingHeuristic();
      double y_ordering_heur = pt_->getNode(y->getVar())
        ->getOrderingHeuristic();
      // use decreasing heuristic of original if tied
      if (x_ordering_heur == y_ordering_heur) {
        return x->getHeur() > y->getHeur();
      } else {
        return x_ordering_heur > y_ordering_heur;
      }
    }
    Pseudotree* pt_;
  };

  std::function<bool(const BFSearchNode*,
                     const BFSearchNode*)> CompNodeOrderingHeurDescFn2;

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
