#ifndef ANYTIMEAOSTAR_H_
#define ANYTIMEAOSTAR_H_

#include "AOStar.h"
#include "BFSearchSpace.h"

namespace daoopt {
enum BoundType {
  kHeuristic, kFeasible
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

class NodeHeurDesc {
 public:
  bool operator()(const SearchNode* x, const SearchNode* y) const {
    return x->getHeur() > y->getHeur();
  }
};


class AnytimeAOStar : virtual public AOStar {
 public:
  AnytimeAOStar(Problem* p, Pseudotree* pt, SearchSpace* space,
                Heuristic* heur, BoundPropagator* prop, ProgramOptions* po);
//  bool solve(size_t nodeLimit);
  bool printStats() const override;
 protected:
  void InitBFSearchSpace() override;
  bool DoSearch() override;
  void ExpandAndRevise(BFSearchNode* node, BoundType bound_type);

  // for feasible tree
  void ExpandAndReviseFeasible(BFSearchNode* node);
  void ArrangeTipNodesFeasible();
  BFSearchNode* ChooseTipNodeFeasible();

  // upper-lower bound updates per node (assumes maximization)
  void UpdateBound(BFSearchNode* node, BoundType bound_type);

  // propagates upper and lower bound values after finding a solution
  bool Repair();

  bool ReviseFeasible(BFSearchNode* node);

  bool FindFeasiblePartialTree();

  // marks the best feasible child of the given node if it exists.
  bool MarkFeasibleChild(BFSearchNode* node);

  void IncExpNodeCount(BoundType bound_type, int node_type);


 protected: 
  std::vector<SearchNode*> tip_nodes_feasible_;
  std::vector<val_t> assignment_feasible_;
  size_t exp_depth_first_or_;
  size_t exp_best_first_or_;
  size_t exp_depth_first_and_;
  size_t exp_best_first_and_;

  std::function<bool(const SearchNode*, const SearchNode*)>
    comp_node_ordering_gap_fn_;
  
};

inline void AnytimeAOStar::IncExpNodeCount(BoundType bound_type,
                                           int node_type) {
  switch(bound_type) {
    case BoundType::kHeuristic: {
      switch(node_type) {
        case NODE_OR: {
          ++exp_best_first_or_;
          break;
        }
        case NODE_AND: {
          ++exp_best_first_and_;
          break;
        }
      }
      break;
    }
    case BoundType::kFeasible: {
      switch(node_type) {
        case NODE_OR: {
          ++exp_depth_first_or_;
          break;
        }
        case NODE_AND: {
          ++exp_depth_first_and_;
          break;
        }
      }
      break;
    }
  }
}

}  // namespace daoopt

#endif
