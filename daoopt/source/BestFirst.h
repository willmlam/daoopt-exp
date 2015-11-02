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
class BestFirst : public Search {
 public:
  BestFirst(Problem* p, Pseudotree* pt, SearchSpace* space,
      Heuristic* heur, ProgramOptions* po);
  bool solve(size_t nodeLimit);

 protected:
  bool AOStar();
  void ExpandAndRevise(BFSearchNode* node);
  void Expand(BFSearchNode* node);
  bool Revise(BFSearchNode* node);
  bool FindBestPartialTree();

  std::string Context(int var, int val, const std::set<int>& ctxt);
  void reset(SearchNode* p);
  SearchNode* nextNode();

  bool doCompleteProcessing(SearchNode* n);
  bool doExpand(SearchNode* n);

  std::vector<SearchNode*> tip_nodes_;

  bool best_first_limit_reached_;
    
};

/* Inline definitions */

}  // namespace daoopt

#endif
