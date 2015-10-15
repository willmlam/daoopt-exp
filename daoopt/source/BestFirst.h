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
    void reset(SearchNode* p);
    SearchNode* nextNode();

    bool doCompleteProcessing(SearchNode* n);
    bool doExpand(SearchNode* n);

    // Use a priority queue for the tip nodes.
    priority_queue<SearchNode*> tip_nodes_;

    bool best_first_limit_reached_;
      
  };

  /* Inline definitions */
  
  inline SearchNode* BestFirst::nextNode() {
  }
}  // namespace daoopt

#endif
