#ifndef BESTFIRST_H_
#define BESTFIRST_H_

#include "Search.h"
#include "BFSearchSpace.h"

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

   protected:
    void reset(SearchNode* p);
    SearchNode* nextNode();

    bool doExpand(SearchNode* n);

    // Stack of nodes representing the open list when search is in 
    // depth-first mode
    stack<SearchNode*> depth_first_open_list_;

    bool best_first_limit_reached_;
      
  };

  /* Inline definitions */
  
  inline SearchNode* BestFirst::nextNode() {
  }
}  // namespace daoopt

#endif
