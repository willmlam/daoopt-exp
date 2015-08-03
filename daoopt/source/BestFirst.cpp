#include "BestFirst.h"
#include "BFSearchSpace.h"

namespace daoopt {

bool BestFirst::solve(size_t nodeLimit) {
  // for when BestFirst search runs out of memory
  // allocate 64K memory.
  char* emergency_memory = new char[65536]; 

  // Initial best-first search phase
  try {
  } catch (std::bad_alloc& ba_exception) {
    delete [] emergency_memory;
    best_first_limit_reached_ = true;
  } catch (...) {
    std::cout << "Non bad alloc error?" << std::endl;
  }

  // If the problem is not already solved, use AOBB the rest of the way
  // (would have returned earlier)
  // We then need to put all the explicated nodes in depth_first_open_list_

}

void BestFirst::reset(SearchNode* p) {
  
}

}  // namespace daoopt
