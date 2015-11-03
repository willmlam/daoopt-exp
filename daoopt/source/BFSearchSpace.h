#ifndef BFSEARCHSPACE_H_
#define BFSEARCHSPACE_H_

// The best-first search space uses a different type of caching to maintain
// the AND/OR graph.

#include "SearchSpace.h"
#include "BFSearchNode.h"

#include <unordered_map>
#include <functional>

using std::hash;
using std::unordered_map;
using std::string;

namespace daoopt {

struct BFSearchState {
  int type_;
  std::string context_;
  bool operator==(const BFSearchState& state) const {
    return (type == state.type && context == a.context);
  }
  BFSearchState(int type, const string& context)
    : type_(type), context_(context) {}
};

struct StateHasher {
  size_t operator()(const BFSearchState& state) const {
    return hash<string>()(state.context_);
  }
};

using AOGraph = unordered_map<BFSearchState, BFSearchNode*, StateHasher>;

class BFSearchSpace : public SearchSpace {
 public:
  void add_node(const BFSearchState& state, BFSearchNode* node);
  bool find_node(const BFSearchState& state) const;
  BFSearchNode* get_node(const BFSearchState& state) const;
  void erase_node(const BFSearchState& state);
  void clear_all_nodes();

  const AOGraph get_nodes() const { return nodes_; }

  BFSearchSpace(Pseudotree* pt, ProgramOptions* opt);
  virtual ~BFSearchSpace();

 protected:
  // Cache for nodes in explicated search space
  AOGraph nodes_;

 private:
  BFSearchSpace(const BFSearchSpace&);
};

inline BFSearchSpace::BFSearchSpace(Pseudotree* pt, ProgramOptions* opt) 
  : SearchSpace(pt, opt) {
}

inline void BFSearchSpace::add_node(const BFSearchState& state,
    BFSearchNode* node) {
  assert(node);
  nodes_.insert(make_pair(state, node));
}

inline bool BFSearchSpace::find_node(const BFSearchState& state) const {
  return ContainsKey(nodes_, state);
}

inline BFSearchNode* BFSearchSpace::get_node(const BFSearchState& state) const {
  AOGraph::iterator it = nodes_.find(state);
  return it != nodes_.end() ? it->second : nullptr;
}

inline void BFSearchSpace::erase(const BFSearchState& state) {
  AOGraph::iterator it = nodes_.find(state);
  if (it != nodes_.end()) {
    BFSearchNode* node = it->second;
    nodes_.erase(it);
    delete node;
  }
}

inline void BFSearchSpace::incNodesExpanded(int nodeType) {
  if (nodeType == NODE_AND) {
    ++stats.numExpAND;
  } else {
    ++stats.numExpOR;
  }
}

}  // namespace daoopt
#endif
