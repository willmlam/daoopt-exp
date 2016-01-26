#ifndef BFSEARCHSPACE_H_
#define BFSEARCHSPACE_H_

// The best-first search space uses a different type of caching to maintain
// the AND/OR graph.

#include "SearchSpace.h"
#include "BFSearchNode.h"
#include "hash_murmur.h"

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
    return (type_ == state.type_ && context_ == state.context_);
  }
  BFSearchState(int type, const string& context)
    : type_(type), context_(context) {}
};

struct StateHasher {
  size_t operator()(const BFSearchState& state) const {
    size_t hash = 0;
    size_t seed = 0x9e3779b9;
    int len = (int)state.context_.size();
    register unsigned char* key = (unsigned char*)state.context_.c_str();
    MurmurHash3_x64_64(key, len, seed, &hash);
    return hash;
  }
};

using AOGraph = boost::unordered_map<BFSearchState, BFSearchNode*, StateHasher>;

class BFSearchSpace : virtual public SearchSpace {
 public:
  void add_node(const BFSearchState& state, BFSearchNode* node);
  bool find_node(int var, const BFSearchState& state) const;
  BFSearchNode* get_node(int var, const BFSearchState& state) const;
  void erase_node(int var, const BFSearchState& state);
  void clear_all_nodes();

  const std::vector<AOGraph>& get_nodes() const { return nodes_; }
  std::vector<AOGraph>& get_nodes() { return nodes_; }

  void IncNodesExpanded(int node_type);

  BFSearchSpace(Pseudotree* pt, ProgramOptions* opt, int size);
  virtual ~BFSearchSpace();

 protected:
  // Cache for nodes in explicated search space
  std::vector<AOGraph> nodes_;

 private:
  BFSearchSpace(const BFSearchSpace&);
};

inline BFSearchSpace::BFSearchSpace(Pseudotree* pt, ProgramOptions* opt,
                                    int size) 
  : SearchSpace(pt, opt) {
  nodes_.resize(size + 1);
}

inline BFSearchSpace::~BFSearchSpace() {
  clear_all_nodes();
  root = nullptr;
}

inline void BFSearchSpace::add_node(const BFSearchState& state,
    BFSearchNode* node) {
  assert(node);
  nodes_[node->getVar()].insert(std::make_pair(state, node));
}

inline bool BFSearchSpace::find_node(int var,
                                     const BFSearchState& state) const {
  return ContainsKey(nodes_[var], state);
}

inline BFSearchNode* BFSearchSpace::get_node(int var, const BFSearchState& state) const {
  AOGraph::const_iterator it = nodes_[var].find(state);
  return it != nodes_[var].end() ? it->second : nullptr;
}

inline void BFSearchSpace::erase_node(int var, const BFSearchState& state) {
  AOGraph::iterator it = nodes_[var].find(state);
  if (it != nodes_[var].end()) {
    BFSearchNode* node = it->second;
    nodes_[var].erase(it);
    delete node;
  }
}

inline void BFSearchSpace::clear_all_nodes() {
  for (auto& var_nodes : nodes_) {
    for (auto& kv : var_nodes) {
      if (kv.second) {
        delete kv.second;
      }
    }
  }
  nodes_.clear();
}

inline void BFSearchSpace::IncNodesExpanded(int node_type) {
  if (node_type == NODE_AND) {
    ++stats.numExpAND;
  } else if (node_type == NODE_OR) {
    ++stats.numExpOR;
  }
}

}  // namespace daoopt
#endif
