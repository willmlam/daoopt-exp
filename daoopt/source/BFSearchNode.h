#ifndef BFSEARCHNODE_H_
#define BFSEARCHNODE_H_

#include "SearchNode.h"

namespace daoopt {

class BFSearchNode;

class BFSearchNode : public SearchNode {
 public:
  inline double get_heur_updated() const {
    return heur_updated_;
  }
  inline void set_heur_updated(double f) {
    heur_updated_ = f;
  }

  inline BFSearchNode* get_current_parent() const {
    return current_parent_;
  }
  inline void set_current_parent(BFSearchNode* node) {
    current_parent_ = node;
  }

  inline int get_index() const {
    return index_;
  }
  inline void set_index(int index) {
    _index = index;
  }
  inline void increment_index() {
    ++_index;
  }
  inline void decrement_index() {
    --_index;
  }

  inline std::vector<BFSearchNode*>& get_parents() const {
    return parents_;
  }
  inline void add_parent(BFSearchNode* node) {
    parents_.push_back(node);
  }
  inline void erase_parent(BFSearchNode* node) {
    parents_.remove(node);  
  }
  inline bool has_parent(BFSearchNode* node) {
    for (auto parent : parents_) {
      if (parent == node) {
        return true;
      }
    }
    return false;
  }

  inline std::vector<BFSearchNode*>& get_children() const {
    return children_;
  }
  inline void add_child(BFSearchNode* node) {
    children_.push_back(node);
  }
  inline void erase_child(BFSearchNode* node) {
    children_.remove(node);
  }
  inline bool has_child(BFSearchNode* node) {
    for (auto child : children_) {
      if (child == node) {
        return true;
      }
    }
    return false;
  }

  virtual void set_best_child(BFSearchNode* node) = 0;
  virtual BFSearchNode* get_best_child() = 0;

  inline size_t get_hash_key() const {
    return hash_key_;
  }
  inline void set_hash_key(size_t hash_key) {
    hash_key_ = hash_key;
  }

  inline bool is_fringe() const {
    return is_fringe_;
  }
  inline void set_fringe(bool flag) {
    is_fringe_ = flag;
  }

  inline bool is_fringe() const {
    return is_fringe_;
  }
  inline void set_fringe(bool flag) {
    is_fringe_ = flag;
  }

  inline bool is_solved() const {
    return is_solved_;
  }
  inline void set_solved(bool flag) {
    is_solved_ = flag;
  }

  inline bool is_lookahead() const {
    return is_lookahead_;
  }
  inline void set_lookahead(bool flag) {
    is_lookahead_ = flag;
  }

  inline bool is_terminal() const {
    return is_terminal_;
  }
  inline void set_terminal(bool flag) {
    is_terminal_ = flag;
  }

  inline bool is_expanded() const {
    return is_expanded_;
  }
  inline void set_expanded(bool flag) {
    is_expanded_ = flag;
  }

  inline bool is_visited() const {
    return is_visited_;
  }
  inline void set_visited(bool flag) {
    is_visited_ = flag;
  }

  inline bool is_deadend() const {
    return is_deadend_;
  }
  inline void set_deadend(bool flag) {
    is_deadend_ = flag;
  }
 protected:
  // Various flags to indicate the status of the node.
  bool is_fringe_;
  bool is_solved_;
  bool is_lookahead_;
  bool is_terminal_;
  bool is_expanded_;
  bool is_visited_;
  bool is_deadend_;

  // Tightened estimate via lookahead
  double heur_updated_;

  std::list<BFSearchNode*> parents_;
  std::list<BFSearchNode*> children_;

  int index_;
  BFSearchNode* current_parent_;

  size_t hash_key_;
};

class BFSearchNodeAND : public BFSearchNode {
 public:
  BFSearchNodeAND(int var, int val, int depth);
  ~BFSearchNodeAND();

  inline int getType() const {
    return NODE_AND;
  }
  inline int getVar() const {
    assert(current_parent_);
    return current_parent_->getVar();
  }
  inline val_t getVal() const {
    return val_;
  }
  inline int getValue() const {
    return m_nodeValue;
  }
  inline int setValue(double value) {
    m_nodeValue = value;
  }

  inline int getDepth() const {
    assert(current_parent_);
    return current_parent_->getDepth();
  }

  /* empty, but required implementations */
  void setCacheContext(const context_t& c) { }
  const context_t& getCacheContext() const { return empty_context_; }
  void setCacheInst(size_t i) { }
  size_t getCacheInst() const { return 0; }
  void getPST(vector<double>& v) const { }
  void setHeurCache(double* d) { }
  double* getHeurCache() const { return nullptr; }
  void clearHeurCache() { }

  void set_best_child(BFSearchNode* node) { }
  BFSearchNode* get_best_child(BFSearchNode* node) { }

 protected:
  val_t val_;
  static context_t empty_context_;
};

class BFSearchNodeOR : public BFSearchNode {
 public:
  BFSearchNodeOR(int var, int depth);
  ~BFSearchNodeOR();
  int getType() const {
    return NODE_OR;
  }
  int getVar() const {
    return var_;
  }
  val_t getVal() const { return NONE; } // "empty" implementation

  inline int getValue() const {
    return m_nodeValue;
  }
  inline int setValue(double value) {
    m_nodeValue = value;
  }

  inline double getLabel(int val) const {
    assert(heur_cache_);
    return heur_cache_[2*val+1];
  }
  inline void setLabel(int val, double w) {
    assert(heur_cache_);
    heur_cache_[2*val+1] = w;
  }

  inline double* getHeurCache() const {
    return heur_cache_;
  }
  inline void setHeurCache(double* heur_cache) const {
    heur_cache_ = heur_cache;
  }
  inline void clearHeurCache() {
    if (heur_cache_) delete[] heur_cache_;
  }

  inline BFSearchNode* get_best_child() {
    return best_child_;
  }
  inline void set_best_child(BFSearchNode* node) {
    best_child_ = node;
  }

  int getDepth() const {
    return m_depth;
  }
 protected:
  int var_;

  double* heur_cache_;
  BFSearchNode* best_child_;
};

inline BFSearchNode::BFSearchNode() 
  : m_flags(0),
  m_depth(NONE),
  m_nodeValue(ELEM_NAN),
  m_heurValue(ELEM_NAN),
  heur_updated_(ELEM_NAN),
  var_(NONE),
  index_(NONE),
  current_parent_(nullptr),
  hash_key_(0) {
}

inline BFSearchNode::~BFSearchNode() {
  parents_.clear();
  children_.clear();
}

inline BFSearchNodeAND::BFSearchNodeAND(int var, int val)
  : BFSearchNode(), var_(var), val_(val) {
}

inline BFSearchNodeAND::~BFSearchNodeAND() {
}

inline BFSearchNodeOR::BFSearchNodeOR(int var, int depth)
  : BFSearchNode(), var_(var), m_depth(depth) {
}

inline BFSearchNodeOR::~BFSearchNodeOR() {
  clearHeurCache();
}

}  // namespace daoopt


#endif
