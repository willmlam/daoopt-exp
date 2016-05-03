#ifndef BFSEARCHNODE_H_
#define BFSEARCHNODE_H_

#include "SearchNode.h"

namespace daoopt {

class BFSearchNode;
class BFSearchNodeAND;
class BFSearchNodeOR;

class BFSearchNode : public SearchNode {
 public:
  BFSearchNode(SearchNode* parent);
  ~BFSearchNode();

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
    index_ = index;
  }
  inline void increment_index() {
    ++index_;
  }
  inline void decrement_index() {
    --index_;
  }

  inline const std::list<BFSearchNode*>& get_parents() const {
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

  inline const std::list<BFSearchNode*>& get_children() const {
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

  inline bool is_solved() const {
    return is_solved_;
  }
  inline void set_solved(bool flag) {
    is_solved_ = flag;
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

  virtual void setWeight(int val, double w) = 0;
  virtual double getWeight(int val) const = 0;

  virtual int getVar() const { return var_; }
  virtual int getDepth() const { return depth_; }

  virtual std::string ToString() = 0;
 protected:
  // Various flags to indicate the status of the node.
  bool is_fringe_;
  bool is_solved_;
  bool is_terminal_;
  bool is_expanded_;
  bool is_visited_;
  bool is_deadend_;

  // Tightened estimate via lookahead
  double heur_updated_;

  std::list<BFSearchNode*> parents_;
  std::list<BFSearchNode*> children_;

  int var_;

  int depth_;

  int index_;
  BFSearchNode* current_parent_;

  size_t hash_key_;
};

class BFSearchNodeAND : public BFSearchNode {
 public:
  BFSearchNodeAND(SearchNode* parent, int var, int val, int depth);
  ~BFSearchNodeAND();

  inline int getType() const {
    return NODE_AND;
  }

  inline val_t getVal() const {
    return val_;
  }
  inline double getValue() const {
    return m_nodeValue;
  }
  inline void setValue(double value) {
    m_nodeValue = value;
  }

  std::string ToString();

  /* empty, but required implementations */
  void setWeight(int val, double w) { assert(false); }
  double getWeight(int val) const { assert(false); return 0; }
  double getLabel() const { 
    assert(false);
  }
  void addSubSolved(double d) {
    sub_solved_ OP_TIMESEQ d;
  }
  double getSubSolved() const {
    return sub_solved_;
  }
  void setInitialBound(double d) { assert(false); }
  double getInitialBound() const { assert(false); return 0; }
  void setComplexityEstimate(double d) { assert(false); }
  double getComplexityEstimate() const { assert(false); return 0; }

  void setCacheContext(const context_t& c) { }
  const context_t& getCacheContext() const { return SearchNode::emptyCtxt; }
  void setCacheInst(size_t i) { }
  size_t getCacheInst() const { return 0; }
  void getPST(vector<double>& v) const { assert(false); }
  void setHeurCache(double* d) { }
  double* getHeurCache() const { return nullptr; }
  void clearHeurCache() { }
  void setOrderingHeurCache(double* d) { }
  double* getOrderingHeurCache() const { return nullptr; }
  void clearOrderingHeurCache() { }

  void set_best_child(BFSearchNode* node) { }
  BFSearchNode* get_best_child() { }

 protected:
  val_t val_;
  double sub_solved_;
};

class BFSearchNodeOR : public BFSearchNode {
 public:
  BFSearchNodeOR(SearchNode* parent, int var, int depth);
  ~BFSearchNodeOR();

  double getLabel() const { assert(false); return 0; } // no label for OR nodes!
  void addSubSolved(double d) { assert(false); } // not applicable for OR nodes
  double getSubSolved() const { assert(false); return 0; } // not applicable for OR nodes

  void setCacheContext(const context_t& t) { assert(false); }
  const context_t& getCacheContext() const {
    assert(false); return emptyCtxt;
  }
  void setCacheInst(size_t i) { assert(false); }
  size_t getCacheInst() const { assert(false); return 0; }
  void setInitialBound(double d) { assert(false); }
  double getInitialBound() const { assert(false); return 0; }

  void setComplexityEstimate(double d) { assert(false); }
  double getComplexityEstimate() const { assert(false); return 0; }

  void getPST(vector<double>& pst) const;

  int getType() const {
    return NODE_OR;
  }

  val_t getVal() const { return NONE; } // "empty" implementation

  std::string ToString();

  inline double getValue() const {
    return m_nodeValue;
  }
  inline void setValue(double value) {
    m_nodeValue = value;
  }

  inline double getWeight(int val) const {
    assert(heur_cache_);
    return heur_cache_[2*val+1];
  }
  inline void setWeight(int val, double w) {
    assert(heur_cache_);
    heur_cache_[2*val+1] = w;
  }

  inline double* getHeurCache() const {
    return heur_cache_;
  }
  inline void setHeurCache(double* heur_cache) {
    heur_cache_ = heur_cache;
  }
  inline void clearHeurCache() {
    if (heur_cache_) delete[] heur_cache_;
  }

  inline double* getOrderingHeurCache() const {
    return ordering_heur_cache_;
  }
  inline void setOrderingHeurCache(double* ordering_heur_cache) {
    ordering_heur_cache_ = ordering_heur_cache;
  }
  inline void clearOrderingHeurCache() {
    if (ordering_heur_cache_) delete[] ordering_heur_cache_;
  }

  inline BFSearchNode* get_best_child() {
    return best_child_;
  }
  inline void set_best_child(BFSearchNode* node) {
    best_child_ = node;
  }

 protected:

  double* heur_cache_;
  double* ordering_heur_cache_;
  BFSearchNode* best_child_;


};

inline BFSearchNode::BFSearchNode(SearchNode* parent) 
  : SearchNode(parent),
    is_fringe_(false),
    is_solved_(false),
    is_terminal_(false),
    is_expanded_(false),
    is_visited_(false),
    is_deadend_(false),
    heur_updated_(ELEM_NAN),
    var_(NONE),
    index_(0),
    current_parent_(nullptr),
    hash_key_(0) {
}

inline void BFSearchNodeOR::getPST(vector<double>& pst) const {
  const SearchNode* curAND = nullptr;
  const SearchNode* curOR = this;

  pst.clear();

  while (curOR->getParent()) {
    curAND = curOR->getParent();
    val_t val = curAND->getVal();
    double label = dynamic_cast<BFSearchNodeOR*>(curAND->getParent())->
        getWeight(val);
    label OP_TIMESEQ curAND->getSubSolved();
    NodeP* children = curAND->getChildren();
    for (size_t i = 0; i< curAND->getChildCountFull(); ++i) {
      if (children[i] && children[i] != curOR) {
        label OP_TIMESEQ children[i]->getHeur();
      }
    }
    pst.push_back(label);
    
    curOR = curAND->getParent();
    pst.push_back(curOR->getValue());
  }
}

inline std::string BFSearchNodeAND::ToString() {
  std::ostringstream oss;
  oss << "AND node: (x" << getVar() << "," << getVal() << ")"
      << ", h = " << (getHeur() == 0 ? 0 : -getHeur())
      << ", q = " << (getValue() == 0 ? 0 : -getValue())
      << ", oh = " << getOrderingHeur()
      << ", ub = inf" 
      << ", depth = " << getDepth();
  oss << ", children { ";
  for (BFSearchNode* c : get_children()) {
    int var = (int) c->getVar();
    oss << var << ":" << (c->getHeur() == 0 ? 0 : -c->getHeur()) << " ";
  }

  oss << "}";
  oss << ", expanded = " << (is_expanded() ? "YES" : "NO");
  oss << ", solved = " << (is_solved() ? "YES" : "NO");

  return oss.str();
}

inline std::string BFSearchNodeOR::ToString() {
  std::ostringstream oss;
  oss << "OR node: (x" << getVar() << ")"
      << ", h = " << (getHeur() == 0 ? 0 : -getHeur())
      << ", q = " << (getValue() == 0 ? 0 : -getValue())
      << ", oh = " << getOrderingHeur()
      << ", ub = inf" 
      << ", depth = " << getDepth()
      << ", weights { ";
  for (BFSearchNode* c : get_children()) {
    int val = (int) c->getVal();
    oss << val << ":" << -getWeight(val) << " ";
  }

  oss << "}, children { ";
  for (BFSearchNode* c : get_children()) {
    int val = (int) c->getVal();
    oss << val << ":" << (c->getHeur() == 0 ? 0 : -c->getHeur()) << " ";
  }

  oss << "}";
  oss << ", expanded = " << (is_expanded() ? "YES" : "NO");
  oss << ", solved = " << (is_solved() ? "YES" : "NO");

  return oss.str();
}

inline BFSearchNode::~BFSearchNode() {
  parents_.clear();
  children_.clear();
}

inline BFSearchNodeAND::BFSearchNodeAND(SearchNode* parent, int var, int val,
                                        int depth)
  : BFSearchNode(parent) {
  var_ = var;
  val_ = val;
  depth_ = depth;
  sub_solved_ = ELEM_ONE;
}

inline BFSearchNodeAND::~BFSearchNodeAND() {
}


inline BFSearchNodeOR::BFSearchNodeOR(SearchNode* parent, int var, int depth)
  : BFSearchNode(parent) {
  var_ = var;
  depth_ = depth;
  best_child_ = nullptr;
  heur_cache_ = nullptr;
  ordering_heur_cache_ = nullptr;
}

inline BFSearchNodeOR::~BFSearchNodeOR() {
  clearHeurCache();
}

}  // namespace daoopt


#endif
