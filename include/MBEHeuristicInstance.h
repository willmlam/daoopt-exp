// Data structure to compile a heuristic that can share messages with
// other compiled heuristics

#ifndef MBEHEURISTICINSTANCE_H_
#define MBEHEURISTICINSTANCE_H_

#include "Function.h"
#include "utils.h"

class MBEHeuristicInstance;
class SearchNode;

// Class to store the outgoing messages from a bucket
class ConditionedMessages {
    // The messages that are sent from this bucket
    vector<Function*> m_functions;

    // Contains the path that the corresponding message 
    vector<vector<int> *> m_paths;

    // Track which compilation instance computed these messages
    MBEHeuristicInstance *m_owner;

    // Accurate: if the messages contain the same information as 
    // unbounded BE
    bool m_isAccurate;


    public:
    ConditionedMessages(MBEHeuristicInstance *owner, bool isAccurate) 
        : m_owner(owner), m_isAccurate(isAccurate) {
    }

    const vector<Function*> &getFunctions() const {
        return m_functions;
    }

    const vector<vector<int> *> &getPaths() const {
        return m_paths;
    }

    // Return the bucket that message i goes into
    const int getTargetVar(int i) const {
        return m_paths[i]->back();
    }

    bool isAccurate() const {
        return m_isAccurate;
    }

    void addFunction(Function *f, vector<int> *path) {
        m_functions.push_back(f);
        m_paths.push_back(path);
    }

    MBEHeuristicInstance *getOwner() const {
        return m_owner;
    }

    ~ConditionedMessages() {
        for (unsigned i = 0; i < m_functions.size(); ++i) {
            delete m_functions[i];
        }
        m_functions.clear();
        for (unsigned i = 0; i < m_paths.size(); ++i) {
            delete m_paths[i];
        }
        m_paths.clear();
    }
};

class MBEHeuristicInstance {
    static int maxNumActive;
    static int currentNumActive;
    int m_var;
    int m_depth;
    vector<vector<Function*> > m_augmented;
    vector<vector<Function*> > m_intermediate;
    vector<ConditionedMessages*> m_computedMessages;

    vector<bool> m_accurateHeurIn;

    //map<int,val_t> m_assignment;

    // Keep track of the original node using this heuristic
    SearchNode *m_owner;

    // Parent heursitic
    MBEHeuristicInstance *m_parent;

    public:
    MBEHeuristicInstance(int nVars, int var, MBEHeuristicInstance *parent);

    inline static int getCurrentNumActive() {
        return currentNumActive;
    }
    inline static int getMaxNumActive() {
        return maxNumActive;
    }

    inline int getVar() const {
        return m_var;
    }

    inline int getDepth() const {
        return m_depth;
    }

    inline void setDepth(int depth) {
        m_depth = depth;
    }

    inline vector< vector<Function*> > &getAugmented() {
        return m_augmented;
    }

    inline vector< vector<Function*> > &getIntermediate() {
        return m_intermediate;
    }

    inline vector<ConditionedMessages*> &getComputedMessages() {
        return m_computedMessages;
    }

    inline vector<bool> &getAccurateHeurIn() {
        return m_accurateHeurIn;
    }

    inline SearchNode *getOwner() const {
        return m_owner;
    }

    inline MBEHeuristicInstance *getParent() const {
        return m_parent;
    }

    inline void setOwner(SearchNode *owner) {
        m_owner = owner;
    }

    /*
    inline const map<int,val_t> &getAssignment() const {
        return m_assignment;
    }

    inline void setAssignment(const map<int,val_t> &assignment) {
        m_assignment = assignment;
    }

    bool isConsistent(const map<int,val_t> &assignment) const {
        map<int,val_t>::const_iterator it = m_assignment.begin();
        for (; it != m_assignment.end(); ++it) {
            map<int,val_t>::const_iterator fit = assignment.find(it->first);
            if (fit != assignment.end() && fit->second != it->second)
                return false;
        }
        return true;
    }
    */


    void populateMessages(int var, vector<bool> &visited);

    ~MBEHeuristicInstance();
};

inline MBEHeuristicInstance::MBEHeuristicInstance(int nVars, int var, MBEHeuristicInstance *parent) 
: m_var(var), m_depth(-1), m_accurateHeurIn(nVars, true), m_parent(parent) {
    m_augmented.resize(nVars);
    m_intermediate.resize(nVars);
    m_computedMessages.resize(nVars, NULL);
    ++currentNumActive;
    if (currentNumActive > maxNumActive) maxNumActive = currentNumActive;
}
inline MBEHeuristicInstance::~MBEHeuristicInstance() {
    for (unsigned i = 0; i < m_computedMessages.size(); ++i) {
        if (m_computedMessages[i] && this == m_computedMessages[i]->getOwner())
            delete m_computedMessages[i];
    }
    --currentNumActive;
}


#endif

