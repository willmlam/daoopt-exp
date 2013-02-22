#include "MBEHeuristicInstance.h"

void MBEHeuristicInstance::populateMessages(int var, vector<bool> &visited) {
    ConditionedMessages *cm = m_computedMessages[var];
    assert(cm);

    const vector< vector<int> *> paths = cm->getPaths();
    unsigned i, j;
    int v;
    for (i = 0; i < paths.size(); ++i) {
        for (j = 0; j < paths[i]->size() - 1; ++j) {
            v = paths[i]->at(j);
            m_intermediate[v].push_back(cm->getFunctions()[i]);
            //if (v == m_var) return;
        }
        v = paths[i]->at(j);
        m_augmented[v].push_back(cm->getFunctions()[i]);
        m_accurateHeurIn[v] = 
            cm->isAccurate() &&
            (!visited[v] || m_accurateHeurIn[v]);
        visited[v] = true;
    }
}
