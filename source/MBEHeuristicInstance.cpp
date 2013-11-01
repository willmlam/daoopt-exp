#include "MBEHeuristicInstance.h"

int MBEHeuristicInstance::currentNumActive = 0;
int MBEHeuristicInstance::maxNumActive = 0;
double MBEHeuristicInstance::currentMemory = 0.0;
double MBEHeuristicInstance::maxMemory = 0.0;

void MBEHeuristicInstance::populateMessages(int var, vector<bool> &visited) {
    ConditionedMessages *cm = m_computedMessages[var];
    assert(cm);

    const vector< vector<int> *> paths = cm->getPaths();
    unsigned i, j;
    int v;
    //cout << endl;
    for (i = 0; i < paths.size(); ++i) {
        for (j = 0; j < paths[i]->size() - 1; ++j) {
            v = paths[i]->at(j);
            //cout << " " << v;
            m_intermediate[v].push_back(cm->getFunctions()[i]);
            //if (v == m_var) return;
        }
        v = paths[i]->at(j);
        //cout << " " << v << endl;
        m_augmented[v].push_back(cm->getFunctions()[i]);
        m_accurateHeurIn[v] = 
            cm->isAccurate() &&
            (!visited[v] || m_accurateHeurIn[v]);
        visited[v] = true;
    }
}
