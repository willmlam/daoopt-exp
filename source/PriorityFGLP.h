#ifndef PRIORITYFGLP_H_
#define PRIORITYFGLP_H_

#include "FGLP.h"
#include "mex/indexedHeap.h"

class PriorityFGLP : public FGLP {
public:
    PriorityFGLP(Problem *p, bool use_nullary_shift = false);
    PriorityFGLP(PriorityFGLP *parent_fglp, const map<int,val_t> &assignment);

    virtual void Run(int max_updates, double max_time, 
            double tolerance=DEFAULT_TOLERANCE);

    inline const mex::indexedHeap &var_priority() const { return var_priority_; }
    inline void set_var_priority(const mex::indexedHeap &vp) { var_priority_ = vp; }
    inline ~PriorityFGLP();
protected:
    virtual void Condition(const vector<Function*> &fns,
            const map<int,val_t> &assignment);

    double MessageDist(double *m1, double *m2, double ns, int var);

    mex::indexedHeap var_priority_;

    
};

inline PriorityFGLP::~PriorityFGLP() {
}

#endif
