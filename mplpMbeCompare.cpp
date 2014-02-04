#include "_base.h"

#include "Problem.h"
#include "MiniBucketElim.h"
#include "FGLP.h"
#include "mex/mplp.h"
#include "ProgramOptions.h"
#include <random>
using namespace std;

// A small driver to only evaluate the heuristics used in daoopt
//

random_device rd;
mt19937 gen(1234567890);
uniform_real_distribution<> dis(0,1);


// Changes the assignments of already assigned values
void shuffleAssignments(const vector<val_t> &domains, vector<val_t> &assignment) {
    for (unsigned int i=0; i<domains.size(); ++i) {
        if (assignment[i] != NONE) {
            assignment[i] = int(dis(gen)*domains[i]);
        }
    }
}

double evaluateMBE(const shared_ptr<MiniBucketElim> &mbe, int var, const vector<val_t> &assignment) {
    // Evaluate heuristic
    const auto &mesAug = mbe->getRootInstance()->getAugmented(); 
    const auto &mesInt = mbe->getRootInstance()->getIntermediate();

    double h = ELEM_ONE;
    for (Function *f : mesAug[var]) {
        h OP_TIMESEQ f->getValue(assignment);
    }
    for (Function *f : mesInt[var]) {
        h OP_TIMESEQ f->getValue(assignment);
    }

    return h;
}

double evaluateMPLP(const shared_ptr<Problem> &p, const shared_ptr<ProgramOptions> &po, const vector<int> &ordering, const vector<val_t> &assignment, double &runtime, int &iters) {
    map<int,val_t> mAssn;
    for (unsigned int i=0; i<assignment.size(); ++i) {
        if (assignment[i] != NONE) mAssn[i] = assignment[i];
    }
    cout << mAssn << endl;

    double globalConstant = ELEM_ONE;
    int arityZeroFnIdx = NONE;
    const vector<Function*> &fns = p->getFunctions();
    vector<Function*> newFunctions;
    for (size_t i=0; i<fns.size();++i) {
        if (fns[i]->getArity() == 0) {
            globalConstant OP_TIMESEQ fns[i]->getTable()[0];
            arityZeroFnIdx = i;
            continue;
        }
        Function *new_fn = fns[i]->substitute(mAssn);
        if (new_fn->isConstant()) {
            globalConstant OP_TIMESEQ new_fn->getTable()[0];
            delete new_fn;
        } else {
            newFunctions.push_back(new_fn);
        }
    }
    Function *constFun = fns[arityZeroFnIdx]->clone();
    constFun->getTable()[0] = globalConstant;
    newFunctions.push_back(constFun);

    mex::vector<mex::Factor> fs(newFunctions.size());
    for (unsigned int i=0;i<newFunctions.size();++i) 
        fs[i] = newFunctions[i]->asFactor().exp();

    mex::mplp _mplp(fs);
    _mplp.setProperties("Schedule=Fixed,Update=Var,StopIter=100,StopObj=-1,StopMsg=-1,StopMsg=-1,StopTime=-1");
    _mplp.init();
    cout << "Initial UB: " << _mplp.ub() << endl;
    char opt[50];
    if (po->ndfglp > 0) { 
        sprintf(opt,"StopIter=%d",po->ndfglp); _mplp.setProperties(opt); 
    }
    if (po->ndfglps > 0) { 
        sprintf(opt,"StopTime=%f",po->ndfglps); _mplp.setProperties(opt); 
    }
    if (po->ndfglpt > 0) { 
        sprintf(opt,"StopObj=%f",po->ndfglpt); _mplp.setProperties(opt); 
    }

    _mplp.run();

    runtime = 0;
    iters = 0;
    return _mplp.ub();
        
}



int main(int argc, char **argv) {
    shared_ptr<ProgramOptions> po;
    for (int i=0; i<argc; ++i)
        cout << argv[i] << ' ';
    cout << endl;
    // parse command line
    po.reset(parseCommandLine(argc, argv));
    if (!po) {
        err_txt("Error parsing command line.");
        return false;
    }

    if (po->seed == NONE)
        po->seed = time(0);
    rand::seed(po->seed);

    // Load problem
    shared_ptr<Problem> p(new Problem());
    if (!p->parseUAI(po->in_problemFile, po->in_evidenceFile))
        exit(0);
    cout << "Created problem with " << p->getN()
        << " variables and " << p->getC() << " functions." << endl;
    p->removeEvidence();
    cout << "Removed evidence, now " << p->getN()
        << " variables and " << p->getC() << " functions." << endl;


    Graph g(p->getN());
    for(Function *f : p->getFunctions()) {
        g.addClique(f->getScopeVec());
    }
    cout << "Graph with " << g.getStatNodes() << " nodes and "
        << g.getStatEdges() << " edges created." << endl;

    // Find variable ordering
    vector<int> elim;
    int w = numeric_limits<int>::max();
    bool orderFromFile = false;
    if (!po->in_orderingFile.empty()) {
        orderFromFile = p->parseOrdering(po->in_orderingFile, elim);
    }

    shared_ptr<Pseudotree> pt(new Pseudotree(p.get(),po->subprobOrder));

    if (orderFromFile) {
        pt->build(g, elim, po->cbound);
        w = pt->getWidth();
        cout << "Read elimination ordering from file " << po->in_orderingFile
            << " (" << w << '/' << pt->getHeight() << ")." << endl;
    }
    else {
        cout << "Ordering required." << endl;
        return 0;
    }

    cout << "Rebuilding pseudo tree as chain." << endl;
    pt->buildChain(g,elim,po->cbound);

    p->addDummy();
    pt->addFunctionInfo(p->getFunctions());
    pt->addDomainInfo(p->getDomains());


    // Build MBE first once.
    shared_ptr<MiniBucketElim> mbe(new MiniBucketElim(p.get(),pt.get(),po.get(),po->ibound));
    // test the memory usage first
    size_t sz = mbe->build(NULL,false);
    double szMB = (sz/(1024*1024.0)) * sizeof(double);
    cout << "\tMB size: " << szMB << " MBytes" << endl;
    if (po->memlimit != NONE && szMB > po->memlimit) {
        cout << "i-bound exceeds memory limits!";
        return 0;
    }
    mbe->build();



    // Evaluation code

    int nvAssign = po->cutoff_depth; // ...Borrowing the option variable

    vector<val_t> assignment(p->getN(),NONE);
    PseudotreeNode *node = pt->getRoot();
    for (int i=0; i<nvAssign; ++i) {
        node=node->getChildren()[0];
        assignment[node->getVar()] = 0;
    }

    cout << endl;
    cout << "<MBEmemory>=" << szMB << endl << endl;
    cout << "MBEValue,MPLPValue,difference,MPLPiters,MPLPtime" << endl;
    for (unsigned int k=0;k<1;++k) {
        shuffleAssignments(p->getDomains(), assignment);
        double mplpRunTime = 0;
        int mplpIters = 0;

        double mbeValue = evaluateMBE(mbe,node->getVar(),assignment);
        double mplpValue = evaluateMPLP(p,po,elim,assignment,mplpRunTime,mplpIters);
        double diff = mplpValue - mbeValue;
        if(std::isnan(diff)) diff = 0;


        cout << mbeValue << "," << mplpValue << "," << diff 
            << "," << mplpIters << "," << mplpRunTime << endl;
    }

    return 0;
}
