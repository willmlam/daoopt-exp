#include "_base.h"

#include "Problem.h"
#include "MiniBucketElim.h"
#include "FGLP.h"
#include "ProgramOptions.h"
#include <random>
using namespace std;

// A small driver to only evaluate the heuristics used in daoopt
//

random_device rd;
mt19937 gen(1234567890);
uniform_real_distribution<> dis(0,1);

shared_ptr<Problem> runFGLP(const shared_ptr<Problem> &p, const shared_ptr<ProgramOptions> &po, const vector<int> &ordering, int var) {
    map<int,val_t> mAssn;
    mAssn[var] = 0;

    shared_ptr<FGLP> fglp(new FGLP(p->getN(),p->getDomains(),p->getFunctions(),ordering,mAssn));
    fglp->run(po->ndfglp,po->ndfglps,po->ndfglpt);

    shared_ptr<Problem> pNew(new Problem(p.get()));
    pNew->replaceFunctions(fglp->getFactors());
    return pNew;
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
    if (!p->parseUAI(po->in_problemFile, po->in_evidenceFile, po->collapse))
        exit(0);
    cout << "Created problem with " << p->getN()
        << " variables and " << p->getC() << " functions." << endl;
    p->removeEvidence();
    cout << "Removed evidence, now " << p->getN()
        << " variables and " << p->getC() << " functions." << endl;

    if (po->perturb > 0) {
        p->perturbDeterminism(po->perturb); 
    }


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


    // Evaluation code

    int nvAssign = po->cutoff_depth; // ...Borrowing the option variable

    vector<val_t> assignment(p->getN(),NONE);
    PseudotreeNode *node = pt->getRoot();
    for (int i=0; i<nvAssign; ++i) {
        node=node->getChildren()[0];
        assignment[node->getVar()] = 0;
    }

    vector <int> linearOrdering;
    for (int i=0; i<elim.size();++i) linearOrdering.push_back(i);

    // Preprocess problem first to convergence
    shared_ptr<FGLP> fglp(new FGLP(p->getN(),p->getDomains(),p->getFunctions(),linearOrdering));
    fglp->run(po->ndfglp,po->ndfglps,po->ndfglpt);

    p->replaceFunctions(fglp->getFactors());

    map<int,val_t> mAssn;
    mAssn[linearOrdering[0]] = 0;

    shared_ptr<FGLP> fglp2(new FGLP(p->getN(),p->getDomains(),p->getFunctions(),linearOrdering,mAssn));
    fglp2->run(po->ndfglp,po->ndfglps,po->ndfglpt);
    fglp2->setVerbose(true);

    p->condition(mAssn);

    const vector<Function*> &fOrig = p->getFunctions();
    const vector<Function*> &fShifted = fglp2->getFactors();

    assert(fOrig.size() == fShifted.size());

    vector<double> changes(fOrig.size(),0);

    for (size_t i = 0; i < fOrig.size(); ++i) {
        size_t tSize = fOrig[i]->getTableSize();
        double oMax = -std::numeric_limits<double>::infinity();
        double sMax = -std::numeric_limits<double>::infinity();
        for (size_t j = 0; j < tSize; ++j) {
            oMax = max(oMax,fOrig[i]->getTable()[j]);
            sMax = max(sMax,fShifted[i]->getTable()[j]);
        }
        changes[i] += sMax - oMax;
    }

    for (size_t i = 0; i < fOrig.size(); ++i) {
        cout << *fOrig[i] << ", change: " << changes[i] << endl;
    }




    return 0;
}
