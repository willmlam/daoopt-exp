#include "_base.h"

#include "Problem.h"
#include "MiniBucketElim.h"
#include "FGLP.h"
#include "ProgramOptions.h"
#include "Graph.h"
#include <random>
using namespace std;

// A small driver to only evaluate the heuristics used in daoopt
//

random_device rd;
mt19937 gen(1234567890);
uniform_real_distribution<> dis(0,1);

vector<int> bfsOrdering(const shared_ptr<Problem> &p, const vector<int> &ordering, int var) {
    Graph g(p->getN());
    for (Function *f : p->getFunctions()) {
        g.addClique(f->getScopeVec());
    }

    vector<int> retOrdering;
    SETCLASS<int> nodes(++ordering.begin(), ordering.end());

    queue<int> bfs;
    bfs.push(var);
    nodes.erase(bfs.front());
    while (!bfs.empty()) {
        int currentVar = bfs.front();
        cout << currentVar << endl;
        retOrdering.push_back(currentVar); 
        bfs.pop();
        const set<int> &nb = g.getNeighbors(currentVar);
        for (int vn : nb) {
            if (nodes.find(vn) != nodes.end()) {
                bfs.push(vn); nodes.erase(vn);
            }
        }
    }
    return retOrdering;
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

    vector <int> linearOrdering;
    for (size_t i=0; i<elim.size();++i) linearOrdering.push_back(i);

    // Preprocess problem first to convergence
    shared_ptr<FGLP> fglp(new FGLP(p->getN(),p->getDomains(),p->getFunctions(),linearOrdering,po->useNullaryShift));
    fglp->run(1000,po->ndfglps,po->ndfglpt);


    /*
    for (size_t fid = 0; fid < p->getFunctions().size(); ++fid) {
        Function *f = p->getFunctions()[fid];
        Function *fs = fglp->getFactors()[fid];
        cout << *f << endl;
        for (size_t k = 0; k < f->getTableSize(); ++k) {
            cout << f->getTable()[k] << " " << fs->getTable()[k] << endl;
        }
        cout << endl;
    }
    */

    p->replaceFunctions(fglp->getFactors(),true);


    vector<double> fMax(p->getFunctions().size(),-std::numeric_limits<double>::infinity());
    for (size_t i = 0; i < p->getFunctions().size(); ++i) {
        Function *f = p->getFunctions()[i];
        for (size_t j = 0; j < f->getTableSize(); ++j) {
            fMax[f->getId()] = max(fMax[f->getId()],f->getTable()[j]);
        }
    }

    map<int,val_t> mAssn;
    int vA = 0;
    mAssn[vA] = 0;

    vector<int> bOrdering = bfsOrdering(p,linearOrdering,vA);

    shared_ptr<FGLP> fglp2(new FGLP(p->getN(),p->getDomains(),p->getFunctions(),linearOrdering,mAssn,po->useNullaryShift));
//fglp2->run(10,po->ndfglps,po->ndfglpt);
    fglp2->setVerbose(true);
    fglp2->run(po->ndfglp,po->ndfglps,po->ndfglpt);



    p->condition(mAssn);
    vector<double> fMaxAfter(p->getFunctions().size(),-std::numeric_limits<double>::infinity());
    for (size_t i = 0; i < p->getFunctions().size(); ++i) {
        Function *f = p->getFunctions()[i];
        for (size_t j = 0; j < f->getTableSize(); ++j) {
            fMaxAfter[f->getId()] = max(fMaxAfter[f->getId()],f->getTable()[j]);
        }
    }

    const vector<Function*> &fOrig = p->getFunctions();
    const vector<Function*> &fShifted = fglp2->getFactors();

    assert(fOrig.size() == fShifted.size());

    vector<double> changes(fOrig.size(),0);
    vector<double> changesMin(fOrig.size(),0);

    double totalPositiveChange = 0;
    double totalNegativeChange = 0;


    for (size_t i = 0; i < fOrig.size(); ++i) {
        size_t tSize = fOrig[i]->getTableSize();
        double oMax = -std::numeric_limits<double>::infinity();
        double sMax = -std::numeric_limits<double>::infinity();
        double oMin = std::numeric_limits<double>::infinity();
        double sMin = std::numeric_limits<double>::infinity();
        for (size_t j = 0; j < tSize; ++j) {
            oMax = max(oMax,fOrig[i]->getTable()[j]);
            sMax = max(sMax,fShifted[i]->getTable()[j]);
            oMin = min(oMin,fOrig[i]->getTable()[j]);
            sMin = min(sMin,fShifted[i]->getTable()[j]);
        }
        changes[i] = sMax - oMax;
        changesMin[i] = sMin - oMin;
        if (changes[i] > 0) totalPositiveChange += changes[i];
        else totalNegativeChange += changes[i];
    }

    for (size_t i = 0; i < fMax.size(); ++i) {
        double change = fMaxAfter[i] - fMax[i];
        if (change == 0.0) continue;
        cout << i << " : " << fMax[i] << " " << fMaxAfter[i] << " (change: " << change << ")" << endl;
    }


    /*
    vector<int> fOrdering;
    for (int v : bOrdering) {
        for (Function *f : fOrig) {
            if (f->hasInScope(v)) {
                if (find(fOrdering.begin(),fOrdering.end(),f->getId()) == fOrdering.end()) {
                    fOrdering.push_back(f->getId());
                }
            }
        }
    }
    assert(fOrdering.size() <= fOrig.size());
    */
//    for (int fid : fOrdering) {
    for (size_t fid = 0; fid < fOrig.size(); ++fid) {
        if (changes[fid] == 0.0) continue;
        cout << *fOrig[fid] << endl;
        cout << "change: " << changes[fid] /*<< " | " << changesMin[fid]*/ << endl;
    }

    /*
    for (size_t fid = 0; fid < fOrig.size(); ++fid) {
        Function *f = fOrig[fid];
        Function *fs = fShifted[fid];
        cout << *f << endl;
        for (size_t k = 0; k < f->getTableSize(); ++k) {
            cout << f->getTable()[k] << " " << fs->getTable()[k] << endl;
        }
        cout << endl;
    }
    */

    cout << "total positive change: " << totalPositiveChange << endl;
    cout << "total negative change: " << totalNegativeChange << endl;
    cout << "net change: " << totalPositiveChange + totalNegativeChange << endl;

    vector<double> fMaxFinal(fShifted.size(),-std::numeric_limits<double>::infinity());
    vector<double> fMinFinal(fShifted.size(),std::numeric_limits<double>::infinity());
    for (size_t i = 0; i < fShifted.size(); ++i) {
        Function *f = fShifted[i];
        for (size_t j = 0; j < f->getTableSize(); ++j) {
            fMaxFinal[f->getId()] = max(fMaxFinal[f->getId()],f->getTable()[j]);
            fMinFinal[f->getId()] = min(fMinFinal[f->getId()],f->getTable()[j]);
        }
    }

    for (size_t i = 0; i > fMaxFinal.size(); ++i) {
        cout << i << " : " << fMaxFinal[i] /*<< " , " << fMinFinal[i] */<< endl;
    }



    return 0;
}
