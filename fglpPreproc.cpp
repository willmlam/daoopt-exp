#include "Problem.h"
#include "MiniBucketElim.h"
#include "FGLP.h"
#include "ProgramOptions.h"
using namespace std;

// A small driver to only evaluate the heuristics used in daoopt


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

    vector<int> elim;
    bool orderFromFile = false;
    if (!po->in_orderingFile.empty()) {
        orderFromFile = p->parseOrdering(po->in_orderingFile, elim);
    }

    // do FGLP preprocessing here
    shared_ptr<FGLP> fglp(new FGLP(p->getN(),p->getDomains(),p->getFunctions(),elim));
    fglp->run(po->mplp, po->mplps, po->mplpt);
            
    p->replaceFunctions(fglp->getFactors(),true);

    // Output reduced network?
    if (!po->out_reducedFile.empty()) {
        cout << "Writing reduced network to " << po->out_reducedFile << endl;
        p->writeUAI(po->out_reducedFile);
    }

    return 0;
}
