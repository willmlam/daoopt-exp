#include "FGLP.h"
#include "MiniBucketElim.h"
#include "Problem.h"
#include "ProgramOptions.h"
#include "UAI2012.h"
using namespace std;

string UAI2012::filename = "";
string out_bound_file = "";
time_t _time_start, _time_pre;

int main(int argc, char** argv) {
  unique_ptr<ProgramOptions> po;
  for (int i = 0; i < argc; ++i) cout << argv[i] << ' ';
  cout << endl;
  // parse command line
  po.reset(parseCommandLine(argc, argv));
  if (!po) {
    err_txt("Error parsing command line.");
    return false;
  }

  if (po->seed == NONE) po->seed = time(0);
  rand::seed(po->seed);

  // Load problem
  unique_ptr<Problem> p(new Problem());
  if (!p->parseUAI(po->in_problemFile, po->in_evidenceFile, po->collapse))
    exit(0);
  cout << "Created problem with " << p->getN() << " variables and " << p->getC()
       << " functions." << endl;

  p->removeEvidence();
  cout << "Removed evidence, now " << p->getN() << " variables and " 
       << p->getC() << " functions." << endl;

  vector<int> elim;
  bool orderFromFile = false;
  if (!po->in_orderingFile.empty()) {
    orderFromFile = p->parseOrdering(po->in_orderingFile, elim);
  }

  if (!orderFromFile) return 1;

  // do FGLP preprocessing here
  unique_ptr<FGLP> fglp(new FGLP(p.get(), po->useNullaryShift));
  fglp->Run(po->mplp, po->mplps, po->mplpt);

  p->replaceFunctions(fglp->factors(), true);

  // Output reduced network?
  if (!po->out_reducedFile.empty()) {
    cout << "Writing reduced network to " << po->out_reducedFile << endl;
    p->writeUAI(po->out_reducedFile);
  }

  return 0;
}
