/*
 * daoopt.cpp
 *
 *  Copyright (C) 2008-2012 Lars Otten
 *  This file is part of DAOOPT.
 *
 *  DAOOPT is free software: you can redistribute it and/or modify it
 *  under the terms of the GNU General Public License as published by
 *  the Free Software Foundation, either version 3 of the License, or
 *  (at your option) any later version.
 *
 *  DAOOPT is distributed in the hope that it will be useful, but
 *  WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU General Public License for more details.
 *
 *  You should have received a copy of the GNU General Public License
 *  along with DAOOPT.  If not, see <http://www.gnu.org/licenses/>.
 *
 *  Created on: Oct 13, 2008
 *      Author: Lars Otten <lotten@ics.uci.edu>
 */

#include "Main.h"
#include "MiniBucket.h"
#include "Function.h"
#include "Problem.h"


using namespace std;
using namespace daoopt;

/* define to enable diagnostic output of memory stats */
//#define MEMDEBUG

void printFunction(Function *f) {
    cout << "Scope: ";
    for (auto k : f->getScopeVec()) {
        cout << k << " ";
    }
    cout << endl << endl;
    for (unsigned int k = 0; k < f->getTableSize(); ++k) {
        cout << ELEM_DECODE(f->getTable()[k]) << endl;
    }
    cout << endl << endl;
}

int main(int argc, char** argv) {
    /*
    string uaifile = "test.uai";
    string evidfile = "";
    Problem *p = new Problem();
    p->parseUAI(uaifile, evidfile);
    const vector<Function*> &functions = p->getFunctions();
    for (int i = 0; i < int(functions.size()); ++i) {
        printFunction(functions[i]);
    }
    cout << endl << endl << endl;

    map<int,val_t> cond;
    cond[2] = 0;

    MiniBucket mb(0,10,p);
    MiniBucket mb2(0,10,p);
    mb.addFunction(functions[4]);
    mb2.addFunction(functions[5]);

    set<int> sepScope(intersection(mb.getJointScope(),mb2.getJointScope()));
    set<int> elimVars1(setminus(mb.getJointScope(), sepScope));
    set<int> elimVars2(setminus(mb2.getJointScope(), sepScope));
    
//    Function *f1mm = mb.eliminate(true,elimVars1);
    Function *f1mm = mb.conditionEliminate(true,cond,elimVars1);
    printFunction(f1mm);
    Function *f2mm = mb2.eliminate(true,elimVars2);
    printFunction(f2mm);

    size_t tablesize = 1;
    for (set<int>::iterator sit=sepScope.begin(); sit!=sepScope.end(); ++sit)
        tablesize *= p->getDomainSize(*sit);

    double *avgMMTable = new double[tablesize];
    for (unsigned i = 0; i < tablesize; ++i) {
        avgMMTable[i] = ELEM_ONE;
        avgMMTable[i] OP_TIMESEQ f1mm->getTable()[i];
        avgMMTable[i] OP_TIMESEQ f2mm->getTable()[i];
        cout << ELEM_DECODE(avgMMTable[i]) << endl;
        avgMMTable[i] = OP_ROOT(avgMMTable[i],2);
    }
    int ammid = 0;
    Function *avgMMf = new FunctionBayes(ammid,p,sepScope,avgMMTable,tablesize);
    printFunction(avgMMf);

    Function *msg1 = mb.conditionEliminateMM(true,cond,f1mm,avgMMf);
    printFunction(msg1);

    Function *msg2 = mb2.eliminateMM(true,f2mm,avgMMf);
    printFunction(msg2);
    
     
    return 0;
    */

  Main main;

  if (!main.start())
    exit(1);
  if (!main.parseOptions(argc, argv))
    exit(1);
  if (!main.outputInfo())
    exit(1);
  if (!main.loadProblem())
    exit(1);
  if (!main.runSLS())
    exit(1);
  if (!main.findOrLoadOrdering())
    exit(1);
  if (!main.initDataStructs())
    exit(1);
  if (!main.compileHeuristic())
    exit(1);
  if (!main.runLDS())
    exit(1);
  if (!main.finishPreproc())
    exit(1);
  if (!main.runSearch())
    exit(1);
  if (!main.outputStats())
    exit(1);

  return 0;

}

