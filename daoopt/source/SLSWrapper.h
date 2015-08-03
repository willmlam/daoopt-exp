/*
 * SLSWrapper.h
 *
 *  Wrapper around the SLS4MPE code by Frank Hutter, so it can
 *  be used as a library.
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
 *  Created on: Nov 17, 2011
 *      Author: Lars Otten <lotten@ics.uci.edu>
 */

#ifndef SLSWRAPPER_H_
#define SLSWRAPPER_H_

#include "Problem.h"

#ifdef ENABLE_SLS

#include "sls4mpe/main_algo.h"
#include "sls4mpe/timer.h"
#include "sls4mpe/global.h"
#include "sls4mpe/ProblemReader.h"

namespace daoopt {

class SLSWrapper {
protected:
  double m_likelihood;
  int* m_assignment;
  Problem* m_problem;

public:
  bool init(Problem* prob, int iter, int time);
  virtual bool run();
  double getSolution(vector<val_t>* tuple = NULL) const;
  void reportSolution(double cost, int num_vars, int* assignment);

public:
  SLSWrapper();
  virtual ~SLSWrapper();
};


/* Inline definitions */
inline SLSWrapper::SLSWrapper() :
    m_likelihood(0.0), m_assignment(NULL) {
  /* nothing here */
}

inline SLSWrapper::~SLSWrapper() {
  if (m_assignment)
    delete[] m_assignment;
}

class SLSWrapperInc : public SLSWrapper {
protected:
	bool m_solution_found;
	//bool m_solution_optimal; // sls4mpe has Hybrid of MB*, ILS, GLS+, proving optimality by UB==LB with restarts! could use this...
	int m_num_trial;
	double m_time_elapsed;
	double m_converge_rate; // maximum time allowed for convergence, from program option

public:
	bool run();
	SLSWrapperInc(double converge_rate = 2.0) : SLSWrapper(), 
		m_solution_found(false), m_num_trial(0), m_time_elapsed(0.0), m_converge_rate(converge_rate) {}
	virtual ~SLSWrapperInc() {}
};

class SLSWrapperHybrid : public SLSWrapper {
	/* Wrapper for Hybrid SLS, it proves optimality! */

public:
	bool run();
	SLSWrapperHybrid() : SLSWrapper() {}
	virtual ~SLSWrapperHybrid() {}
};


}  // namespace daoopt

#endif /* ENABLE_SLS */

#endif /* SLSWRAPPER_H_ */
