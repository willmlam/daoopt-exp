/*
 * DaooptInterface.cpp
 *
 *  Created on: May 24, 2014
 *      Author: lars
 */

#include "DaooptInterface.h"

namespace daoopt {

DaooptInterface::DaooptInterface()
    : m_initialized(false), m_preprocessed(false) {}

bool DaooptInterface::initialize() {
  m_main.reset(new daoopt::Main());
  m_initialized = true;
  m_preprocessed = false;
  return true;
}

bool DaooptInterface::preprocess(int argc, char** argv) {
  if (!m_initialized || !m_main.get())
    return false;
  if (!m_main->start())
    return false;
  if (!m_main->parseOptions(argc, argv))
    return false;
  return preprocessCommon();
}

bool DaooptInterface::preprocess(const ProgramOptions& options) {
  if (!m_initialized || !m_main.get())
    return false;
  if (!m_main->start())
    return false;
  if (!m_main->setOptions(options))
    return false;
  return preprocessCommon();
}

bool DaooptInterface::preprocessCommon() {
  if (!m_main->outputInfo())
    return false;
  if (!m_main->loadProblem())
    return false;
  if (!m_main->preprocessHeuristic())
    return false;
  if (!m_main->runSLS())
    return false;
  if (!m_main->findOrLoadOrdering())
    return false;
  if (!m_main->initDataStructs())
    return false;
  if (!m_main->compileHeuristic())
    return false;
  if (!m_main->runLDS())
    return false;
  if (!m_main->finishPreproc())
    return false;
  m_preprocessed = true;
  return true;
}

int DaooptInterface::stopSLS(void)
{
	if (! m_initialized || ! m_preprocessed) 
		return 1 ;
	m_main.get()->stopSLS() ;
	return 0 ;
}

int DaooptInterface::runSLS(int slsIter_, int slsTimePerIter_)
{
	if (! m_initialized || ! m_preprocessed) 
		return 1 ;
	m_main.get()->setSLSOptions(slsIter_, slsTimePerIter_) ;
	if (! m_main->runSLS())
		return 2 ;
	return 0 ;
}

double DaooptInterface::estimate(size_t sampleSize) const {
  assert(m_initialized && m_preprocessed);
  return m_main->runEstimation(sampleSize);
}

bool DaooptInterface::solve(size_t nodeLimit) {
  assert(m_initialized && m_preprocessed);
  if (!m_main->runSearch(nodeLimit))
    return false;
  if (!m_main->outputStats())
    return false;
  return true;
}

double DaooptInterface::getSolution(vector<int>* assignment) const {
  double cost = m_main->getSolution();
  if (assignment != NULL) {
//    const vector<val_t>& assg = m_main->getSolutionAssg();
    vector<val_t> assg ;
	m_main->getSolutionAssgOrg(assg);
    assignment->resize(0);  // Kill existing contents
    for (vector<val_t>::const_iterator it = assg.begin(); it != assg.end(); ++it) {
      assignment->push_back((int) *it);
    }
  }
  return cost;
}

int DaooptInterface::outputStatistics(void)
{
	return m_main->outputStatsToFile() ;
}

}  // namespace daoopt
