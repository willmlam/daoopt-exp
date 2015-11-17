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

bool DaooptInterface::preprocessAndSLS(int argc, char** argv) {
	if (!m_initialized || !m_main.get())
		return false;
	if (!m_main->start())
		return false;
	if (!m_main->parseOptions(argc, argv))
		return false;

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
	if (!m_main->compileHeuristicBeforeBuild())
		return false;

	return true;
}

bool DaooptInterface::preprocessAndSLS(const ProgramOptions& options) {
	if (!m_initialized || !m_main.get())
		return false;
	if (!m_main->start())
		return false;
	if (!m_main->setOptions(options))
		return false;

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
	if (!m_main->compileHeuristicBeforeBuild())
		return false;
	return true;
}

bool DaooptInterface::preprocessNoSLS(const ProgramOptions& options, int & ErrorCode) {
	ErrorCode = 0 ;
	if (!m_initialized || !m_main.get())
		{ ErrorCode = 1 ; return false; }
	if (!m_main->start())
		{ ErrorCode = 2 ; return false; }
	if (!m_main->setOptions(options))
		{ ErrorCode = 3 ; return false; }
	if (!m_main->outputInfo())
		{ ErrorCode = 4 ; return false; }
	if (!m_main->loadProblem())
		{ ErrorCode = 5 ; return false; }
	if (!m_main->preprocessHeuristic())
		{ ErrorCode = 6 ; return false; }
	if (!m_main->findOrLoadOrdering())
		{ ErrorCode = 7 ; return false; }
	if (!m_main->initDataStructs())
		{ ErrorCode = 8 ; return false; }
	//if (!m_main->compileHeuristicBeforeBuild())
	//	return false;
	return true;
}

bool DaooptInterface::execSLS(int nSlsIter, int slsTimePerIter) {
	m_main.get()->setSLSOptions(nSlsIter, slsTimePerIter);
	if (!m_main->runSLS())
		return false;
	return true;
}

bool DaooptInterface::exportProblemStats(int& fN, int& fF, int& fK, int& fS, int& fW, int& fH) {
	m_main->exportProblemStats(fN, fF, fK, fS, fW, fH);
	return true;
}

bool DaooptInterface::compileFGLP(int mplp, int mplps) {
	//if (!m_main->compileHeuristicBeforeBuild())
	//	return false;
	m_main.get()->setFGLPOptions(mplp, mplps);
	if (!m_main->compileHeuristic(0))
		return false;
	m_main.get()->setFGLPOptions(0, 0);
	return true;
}

bool DaooptInterface::compileJGLP(int jglp, int jglps, int ibound) {
	m_main.get()->setIboundOptions(ibound); // need to pass BB ibound to daoopt m_options ibound
	m_main.get()->setJGLPOptions(jglp, jglps);

	if (!m_main->compileHeuristicBeforeBuild())
		return false;
	if (!m_main->compileHeuristic(1))
		return false;
	m_main.get()->setJGLPOptions(0, 0);
	return true;
}

bool DaooptInterface::compileMbeAndFinishPreprocessing() {
	if (!m_main->compileHeuristic(2))
		return false;
	if (!m_main->compileHeuristicAfterBuild())
		return false;

	if (!m_main->runLDS())
		return false;
	if (!m_main->finishPreproc())
		return false;
	m_preprocessed = true;
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
#ifndef NO_ASSIGNMENT
  if (assignment != NULL) {
//    const vector<val_t>& assg = m_main->getSolutionAssg();
    vector<val_t> assg ;
	m_main->getSolutionAssgOrg(assg);
    assignment->resize(0);  // Kill existing contents
    for (vector<val_t>::const_iterator it = assg.begin(); it != assg.end(); ++it) {
      assignment->push_back((int) *it);
    }
  }
#endif
  return cost;
}

int DaooptInterface::outputStatistics(void)
{
	return m_main->outputStatsToFile() ;
}

int64 DaooptInterface::getNumORnodesExpanded(void)
{
	return m_main->getNumORnodesExpanded() ;
}

}  // namespace daoopt
