#ifndef MINIBUCKETELIMINC_H_
#define MINIBUCKETELIMINC_H_

#include "Heuristic.h"
#include "Function.h"
#include "Problem.h"
#include "ProgramOptions.h"
#include "Pseudotree.h"
#include "utils.h"

#include "MiniBucket.h"

#include "mex/include/mbe.h"
#include "FGLP.h"
#include "PriorityFGLP.h"

#include "MiniBucketElim.h"

namespace daoopt {

	class MiniBucketElimInc : public MiniBucketElim {
		friend class MiniBucket;

	protected: // additional members to control incremental behavior (fglp runs incrementally)
		bool m_jglp_converged; // test convergence for incremental JGLP
		double m_jglps_inc;  // m_options->jglps is global time; jglps_inc -= time_elapsed while JGLP(inc)
		double m_jglp_delta; // tolerance parameter for escaping single JGLP run
		double m_jglp_converge_delta; // converge paramter for escaping incremental JGLP run
		bool m_fglp_converged; // same purpose as jglp_converged

	public: // methods to override same signature with slightly different implementation

		virtual size_t build(int task, const vector<val_t>* assignment = NULL, bool computeTables = true);
		
		bool DoJGLP();

		MiniBucketElimInc(Problem* p, Pseudotree* pt, ProgramOptions* po, int ib);
		virtual ~MiniBucketElimInc();
	};

	inline MiniBucketElimInc::MiniBucketElimInc(Problem* p, Pseudotree* pt, ProgramOptions* po, int ib)
		: MiniBucketElim(p, pt, po, ib), m_jglp_converged(false), m_fglp_converged(false), 
		m_jglps_inc(po->jglps), m_jglp_delta(po->jglpt), m_jglp_converge_delta(po->jglpc) {}

	inline MiniBucketElimInc::~MiniBucketElimInc() { reset(); }


}
#endif /* MINIBUCKETELIMINC_H_*/
