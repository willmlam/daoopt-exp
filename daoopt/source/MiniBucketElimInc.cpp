#include "MiniBucketElimInc.h"

namespace daoopt {

	size_t MiniBucketElimInc::build(int task, const vector<val_t>* assignment, bool computeTables) {

#ifdef DEBUG
		cout << "$ Building MBE(" << m_ibound << ")" << endl;
#endif

		this->reset();
		// keep track of total memory consumption for MBEMM
		size_t memSize = 0;

		switch (task) {
			case 0:
			{
				/*
				FGLP
				*/
				if (computeTables && (m_options->mplp > 0 || m_options->mplps > 0)) {
					cout << "Running FGLP" << endl;
					DoFGLP();
					m_pseudotree->resetFunctionInfo(m_problem->getFunctions());
				}
				break;
			}
			case 1:
			{
				/*
				JGLP
				*/

				int max_ibound_mbemm = -1;
				if (computeTables && (m_options->jglp > 0 || m_options->jglps > 0)) {
					if (m_options->jglpi > 0 && m_options->memlimit != NONE) {
						LimitJGLPIBound(m_options->memlimit, NULL); // set m_options->jglpi to be max. available we pass this since jglpi = 0 we don't set this parameter
						max_ibound_mbemm = m_options->jglpi;
						m_options->jglpi /= 2;
						cout << "Adjusted JGLP i-bound to maximum i-bound / 2" << endl;
					}
					else {
						max_ibound_mbemm = m_ibound;
						m_options->jglpi = m_ibound / 2; // we set m_ibound
						cout << "Setting JGLP i-bound to the current i-bound / 2" << endl;
					}

					cout << "Running JGLP with incremental i-bounds " << endl;
					cout << "Max. initial i-bound " << m_options->jglpi << endl;

					/*m_jglp_converged = false;
					m_jglps_inc = m_options->jglps; // get global time constraint for JGLP
					m_jglp_delta = m_options->jglpt;
					m_jglp_converge_delta = m_options->jglpc;*/

					/*
					d is log(Z1) - log(Z2), Z2=Z1/10^(d) -->
					d ~ ratio: (1e-1, 0.79432823472),  (1e-2, 0.97723722095), (1e-3, 0.99770006382), (1e-4, 0.99976976799), (1e-10, 0.99999999977)
					when d = 1e-4, the next Z value is about 99.977 % of the previous one, not a big decrease!
					differece will be very small value because Z1, Z2 is also very small, diff = Z2*(1-ratio) << small
					*/

					bool ranJglp = false; // test if we ran JGLP at least once
					int jglpi_max = m_options->jglpi; // initial maximum
					int i_iter = m_problem->getR(); // start from max function arity, e.g. 3 or 4
					while (true)
					{
						m_options->jglpi = i_iter;
						cout << "Running JGLP with i-bound " << i_iter << endl;
						ranJglp = DoJGLP(); // return true if cost shifted
						m_pseudotree->resetFunctionInfo(m_problem->getFunctions()); // inefficiency [copy update copy..]
						if ((i_iter + 2 > jglpi_max) || (m_jglp_converged) || (m_jglps_inc <= 0.0)) {
							break;
						}
						else
							i_iter += 2;
					}
					m_options->ibound = max_ibound_mbemm;
					//cout << "Readjusting minibucket i-bound" << endl;
					limitSize(m_options->memlimit, NULL); // this will give Max ibound available
					cout << "Finish running JGLP with incremental i-bounds" << endl;
				}
				break;
			}
			case 2:
			{
					/*
					MBE-MM
					*/
					vector<int> elimOrder; // will hold dfs order
					findDfsOrder(elimOrder); // computes dfs ordering of relevant subtree

					m_augmented.resize(m_problem->getN());
					m_intermediate.resize(m_problem->getN());
					_MiniBuckets.resize(m_problem->getN());

					// ITERATES OVER BUCKETS, FROM LEAVES TO ROOT
					for (vector<int>::reverse_iterator itV = elimOrder.rbegin(); itV != elimOrder.rend(); ++itV) {
						int v = *itV;  // this is the variable being eliminated
						// partition functions into minibuckets
						vector<MiniBucket>& minibuckets = _MiniBuckets[v];
						minibuckets.clear();

#ifdef DEBUG
						cout << "$ Bucket for variable " << *itV << endl;
#endif

						// collect relevant functions in funs
						vector<Function*> funs;
						const vector<Function*>& fnlist = m_pseudotree->getFunctions(*itV);
						funs.insert(funs.end(), fnlist.begin(), fnlist.end());
						funs.insert(funs.end(), m_augmented[*itV].begin(), m_augmented[*itV].end());
#ifdef DEBUG
						for (vector<Function*>::iterator itF = funs.begin(); itF != funs.end(); ++itF)
							cout << ' ' << (**itF);
						cout << endl;
#endif

						// compute global upper bound for root (dummy) bucket
						if (*itV == elimOrder[0]) {// variable is dummy root variable
							if (computeTables && assignment) { // compute upper bound if assignment is given
								m_globalUB = ELEM_ONE;
								for (vector<Function*>::iterator itF = funs.begin(); itF != funs.end(); ++itF)
									m_globalUB OP_TIMESEQ(*itF)->getValue(*assignment);
								cout << "    MBE-ALL  = " << SCALE_LOG(m_globalUB) << " (" << SCALE_NORM(m_globalUB) << ")" << endl;
								m_globalUB OP_DIVIDEEQ m_problem->globalConstInfo();  // for backwards compatibility of output
								cout << "    MBE-ROOT = " << SCALE_LOG(m_globalUB) << " (" << SCALE_NORM(m_globalUB) << ")" << endl;
							}
							continue; // skip the dummy variable's bucket
						}

						// sort functions by decreasing scope size
						sort(funs.begin(), funs.end(), scopeIsLarger);

						// partition functions into minibuckets
						//    vector<Function*>::iterator itF; bool placed;
						for (vector<Function*>::iterator itF = funs.begin(); itF != funs.end();
							++itF) {
							bool placed = false;
							for (vector<MiniBucket>::iterator itB = minibuckets.begin();
								!placed && itB != minibuckets.end(); ++itB)
							{
								if (itB->allowsFunction(*itF)) { // checks if function fits into bucket
									itB->addFunction(*itF);
									placed = true;
								}
							}
							if (!placed) { // no fit, need to create new bucket
								MiniBucket mb(*itV, m_ibound, m_problem);
								mb.addFunction(*itF);
								minibuckets.push_back(mb);
							}
						}
						// Moment-matching step, performed only if we have partitioning.
						vector<Function*> max_marginals;
						std::unique_ptr<Function> average_mm_function;
						if (computeTables && m_options->match && minibuckets.size() > 1) {

							set<int> scope_intersection;
							bool first_mini_bucket = true;
							for (const MiniBucket& mini_bucket : minibuckets) {
								if (first_mini_bucket) {
									scope_intersection = mini_bucket.getJointScope();
									first_mini_bucket = false;
								}
								else {
									scope_intersection =
										intersection(scope_intersection, mini_bucket.getJointScope());
								}
							}
							for (MiniBucket& mini_bucket : minibuckets) {
								set<int> elim_vars =
									setminus(mini_bucket.getJointScope(), scope_intersection);
								max_marginals.push_back(
									mini_bucket.eliminate(computeTables, elim_vars));
							}

							// Find average max-marginals (geometric mean)
							size_t table_size = 1;
							for (const int var : scope_intersection) {
								table_size *= m_problem->getDomainSize(var);
							}

							double* average_mm_table = new double[table_size];
							for (size_t i = 0; i < table_size; ++i) {
								average_mm_table[i] = ELEM_ONE;
							}
							for (const Function* max_marginal : max_marginals) {
								for (size_t i = 0; i < table_size; ++i) {
									average_mm_table[i] OP_TIMESEQ max_marginal->getTable()[i];
								}
							}
							for (size_t i = 0; i < table_size; ++i) {
								average_mm_table[i] = OP_ROOT(average_mm_table[i], minibuckets.size());
							}
							int dummy_id = 0;
							average_mm_function.reset(
								new FunctionBayes(dummy_id, m_problem, scope_intersection,
								average_mm_table, table_size));
						}

						// minibuckets for current bucket are now ready, process each
						// and place resulting function
						int bucket_idx = 0;
						for (MiniBucket& mini_bucket : minibuckets) {
							Function* new_function;
							if (!computeTables || !m_options->match || minibuckets.size() <= 1) {
								new_function = mini_bucket.eliminate(computeTables);
							}
							else {
								new_function = mini_bucket.eliminateMM(computeTables,
									max_marginals[bucket_idx++],
									average_mm_function.get());
							}

							const set<int>& new_scope = new_function->getScopeSet();
							memSize += new_function->getTableSize();
							// go up in tree to find target bucket
							PseudotreeNode* n = m_pseudotree->getNode(*itV)->getParent();
							while (new_scope.find(n->getVar()) == new_scope.end() &&
								n != m_pseudotree->getRoot()) {
								m_intermediate[n->getVar()].push_back(new_function);
								n = n->getParent();
							}
							// matching bucket found OR root of pseudo tree reached
							m_augmented[n->getVar()].push_back(new_function);
						}
						// all minibuckets processed and resulting functions placed
					}

#ifdef DEBUG
/*					// output augmented and intermediate buckets
					if (computeTables)
					for (int i = 0; i<m_problem->getN(); ++i) {
						cout << "$ AUG" << i << ": " << m_augmented[i] << " + " << m_intermediate[i] << endl;
					}*/
#endif

					// clean up for estimation mode
					if (!computeTables) {
						this->reset();
						/*    for (vector<vector<Function*> >::iterator itA = m_augmented.begin(); itA!=m_augmented.end(); ++itA)
						for (vector<Function*>::iterator itB = itA->begin(); itB!=itA->end(); ++itB)
						delete *itB;
						m_augmented.clear();
						m_intermediate.clear();*/
					}

			}
			default:
				break;
		}
		return memSize;
	}



	bool MiniBucketElimInc::DoJGLP() {
		assert(m_pseudotree);
		bool changed_functions = false;
		//DEBUG
		//m_options->jglp = 2; // DBG just see how incremental i-bound works

		//if (m_options && (m_options->jglp > 0 || m_options->jglps >0))  {
		//if (m_options && (m_options->jglp > 0 || jglps_inc > 0))  {
		if (m_options && (m_jglps_inc > 0.0) && (!m_jglp_converged))  {
			//cout << "JGLP inc @1 :: " << jglps_inc << endl;
			mex::mbe _jglp(CopyFactors());
			mex::VarOrder var_order(m_pseudotree->getElimOrder().begin(),
				--m_pseudotree->getElimOrder().end());
			_jglp.setOrder(var_order);

			mex::VarOrder parents(m_problem->getN() - 1); // copy pseudotree information
			for (int i = 0; i < m_problem->getN() - 1; ++i) {
				int parent_var = m_pseudotree->getNode(i)->getParent()->getVar();
				parents[i] = parent_var == m_pseudotree->getRoot()->getVar()
					? -1 : parent_var;
			}
			_jglp.setPseudotree(parents);
			_jglp.setIBound(m_options->jglpi);
			_jglp.setProperties("DoMatch=1,DoFill=1,DoJG=1,DoMplp=0");
			_jglp.init();
			
			//_jglp.tighten(m_options->jglp > 0 ? m_options->jglp : 100, m_options->jglps);
			time_t before_tighten, after_tighten;
			time(&before_tighten);
// 2015-04 changes to make it compile
//			m_jglp_converged = _jglp.tighten(m_jglp_converge_delta, m_options->jglp > 0 ? m_options->jglp : 100, m_jglps_inc, m_jglp_delta); // get dObj directly in future?
			_jglp.tighten(m_options->jglp > 0 ? m_options->jglp : 100, m_jglps_inc, m_jglp_delta); // get dObj directly in future?
			time(&after_tighten);
			m_jglps_inc = m_jglps_inc - difftime(after_tighten, before_tighten);
			//cout << "JGLP inc @2 :: " << jglps_inc << endl;
			RewriteFactors(_jglp.factors(), _jglp.logZ());
			changed_functions = true;
		}

		// reduce jglps
		return changed_functions;
	}

}
