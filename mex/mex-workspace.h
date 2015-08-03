// Wrapper for UAI competition code

#ifndef MEX_WORKSPACE_HXX_INCLUDED
#define MEX_WORKSPACE_HXX_INCLUDED

#include <cstdio>
#include <iostream>
#include <sstream>
#include <fstream>

#include "factorgraph.h"

#include "mplp.h"
#include "mbe.h"

#include "lbp.h"
#include "gbp.h"
#include "enum.h"

namespace mex
{

MEX_ENUM(Task, PR, MAR, MPE) ;

class Workspace
{

public :

	// parameters
	double _MemLimit ;

	mex::graphModel *_gmo ;
	mex::vector<Factor> _bel ;
	mex::vector<size_t> _evid ;
	VarSet _evVar ;
	Task _task ;

	std::vector<int> _VarElimOrder ;
	int _InducedWidth ;

	double _lnZtot ;

public :

	int SetTask(const char *task) 
	{
		try 
			{ _task = Task(task) ; }
		catch (...) 
			{ return 1 ; }
		return 0 ;
	}

	int LoadUAI(const char *buf, int L)
	{
		std::istringstream is(buf, L) ;
		mex::vector<Factor> forig = Factor::readUai10(is) ;
		if (0 == forig.size()) 
			return 1 ;
		if (NULL != _gmo) 
			{ delete _gmo ; _gmo = NULL ; }
		_gmo = new mex::graphModel(forig) ;
		if (NULL == _gmo) 
			return 1 ;
		size_t nvar = _gmo->nvar() ;
		_bel.resize(nvar) ;
		_evid.resize(nvar, 0) ;
		return 0 ;
	}

	int LoadUAIEvidence(const char *buf, int L)
	{
		if (NULL == _gmo) 
			return 1 ;
		size_t nvar = _gmo->nvar() ;

		std::istringstream is2(buf, L) ;
		int nEvidVar; is2 >> nEvidVar ;
		for (size_t i = 0 ; i < nEvidVar ; i++) {
			uint32_t vid ; size_t vval ; is2 >> vid >> vval ; 
			_evid[vid] = vval ;
			_evVar |= _gmo->var(vid) ;
			_bel[vid] = Factor::delta(_gmo->var(vid), vval) ;
			}
//		std::cout<<"Evidence on variables "<<evVar<<"\n";
		_gmo->condition(_evVar,_evid) ;

		return 0 ;
	}

	int SetVarElimOrdering(std::vector<int> & Order, int induced_width)
	{
		if (Order.size() != _gmo->nvar()) 
			return 1 ;
		_VarElimOrder = Order ;
		_InducedWidth = induced_width ;
		return 0 ;
	}

	int Solve(void) ;

public :

	void Delete(void)
	{
		if (NULL != _gmo) 
			{ delete _gmo ; _gmo = NULL ; }
		_bel.clear() ;
		_evid.clear() ;
		_evVar.clear() ;

		_VarElimOrder.clear() ;
		_InducedWidth = -1 ;

		_lnZtot = DBL_MAX ;
	}

	Workspace(void) 
		:
		_MemLimit(1073741824), // 1GB
		_gmo(NULL), 
		_InducedWidth(-1), 
		_lnZtot(DBL_MAX)
	{
	}

	virtual ~Workspace(void)
	{
		Delete() ;
	}

} ;

} // namespace mex

#endif // MEX_WORKSPACE_HXX_INCLUDED
