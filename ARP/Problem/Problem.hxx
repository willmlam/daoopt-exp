#ifndef ARP_HXX_INCLUDED
#define ARP_HXX_INCLUDED

#include <stdlib.h>
#include <string>

#include "Function.hxx"
//#include "ProblemGraphNode.hxx"
#include <vector>
using namespace std ;

#if defined LINUX
#define stricmp strcasecmp
#endif

#ifndef ARP_nInfinity
#define ARP_nInfinity (-std::numeric_limits<double>::infinity())
#endif // ARP_nInfinity
#ifndef ARP_pInfinity
#define ARP_pInfinity (std::numeric_limits<double>::infinity())
#endif // ARP_pInfinity

namespace ARE // Automated Reasoning Engine
{

class ARP // Automated Reasoning Problem
{
protected :
	std::string _Name ;
public :
	inline const std::string & Name(void) const { return _Name ; }
	inline void SetName(const std::string & Name) { _Name = Name ; }

protected :
	std::string _FileName ;
	std::string _EvidenceFileName ;
public :
	inline const std::string & FileName(void) const { return _FileName ; }
	inline const std::string & EvidenceFileName(void) const { return _EvidenceFileName ; }
	int32_t GetFilename(const std::string & Dir, std::string & fn) ; // filename will not have an extension
	inline void SetFileName(const std::string & FileName) { _FileName = FileName ; }

protected :
	FILE *_fpLOG ;
public :
	inline FILE * & fpLOG(void) { return _fpLOG ; }

protected :
	int32_t _nVars ; // number of variables; indeces run [0,_nVars-1]
	int32_t *_K ; // domain size of each variable; size of this array is _nVars.
	int32_t *_Value ; // current value of each variable; for var i, legal values run [0, _K[i]-1].
	int32_t _nSingletonDomainVariables ; // number of singleton-domain variables in the problem
	ARE_Function_TableType _AssignmentValue ; // value of the current assignment ;
public :
	inline int32_t N(void) const { return _nVars ; }
	inline int32_t K(int32_t IDX) const { return _K[IDX] ; }
	inline const int32_t *K(void) const { return _K ; }
	inline int32_t Value(int32_t IDX) const { return _Value[IDX] ; }
	inline int32_t *ValueArray(void) { return _Value ; }
	inline int32_t SetValue(int32_t IDX, int32_t v) { _Value[IDX] = v ; }
	inline int32_t nSingletonDomainVariables(void) const { return _nSingletonDomainVariables ; }
	inline ARE_Function_TableType & AssignmentValue(void) { return _AssignmentValue ; }

	int32_t ComputeSumOfDomainSizes(bool IgnoreSingletons)
	{
		int32_t sum = 0 ;
		for (int32_t i = _nVars-1 ; i >= 0 ; i--) {
			if (IgnoreSingletons && _K[i] < 2) continue ;
			if (_K[i] > 0) 
				sum += _K[i] ;
			}
		return sum ;
	}

	inline int32_t SetN(int32_t N)
	{
		if (_nVars < 0) 
			return 1 ;
		if (_nVars == N) 
			return 0 ;
		_nVars = N ;
		if (NULL != _K) { delete [] _K ; _K = NULL ; }
		if (NULL != _Value) { delete [] _Value ; _Value = NULL ; }
		if (_nVars > 0) {
			_K = new int32_t[_nVars] ;
			_Value = new int32_t[_nVars] ;
			if (NULL == _K || NULL == _Value) {
				Destroy() ;
				return 1 ;
				}
			for (int32_t i = 0 ; i < _nVars ; i++) {
				_K[i] = -1 ;
				_Value[i] = -1 ;
				}
			}
		return 0 ;
	}
	inline int32_t SetK(int32_t *K)
	{
		if (NULL == _K) 
			return 1 ;
		for (int32_t i = 0 ; i < _nVars ; i++) 
			_K[i] = K[i] ;
		return 0 ;
	}
	inline int32_t SetK(int32_t K)
	{
		if (NULL == _K || K < 0)
			return 1;
		for (int32_t i = 0; i < _nVars; i++)
			_K[i] = K;
		return 0;
	}

	// ***************************************************************************************************
	// Functions of the problem
	// ***************************************************************************************************

protected :

	bool _FunctionsAreConvertedToLogScale ;
	int32_t _nFunctions ; // number of functions in the  problem.
	Function **_Functions ; // a list of  functions.
	int64_t _FunctionsSpace ; // space, in bytes, taken up by all functions.
	int32_t _nFunctionsIrrelevant ; // number of functions in the problem not relevent to the query.

public :

	inline bool FunctionsAreConvertedToLogScale(void) const { return _FunctionsAreConvertedToLogScale ; }
	inline int32_t nFunctions(void) const { return _nFunctions ; }
	inline ARE::Function *getFunction(int32_t IDX) const { return _Functions[IDX] ; }
	inline int64_t FunctionsSpace(void) const { return _FunctionsSpace ; }
	inline int32_t nFunctionsIrrelevant(void) const { return _nFunctionsIrrelevant ; }

	// in bytes, the space required
	int64_t ComputeFunctionSpace(void) ;

	int32_t ConvertFunctionsToLogScale(void) ;
	
	// return non-0 iff something is wrong with any of the functions
	int32_t CheckFunctions(void) ;

	int32_t DeleteDuplicateFunctions(void) ;

	// find Bayesian CPT for the given variable
	inline ARE::Function *GetCPT(int32_t V)
	{
		for (int32_t i = 0 ; i < _nFunctions ; i++) {
			ARE::Function *f = _Functions[i] ;
			if (NULL == f) continue ;
			if (V != f->BayesianCPTChildVariable()) 
				continue ;
			return f ;
			}
		return NULL ;
	}

	// given a variable, it will compute all ancestors, assuming functions are Bayesian
	int32_t ComputeBayesianAncestors(int32_t V, int32_t *AncestorFlags, int32_t *Workspace) ;

	// this function checks if all functions have tables; if not, it will fill in the table, checking the type. 
	// e.g. if the function is a Bayesian CPT, it will fill it in appropriately.
	int32_t FillInFunctionTables(void) ;

	double ComputeAvgFn0Tightness(void)
	{
		double n = 0.0, sum = 0.0 ;
		for (int32_t i = 0 ; i < _nFunctions ; i++) {
			ARE::Function *f = _Functions[i] ;
			if (NULL == f) continue ;
			double zeroT = f->Compute0Tightness() ;
			if (zeroT >= 0.0) 
				{ n += 1.0 ; sum += zeroT ; }
			}
		return n > 0.0 ? sum/n : -1.0 ;
	}

	// ***************************************************************************************************
	// for each variable, functions adjacent to it
	// ***************************************************************************************************

protected :
	Function **_StaticAdjFnTotalList ; // helper list for maitaining adj fn lists for each var; we allocate a static list enough for all vars.
	int32_t _StaticAdjFnTotalListSize ; // _StaticAdjFnTotalList allocated length.
	int32_t *_nAdjFunctions ; // for each variable, the number of functions it participates in
	int32_t *_AdjFunctions ; // for each variable, idx into _StaticAdjFnTotalList[] array where Fn ptrs list (FNs that this var participates in) for that var start; length of list is _nAdjFunctions[].
public :
	inline int32_t nAdjFunctions(int32_t idxVar) { return _nAdjFunctions[idxVar] ; }
	inline Function *AdjFunction(int32_t idxVar, int32_t idxFn) { return _StaticAdjFnTotalList[_AdjFunctions[idxVar] + idxFn] ; }
	inline int32_t nAdjacentFunctions_QueryReleventFunctionsOnly(int32_t idxVar) 
	{
		int32_t n = 0 ;
		for (int32_t i = 0 ; i < _nAdjFunctions[idxVar] ; i++) {
			Function *f = AdjFunction(idxVar, i) ;
			if (f->IsQueryIrrelevant()) 
				continue ;
			++n ;
			}
		return n ;
	}
	inline Function *AdjacentFunction_QueryRelevantFunctionsOnly(int32_t idxVar, int32_t idxFn) 
	{
		int32_t n = 0 ;
		for (int32_t i = 0 ; i < _nAdjFunctions[idxVar] ; i++) {
			Function *f = AdjFunction(idxVar, i) ;
			if (f->IsQueryIrrelevant()) 
				continue ;
			if (n == idxFn) 
				return f ;
			++n ;
			}
		return NULL ;
	}
	int32_t DestroyAdjFnList(void) ;
	int32_t ComputeAdjFnList(bool IgnoreIrrelevantFunctions) ;

	// ***************************************************************************************************
	// problem graph
	// ***************************************************************************************************

protected :
	int32_t *_StaticVarTotalList ; // helper list for maitaining adj var for each var; we allocate a static list enough for all vars.
	int32_t _StaticVarTotalListSize ; // _StaticVarTotalList allocated length.
	int32_t *_Degree ; // for each variable, its degree in the graph
	int32_t *_AdjVars ; // for each variable, idx into _StaticVarTotalList[] array where AdjVar list (vars that this var is adj to) for that var start; length of list is _Degree[].
	int32_t _nSingletonVariables ;
public :
	inline int32_t Degree(int32_t idxVar) { return _Degree[idxVar] ; }
	inline int32_t AdjVar(int32_t idxVar, int32_t idxAdjVar) { return _StaticVarTotalList[_AdjVars[idxVar] + idxAdjVar] ; }
	int32_t DestroyAdjVarList(void) ;
	int32_t ComputeAdjVarList(void) ; // when computing this, we assume AdjFn[] is computed.

	// assuming variable V is adjacent to U, this fn will remove V from the list of variables adjacent to U.
	int32_t RemoveAdjVar(int32_t U, int32_t V) ;

	// ***************************************************************************************************
	// connected components of the problem (graph)
	// ***************************************************************************************************

protected :
	int32_t _nConnectedComponents ;
public :
	inline int32_t nConnectedComponents(void) const { return _nConnectedComponents ; }
	int32_t ComputeConnectedComponents(void) ; // when computing this, we assume _AdjVars[] is computed.

	// ***************************************************************************************************
	// variable ordering
	// ***************************************************************************************************

protected :
	int32_t _VarOrdering_InducedWidth ; // width of the ordering.
	int32_t *_VarOrdering_VarList ; // a list of variables, in the order, [0] is root of bucket tree, for any i, descendants of i must be at indeces >i.
	int32_t *_VarOrdering_VarPos ; // for each variable, its position in the order.
public :
	inline int32_t VarOrdering_InducedWidth(void) const { return _VarOrdering_InducedWidth ; }
	inline const int32_t *VarOrdering_VarList(void) const { return _VarOrdering_VarList ; }
	inline const int32_t *VarOrdering_VarPos(void) const { return _VarOrdering_VarPos ; }
	void DestroyVarOrdering(void) ;
	inline bool HasVarOrdering(void) const 
		// 2015-05-18 KK : if variable order is provided by the user, we may not know the induced width.
		{ return /* _VarOrdering_InducedWidth >= 0 && _VarOrdering_InducedWidth < INT_MAX && */ NULL != _VarOrdering_VarList && NULL != _VarOrdering_VarPos ; }

	// get variable order in various formats
	int32_t GetVarElimOrdering(std::vector<int32_t> & Order, int32_t & induced_width) ;

	// set variable elimination order; we assume the order in the input array is elimination order, i.e. var[0] is the first var to be eliminated, etc.
	int32_t SetVarElimOrdering(const int32_t *VarListInElimOrder, int32_t induced_width) ;
	int32_t SetVarBTOrdering(const int32_t *VarListInBTOrder, int32_t induced_width) ;

	// load variable order from a string buffer; we assume it is in the format "{var0;var1;...;varN}".
	// 0=OrderType means variable elimination order.
	// 1=OrderType means variable bucket-tree order.
	int32_t LoadVariableOrderingFromBuffer(int32_t OrderType, const char *SerializedVarListInElimOrder) ;

	// in the following, VarList is order such that [0] is root of the bucket tree, and for any [i], descendants of [i] must be at indeces >i.
	int32_t ComputeInducedWidth(const int32_t *VarList, const int32_t *Var2PosMap, int32_t & InducedWidth) ;
	int32_t TestVariableOrdering(const int32_t *VarList, const int32_t *Var2PosMap) ;

	// ***************************************************************************************************
	// query
	// ***************************************************************************************************

protected :

	int32_t _QueryVariable ; // if this is set, we want to compute a distribution on this var
	int32_t _FnCombinationType ; // 0 = undef, 1 = product, 2 = sum
	int32_t _VarElimType ; // 0 = undef, 1 = sum, 2 = max, 3 = min

public :

	inline int32_t FnCombinationType(void) { return _FnCombinationType ; }
	inline int32_t VarEliminationType(void) { return _VarElimType ; }
	inline int32_t QueryVariable(void) { return _QueryVariable ; }
	inline bool IsOptimizationProblem(void) const { return VAR_ELIMINATION_TYPE_MAX == _VarElimType || VAR_ELIMINATION_TYPE_MIN == _VarElimType ; }

	inline void ApplyFnCombinationOperator(ARE_Function_TableType & V, ARE_Function_TableType v)
	{
		if (FN_COBINATION_TYPE_PROD == _FnCombinationType) {
			if (_FunctionsAreConvertedToLogScale) 
				V += v ;
			else 
				V *= v ;
			}
		else if (FN_COBINATION_TYPE_SUM == _FnCombinationType) {
			if (_FunctionsAreConvertedToLogScale) 
				LOG_OF_SUM_OF_TWO_NUMBERS_GIVEN_AS_LOGS(V, V, v)
			else 
				V += v ;
			}
	}

	inline void SetQueryVariable(int32_t v) { _QueryVariable = v ; }
	inline void SetOperators(int32_t FnCombinationType, int32_t VarElimType) { _FnCombinationType = FnCombinationType ; _VarElimType = VarElimType ; }
	int32_t SetOperators(const char *query)
	{
		if (0 == stricmp("product-sum", query)) 
			{ _FnCombinationType = FN_COBINATION_TYPE_PROD ; _VarElimType = VAR_ELIMINATION_TYPE_SUM ; return 0 ; }
		else if (0 == stricmp("product-max", query)) 
			{ _FnCombinationType = FN_COBINATION_TYPE_PROD ; _VarElimType = VAR_ELIMINATION_TYPE_MAX ; return 0 ; }
		else if (0 == stricmp("product-min", query)) 
			{ _FnCombinationType = FN_COBINATION_TYPE_PROD ; _VarElimType = VAR_ELIMINATION_TYPE_MIN ; return 0 ; }
		else if (0 == stricmp("sum-sum", query)) 
			{ _FnCombinationType = FN_COBINATION_TYPE_SUM ; _VarElimType = VAR_ELIMINATION_TYPE_SUM ; return 0 ; }
		else if (0 == stricmp("sum-max", query)) 
			{ _FnCombinationType = FN_COBINATION_TYPE_SUM ; _VarElimType = VAR_ELIMINATION_TYPE_MAX ; return 0 ; }
		else if (0 == stricmp("sum-min", query)) 
			{ _FnCombinationType = FN_COBINATION_TYPE_SUM ; _VarElimType = VAR_ELIMINATION_TYPE_MIN ; return 0 ; }
		return 1 ;
	}
	int32_t SetOperators(const char *combop, const char *elimop)
	{
		_FnCombinationType = FN_COBINATION_TYPE_NONE ;
		_VarElimType = VAR_ELIMINATION_TYPE_NONE ;
		if (0 == stricmp("product", combop)) 
			_FnCombinationType = FN_COBINATION_TYPE_PROD ;
		else if (0 == stricmp("sum", combop)) 
			_FnCombinationType = FN_COBINATION_TYPE_SUM ;
		if (0 == stricmp("sum", elimop)) 
			_VarElimType = VAR_ELIMINATION_TYPE_SUM ;
		else if (0 == stricmp("max", elimop)) 
			_VarElimType = VAR_ELIMINATION_TYPE_MAX ;
		else if (0 == stricmp("min", elimop)) 
			_VarElimType = VAR_ELIMINATION_TYPE_MIN ;
		return 0 ;
	}

	// compute which functions are relevant/irrelevant to the query; we process recursively from leaves up, checking variables that are in a single function.
	// VarElimOp=0 means sum, VarElimOp=1 means max.
	// CombinationOp=0 means product, CombinationOp=1 means sum.
	// we assume Factor is initialized.
	int32_t ComputeQueryRelevance_VarElimination(ARE_Function_TableType & Factor, char VarElimOp, char CombinationOp) ;

	// ***************************************************************************************************
	// misc
	// ***************************************************************************************************

public :

	// assuming the problem was just constructed/loaded, perform initial analysis.
	// as input we assume just variables, domains, functions.
	// e.g. this fn does not compute variable order. 
	// this fn should run fast.
	int32_t PerformPostConstructionAnalysis(void) ;

	// aliminate all evidence from functions. we assume evidence is set as _Value[var] = value.
	int32_t EliminateEvidence(void) ;

	// given evidence Var=Val, this function will eliminate the evidence variable from all functions, by assigning Var=Val.
	// If as a result, a fn becomes a const fn, it will add it to 
	int32_t EliminateEvidenceVariable(int32_t Var, int32_t Val) ;

	// eliminate all variables, whose domain has just one value, from all functions.
	// all functions that have no arguments left, as the result, will be turned into const functions.
	int32_t EliminateSingletonDomainVariables(void) ;

	// singleton consistency; for probabilistic network, check each domain value against 0-probabilities.
	// return value  0 means ok.
	// return value -1 means domain of some variable became empty.
	int32_t ComputeSingletonConsistency(int32_t & nNewSingletonDomainVariables) ;
	// The following function returns is_consistent[i][j] where i indexes variables and j indexes the domains of variables
	// is_consistent[i][j] is true if the variable value combination is consistent and false otherwise
	void SingletonConsistencyHelper(vector<vector<bool>> & is_consistent) ;

	// this function generates a random Bayesian network structure, that is it does not generate actual tables, 
	// just the variables (arguments) for each function (CPT).
	// tables can be generated by calling FillInFunctionTables().
	int32_t GenerateRandomUniformBayesianNetworkStructure(int32_t N, int32_t K, int32_t P, int32_t C, int32_t ProblemCharacteristic) ;

public :

	virtual int32_t SaveXML(const std::string & Dir) ;
	virtual int32_t SaveUAI08(const std::string & Dir) ;

	// load from file/buffer.
	int32_t LoadFromFile(const std::string & FileName) ;
	int32_t LoadFromFile_Evidence(const std::string & FileName, int32_t & nEvidenceVars) ;

	// this function loads problem from UAI format file; we assume data is already in buf.
	int32_t LoadFromBuffer(const char *format, const char *buf, int32_t L) ;
	int32_t LoadFromBuffer_Evidence(const char *format, const char *buf, int32_t L, int32_t & nEvidenceVars) ;

	int32_t LoadUAIFormat(const char *buf, int32_t L) ;
	int32_t LoadUAIFormat_Evidence(const char *buf, int32_t L, int32_t & nEvidenceVars) ;

public :

	void Destroy(void)
	{
		DestroyAdjVarList() ;
		DestroyAdjFnList() ;
		DestroyVarOrdering() ;
		if (NULL != _Functions) {
			for (int32_t i = 0 ; i < _nFunctions ; i++) {
				if (NULL != _Functions[i]) 
					delete _Functions[i] ;
				}
			delete [] _Functions ;
			_Functions = NULL ;
			}
		_FunctionsAreConvertedToLogScale = false ;
		_nFunctions = 0 ;
		_nFunctionsIrrelevant = -1 ;
		if (NULL != _K) {
			delete [] _K ;
			_K = NULL ;
			}
		if (NULL != _Value) {
			delete [] _Value ;
			_Value = NULL ;
			}
		_nVars = 0 ;
		_nSingletonDomainVariables = -1 ;
		_nConnectedComponents = -1 ;
		_FileName.erase() ;
		_QueryVariable = -1 ;
	}
	ARP(const char *Name = NULL)
		:
		_fpLOG(NULL), 
		_nVars(0), 
		_K(NULL), 
		_Value(NULL), 
		_nSingletonDomainVariables(-1), 
		_FunctionsAreConvertedToLogScale(false), 
		_nFunctions(0),
		_Functions(NULL), 
		_FunctionsSpace(0), 
		_nFunctionsIrrelevant(-1), 
		_StaticAdjFnTotalList(NULL), 
		_StaticAdjFnTotalListSize(0), 
		_nAdjFunctions(NULL), 
		_AdjFunctions(NULL), 
		_StaticVarTotalList(NULL), 
		_StaticVarTotalListSize(0), 
		_Degree(NULL), 
		_AdjVars(NULL), 
		_nSingletonVariables(-1), 
		_nConnectedComponents(-1), 
		_VarOrdering_InducedWidth(-1), 
		_VarOrdering_VarList(NULL), 
		_VarOrdering_VarPos(NULL),
		_QueryVariable(-1), 
		_FnCombinationType(0), 
		_VarElimType(0)
	{
		if (NULL != Name) 
			_Name = Name ;
	}
	virtual ~ARP(void)
	{
		Destroy() ;
	}
} ;

int32_t fileload_getnexttoken(const char * & buf, int32_t & L, const char * & B, int32_t & l, bool IncludeSpecialSymbols) ;

} // namespace ARE

#endif // ARP_HXX_INCLUDED
