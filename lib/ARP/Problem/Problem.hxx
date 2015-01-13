#ifndef ARP_HXX_INCLUDED
#define ARP_HXX_INCLUDED

#include <stdlib.h>
#include <string>

#include <Function.hxx>
//#include <ProblemGraphNode.hxx>
#include <vector>
using namespace std ;

#if defined LINUX
#define stricmp strcasecmp
#endif

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
	int GetFilename(const std::string & Dir, std::string & fn) ; // filename will not have an extension
	inline void SetFileName(const std::string & FileName) { _FileName = FileName ; }

protected :
	FILE *_fpLOG ;
public :
	inline FILE * & fpLOG(void) { return _fpLOG ; }

protected :
	int _nVars ; // number of variables; indeces run [0,_nVars-1]
	int *_K ; // domain size of each variable; size of this array is _nVars.
	int *_Value ; // current value of each variable; for var i, legal values run [0, _K[i]-1].
	int _nSingletonDomainVariables ; // number of singleton-domain variables in the problem
public :
	inline int N(void) const { return _nVars ; }
	inline int K(int IDX) const { return _K[IDX] ; }
	inline const int *K(void) const { return _K ; }
	inline int Value(int IDX) const { return _Value[IDX] ; }
	inline int nSingletonDomainVariables(void) const { return _nSingletonDomainVariables ; }

	inline int SetN(int N)
	{
		if (_nVars < 0) 
			return 1 ;
		if (_nVars == N) 
			return 0 ;
		_nVars = N ;
		if (NULL != _K) { delete [] _K ; _K = NULL ; }
		if (NULL != _Value) { delete [] _Value ; _Value = NULL ; }
		if (_nVars > 0) {
			_K = new int[_nVars] ;
			_Value = new int[_nVars] ;
			if (NULL == _K || NULL == _Value) {
				Destroy() ;
				return 1 ;
				}
			for (int i = 0 ; i < _nVars ; i++) {
				_K[i] = -1 ;
				_Value[i] = -1 ;
				}
			}
		return 0 ;
	}
	inline int SetK(int *K)
	{
		if (NULL == _K) 
			return 1 ;
		for (int i = 0 ; i < _nVars ; i++) 
			_K[i] = K[i] ;
		return 0 ;
	}

	// ***************************************************************************************************
	// Functions of the problem
	// ***************************************************************************************************

protected :

	bool _FunctionsAreConvertedToLogScale ;
	int _nFunctions ; // number of functions in the  problem.
	Function **_Functions ; // a list of  functions.
	__int64 _FunctionsSpace ; // space, in bytes, taken up by all functions.
	int _nFunctionsIrrelevant ; // number of functions in the problem not relevent to the query.

public :

	inline bool FunctionsAreConvertedToLogScale(void) const { return _FunctionsAreConvertedToLogScale ; }
	inline int nFunctions(void) const { return _nFunctions ; }
	inline ARE::Function *getFunction(int IDX) const { return _Functions[IDX] ; }
	inline __int64 FunctionsSpace(void) const { return _FunctionsSpace ; }
	inline int nFunctionsIrrelevant(void) const { return _nFunctionsIrrelevant ; }

	// in bytes, the space required
	__int64 ComputeFunctionSpace(void) ;

	int ConvertFunctionsToLogScale(void) ;
	
	// return non-0 iff something is wrong with any of the functions
	int CheckFunctions(void) ;

	// find Bayesian CPT for the given variable
	inline ARE::Function *GetCPT(int V)
	{
		for (int i = 0 ; i < _nFunctions ; i++) {
			ARE::Function *f = _Functions[i] ;
			if (NULL == f) continue ;
			if (V != f->BayesianCPTChildVariable()) 
				continue ;
			return f ;
			}
		return NULL ;
	}

	// given a variable, it will compute all ancestors, assuming functions are Bayesian
	int ComputeBayesianAncestors(int V, int *AncestorFlags, int *Workspace) ;

	// this function checks if all functions have tables; if not, it will fill in the table, checking the type. 
	// e.g. if the function is a Bayesian CPT, it will fill it in appropriately.
	int FillInFunctionTables(void) ;

	// ***************************************************************************************************
	// for each variable, functions adjacent to it
	// ***************************************************************************************************

protected :
	Function **_StaticAdjFnTotalList ; // helper list for maitaining adj fn lists for each var; we allocate a static list enough for all vars.
	int _StaticAdjFnTotalListSize ; // _StaticAdjFnTotalList allocated length.
	int *_nAdjFunctions ; // for each variable, the number of functions it participates in
	int *_AdjFunctions ; // for each variable, idx into _StaticAdjFnTotalList[] array where Fn ptrs list (FNs that this var participates in) for that var start; length of list is _nAdjFunctions[].
public :
	inline int nAdjFunctions(int idxVar) { return _nAdjFunctions[idxVar] ; }
	inline Function *AdjFunction(int idxVar, int idxFn) { return _StaticAdjFnTotalList[_AdjFunctions[idxVar] + idxFn] ; }
	inline int nAdjacentFunctions_QueryReleventFunctionsOnly(int idxVar) 
	{
		int n = 0 ;
		for (int i = 0 ; i < _nAdjFunctions[idxVar] ; i++) {
			Function *f = AdjFunction(idxVar, i) ;
			if (f->IsQueryIrrelevant()) 
				continue ;
			++n ;
			}
		return n ;
	}
	inline Function *AdjacentFunction_QueryRelevantFunctionsOnly(int idxVar, int idxFn) 
	{
		int n = 0 ;
		for (int i = 0 ; i < _nAdjFunctions[idxVar] ; i++) {
			Function *f = AdjFunction(idxVar, i) ;
			if (f->IsQueryIrrelevant()) 
				continue ;
			if (n == idxFn) 
				return f ;
			++n ;
			}
		return NULL ;
	}
	int DestroyAdjFnList(void) ;
	int ComputeAdjFnList(bool IgnoreIrrelevantFunctions) ;

	// ***************************************************************************************************
	// problem graph
	// ***************************************************************************************************

protected :
	int *_StaticVarTotalList ; // helper list for maitaining adj var for each var; we allocate a static list enough for all vars.
	int _StaticVarTotalListSize ; // _StaticVarTotalList allocated length.
	int *_Degree ; // for each variable, its degree in the graph
	int *_AdjVars ; // for each variable, idx into _StaticVarTotalList[] array where AdjVar list (vars that this var is adj to) for that var start; length of list is _Degree[].
	int _nSingletonVariables ;
public :
	inline int Degree(int idxVar) { return _Degree[idxVar] ; }
	inline int AdjVar(int idxVar, int idxAdjVar) { return _StaticVarTotalList[_AdjVars[idxVar] + idxAdjVar] ; }
	int DestroyAdjVarList(void) ;
	int ComputeAdjVarList(void) ; // when computing this, we assume AdjFn[] is computed.

	// assuming variable V is adjacent to U, this fn will remove V from the list of variables adjacent to U.
	int RemoveAdjVar(int U, int V) ;

	// ***************************************************************************************************
	// connected components of the problem (graph)
	// ***************************************************************************************************

protected :
	int _nConnectedComponents ;
public :
	int ComputeConnectedComponents(void) ; // when computing this, we assume _AdjVars[] is computed.

	// ***************************************************************************************************
	// variable ordering
	// ***************************************************************************************************

protected :
	int _VarOrdering_InducedWidth ; // width of the ordering.
	int *_VarOrdering_VarList ; // a list of variables, in the order, [0] is root of bucket tree, for any i, descendants of i must be at indeces >i.
	int *_VarOrdering_VarPos ; // for each variable, its position in the order.
public :
	inline int VarOrdering_InducedWidth(void) const { return _VarOrdering_InducedWidth ; }
	inline const int *VarOrdering_VarList(void) const { return _VarOrdering_VarList ; }
	inline const int *VarOrdering_VarPos(void) const { return _VarOrdering_VarPos ; }
	void DestroyVarOrdering(void) ;
	inline int HasVarOrdering(void) const { return NULL != _VarOrdering_VarList && NULL != _VarOrdering_VarPos ; }

	// get variable order in various formats
	int GetVarElimOrdering(std::vector<int> & Order, int & induced_width) ;

	// set variable elimination order; we assume the order in the input array is elimination order, i.e. var[0] is the first var to be eliminated, etc.
	int SetVarElimOrdering(const int *VarListInElimOrder, int induced_width) ;
	int SetVarBTOrdering(const int *VarListInBTOrder, int induced_width) ;

	// load variable order from a string buffer; we assume it is in the format "{var0;var1;...;varN}".
	// 0=OrderType means variable elimination order.
	// 1=OrderType means variable bucket-tree order.
	int LoadVariableOrderingFromBuffer(int OrderType, const char *SerializedVarListInElimOrder) ;

	// in the following, VarList is order such that [0] is root of the bucket tree, and for any [i], descendants of [i] must be at indeces >i.
	int ComputeInducedWidth(const int *VarList, const int *Var2PosMap, int & InducedWidth) ;
	int TestVariableOrdering(const int *VarList, const int *Var2PosMap) ;

	// ***************************************************************************************************
	// query
	// ***************************************************************************************************

protected :

	int _QueryVariable ; // if this is set, we want to compute a distribution on this var
	int _FnCombinationType ; // 0 = undef, 1 = product, 2 = sum
	int _VarElimType ; // 0 = undef, 1 = sum, 2 = max, 3 = min

public :

	inline int FnCombinationType(void) { return _FnCombinationType ; }
	inline int VarEliminationType(void) { return _VarElimType ; }
	inline int QueryVariable(void) { return _QueryVariable ; }

	inline void SetQueryVariable(int v) { _QueryVariable = v ; }
	inline void SetOperators(int FnCombinationType, int VarElimType) { _FnCombinationType = FnCombinationType ; _VarElimType = VarElimType ; }
	int SetOperators(const char *query)
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
	int SetOperators(const char *combop, const char *elimop)
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
	int ComputeQueryRelevance_VarElimination(ARE_Function_TableType & Factor, char VarElimOp, char CombinationOp) ;

	// ***************************************************************************************************
	// misc
	// ***************************************************************************************************

public :

	// assuming the problem was just constructed/loaded, perform initial analysis.
	// as input we assume just variables, domains, functions.
	// e.g. this fn does not compute variable order. 
	// this fn should run fast.
	int PerformPostConstructionAnalysis(void) ;

	// aliminate all evidence from functions. we assume evidence is set as _Value[var] = value.
	int EliminateEvidence(void) ;

	// given evidence Var=Val, this function will eliminate the evidence variable from all functions, by assigning Var=Val.
	// If as a result, a fn becomes a const fn, it will add it to 
	int EliminateEvidenceVariable(int Var, int Val) ;

	// eliminate all variables, whose domain has just one value, from all functions.
	// all functions that have no arguments left, as the result, will be turned into const functions.
	int EliminateSingletonDomainVariables(void) ;

	// singleton consistency; for probabilistic network, check each domain value against 0-probabilities.
	// return value  0 means ok.
	// return value -1 means domain of some variable became empty.
	int ComputeSingletonConsistency(int & nNewSingletonDomainVariables) ;
	// The following function returns is_consistent[i][j] where i indexes variables and j indexes the domains of variables
	// is_consistent[i][j] is true if the variable value combination is consistent and false otherwise
	void SingletonConsistencyHelper(vector<vector<bool>> & is_consistent) ;

	// this function generates a random Bayesian network structure, that is it does not generate actual tables, 
	// just the variables (arguments) for each function (CPT).
	// tables can be generated by calling FillInFunctionTables().
	int GenerateRandomUniformBayesianNetworkStructure(int N, int K, int P, int C, int ProblemCharacteristic) ;

public :

	virtual int SaveXML(const std::string & Dir) ;
	virtual int SaveUAI08(const std::string & Dir) ;

	// load from file/buffer.
	int LoadFromFile(const std::string & FileName) ;
	int LoadFromFile_Evidence(const std::string & FileName, int & nEvidenceVars) ;

	// this function loads problem from UAI format file; we assume data is already in buf.
	int LoadFromBuffer(const char *format, const char *buf, int L) ;
	int LoadFromBuffer_Evidence(const char *format, const char *buf, int L, int & nEvidenceVars) ;

	int LoadUAIFormat(const char *buf, int L) ;
	int LoadUAIFormat_Evidence(const char *buf, int L, int & nEvidenceVars) ;

public :

	void Destroy(void)
	{
		DestroyAdjVarList() ;
		DestroyAdjFnList() ;
		DestroyVarOrdering() ;
		if (NULL != _Functions) {
			for (int i = 0 ; i < _nFunctions ; i++) {
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

int fileload_getnexttoken(const char * & buf, int & L, const char * & B, int & l) ;

} // namespace ARE

#endif // ARP_HXX_INCLUDED
