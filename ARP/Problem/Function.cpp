#include <stdlib.h>
#include "float.h"

#include "Globals.hxx"

#include "Utils/Sort.hxx"
#include "Utils/MersenneTwister.h"
#include "Problem.hxx"
#include "Function.hxx"
#include "Bucket.hxx"

static MTRand RNG ;

int32_t ARE::Function::SetArguments(int32_t N, const int32_t *Arguments, int32_t ExcludeThisVar)
{
	_BayesianCPTChildVariable = -1 ;

	int32_t i, j, n = N ;

	for (i = 0 ; i < N ; i++) {
		if (Arguments[i] < 0) // neg var not allowed
			return ERRORCODE_BadVarInFnArgumentList ;
		if (Arguments[i] == ExcludeThisVar) 
			{ --n ; continue ; }
		for (j = i+1 ; j < N ; j++) {
			if (Arguments[i] == Arguments[j]) 
				return ERRORCODE_DuplicateVarInFnArgumentList ;
			}
		}

	if (n != _nArgs) {
		if (NULL != _Arguments) {
			delete [] _Arguments ;
			_Arguments = _ArgumentsPermutationList = NULL ;
			}
		_nArgs = 0 ;
		}

	if (NULL != _SortedArgumentsList) { delete [] _SortedArgumentsList ; _SortedArgumentsList = NULL ; }

	if (n <= 0) 
		return 0 ;

	if (NULL == _Arguments) {
		_Arguments = new int32_t[2*n] ;
		if (NULL == _Arguments) 
			return 1 ;
		_nArgs = n ;
		_ArgumentsPermutationList = _Arguments + _nArgs ;
		}

	if (ExcludeThisVar < 0) 
		memcpy(_Arguments, Arguments, _nArgs*sizeof(int32_t)) ;
	else {
		n = 0 ;
		for (i = 0 ; i < N ; i++) {
			if (Arguments[i] != ExcludeThisVar) 
				_Arguments[n++] = Arguments[i] ;
			}
		}
	if (ARE_Function_Type_BayesianCPT == _Type) 
		_BayesianCPTChildVariable = _Arguments[n-1] ;

	return 0 ;
}


int32_t ARE::Function::ConvertTableToLogScale(void)
{
	if (_ConstValue >= 0.0) {
		_ConstValue = log10(_ConstValue) ;
/*		if (_ConstValue > 0.0 && NULL != ARE::fpLOG) {
			fprintf(ARE::fpLOG, "\nFunction::ConvertTableToLog(): ERROR ConstValue=%g (should be <= 0)", (double) _ConstValue) ;
			fflush(ARE::fpLOG) ;
			}*/
		}
	if (NULL == _TableData) 
		return 0 ;
	for (int64_t i = _TableSize-1 ; i >= 0 ; i--) {
		_TableData[i] = log10(_TableData[i]) ;
/*		if (_TableData[i] > 0.0 && NULL != ARE::fpLOG) {
			fprintf(ARE::fpLOG, "\nFunction::ConvertTableToLog(): ERROR data[%lld]=%g (should be <= 0)", i, (double) _TableData[i]) ;
			fflush(ARE::fpLOG) ;
			}*/
		}
	return 0 ;
}


int64_t ARE::Function::ComputeTableSize(void) 
{
	if (_nArgs < 1) {
		_TableSize = 0 ;
		return 0 ;
		}
	_TableSize = 1 ;
	for (int32_t i = 0 ; i < _nArgs ; i++) {
		_TableSize *= _Problem->K(_Arguments[i]) ;
		if (_TableSize < 0) {
			// overflow ?
			_TableSize = -1 ;
			return _TableSize ;
			}
		}
	return _TableSize ;
}


double ARE::Function::GetTableSize_Log10(void) 
{
	if (_nArgs <= 0) 
		return -1.0 ;
	double size = 0.0 ;
	for (int32_t i = 0 ; i < _nArgs ; i++) {
		int32_t k = _Problem->K(_Arguments[i]) ;
		if (k <= 0) 
			return -1.0 ;
		size += log10((double) k) ;
		}
	return size ;
}


int64_t ARE::Function::ComputeTableSpace(void) 
{
	ComputeTableSize() ;
	if (_TableSize <= 0) 
		return 0 ;
	return sizeof(ARE_Function_TableType)*_TableSize ;
}

double ARE::Function::GetTableSpace_Log10(void) 
{
	double x = GetTableSize_Log10() ;
	if (x < 0.0) 
		return -1.0 ;
	return log10(sizeof(ARE_Function_TableType)) + x ;
}

/*
int32_t ARE::Function::GenerateRandomBayesianSignature(int32_t N, int32_t ChildIDX, int32_t & nFVL, int32_t *FreeVariableList)
{
	return 0 ;
}
*/

int32_t ARE::Function::GenerateRandomBayesianSignature(int32_t N, int32_t ChildIDX)
{
	if (ARE_Function_Type_BayesianCPT != _Type) 
		return 1 ;
	if (N < 1 || ChildIDX < 0) 
		return 1 ;
	if (ChildIDX < N-1) 
		// cannot pick N-1 parents for the CPT
		return 1 ;

	// indeces run [0,(Problem->N)-1]; 
	// we assume that children of CPTs are from [P,(Problem->N)-1] and the parents of a CPT are [0,ChildIDX-1].

	int32_t i, j ;

	Destroy() ;
	_nArgs = N ;

	_Arguments = new int32_t[2*_nArgs] ;
	if (NULL == _Arguments) goto failed ;
	for (i = 0 ; i < _nArgs ; i++) 
		_Arguments[i] = -1 ;
	_ArgumentsPermutationList = _Arguments + _nArgs ;

	// place the child as the last argument
	_Arguments[_nArgs-1] = ChildIDX ;
	_BayesianCPTChildVariable = ChildIDX ;

	if (_nArgs > 1) {
		int32_t number_of_parents = _nArgs-1 ;
		if (number_of_parents == ChildIDX) {
			for (i = 0 ; i < number_of_parents ; i++) 
				_Arguments[i] = i ;
			}
		else {
			for (i = 0 ; i < number_of_parents ;) {
				_Arguments[i] = RNG.randInt(ChildIDX-1) ;
				// check that this parent is not already generated
				for (j = 0 ; j < i ; j++) {
					if (_Arguments[j] == _Arguments[i]) break ;
					}
				if (j == i) i++ ;
				}
			}
		}

	ComputeTableSize() ;

	return 0 ;
failed :
	Destroy() ;
	return 1 ;
}


int32_t ARE::Function::FillInRandomBayesianTable(void)
{
	if (ARE_Function_Type_BayesianCPT != _Type) 
		return 1 ;
	if (NULL == _Arguments || _nArgs < 1) 
		return 1 ;

	// indeces run [0,(Problem->N)-1]; 
	// we assume that children of CPTs are from [P,(Problem->N)-1] and the parents of a CPT are [0,ChildIDX-1].

	int32_t i, j ;

	// compute table size; allocate table
	ComputeTableSize() ;
//	DestroyTable() ;
	if (0 != AllocateInMemoryAsSingleTableBlock()) 
		{ Destroy() ; return 1 ; }
	ARE_Function_TableType *table = TableData() ;
	if (NULL == table) 
		{ Destroy() ; return 1 ; }
	int32_t number_of_parents = _nArgs - 1 ;
	int32_t ChildIDX = _Arguments[number_of_parents] ;

	// generate table contents
	int32_t parent_value_combination_array[MAX_NUM_ARGUMENTS_PER_FUNCTION] ; // up to 'MAX_NUM_ARGUMENTS_PER_FUNCTION' parents allowed
	if (number_of_parents > 0) {
		for (j = 0 ; j < number_of_parents - 1 ; j++) parent_value_combination_array[j] = 0 ;
		parent_value_combination_array[number_of_parents - 1] = -1 ;
		}
	// process all parent value combinations
	int32_t k = _Problem->K(ChildIDX) ;
	for (j = 0 ; j < _TableSize ; j += k) {
		// enumerate parent value combination. find next combination from the current.
		for (i = number_of_parents - 1 ; i >= 0 ; i--) {
			if (++parent_value_combination_array[i] < _Problem->K(_Arguments[i])) break ;
			parent_value_combination_array[i] = 0 ;
			}

		ARE_Function_TableType x = (ARE_Function_TableType) 0.0 ; // normalization constant
		ARE_Function_TableType p_of_x ;
		for (i = 0 ; i < k ; i++) {
			p_of_x = 0.0001 + RNG.rand() ;
			x += table[j + i] = (ARE_Function_TableType) p_of_x ;
			}
		// check if all zeros. add some constant to x and to the prob-matrix.
		if (x < 1.0e-15) {
			x += (ARE_Function_TableType) 0.1 ;
			table[RNG.randInt(k-1)] += (ARE_Function_TableType) 0.1 ;
			}
		for (i = 0 ; i < k ; i++) table[j + i] /= x ;
		}

	return 0 ;
failed :
	Destroy() ;
	return 1 ;
}


int32_t ARE::Function::CheckIntegrity(void)
{
	if (NULL == _Problem) 
		return 1 ;

	int32_t i, n = _Problem->N() ;
	for (i = 0 ; i < _nArgs ; i++) {
		if (_Arguments[i] < 0 || _Arguments[i] >= n) 
			return 2 ;
		}

	return CheckData() ;
}


int32_t ARE::Function::AllocateInMemoryAsSingleTableBlock(void)
{
	DestroyTableData() ;

	// make sure table size is known
	_TableSize = ComputeTableSize() ;

	int32_t res = AllocateTableData() ;
	return res ;
}


int32_t ARE::Function::ReOrderArguments(int32_t nAF, const int32_t *AF, int32_t nAB, const int32_t *AB)
{
	if (_nArgs < 1 || NULL == _Arguments) 
		return 0 ;
	if (_nArgs > MAX_NUM_VARIABLES_PER_BUCKET) 
		return 1 ;

	// compute new order
//	int32_t *neworder = new int32_t[_nArgs << 1] ;
//	if (NULL == neworder) 
//		return 1 ;
	int32_t neworder[MAX_NUM_VARIABLES_PER_BUCKET] ;
//	int32_t *new2old_mapping = neworder + _nArgs ;
	int32_t new2old_mapping[MAX_NUM_VARIABLES_PER_BUCKET] ;
	int32_t i, j, b = 0, e = _nArgs ;
	// separate variables of _Arguments into : in AF vs. in AB
	for (i = 0 ; i < _nArgs ; i++) {
		int32_t a = _Arguments[i] ;
		for (j = 0 ; j < nAF ; j++) {
			if (AF[j] == a) 
				break ;
			}
		if (j < nAF) 
			// a is in AF
			neworder[b++] = a ;
		else 
			// a is in AB (or must be in AB)
			neworder[--e] = a ;
		}
	// it should be that b==e; neworder[0,e) are arguments in _Arguments that are in AF, neworder[e,_nArgs) are arguments in _Arguments that are in AB.
	// sort neworder[0,e)/neworder[e,_nArgs) so that it agrees with the variable order in AF/AB.
	b = 0 ;
	for (i = 0 ; i < nAF ; i++) {
		int32_t a = AF[i] ;
		for (j = b ; j < e ; j++) {
			if (neworder[j] == a) 
				break ;
			}
		if (j >= e) 
			// a is not in _Arguments[]
			continue ;
		// place a in the next position (pointed to by 'b') in neworder by swapping
		if (b != j) {
			neworder[j] = neworder[b] ;
			neworder[b] = a ;
			}
		++b ;
		}
	b = e ;
	for (i = 0 ; i < nAB ; i++) {
		int32_t a = AB[i] ;
		for (j = b ; j < _nArgs ; j++) {
			if (neworder[j] == a) 
				break ;
			}
		if (j >= _nArgs) 
			// a is not in _Arguments[]
			continue ;
		// place a in the next position (pointed to by 'b') in neworder by swapping
		if (b != j) {
			neworder[j] = neworder[b] ;
			neworder[b] = a ;
			}
		++b ;
		}

	// check if order has changed
	if (0 == memcmp(_Arguments, neworder, _nArgs*sizeof(int32_t))) {
//		delete [] neworder ;
		return 0 ;
		}

#ifdef _DEBUG
	if (NULL != ARE::fpLOG) {
		char s[256] ;
		if (NULL == _OriginatingBucket) 
			sprintf(s, "\n   o function %d arguments are being reordered (bucket=%d)", (int32_t) _IDX, (int32_t) _Bucket->IDX()) ;
		else 
			sprintf(s, "\n   n function %d arguments are being reordered (bucket=%d)", (int32_t) _OriginatingBucket->IDX(), (int32_t) _Bucket->IDX()) ;
		fwrite(s, 1, strlen(s), ARE::fpLOG) ;
		fflush(ARE::fpLOG) ;
		}
#endif // _DEBUG

	// compute new2old_mapping = for each var in the new ordering, its position in the old (current) ordering
	for (i = 0 ; i < _nArgs ; i++) {
		int32_t a = neworder[i] ;
		for (j = 0 ; j < _nArgs ; j++) {
			if (_Arguments[j] == a) { new2old_mapping[i] = j ; break ; }
			}
		}

	// fix table
	ARE_Function_TableType *newtable = NULL ;
	ARE_Function_TableType *oldtable = TableData() ;
	if (_TableSize > 0 && NULL != oldtable) {
		newtable = new ARE_Function_TableType[_TableSize] ;
		if (NULL == newtable) {
//			delete [] neworder ;
			return 1 ;
			}

		// enumerate all combinations of arguments, using new order.
		// for each combination, compute adr in old table and pull value

		static int32_t value_combination_array_NEW[MAX_NUM_ARGUMENTS_PER_FUNCTION] ; // up to 'MAX_NUM_ARGUMENTS_PER_FUNCTION' arguments allowed
		static int32_t value_combination_array_OLD[MAX_NUM_ARGUMENTS_PER_FUNCTION] ; // up to 'MAX_NUM_ARGUMENTS_PER_FUNCTION' arguments allowed
		for (j = 0 ; j < _nArgs ; j++) value_combination_array_NEW[j] = 0 ;
		value_combination_array_NEW[_nArgs - 1] = -1 ;
		int64_t oldadr ;
		for (j = 0 ; j < _TableSize ; j++) {
			// enumerate value combination. find next combination from the current.
			EnumerateNextArgumentsValueCombination(_nArgs, neworder, value_combination_array_NEW, _Problem->K()) ;
			// convert new combination into old combination
			for (i = 0 ; i < _nArgs ; i++) 
				value_combination_array_OLD[new2old_mapping[i]] = value_combination_array_NEW[i] ;

			// fetch the value of the table cell from the old table and store in the new table
			oldadr = ARE::ComputeFnTableAdr(_nArgs, _Arguments, value_combination_array_OLD, _Problem->K()) ;
//			if (oldadr < 0 || oldadr >= _TableSize) 
//				{ int32_t bad = 1 ; }
			newtable[j] = oldtable[oldadr] ;
			}

		SetTableData(_TableSize, newtable) ;
		}

	// install new arguments list
	memcpy(_Arguments, neworder, _nArgs*sizeof(int32_t)) ;
//	delete [] neworder ;

	return 0 ;
}


int32_t ARE::Function::RemoveVariable(int32_t Var, int32_t Val)
{
	if (_nArgs < 1 || NULL == _Arguments) 
		return 0 ;

	// compute new args array
	int32_t neworder[MAX_NUM_VARIABLES_PER_BUCKET] ;
	int32_t new2old_mapping[MAX_NUM_VARIABLES_PER_BUCKET] ;
	int32_t i, j ;
	int32_t n = 0 ;
	int32_t VarIdx = -1 ;
	for (i = 0 ; i < _nArgs ; i++) {
		int32_t a = _Arguments[i] ;
		if (a == Var) 
			VarIdx = i ;
		else 
			neworder[n++] = a ;
		}
	if (VarIdx < 0) 
		return 1 ;

	if (n < 1) {
		// this fn is a const-value function
		if (NULL != _TableData) 
			_ConstValue = _TableData[Val] ;
		else 
			_ConstValue = 1.0 ;
		_nArgs = 0 ;
		return 0 ;
		}

	// compute new2old_mapping = for each var in the new ordering, its position in the old (current) ordering
	for (i = 0 ; i < n ; i++) {
		int32_t a = neworder[i] ;
		for (j = 0 ; j < _nArgs ; j++) {
			if (_Arguments[j] == a) { new2old_mapping[i] = j ; break ; }
			}
		}

	// dump all table blocks
	int64_t newtablesize = 1 ;
	for (i = 0 ; i < n ; i++) 
		newtablesize *= _Problem->K(neworder[i]) ;

	// fix table
	ARE_Function_TableType *newtable = NULL ;
	ARE_Function_TableType *oldtable = TableData() ;
	if (_TableSize > 0 && NULL != oldtable) {
		newtable = new ARE_Function_TableType[newtablesize] ;
		if (NULL == newtable) 
			return 1 ;

		// enumerate all combinations of arguments, using new order.
		// for each combination, compute adr in old table and pull value

		static int32_t value_combination_array_NEW[MAX_NUM_ARGUMENTS_PER_FUNCTION] ; // up to 'MAX_NUM_ARGUMENTS_PER_FUNCTION' arguments allowed
		static int32_t value_combination_array_OLD[MAX_NUM_ARGUMENTS_PER_FUNCTION] ; // up to 'MAX_NUM_ARGUMENTS_PER_FUNCTION' arguments allowed
		value_combination_array_OLD[VarIdx] = Val ;
		for (j = 0 ; j < n ; j++) value_combination_array_NEW[j] = 0 ;
		value_combination_array_NEW[n - 1] = -1 ;
		int64_t oldadr ;
		for (j = 0 ; j < newtablesize ; j++) {
			// enumerate value combination. find next combination from the current.
			EnumerateNextArgumentsValueCombination(n, neworder, value_combination_array_NEW, _Problem->K()) ;
			// convert new combination into old combination
			for (i = 0 ; i < n ; i++) 
				value_combination_array_OLD[new2old_mapping[i]] = value_combination_array_NEW[i] ;

			// TODO : if var elimination by summation is needed, do it here

			// fetch the value of the table cell from the old table and store in the new table
			oldadr = ARE::ComputeFnTableAdr(_nArgs, _Arguments, value_combination_array_OLD, _Problem->K()) ;
//			if (oldadr < 0 || oldadr >= _TableSize) 
//				{ int32_t bad = 1 ; }
			newtable[j] = oldtable[oldadr] ;
			}

		SetTableData(newtablesize, newtable) ;
		}

	// install new arguments list
	for (i = VarIdx ; i < n ; i++) 
		_Arguments[i] = _Arguments[i+1] ;
	--_nArgs ;

	// delete sorted arguments list
	if (NULL != _SortedArgumentsList) { delete [] _SortedArgumentsList ; _SortedArgumentsList = NULL ; }

	return 0 ;
}


int32_t ARE::Function::RemoveVariableValue(int32_t Var, int32_t Val)
{
	if (_nArgs < 1 || NULL == _Arguments) 
		return 0 ;

	int32_t i, j ;

	// dump all table blocks
	int64_t newtablesize = 1 ;
	int32_t VarIdx = -1 ;
	for (i = 0 ; i < _nArgs ; i++) {
		int32_t a = _Arguments[i] ;
		if (a == Var) {
			VarIdx = i ;
			newtablesize *= _Problem->K(a) - 1 ;
			}
		else 
			newtablesize *= _Problem->K(a) ;
		}
	if (VarIdx < 0) 
		// this should not happen
		return 1 ;
	if (newtablesize < 1) {
		// domain of the variable is empty; function is invalid
		DestroyTableData() ;
		return 0 ;
		}

	// fix table
	ARE_Function_TableType *newtable = NULL ;
	ARE_Function_TableType *oldtable = TableData() ;
	int32_t *newDomainSizes = NULL ;
	if (_TableSize > 0 && NULL != oldtable) {
		newDomainSizes = new int32_t[_Problem->N()] ;
		if (NULL == newDomainSizes) 
			return 1 ;
		newtable = new ARE_Function_TableType[newtablesize] ;
		if (NULL == newtable) {
			delete [] newDomainSizes ;
			return 1 ;
			}

		// enumerate all combinations of arguments; copy all combinations from old table to new table, unless Var=Val.

		static int32_t value_combination_array_OLD[MAX_NUM_ARGUMENTS_PER_FUNCTION] ; // up to 'MAX_NUM_ARGUMENTS_PER_FUNCTION' arguments allowed
		static int32_t value_combination_array_NEW[MAX_NUM_ARGUMENTS_PER_FUNCTION] ; // up to 'MAX_NUM_ARGUMENTS_PER_FUNCTION' arguments allowed
		static int32_t domains_NEW[MAX_NUM_VALUES_PER_VAR_DOMAIN] ;
		for (j = 0 ; j < _nArgs ; j++) value_combination_array_OLD[j] = 0 ;
		value_combination_array_OLD[_nArgs - 1] = -1 ;
		for (j = 0 ; j < _Problem->N() ; j++) 
			newDomainSizes[j] = _Problem->K(j) ;
		newDomainSizes[Var] = _Problem->K(Var) - 1 ;
		int64_t oldadr, newadr ;
		for (j = 0 ; j < _TableSize ; j++) {
			// enumerate (old) value combination. find next combination from the current.
			EnumerateNextArgumentsValueCombination(_nArgs, _Arguments, value_combination_array_OLD, _Problem->K()) ;
			if (Val == value_combination_array_OLD[VarIdx]) 
				// Var=Val; skip this combination.
				continue ;

			// create corresponding new table combination
			for (i = 0 ; i < _nArgs ; i++) 
				value_combination_array_NEW[i] = value_combination_array_OLD[i] ;
			if (value_combination_array_NEW[VarIdx] > Val) 
				value_combination_array_NEW[VarIdx]-- ;

			// fetch the value of the table cell from the old table and store in the new table
			oldadr = ARE::ComputeFnTableAdr(_nArgs, _Arguments, value_combination_array_OLD, _Problem->K()) ;
			newadr = ARE::ComputeFnTableAdr(_nArgs, _Arguments, value_combination_array_NEW, newDomainSizes) ;
			newtable[newadr] = oldtable[oldadr] ;
			}

		SetTableData(newtablesize, newtable) ;
		}

	if (NULL == newtable) {
		delete [] newtable ;
		return 1 ;
		}
	if (NULL == newDomainSizes) {
		delete [] newDomainSizes ;
		return 1 ;
		}

	return 0 ;
}


int32_t ARE::Function::SaveXMLString(const char *prefixspaces, const char *tag, const std::string & Dir, std::string & S)
{
	char s[1024] ;
	int32_t i ;
	if (_IDX >= 0) {
//		GetTableBinaryFilename(Dir, fn) ;
		sprintf(s, "%s<%s IDX=\"%d\" type=\"%s\" nArgs=\"%d\" scope=\"", prefixspaces, tag, _IDX, _Type.c_str(), _nArgs) ;
		}
	else 
		sprintf(s, "%s<%s type=\"%s\" nArgs=\"%d\" scope=\"", prefixspaces, tag, _Type.c_str(), _nArgs) ;
	S += s ;
	for (i = 0 ; i < _nArgs ; i++) {
		sprintf(s, "%d", _Arguments[i]) ;
		if (i > 0) 
			S += ';' ;
		S += s ;
		}
	sprintf(s, "\" tablecelltype=\"%s\"", ARE_Function_TableTypeString) ;
	S += s ;
	if (_TableSize >= 0) {
		sprintf(s, " TableSize=\"%I64d\"", _TableSize) ;
		S += s ;
		}
/*	if (_nTableBlocks >= 0) {
		sprintf(s, " nTableBlocks=\"%I64d\"", _nTableBlocks) ;
		S += s ;
		}
	if (_nTableBlocks >= 0) {
		sprintf(s, " TableBlockSize=\"%I64d\"", _TableBlockSize) ;
		S += s ;
		}*/
	if (_FileName.length() > 0) {
		sprintf(s, " filename=\"%s\"", _FileName.c_str()) ;
		S += s ;
		}
	S += "/>" ;

	return 0 ;
}

