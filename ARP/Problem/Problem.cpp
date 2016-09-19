#include <stdio.h>
#include <stdlib.h>
#include <sstream>
#include <fstream>

// Added by Vibhav for domain pruning implementation
#include "Solver.h"

#include "Utils/MersenneTwister.h"
#include "Utils/Sort.hxx"
#include "Utils/MiscUtils.hxx"
#include "Globals.hxx"
#include "Function.hxx"
#include "Problem.hxx"
//#include "ProblemGraphNode.hxx"

using minisat::Lit;
using minisat::mkLit;
using minisat::Solver;
using minisat::vec;

static MTRand RNG ;

int32_t ARE::ARP::GetFilename(const std::string & Dir, std::string & fn)
{
	if (0 == _Name.length()) 
		return 1 ;
	if (0 == Dir.length()) 
		return 1 ;
	fn = Dir ;
	if ('\\' != fn[fn.length()-1]) 
		fn += '\\' ;
	fn += _Name ;
//	fn += ".xml" ;
	return 0 ;
}


int32_t ARE::ARP::SaveUAI08(const std::string & Dir)
{
	if (NULL == _K) 
		return 1 ;

	std::string fn ;
	if (0 != GetFilename(Dir, fn)) 
		return 1 ;
	if (0 == fn.length()) 
		return 1 ;
	fn += ".uai" ;
	FILE *fp = fopen(fn.c_str(), "w") ;
	if (NULL == fp) 
		return 1 ;

	char s[1024] ;
	int32_t i, j ;

	// save preamble

	sprintf(s, "BAYES\n%d", (int32_t) _nVars) ;
	fwrite(s, 1, strlen(s), fp) ;

	fwrite("\n", 1, 1, fp) ;
	for (i = 0 ; i < _nVars ; i++) {
		if (i > 0) 
			sprintf(s, " %d", (int32_t) _K[i]) ;
		else 
			sprintf(s, "%d", (int32_t) _K[i]) ;
		fwrite(s, 1, strlen(s), fp) ;
		}

	sprintf(s, "\n%d", (int32_t) _nFunctions) ;
	fwrite(s, 1, strlen(s), fp) ;

	for (i = 0 ; i < _nFunctions ; i++) {
		ARE::Function *f = _Functions[i] ;
		if (NULL == f) {
			fwrite("\n0", 1, 2, fp) ;
			continue ;
			}
		sprintf(s, "\n%d", (int32_t) f->N()) ;
		fwrite(s, 1, strlen(s), fp) ;
		for (j = 0 ; j < f->N() ; j++) {
			sprintf(s, " %d", (int32_t) f->Argument(j)) ;
			fwrite(s, 1, strlen(s), fp) ;
			}
		}

	// save function tables

	for (i = 0 ; i < _nFunctions ; i++) {
		ARE::Function *f = _Functions[i] ;
		if (NULL == f) 
			continue ;
		ARE_Function_TableType *td = f->TableData() ;
		if (NULL == td) {
			fwrite("\n\n0", 1, 3, fp) ;
			continue ;
			}
		int64_t ts = f->TableSize() ;
		if (ts < 0) {
			fwrite("\n\n0", 1, 3, fp) ;
			continue ;
			}
		sprintf(s, "\n\n%I64d", ts) ;
		fwrite(s, 1, strlen(s), fp) ;
		int64_t kN = 1, kR = _K[f->Argument(f->N()-1)] ;
		for (j = 0 ; j < f->N() - 1 ; j++) 
			kN *= _K[f->Argument(j)] ;
		// note that ts=kN*kR
		int64_t k, l, K = 0 ;
		for (k = 0 ; k < kN ; k++) {
			fwrite("\n", 1, 1, fp) ;
			for (l = 0 ; l < kR ; l++, K++) {
				sprintf(s, " %g", (double) f->TableEntry(K)) ;
				fwrite(s, 1, strlen(s), fp) ;
				}
			}
		}

	fclose(fp) ;
	return 0 ;
}


int32_t ARE::ARP::SaveXML(const std::string & Dir)
{
	std::string fn ;
	if (0 != GetFilename(Dir, fn)) 
		return 1 ;
	if (0 == fn.length()) 
		return 1 ;
	fn += ".xml" ;
	FILE *fp = fopen(fn.c_str(), "w") ;
	if (NULL == fp) 
		return 1 ;

	char s[1024] ;
	int32_t i ;

	sprintf(s, "<problem name=\"%s\" N=\"%d\" nF=\"%d\">", _Name.c_str(), (int32_t) _nVars, (int32_t) _nFunctions) ;
	fwrite(s, 1, strlen(s), fp) ;

	std::string temp ;
	if (NULL != _K) {
		temp = "\n <variables" ;
		// check if all domains are the same size
		int32_t k = _K[0] ;
		for (i = 1 ; i < _nVars ; i++) {
			if (k != _K[i]) 
				break ;
			}
		if (i >= _nVars) {
			temp += " domainsize=\"" ;
			sprintf(s, "%d", k) ;
			temp += s ;
			temp += "\"" ;
			}
		else {
			temp += " domainsizelist=\"" ;
			for (i = 0 ; i < _nVars ; i++) {
				sprintf(s, "%d", _K[i]) ;
				if (i > 0) 
					temp += ';' ;
				temp += s ;
				}
			}
		temp += "/>" ;
		fwrite(temp.c_str(), 1, temp.length(), fp) ;
		}

	// save functions
	if (_nFunctions > 0 && NULL != _Functions) {
		std::string fFN ;
		for (i = 0 ; i < _nFunctions ; i++) {
			ARE::Function *f = _Functions[i] ;
			if (NULL == f) continue ;
			temp = "\n" ;
			if (0 == f->SaveXMLString(" ", "function", Dir, temp)) 
				fwrite(temp.c_str(), 1, temp.length(), fp) ;
/*
if (CheckFunctions()) {
int32_t error = 1 ;
}
			f->SaveTableBinary(Dir) ;
if (CheckFunctions()) {
int32_t error = 1 ;
}
*/
			}
		}

	// save variable orderings
	if (NULL != _VarOrdering_VarList) {
		fwrite("\n <variableordering>", 1, 20, fp) ;
		if (NULL != _VarOrdering_VarList) {
			temp = "\n  <ordering width=\"" ;
			sprintf(s, "%d", _VarOrdering_InducedWidth) ;
			temp += s ;
			temp += "\" list=\"" ;
			for (i = 0 ; i < _nVars ; i++) {
				sprintf(s, "%d", _VarOrdering_VarList[i]) ;
				if (i > 0) 
					temp += ';' ;
				temp += s ;
				}
			temp += "\"/>" ;
			fwrite(temp.c_str(), 1, temp.length(), fp) ;
			}
		fwrite("\n </variableordering>", 1, 21, fp) ;
		}

	fwrite("\n</problem>", 1, 11, fp) ;
	fclose(fp) ;

/*	// do the actual save of functions as binary files
	for (i = 0 ; i < _nFunctions ; i++) {
		ARE::Function *f = _Functions[i] ;
		if (NULL == f) continue ;
		temp = "\n" ;
		f->SaveTableBinary(Dir) ;
		}*/

//if (CheckFunctions()) {
//int32_t error = 1 ;
//}
	return 0 ;
}


int32_t ARE::fileload_getnexttoken(const char * & buf, int32_t & L, const char * & B, int32_t & l, bool IncludeSpecialSymbols)
{
	B = NULL ;
	l = 0 ;
	// find beginning
	int32_t i ;
	for (i = 0 ; i < L ; i++) { if ('\r' != buf[i] && '\n' != buf[i] && ' ' != buf[i] && '\t' != buf[i]) break ; }
	buf += i ;
	L -= i ;
	if (L <= 0) 
		return 1 ;
	if (IncludeSpecialSymbols) {
		bool token_is_special_symbol = '(' == *buf || ':' == *buf || ')' == *buf ;
		if (token_is_special_symbol) {
			B = buf ;
			l = 1 ;
			buf += l ;
			L -= l ;
			return 0 ;
			}
		}
	// find end
	B = buf ;
	for (l = 1 ; l < L ; l++) {
		if ('\r' == buf[l] || '\n' == buf[l] || ' ' == buf[l] || '\t' == buf[l]) break ;
		if (IncludeSpecialSymbols) {
			bool token_is_special_symbol = '(' == buf[l] || ':' == buf[l] || ')' == buf[l] ;
			if (token_is_special_symbol) 
				break ;
			}
		}
	buf += l ;
	L -= l ;
	return 0 ;
}


int32_t ARE::ARP::LoadFromBuffer(const char *format, const char *buf, int32_t L)
{
	int32_t i ;

	if (NULL != _fpLOG) {
		int64_t tNOW = ARE::GetTimeInMilliseconds() ;
		fprintf(_fpLOG, "\n%I64d ARE::ARP::LoadFromBuffer(); L=%d ...", tNOW, L) ;
		fflush(_fpLOG) ;
		}

	if (0 == stricmp("UAI", format)) {
		if (NULL != _fpLOG) {
			int64_t tNOW = ARE::GetTimeInMilliseconds() ;
			fprintf(_fpLOG, "\n%I64d ARE::ARP::LoadFromBuffer(); format is UAI ...", tNOW) ;
			fflush(_fpLOG) ;
			}
		// figure our type of data : BAYES, MARKOV
		for (i = 0 ; i < L && '\n' != buf[i] && '\r' != buf[i] ; i++) ;
		if (5 == i ? 0 == memcmp(buf, "BAYES", i) : false) {
			if (NULL != _fpLOG) {
				int64_t tNOW = ARE::GetTimeInMilliseconds() ;
				fprintf(_fpLOG, "\n%I64d ARE::ARP::LoadFromBuffer(); UAI type is BAYES ...", tNOW) ;
				fflush(_fpLOG) ;
				}
			i = LoadUAIFormat(buf+6, L-6) ;
			}
		else if (6 == i ? 0 == memcmp(buf, "MARKOV", i) : false) {
			if (NULL != _fpLOG) {
				int64_t tNOW = ARE::GetTimeInMilliseconds() ;
				fprintf(_fpLOG, "\n%I64d ARE::ARP::LoadFromBuffer(); UAI type is MARKOV ...", tNOW) ;
				fflush(_fpLOG) ;
				}
			i = LoadUAIFormat(buf+7, L-7) ;
			}
		else if (6 == i ? 0 == memcmp(buf, "SPARSE", i) : false) {
			if (NULL != _fpLOG) {
				int64_t tNOW = ARE::GetTimeInMilliseconds() ;
				fprintf(_fpLOG, "\n%I64d ARE::ARP::LoadFromBuffer(); UAI type is SPARSE ...", tNOW) ;
				fflush(_fpLOG) ;
				}
			i = LoadUAIFormat(buf+7, L-7) ;
			}
		else {
			if (NULL != _fpLOG) {
				int64_t tNOW = ARE::GetTimeInMilliseconds() ;
				fprintf(_fpLOG, "\n%I64d ARE::ARP::LoadFromBuffer(); UAI type is unknown ...", tNOW) ;
				fflush(_fpLOG) ;
				}
			fprintf(ARE::fpLOG, "\nfile type unknown; will quit ...") ;
			return ERRORCODE_problem_type_unknown ;
			}
		if (0 != i) {
			if (NULL != _fpLOG)
				fprintf(ARE::fpLOG, "\nLoad failed; res = %d", i) ;
			if (i < ERRORCODE_generic)
				i = ERRORCODE_generic ;
			return i ;
			}
		}
	else {
		if (NULL != _fpLOG) {
			int64_t tNOW = ARE::GetTimeInMilliseconds() ;
			fprintf(_fpLOG, "\n%I64d ARE::ARP::LoadFromBuffer(); format is uknown ...", tNOW) ;
			fflush(_fpLOG) ;
			}
		return ERRORCODE_problem_type_unknown ;
		}

	return 0 ;
}


int32_t ARE::ARP::LoadFromBuffer_Evidence(const char *format, const char *buf, int32_t L, int32_t & nEvidenceVars)
{
	nEvidenceVars = 0 ;

	int32_t i ;

	if (0 == stricmp("UAI", format)) {
		i = LoadUAIFormat_Evidence(buf, L, nEvidenceVars) ;
		if (0 != i) {
			if (NULL != ARE::fpLOG) 
				fprintf(ARE::fpLOG, "\nLoad evidence failed; res = %d", i) ;
			return ERRORCODE_generic ;
			}
		}
	else 
		return ERRORCODE_problem_type_unknown ;

	return 0 ;
}


int32_t ARE::ARP::DeleteDuplicateFunctions(void)
{
	// sort functions by scope
	int32_t left[32], right[32];
	QuickSort((void**)_Functions, _nFunctions, left, right, FunctionComparisonByScope_Greater);

	// delete duplicates; [j] is last good ptr. check [i] against [j]; if equal, delete [j]; otherwise set j=i.
	int32_t j = 0, i, k ;
	int32_t nDuplicateFunctions = 0 ;
	for (i = 1; i < _nFunctions; i++) {
		ARE::Function *fi = _Functions[i];
		// compare fi against fj
		if (FunctionComparisonByScope_Equal(_Functions[j], _Functions[i])) {
			ARE::Function *fj = _Functions[j];
			// merge [i] into [j]
			if (0 == fi->N()) {
				ApplyFnCombinationOperator(fj->ConstValue(), fi->ConstValue());
				}
			else {
				for (k = 0; k < fi->N(); k++)
					_Value[fi->Argument(k)] = 0;
				for (int64_t idx = 0; idx < fi->TableSize(); idx++) {
					int64_t adr_i = fi->ComputeFnTableAdr(_Value, _K);
					int64_t adr_j = fj->ComputeFnTableAdr(_Value, _K);
					ApplyFnCombinationOperator((fj->TableData())[adr_i], (fi->TableData())[adr_j]);
					ARE::EnumerateNextArgumentsValueCombinationEx(fi->N(), fi->Arguments(), _Value, _K);
					}
				// reset value so that we won't later think these are evidence variables
				for (k = 0; k < fi->N(); k++)
					_Value[fi->Argument(k)] = -1;
				}
			delete _Functions[i]; _Functions[i] = NULL; nDuplicateFunctions++;
			}
		else
			j = i;
		}
	j = 0;
	for (i = 0; i < _nFunctions; i++) {
		ARE::Function *f = _Functions[i];
		if (NULL != f) {
			_Functions[j] = f;
			f->SetIDX(j++);
			}
		}
	_nFunctions = j;

	return 0 ;
}


int32_t ARE::ARP::LoadFromFile(const std::string & FileName)
{
	int32_t res = -1 ;

	// fn may have dir in it; extract filename.
	int32_t i = FileName.length() - 1 ;
	for (; i >= 0 ; i--) {
#ifdef LINUX
		if ('\\' == FileName[i] || '/' == FileName[i]) 
#else
		if ('\\' == FileName[i] || '//' == FileName[i]) 
#endif
			break ;
		}
	std::string fn(FileName.substr(i+1)), fn_pure, fn_ext ;

	// extract extension
	std::size_t pos_dot = fn.rfind('.');
	if (std::string::npos != pos_dot) 
		{ fn_pure = fn.substr(0, pos_dot) ; fn_ext = fn.substr(pos_dot+1) ; }
	else 
		fn_pure = fn ;

	if (NULL != ARE::fpLOG) 
		fprintf(ARE::fpLOG, "\nWill load problem from %s", FileName.c_str()) ;
	if ("gr" == fn_ext) {
		std::istream *input = NULL ;
		bool delete_istream = false ;
		if ("cin" == fn_pure) {
			input = &(std::cin) ;
			}
		else {
			input = new ifstream(FileName.c_str()) ;
			}

		Destroy();

		std::string s, s1, s2; int32_t nV = -1, nE = -1, nEdgesRead = 0;
		bool have_1stline = false;
		while (true) {
			getline(*input, s);
			if (s.length() <= 0) continue;
			if ('c' == s[0] || 'C' == s[0]) continue;
			std::istringstream line_stream(s);
			if (!have_1stline) {
				line_stream >> s1; line_stream >> s2; line_stream >> nV; line_stream >> nE;
				if (0 == s1.length() || 0 == s2.length()) { res = 1 ; goto done ; }
				if ('p' != s1[0] || "tw" != s2) { res = 1; goto done; }
				have_1stline = true; break;
				}
			if (++nEdgesRead >= nE)
				break;
			}

		SetN(nV) ;
		SetK(1) ;

		_nFunctions = nE;
		if (_nFunctions < 1)
			{ _nFunctions = 0; res = 0; goto done; }
		_Functions = new ARE::Function*[_nFunctions];
		if (NULL == _Functions)
			{ res = 1 ; goto done ; }
		for (i = 0; i < _nFunctions; i++)
			_Functions[i] = NULL;

		// read in function signatures
		int32_t A[2] ;
		while (true) {
			getline(*input, s);
			if (s.length() <= 0) continue;
			if ('c' == s[0] || 'C' == s[0]) continue;
			std::istringstream line_stream(s);
			int32_t v1, v2 ; line_stream >> v1 ; line_stream >> v2 ;
			_Functions[nEdgesRead] = new ARE::Function(NULL, this, nEdgesRead);
			if (NULL == _Functions[nEdgesRead])
				{ res = 1; goto done; }
			ARE::Function *f = _Functions[nEdgesRead];
			f->SetType(ARE_Function_Type_Const);
			A[0] = v1-1 ; A[1] = v2-1 ; // var indeces in input range [1,n]
			if (A[0] < 0 || A[0] >= _nVars || A[1] < 0 || A[1] >= _nVars)
				{ res = 1; goto done; }
			if (0 != f->SetArguments(2, A, -1))
				goto done;
			// sort scope, so that we can sort the function list easily
			int32_t *sorted_scope = f->SortedArgumentsList(true);
			if (NULL == sorted_scope)
				{ res = 1; goto done; }
			int64_t ts = f->ComputeTableSize();
			f->ConstValue() = 0.0 ;
			if (++nEdgesRead >= nE)
				break;
			}

		DeleteDuplicateFunctions() ;

		if (delete_istream) 
			{ delete input ; input = NULL ; }
		res = 0 ; goto done ;
		}
	else {
		FILE *fp = fopen(FileName.c_str(), "rb") ;
		if (NULL == fp) {
			if (NULL != ARE::fpLOG)
				fprintf(ARE::fpLOG, "\nfailed to open file; will quit ...") ;
			return ERRORCODE_cannot_open_file ;
			}
		// get file size
		fseek(fp, 0, SEEK_END) ;
		int32_t filesize = ftell(fp) ;
		fseek(fp, 0, SEEK_SET) ;
//#define bufsize	 1048576
//#define bufsize	16777216
//#define bufsize	33554432
		char *buf = new char[filesize] ;
		if (NULL == buf) {
			fclose(fp) ;
			if (NULL != ARE::fpLOG)
				fprintf(ARE::fpLOG, "\nfailed to allocate buffer for file data; will quit ...") ;
			return ERRORCODE_cannot_allocate_buffer_size_too_large ;
			}
		int32_t L = fread(buf, 1, filesize, fp) ;
		fclose(fp) ;
		if (filesize != L) {
			// we probably did not read in the whole file; buf is too small
			// todo : can use foef() or ferror.
			if (NULL != ARE::fpLOG) 
				fprintf(ARE::fpLOG, "\nfailed to load file; will quit ...") ;
			delete [] buf ;
			return ERRORCODE_file_too_large ;
			}

		res = LoadFromBuffer("UAI", buf, L) ;
		delete [] buf ;
		}

done :
	return res ;
}


int32_t ARE::ARP::LoadUAIFormat(const char *buf, int32_t L)
{
	if (NULL == buf) 
		return 1 ;

	int32_t ret = 1, i, j ;
	char *BUF = NULL ;
	int32_t *K = NULL ;

	Destroy() ;

	{

	const char *token = NULL ;
	int32_t l = 0 ;
/*
	if (0 != fileload_getnexttoken(buf, L, token, l)) 
		goto done ;
	if (5 != l ? true : 0 != memcmp(token, "BAYES", l)) 
		goto done ;
*/
	// get # of variables
	if (0 != fileload_getnexttoken(buf, L, token, l, false)) 
		goto done ;
	int32_t N = atoi(token) ;

	if (N > 0) {
		// get domain size of each variable
		K = new int32_t[N] ;
		if (NULL == K) goto done ;
		for (i = 0 ; i < N ; i++) {
			if (0 != fileload_getnexttoken(buf, L, token, l, false)) 
				goto done ;
			K[i] = atoi(token) ;
			if (K[i] < 0 || K[i] > MAX_NUM_VALUES_PER_VAR_DOMAIN) 
				{ ret = ERRORCODE_VarDomainSizeTooLarge ; goto done ; }
			}
		if (0 != SetN(N)) 
			goto done ;
		if (0 != SetK(K)) 
			goto done ;
		}

	// get # of functions
	if (0 != fileload_getnexttoken(buf, L, token, l, false)) 
		goto done ;
	_nFunctions = atoi(token) ;
	if (_nFunctions < 1) 
		{ _nFunctions = 0 ; ret = 0 ; goto done ; }
	_Functions = new ARE::Function*[_nFunctions] ;
	if (NULL == _Functions) 
		goto done ;
	for (i = 0 ; i < _nFunctions ; i++) 
		_Functions[i] = NULL ;

	// read in function signatures
	int32_t A[MAX_NUM_ARGUMENTS_PER_FUNCTION] ;
	for (i = 0 ; i < _nFunctions ; i++) {
		if (0 != fileload_getnexttoken(buf, L, token, l, false)) 
			goto done ;
		int32_t nA = atoi(token) ;
// 2014-11-09 KK : allow const functions; a function with nA==0 is a const fn.
//		if (0 == nA) 
//			// this means function is missing essentially
//			continue ;
		if (nA < 0) 
			goto done ;
		_Functions[i] = new ARE::Function(NULL, this, i) ;
		if (NULL == _Functions[i]) 
			goto done ;
		ARE::Function *f = _Functions[i] ;
		f->SetType(ARE_Function_Type_RealCost) ;
		for (j = 0 ; j < nA ; j++) {
			if (0 != fileload_getnexttoken(buf, L, token, l, false)) 
				goto done ;
			A[j] = atoi(token) ;
			if (A[j] < 0 || A[j] >= _nVars) 
				goto done ;
			}
		if (0 != f->SetArguments(nA, A, -1)) 
			goto done ;
		// sort scope, so that we can sort the function list easily
		if (nA > 0) {
			int32_t *sorted_scope = f->SortedArgumentsList(true) ;
			if (NULL == sorted_scope) 
				goto done ;
			}
		int64_t ts = f->ComputeTableSize() ;
		}

	// create function tables
	for (i = 0 ; i < _nFunctions ; i++) {
		ARE::Function *f = _Functions[i] ;
		if (NULL == f) 
			goto done ;
		f->ComputeTableSize() ;
		f->AllocateInMemoryAsSingleTableBlock() ;
		ARE_Function_TableType *td = f->TableData() ;
		if (NULL == td) {
			if (f->TableSize() > 0) 
				goto done ;
			else 
				continue ;
			}
		}

	// load function tables
	for (i = 0 ; i < _nFunctions ; i++) {
		ARE::Function *f = _Functions[i] ;
		if (NULL == f) 
			goto done ;
		if (f->TableSize() < 0) 
			continue ; // error
		ARE_Function_TableType *data = f->TableData() ; // this may be NULL if const fn
		// load # of entries
		if (0 != fileload_getnexttoken(buf, L, token, l, false)) 
			goto done ;
		int32_t n = atoi(token) ;
		if (n <= 0) 
			goto done ;
		if (f->N() > 0 ? n != f->TableSize() : 1 != n) 
			goto done ;
		for (j = 0 ; j < n ;) {
			if (0 != fileload_getnexttoken(buf, L, token, l, true)) 
				goto done ;
			// check sparse-factor encoding "(x:n)" where x is real number and n is an int32_t. meaning of this is, next entries in the table equal x.
			if (1 == l && '(' == *token) {
				// next token is x
				if (0 != fileload_getnexttoken(buf, L, token, l, true)) 
					goto done ;
				double x = atof(token) ;
				// next token is ":"
				if (0 != fileload_getnexttoken(buf, L, token, l, true)) 
					goto done ;
				if (1 != l || ':' != *token) 
					goto done ;
				// next token is n
				if (0 != fileload_getnexttoken(buf, L, token, l, true)) 
					goto done ;
				int32_t nEntriesLeft = n - j ;
				int32_t nEntries = atoi(token) ;
				if (nEntries <= 0 || nEntries > nEntriesLeft) 
					goto done ;
				// next token is ")"
				if (0 != fileload_getnexttoken(buf, L, token, l, true)) 
					goto done ;
				if (1 != l || ')' != *token) 
					goto done ;
				if (f->N() <= 0) {
					f->ConstValue() = x ;
					j += nEntries ;
					}
				else {
					for (int32_t ni = 0 ; ni < nEntries ; ni++) 
						data[j++] = x ;
					}
				}
			else {
				// token is a table entry
				double x = atof(token) ;
				if (f->N() <= 0) 
					{ f->ConstValue() = x ; j++ ; }
				else 
					data[j++] = x ;
				}
			}
		}

	DeleteDuplicateFunctions() ;

	ret = 0 ;

	}

done :

	if (NULL != BUF) 
		delete [] BUF ;
	if (NULL != K) 
		delete [] K ;
	if (0 != ret) 
		Destroy() ;

	return ret ;
}


int32_t ARE::ARP::LoadFromFile_Evidence(const std::string & FileName, int32_t & nEvidenceVars)
{
	nEvidenceVars = 0 ;

	int32_t ret = 1 ;
	const char *buf = NULL ;
	char *BUF = NULL ;
	int32_t *K = NULL ;

	FILE *fp = fopen(FileName.c_str(), "r") ;
	if (NULL == fp) 
		return 1 ;

	// get file size
	fseek(fp, 0, SEEK_END) ;
	int32_t filesize = ftell(fp) ;
	if (filesize <= 0) {
		printf("\nevidencefile empty; will quit ...") ;
		return 1 ;
		}
	fseek(fp, 0, SEEK_SET) ;

	BUF = new char[filesize] ;
	if (NULL == BUF) {
		printf("\nfailed to allocate memory for evidencefile; will quit ...") ;
		return 1 ;
		}
	int32_t L = fread(BUF, 1, filesize, fp) ;
	fclose(fp) ;
	if (filesize == L) {
		delete [] BUF ;
		printf("\nfailed to load evidence file ...") ;
		return 1 ;
		}

	int32_t res = LoadUAIFormat_Evidence(BUF, L, nEvidenceVars) ;
	if (0 == res) 
		_EvidenceFileName = FileName ;
	delete [] BUF ;
	return res ;
}


int32_t ARE::ARP::LoadUAIFormat_Evidence(const char *BUF, int32_t L, int32_t & nEvidenceVars)
{
	// file format is : <nEvidence> followed by nEvidence times <evidVar> <evidValue>.
	// it will set evidence as the current value of the variable.

	nEvidenceVars = 0 ;

	int32_t ret = 1, i ;
	const char *buf = NULL ;
	int32_t *K = NULL ;

	const char *token = NULL ;
	int32_t l = 0 ;

	int32_t *var_loaded = NULL ;
	int32_t *val_loaded = NULL ;

	{
	// get # of evidence
	buf = BUF ;
	if (0 != fileload_getnexttoken(buf, L, token, l, false)) 
		goto done ;
	int32_t N = atoi(token) ;
	if (N < 1) 
		goto doneok ;

	// need _K and _Value arrays to exist
	if (NULL == _K || NULL == _Value) 
		goto done ;

	var_loaded = new int32_t[N] ;
	val_loaded = new int32_t[N] ;
	if (NULL == var_loaded || NULL == val_loaded) 
		goto done ;

	// read in evidence
	int32_t n = 0 ;
	for (i = 0 ; i < N ; i++) {
		if (0 != fileload_getnexttoken(buf, L, token, l, false)) 
			goto done ;
		int32_t var = atoi(token) ;
		if (0 != fileload_getnexttoken(buf, L, token, l, false)) 
			goto done ;
		int32_t val = atoi(token) ;
		if (var < 0 || var >= _nVars) 
			goto done ;
		if (val < 0 || val >= _K[var]) 
			goto done ;
		var_loaded[n] = var ;
		val_loaded[n] = val ;
		++n ;
		}
	for (i = 0 ; i < n ; i++) {
		_Value[var_loaded[i]] = val_loaded[i] ;
		++nEvidenceVars ;
		}
	}

doneok : 
	ret = 0 ;
done :
	if (NULL != val_loaded) 
		delete [] val_loaded ;
	if (NULL != var_loaded) 
		delete [] var_loaded ;
	return ret ;
}


int32_t ARE::ARP::PerformPostConstructionAnalysis(void) 
{
	int32_t i = 0 ;
	if (0 == i) 
		i = ComputeAdjFnList(false) ;
	if (0 == i) 
		i = ComputeAdjVarList() ;
	if (0 == i) 
		i = ComputeConnectedComponents() ;
	return i ;
}


int64_t ARE::ARP::ComputeFunctionSpace(void)
{
	int64_t n ;
	_FunctionsSpace = 0 ;

	for (int32_t i = 0 ; i < _nFunctions ; i++) {
		ARE::Function *f = _Functions[i] ;
		if (NULL == f) continue ;
		n = f->ComputeTableSpace() ;
		if (n > 0) 
			_FunctionsSpace += n ;
		}

	return _FunctionsSpace ;
}


int32_t ARE::ARP::ConvertFunctionsToLogScale(void)
{
	if (_FunctionsAreConvertedToLogScale) 
		return 0 ;
	_FunctionsAreConvertedToLogScale = true ;
	for (int32_t i = 0 ; i < nFunctions() ; i++) {
		ARE::Function *f = getFunction(i) ;
		if (NULL != f) 
			f->ConvertTableToLogScale() ;
		}
	return 0 ;
}


int32_t ARE::ARP::CheckFunctions(void)
{
	for (int32_t i = 0 ; i < _nFunctions ; i++) {
		ARE::Function *f = _Functions[i] ;
		if (NULL == f) continue ;
		if (f->CheckData()) 
			return 1 ;
		}
	return 0 ;
}


int32_t ARE::ARP::FillInFunctionTables(void)
{
	for (int32_t i = 0 ; i < _nFunctions ; i++) {
		ARE::Function *f = _Functions[i] ;
		if (NULL == f) continue ;
		if (f->ComputeTableSize() <= 0) 
			// this function does not have a table
			continue ;
		if (0 != f->FillInRandomBayesianTable()) 
			return 2 ;
		}

	return 0 ;
}


int32_t ARE::ARP::ComputeQueryRelevance_VarElimination(ARE_Function_TableType & Factor, char VarElimOp, char CombinationOp)
{
	// we assume Factor is initialled

	int32_t i, j, k ;

	if (0 != VarElimOp && 1 != VarElimOp) 
		return 1 ;

	double neutral_value = (0 == CombinationOp) ? 1.0 : 0.0 ;

	if (NULL != ARE::fpLOG) 
		fprintf(ARE::fpLOG, "\nARP::ComputeQueryRelevance_VarElimination() combination_type=%d VarElimOp=%d neutral_value=%g", (int32_t) CombinationOp, (int32_t) VarElimOp, neutral_value) ;

	_nFunctionsIrrelevant = 0 ;
	for (i = 0 ; i < _nFunctions ; i++) {
		if (NULL != _Functions[i]) 
			_Functions[i]->MarkAsQueryRelevant() ;
		}

	if (_nFunctions < 1 || _nVars < 1) 
		return 0 ;

	int32_t nLeaves = 0 ;
//	ARE::Function **FNs = new ARE::Function*[_nFunctions] ;
	int32_t *Leaves = new int32_t[_nVars] ;
	if (NULL == Leaves) 
		return 1 ;
	char *VarIsInQueue = new char [_nVars] ;
	if (NULL == VarIsInQueue) 
		{ delete [] Leaves ; return 1 ; }
	for (i = 0 ; i < _nVars ; i++) 
		VarIsInQueue[i] = 0 ; // 0=not in queue 1=in queue

	// construct an initial list of variables that participate in 1 fn only
	for (i = 0 ; i < _nVars ; i++) {
		if (1 != nAdjFunctions(i)) 
			continue ;
/*		ARPGraphNode & adj = _GraphAdjacencyMatrix[i] ;
		if (adj.nAdjacentFunctions() > 1) 
			// cannot be a leaf if more than 1 adj function
			continue ;
		if (adj.nAdjacentFunctions() < 1) 
			// this is a singleton variable; ignore it.
			continue ;*/
		Leaves[nLeaves++] = i ;
		VarIsInQueue[i] = 1 ;
		}

	while (nLeaves > 0) {
		int32_t var = Leaves[--nLeaves] ;
		VarIsInQueue[var] = 0 ;
//		ARPGraphNode & adj = _GraphAdjacencyMatrix[var] ;
		if (1 != nAdjacentFunctions_QueryReleventFunctionsOnly(var)) {
			if (NULL != ARE::fpLOG) 
				fprintf(ARE::fpLOG, "\nWARNING : var %d is being analyzed for adj-fn query relevance, and 1 != adj.nAdjacentFunctions_QueryReleventFunctionsOnly()", var) ;
			continue ;
			}

		ARE_Function_TableType elim_value = -1.0 ;

		ARE::Function *f = AdjacentFunction_QueryRelevantFunctionsOnly(var, 0) ;
		ARE_Function_TableType *data = f->TableData() ;
		if (NULL == data) 
			goto dump_table ;

		// check that if we eliminated 'var' from this fn, all entries would be the same (some const value).
		{

		int32_t Val[MAX_NUM_VARIABLES_PER_BUCKET], val[MAX_NUM_VARIABLES_PER_BUCKET] ; // this is the current value combination of the arguments of this function (table block).
		int32_t Var[MAX_NUM_VARIABLES_PER_BUCKET] ; // this is the current value combination of the arguments of this function (table block).
		int32_t iVar = -1, n ;
		int64_t idx, size = 1 ;
		for (i = n = 0 ; i < f->N() ; i++) {
			Val[i] = val[i] = 0 ;
			int32_t u = f->Argument(i) ;
			if (var == u) 
				{ iVar = i ; continue ; }
			Var[n++] = u ;
			size *= _K[u] ;
			}
		bool fn_is_relevant = false ;
		for (idx = 0 ; idx < size ; idx++) {
			for (i = 0 ; i < n ; i++) 
				val[i + (i < iVar ? 0 : 1)] = Val[i] ;
			double v = 0.0 ;
			for (k = 0 ; k < _K[var] ; k++, idx++) {
				val[iVar] = k ;
				idx = ComputeFnTableAdr(f->N(), f->Arguments(), val, _K) ;
				if (0 == VarElimOp) 
					v += f->TableEntry(idx) ;
				else if (1 == VarElimOp) 
					{ if (f->TableEntry(idx) > v) v = f->TableEntry(idx) ; }
				}
			if (elim_value < 0.0) {
				elim_value = v ;
				}
			else if (fabs(v - elim_value) > 0.000001) {
				fn_is_relevant = true ;
				break ;
				}
			ARE::EnumerateNextArgumentsValueCombination(n, Var, Val, _K) ;
			}
/*
		ARE::Function *f = GetCPT(var) ;
		if (NULL == f) {
			// this is an error
			if (NULL != ARE::fpLOG) 
				fprintf(ARE::fpLOG, "\nERROR : var %d is being analyzed for adj-fn query relevance, and it has no CPT", var) ;
			continue ;
			}
		ARE::FunctionTableBlock *ftb = f->Table() ;
		if (NULL == ftb) 
			goto dump_table ;
		int64_t n = f->TableSize() ;
		int32_t cVar = f->BayesianCPTChildVariable() ;
		int64_t idx ;
		bool fn_is_relevant = false ;
		for (idx = 0 ; idx < n ; idx++) {
			double sum = 0.0 ;
			for (k = 0 ; k < _K[var] ; k++, idx++) 
				sum += ftb->Entry(idx) ;
			if (fabs(sum - 1.0) > 0.000001) {
				fn_is_relevant = true ;
				break ;
				}
			}
*/
		if (fn_is_relevant) 
			continue ;
		}

dump_table :
		if (f->IsQueryIrrelevant()) {
			// this is an error
			if (NULL != ARE::fpLOG) 
				fprintf(ARE::fpLOG, "\nERROR : var %d is being analyzed for adj-fn query relevance, and its CPT %d is already marked as irrelevant", var, f->IDX()) ;
			}
		else {
			if (NULL != ARE::fpLOG) 
				fprintf(ARE::fpLOG, "\nARP::ComputeQueryRelevance_VarElimination() fn=%d", (int32_t) f->IDX()) ;
			f->MarkAsQueryIrrelevant() ;
			++_nFunctionsIrrelevant ;
			if (fabs(elim_value - neutral_value) > 0.000001) {
				if (0 == CombinationOp) 
					Factor *= elim_value ;
				else if (1 == CombinationOp) 
					Factor += elim_value ;
				}
			}
		// check if other variables in the fn need checking
		for (k = 0 ; k < f->N() ; k++) {
			int32_t u = f->Argument(k) ;
			if (u == var) 
				continue ;
			if (0 != VarIsInQueue[u]) 
				continue ;
//			ARPGraphNode & u_adj = _GraphAdjacencyMatrix[u] ;
			if (1 != nAdjacentFunctions_QueryReleventFunctionsOnly(u)) 
				continue ;
			Leaves[nLeaves++] = u ;
			VarIsInQueue[u] = 1 ;
			}
		}

/*
static int32_t S = 0 ;
static int32_t E = 0 ;
j = 0 ;
for (i = 0 ; i < _nFunctions ; i++) {
	ARE::Function *f = _Functions[i] ;
	if (! f->IsQueryIrrelevant()) 
		continue ;
	if (j < S || j > E) 
		f->MarkAsQueryRelevant() ;
	++j ;
	}
*/

	delete [] VarIsInQueue ;
	delete [] Leaves ;

	if (NULL != ARE::fpLOG) {
		fprintf(ARE::fpLOG, "\nARP::ComputeQueryRelevance_VarElimination() done ...") ;
		fflush(ARE::fpLOG) ;
		}

	return 0 ;
}


int32_t ARE::ARP::ComputeBayesianAncestors(int32_t V, int32_t *AncestorFlags, int32_t *Workspace)
{
	int32_t i, j ;

	for (i = 0 ; i < _nVars ; i++) 
		AncestorFlags[i] = 0 ;

	Workspace[0] = V ;
	int32_t n = 1 ;
	while (n > 0) {
		int32_t v = Workspace[0] ;
		Workspace[0] = Workspace[--n] ;
		if (0 != AncestorFlags[v]) 
			continue ;
		AncestorFlags[v] = 1 ;
		ARE::Function *f = GetCPT(v) ;
		if (NULL == f) 
			continue ;
		for (i = f->N() - 2 ; i >= 0 ; i--) {
			int32_t u = f->Argument(i) ;
			if (0 != AncestorFlags[u]) 
				continue ;
			for (j = 0 ; j < n ; j++) {
				if (Workspace[j] == u) 
					break ;
				}
			if (j >= n) 
				Workspace[n++] = u ;
			}
		}

	return 0 ;
}


int32_t ARE::ARP::GenerateRandomUniformBayesianNetworkStructure(int32_t N, int32_t K, int32_t P, int32_t C, int32_t ProblemCharacteristic)
{
	Destroy() ;

	int32_t i, j ;

	// check if parameters are valid
	if (N < 1 || K < 1 || K > MAX_NUM_VALUES_PER_VAR_DOMAIN) return 1 ;
	if (C < 0 ||P < 1) return 1 ;
	if (N <= P) return 1 ;
	if (P+1 > MAX_NUM_ARGUMENTS_PER_FUNCTION) return 1 ;

	// safety check.
	// notice, that C cannot be larger than N-P since for last P variables in the ordering,
	// we don't have enough parents to put in the parent set.
	if (C > N - P) C = N - P ;

	_nVars = N ;
	_K = new int32_t[_nVars] ;
	if (NULL == _K) 
		goto failed ;
	_Value = new int32_t[_nVars] ;
	if (NULL == _Value) 
		goto failed ;
	for (i = 0 ; i < _nVars ; i++) {
		_K[i] = K ;
		_Value[i] = -1 ;
		}
	_nFunctions = 0 ;
	_Functions = new ARE::Function*[_nVars] ;
	if (NULL == _Functions) goto failed ;
	for (i = 0 ; i < _nVars ; i++) 
		_Functions[i] = NULL ;

	if (1 == ProblemCharacteristic) {
		int32_t ret = 1 ;
		int32_t *space = new int32_t[4*_nVars] ;
		if (NULL == space) 
			goto failed ;
		// idea : 
		// 1) maintain a set of current-priors ; i.e. variables that are parents of a CPT such that these variables themselves have no CPT
		// 2) recursively, pick a var from the set of current-priors, and then generate parent as follows: 
		//    pick a var, randomly among all variables, and check if current child is an ancestor of this picked var;
		//    if yes, this picked var is no good.
		int32_t *L = space ; // L is the set of leaves
		int32_t *A = space + _nVars ; // ancestor bits
		int32_t *W = space + 2*_nVars ; // workspace bits
		int32_t *PL = space + 3*_nVars ; // potential parent list for child

		L[0] = _nVars-1 ;
		int32_t nL = 1 ;

		while (_nFunctions < C) {
			// pick child
			int32_t child = -1 ;
			if (nL <= 0) goto probchar1_done ;
			else if (1 == nL) { child = L[0] ; nL = 0 ; }
			else { int32_t idx = RNG.randInt(nL-1) ; child = L[idx] ; L[idx] = L[--nL] ; }
			// generate a list of possible parents
			int32_t nPL = 0 ;
			for (i = 0 ; i < _nVars ; i++) {
				if (child == i) continue ;
				ComputeBayesianAncestors(i, A, W) ;
				if (0 == A[child]) 
					PL[nPL++] = i ;
				}
			// PL is now a list of potential parents; generate actual parents
			int32_t args[MAX_NUM_ARGUMENTS_PER_FUNCTION] ;
			if (nPL < P) 
				goto probchar1_done ;
			else if (nPL == P) {
				for (i = 0 ; i < P ; i++) args[i] = PL[i] ;
				nPL = 0 ;
				}
			else {
				for (i = 0 ; i < P ; i++) {
					int32_t idx = RNG.randInt(nPL-1) ;
					args[i] = PL[idx] ;
					PL[idx] = PL[--nPL] ;
					}
				}

			// create function
			ARE::Function *f = new ARE::Function(NULL, this, child) ;
			if (NULL == f) goto probchar1_done ;
			f->SetType(ARE_Function_Type_BayesianCPT) ;
			_Functions[child] = f ;
			++_nFunctions ;
			args[P] = child ;
			if (0 != f->SetArguments(P+1, args)) 
				goto probchar1_done ;

			// add arguments of f to L
			for (i = 0 ; i < P ; i++) {
				int32_t v = args[i] ;
				if (NULL != _Functions[v]) 
					continue ;
				for (j = 0 ; j < nL ; j++) { if (v == L[j]) break ; }
				if (j >= nL) 
					L[nL++] = v ;
				}
			}

		ret = 0 ;
probchar1_done :
		delete [] space ;
		if (0 != ret) 
			goto failed ;
		}
	else {
		// create C conditional probabilities
		while (_nFunctions < C) {
			// we will let the CSP_Bayesian_PT constructor create the parent set and the probability matrix for each conditional probability.
			// indeces run [0,N-1]; we assume that children of CPTs are from [P,N-1] and the parents of a CPT are [0,ChildIDX-1].
			// randomly pick a childIDX.
			int32_t child = P + RNG.randInt(_nVars - P-1) ;
			// check if a conditional probability has been created for this variable
			if (NULL != _Functions[child]) continue ;
			ARE::Function *f = new ARE::Function(NULL, this, child) ;
			if (NULL == f) goto failed ;
			f->SetType(ARE_Function_Type_BayesianCPT) ;
			_Functions[child] = f ;
			++_nFunctions ;
			if (0 != f->GenerateRandomBayesianSignature(P+1, child)) goto failed ;
//			if (0 != f->GenerateRandomBayesianTable(P+1, child)) goto failed ;
			}
		}

	// create prior probability matrices for variables which have no conditional probability matrix
	for (i = 0 ; i < _nVars ; i++) {
		if (NULL != _Functions[i]) continue ;
		ARE::Function *f = new ARE::Function(NULL, this, i) ;
		if (NULL == f) goto failed ;
		f->SetType(ARE_Function_Type_BayesianCPT) ;
		_Functions[i] = f ;
		++_nFunctions ;
		if (0 != f->GenerateRandomBayesianSignature(1, i)) goto failed ;
//		if (0 != f->GenerateRandomBayesianTable(1, i)) goto failed ;
		}

	return 0 ;
failed :
	Destroy() ;
	return 1 ;
}


int32_t ARE::ARP::DestroyAdjFnList(void)
{
	if (NULL != _StaticAdjFnTotalList) {
		delete [] _StaticAdjFnTotalList ;
		_StaticAdjFnTotalList = NULL ;
		}
	_StaticAdjFnTotalListSize = 0 ;
	if (NULL != _nAdjFunctions) {
		delete [] _nAdjFunctions ;
		_nAdjFunctions = NULL ;
		}
	if (NULL != _AdjFunctions) {
		delete [] _AdjFunctions ;
		_AdjFunctions = NULL ;
		}
	return 0 ;
}


int32_t ARE::ARP::ComputeAdjFnList(bool IgnoreIrrelevantFunctions)
{
	if (0 != DestroyAdjFnList()) 
		return 1 ;
	if (_nVars < 1) 
		return 0 ;

	int32_t i, j, n = 0 ;
	for (i = 0 ; i < _nFunctions ; i++) {
		Function *f = _Functions[i] ;
		if (NULL == f) continue ;
		if (IgnoreIrrelevantFunctions && f->IsQueryIrrelevant()) 
			continue ;
		n += f->N() ;
		}
	_nAdjFunctions = new int32_t[_nVars] ;
	_AdjFunctions = new int32_t[_nVars] ;
	_StaticAdjFnTotalList = n > 0 ? new Function*[n] : NULL ;
	if (NULL == _nAdjFunctions || NULL == _AdjFunctions || (n > 0 && NULL == _StaticAdjFnTotalList)) {
		DestroyAdjFnList() ;
		return 1 ;
		}
	_StaticAdjFnTotalListSize = n ;

	for (i = 0 ; i < _nVars ; i++) {
		_nAdjFunctions[i] = 0 ;
		_AdjFunctions[i] = -1 ;
		}

	if (_nFunctions < 1) 
		return 0 ;

	// compute num adj FNs for each var
	for (i = 0 ; i < _nFunctions ; i++) {
		Function *f = _Functions[i] ;
		if (NULL == f) continue ;
		if (IgnoreIrrelevantFunctions && f->IsQueryIrrelevant()) 
			continue ;
		for (j = 0 ; j < f->N() ; j++) {
			int32_t var = f->Argument(j) ;
			++_nAdjFunctions[var] ;
			}
		}
	// prep _AdjFunctions[]
	n = 0 ;
	for (i = 0 ; i < _nVars ; i++) {
		if (_nAdjFunctions[i] > 0) {
			_AdjFunctions[i] = n ;
			n += _nAdjFunctions[i] ;
			_nAdjFunctions[i] = 0 ; // reset to 0; we will compute correct value later.
			}
		}

	// fill in adj FNs ptrs for each var
	for (i = 0 ; i < _nFunctions ; i++) {
		Function *f = _Functions[i] ;
		if (NULL == f) continue ;
		if (IgnoreIrrelevantFunctions && f->IsQueryIrrelevant()) 
			continue ;
		for (j = 0 ; j < f->N() ; j++) {
			int32_t var = f->Argument(j) ;
			_StaticAdjFnTotalList[_AdjFunctions[var] + _nAdjFunctions[var]] = f ;
			++_nAdjFunctions[var] ;
			}
		}

	// DEBUG : check how many variables have functions associated with them
	for (i = 0 ; i < _nVars ; i++) {
		int32_t nFNs = 0 ;
		for (j = 0 ; j < nAdjFunctions(i) ; j++) {
			ARE::Function *f = AdjFunction(i, j) ;
			if (NULL == f) 
				continue ;
			++nFNs ;
			}
		if (0 == nFNs) {
			int32_t strange = 1 ;
			}
		}

	return 0 ;
}


int32_t ARE::ARP::DestroyAdjVarList(void)
{
	_nSingletonVariables = -1 ;
	if (NULL != _StaticVarTotalList) {
		delete [] _StaticVarTotalList ;
		_StaticVarTotalList = NULL ;
		}
	_StaticVarTotalListSize = 0 ;
	if (NULL != _Degree) {
		delete [] _Degree ;
		_Degree = NULL ;
		}
	if (NULL != _AdjVars) {
		delete [] _AdjVars ;
		_AdjVars = NULL ;
		}
	return 0 ;
}


int32_t ARE::ARP::ComputeAdjVarList(void)
{
	if (0 != DestroyAdjVarList()) 
		return 1 ;
	if (_nVars < 1) 
		return 0 ;

	_nSingletonVariables = 0 ;

	// compute upper bound on the size of _StaticVarTotalList
	int32_t i, j, k, UBvarlistsize = 0 ;
	for (i = 0 ; i < _nFunctions ; i++) {
		Function *f = _Functions[i] ;
		if (NULL == f) continue ;
		j = f->N() * (f->N() - 1) ;
		UBvarlistsize += j ;
		}

	_Degree = new int32_t[_nVars] ;
	_AdjVars = new int32_t[_nVars] ;
	_StaticVarTotalList = UBvarlistsize > 0 ? new int32_t[UBvarlistsize] : NULL ;
	if (NULL == _Degree || NULL == _AdjVars || (UBvarlistsize > 0 && NULL == _StaticVarTotalList)) {
		DestroyAdjVarList() ;
		return 1 ;
		}
	_StaticVarTotalListSize = UBvarlistsize ;

	for (i = 0 ; i < _nVars ; i++) {
		_Degree[i] = 0 ;
		_AdjVars[i] = -1 ;
		}

	int32_t *vars = new int32_t[UBvarlistsize] ;
	if (NULL == vars) {
		DestroyAdjVarList() ;
		return 1 ;
		}
	int32_t left[32], right[32] ;
	int32_t idxStaticVarTotalList = 0 ;
	for (i = 0 ; i < _nVars ; i++) {
		int32_t n = 0 ;
		// collect all adj vars from adj FNs; this list may have duplicates.
		for (j = nAdjFunctions(i)-1 ; j >= 0 ; j--) {
			Function *f = AdjFunction(i, j) ;
			if (NULL == f) continue ;
			for (k = 0 ; k < f->N() ; k++) {
				int32_t v = f->Argument(k) ;
				if (v == i) 
					continue ;
				if (n >= UBvarlistsize) 
					{ DestroyAdjVarList() ; return ERRORCODE_VarDegreeTooLarge ; }
				vars[n++] = v ;
				}
			}
		if (0 == n) {
			++_nSingletonVariables ;
			continue ;
			}
		// eliminate duplicates
		QuickSortLong2(vars, n, left, right) ;
		k = 0 ;
		for (j = 1 ; j < n ; j++) {
			if (vars[k] == vars[j]) 
				continue ;
			vars[++k] = vars[j] ;
			}
		// fill in adj var list
		_Degree[i] = ++k ;
		_AdjVars[i] = idxStaticVarTotalList ;
		for (j = 0 ; j < k ; j++, idxStaticVarTotalList++) 
			_StaticVarTotalList[idxStaticVarTotalList] = vars[j] ;
		}
	delete [] vars ;

	return 0 ;
}


int32_t ARE::ARP::ComputeConnectedComponents(void)
{
	if (_nVars <= 0) 
		return 0 ;

	int32_t i, j ;

	_nConnectedComponents = 0 ;
	int32_t *temp = new int32_t[3*_nVars] ;
	if (NULL == temp) 
		return ERRORCODE_memory_allocation_failure ;
	for (i = 0 ; i < _nVars ; i++) temp[i] = 0 ;
	for (i = 0 ; i < _nVars ; i++) {
		if (0 != temp[i]) continue ;
		// variable i belongs to a new component; do DFS of the component starting from i.
		temp[i] = ++_nConnectedComponents ;
		int32_t *nAdj = temp + _nVars ; // for each variable, number of ajacent variables we have traversed
		int32_t *vStack = temp + 2*_nVars ; // variable stack
		vStack[0] = i ;
		nAdj[0] = 0 ;
		j = 0 ; // stacksize - 1
		while (j >= 0) {
			int32_t v = vStack[j] ;
			if (nAdj[j] >= _Degree[v]) { --j ; continue ; }
			int32_t u = AdjVar(v, nAdj[j]++) ;
			if (0 != temp[u]) continue ;
			temp[u] = _nConnectedComponents ;
			vStack[++j] = u ;
			nAdj[j] = 0 ;
			}
		}
	delete [] temp ;

	return 0 ;
}

/*
void ARE::ARP::DestroyGraph(void)
{
	if (NULL != _GraphAdjacencyMatrix) {
		delete [] _GraphAdjacencyMatrix ;
		_GraphAdjacencyMatrix = NULL ;
		}
	_nConnectedComponents = _nSingletonVariables = 0 ;
}


int32_t ARE::ARP::ComputeGraph(void)
{
	DestroyGraph() ;
	if (_nVars < 1) 
		return 0 ;

	int32_t i, j ;

	{

	_GraphAdjacencyMatrix = new ARPGraphNode[_nVars] ;
	if (NULL == _GraphAdjacencyMatrix) goto failed ;
	// process one variable at a time
	for (i = 0 ; i < _nVars ; i++) {
		int32_t n = 0 ;
		// compute degree of i
		for (j = 0 ; j < _nFunctions ; j++) {
			ARE::Function *f = _Functions[j] ;
			if (NULL == f) continue ;
			if (f->ContainsVariable(i)) 
				n++ ;
			}
		// allocate function list
		ARPGraphNode & node = _GraphAdjacencyMatrix[i] ;
		if (0 != node.SetN(n)) 
			goto failed ;
		node.SetV(i) ;
		n = 0 ;
		// fill in function list
		for (j = 0 ; j < _nFunctions ; j++) {
			ARE::Function *f = _Functions[j] ;
			if (NULL == f) continue ;
			if (f->ContainsVariable(i)) 
				node.SetFunction(n++, f) ;
			}
		// finish
		node.ComputeAdjacentVariables() ;
		}

	// run test of nodes
	for (i = 0 ; i < _nVars ; i++) {
		ARPGraphNode & node = _GraphAdjacencyMatrix[i] ;
		if (node.test() > 0) {
			int32_t error = 1 ;
			}
		}

	// compute num of singleton variables
	for (i = 0 ; i < _nVars ; i++) {
		ARPGraphNode & node = _GraphAdjacencyMatrix[i] ;
		if (node.Degree() < 0) {
			int32_t error = 1 ;
			}
		else if (0 == node.Degree()) 
			_nSingletonVariables++ ;
		}

	// compute number of connected components
	int32_t *temp = new int32_t[3*_nVars] ;
	if (NULL == temp) goto failed ;
	for (i = 0 ; i < _nVars ; i++) temp[i] = 0 ;
	for (i = 0 ; i < _nVars ; i++) {
		if (0 != temp[i]) continue ;
		// variable i belongs to a new component; do DFS of the component starting from i.
		temp[i] = ++_nConnectedComponents ;
		int32_t *nAdj = temp+_nVars ; // for each variable, number of ajacent variables we have traversed
		int32_t *vStack = temp+2*_nVars ; // variable stack
		vStack[0] = i ;
		nAdj[0] = 0 ;
		j = 0 ; // stacksize - 1
		while (j >= 0) {
			int32_t v = vStack[j] ;
			ARPGraphNode & node = _GraphAdjacencyMatrix[v] ;
			if (nAdj[j] >= node.Degree()) { --j ; continue ; }
			int32_t u = node.AdjacentVariable(nAdj[j]++) ;
			if (temp[u]) continue ;
			temp[u] = _nConnectedComponents ;
			vStack[++j] = u ;
			nAdj[j] = 0 ;
			}
		}
	delete [] temp ;

	}

	return 0 ;
failed :
	DestroyGraph() ;
	return 1 ;
}
*/

void ARE::ARP::DestroyVarOrdering(void)
{
	if (NULL != _VarOrdering_VarList) {
		delete [] _VarOrdering_VarList ;
		_VarOrdering_VarList = NULL ;
		}
	if (NULL != _VarOrdering_VarPos) {
		delete [] _VarOrdering_VarPos ;
		_VarOrdering_VarPos = NULL ;
		}
	_VarOrdering_InducedWidth = -1 ;
}


int32_t ARE::ARP::GetVarElimOrdering(std::vector<int32_t> & Order, int32_t & induced_width)
{
	Order.clear() ;
	induced_width = -1 ;
	if (NULL == _VarOrdering_VarList) 
		return 0 ;
	for (int32_t i = _nVars - 1 ; i >= 0 ; i--) 
		Order.push_back(_VarOrdering_VarList[i]) ;
	induced_width = _VarOrdering_InducedWidth ;
	return 0 ;
}


int32_t ARE::ARP::SetVarElimOrdering(const int32_t *VarListInElimOrder, int32_t induced_width)
{
	if (_nVars <= 0) {
		DestroyVarOrdering() ;
		_VarOrdering_InducedWidth = 0 ;
		return 0 ;
		}
	if (NULL == _VarOrdering_VarList) 
		_VarOrdering_VarList = new int32_t[_nVars] ;
	if (NULL == _VarOrdering_VarPos) 
		_VarOrdering_VarPos = new int32_t[_nVars] ;
	if (NULL == _VarOrdering_VarList || NULL == _VarOrdering_VarPos) {
		DestroyVarOrdering() ;
		return 1 ;
		}
	for (int32_t i = 0 ; i < _nVars ; i++) {
		_VarOrdering_VarList[i] = _VarOrdering_VarPos[i] = -1 ;
		}
	for (int32_t i = 0, j = _nVars-1 ; i < _nVars ; i++, j--) {
		if (VarListInElimOrder[j] < 0 || VarListInElimOrder[j] >= _nVars) 
			return ERRORCODE_InvalidInputData ;
		_VarOrdering_VarList[i] = VarListInElimOrder[j] ;
		_VarOrdering_VarPos[_VarOrdering_VarList[i]] = i ;
		}
	for (int32_t i = 0 ; i < _nVars ; i++) {
		if (_VarOrdering_VarList[i] < 0 || _VarOrdering_VarList[i] >= _nVars) 
			{ DestroyVarOrdering() ; return ERRORCODE_InvalidInputData ; }
		}
	_VarOrdering_InducedWidth = induced_width ;
	return 0 ;
}


int32_t ARE::ARP::SetVarBTOrdering(const int32_t *VarListInBTOrder, int32_t induced_width)
{
	if (_nVars <= 0) {
		DestroyVarOrdering() ;
		return 0 ;
		}
	if (NULL == _VarOrdering_VarList) 
		_VarOrdering_VarList = new int32_t[_nVars] ;
	if (NULL == _VarOrdering_VarPos) 
		_VarOrdering_VarPos = new int32_t[_nVars] ;
	if (NULL == _VarOrdering_VarList || NULL == _VarOrdering_VarPos) {
		DestroyVarOrdering() ;
		return 1 ;
		}
	for (int32_t i = 0 ; i < _nVars ; i++) {
		_VarOrdering_VarList[i] = _VarOrdering_VarPos[i] = -1 ;
		}
	for (int32_t i = 0 ; i < _nVars ; i++) {
		if (VarListInBTOrder[i] < 0 || VarListInBTOrder[i] >= _nVars) 
			return ERRORCODE_InvalidInputData ;
		_VarOrdering_VarList[i] = VarListInBTOrder[i] ;
		_VarOrdering_VarPos[_VarOrdering_VarList[i]] = i ;
		}
	for (int32_t i = 0 ; i < _nVars ; i++) {
		if (_VarOrdering_VarList[i] < 0 || _VarOrdering_VarList[i] >= _nVars) 
			{ DestroyVarOrdering() ; return ERRORCODE_InvalidInputData ; }
		}
	_VarOrdering_InducedWidth = induced_width ;
	return 0 ;
}


int32_t ARE::ARP::LoadVariableOrderingFromBuffer(int32_t OrderType, const char *SerializedVarListInElimOrder)
{
	if (_nVars <= 0) 
		return 0 ;
	const char *s = SerializedVarListInElimOrder ;
	// find '{'
	for (; 0 != *s && '{' != *s ; s++) ;
	if (0 == *s) 
		return ERRORCODE_InvalidInputData ;
	int32_t n = 0 ;
	int32_t *vars = new int32_t[_nVars] ;
	if (NULL == vars) 
		return ERRORCODE_memory_allocation_failure ;
	const char *sB = ++s ;
	char stemp[32] ;
	while (0 != *s && '}' != *s) {
		if (';' != *s) { ++s ; continue ; }
		int32_t l = s - sB ;
		if (l > 0 && n < _nVars && l < 32) {
			memcpy(stemp, sB, l) ;
			stemp[l] = 0 ;
			vars[n++] = atoi(stemp) ;
			}
		sB = ++s ;
		}
	int32_t l = s - sB ;
	if (l > 0 && n < _nVars && l < 32) {
		memcpy(stemp, sB, l) ;
		stemp[l] = 0 ;
		vars[n++] = atoi(stemp) ;
		}
	if (n != _nVars) 
		{ delete [] vars ; return ERRORCODE_InvalidInputData ; }
	int32_t res = (1 == OrderType) ? SetVarBTOrdering(vars, -1) : SetVarElimOrdering(vars, -1) ;
	delete [] vars ;
	return res ;
}


int32_t ARE::ARP::TestVariableOrdering(const int32_t *VarList, const int32_t *Var2PosMap)
{
	int32_t i, n, m ;

	// check no var/pos is out of bounds
	for (i = 0 ; i < _nVars ; i++) {
		if (VarList[i] < 0 || VarList[i] >= _nVars) 
			return 1 ;
		if (Var2PosMap[i] < 0 || Var2PosMap[i] >= _nVars) 
			return 2 ;
		}

	// check sum of VarList/Var2PosMap; should be (N choose 2) = N*(N-1)/2
	n = m = 0 ;
	for (i = 0 ; i < _nVars ; i++) {
		n += VarList[i] ;
		m += Var2PosMap[i] ;
		}
	int32_t test = (_nVars*(_nVars-1)) >> 1 ;
	if (n != test) 
		return 3 ;
	if (m != test) 
		return 3 ;

	return 0 ;
}


int32_t ARE::ARP::RemoveAdjVar(int32_t U, int32_t V)
{
	int32_t *adjlist = _StaticVarTotalList + _AdjVars[U] ;
	for (int32_t i = _Degree[U]-1 ; i >= 0 ; i--) {
		if (V == adjlist[i]) {
			adjlist[i] = adjlist[--_Degree[U]] ;
			break ;
			}
		}
	return 0 ;
}


int32_t ARE::ARP::EliminateEvidence(void)
{
	int32_t i, nEvidVars = 0 ;
	for (i = 0 ; i < N() ; i++) {
		int32_t val = Value(i) ;
		if (val < 0) 
			continue ;
		EliminateEvidenceVariable(i, val) ;
		++nEvidVars ;
		}
	return 0 ;
}


int32_t ARE::ARP::EliminateEvidenceVariable(int32_t Var, int32_t Val)
{
	int32_t i ;

	// remove this var from all adj FNs
	for (i = nAdjFunctions(Var) - 1 ; i >= 0 ; i--) {
		ARE::Function *f = AdjFunction(Var, i) ;
		if (NULL == f) continue ;
		f->RemoveVariable(Var, Val) ; // note that when this var has no more arguments, this fn will turn it into a const fn
		}
	// this var no longer participates in any FNs
	_nAdjFunctions[Var] = 0 ;
	_AdjFunctions[Var] = -1 ;

	// remove this var from the adj list of any var that this var was adj to
	for (i = Degree(Var) - 1 ; i >= 0 ; i--) {
		int32_t v = AdjVar(Var, i) ;
		RemoveAdjVar(v, Var) ;
		}
	_Degree[Var] = 0 ;
	_AdjVars[Var] = -1 ;

	return 0 ;
}


int32_t ARE::ARP::EliminateSingletonDomainVariables(void)
{
	int32_t i, j ;

	_nSingletonDomainVariables = 0 ;
	for (i = 0 ; i < _nVars ; i++) {
		if (_K[i] > 1) continue ;
		++_nSingletonDomainVariables ;
		for (j = 0 ; j < nAdjFunctions(i) ; j++) {
			ARE::Function *f = AdjFunction(i, j) ;
			if (NULL == f) continue ;
			// since it is singleton domain, we don't have to do anything (e.g. change fn table), 
			// just remove the variable from the scope.
#ifdef _DEBUG
			if (NULL != ARE::fpLOG) {
				fprintf(ARE::fpLOG, "\np.EliminateSingletonDomainVariables(): removing var %d from fn %d ...", i, (int32_t) f->IDX()) ;
				fflush(ARE::fpLOG) ;
				}
#endif
			f->RemoveVariableFromArguments(i) ;
			// if the fn has no arguments, make it a const function
			if (0 == f->N()) {
#ifdef _DEBUG
				if (NULL != ARE::fpLOG) {
					fprintf(ARE::fpLOG, "\np.EliminateSingletonDomainVariables(): fn %d has no arguments left ...", (int32_t) f->IDX()) ;
					fflush(ARE::fpLOG) ;
					}
#endif
				if (NULL != f->TableData()) 
					f->ConstValue() = f->TableEntry(0) ;
				else 
					f->ConstValue() = 1.0 ;
				}
//			f->RemoveEvidenceVariable(i, _K[i]) ;
			}
		// this var no longer participates in any FNs
		_nAdjFunctions[i] = 0 ;
		_AdjFunctions[i] = -1 ;
		// remove this var from the adj list of any var that this var was adj to
		for (j = Degree(i) - 1 ; j >= 0 ; j--) {
			int32_t v = AdjVar(i, j) ;
			RemoveAdjVar(v, i) ;
			}
		_Degree[i] = 0 ;
		_AdjVars[i] = -1 ;
		}

	return 0 ;
}

/*
int32_t ARE::ARP::ComputeMinDegreeOrdering(void)
{
	DestroyMinDegreeOrdering() ;
	if (_nVars < 1) 
		return 0 ;
	if (NULL == _GraphAdjacencyMatrix) 
		return 1 ;

	int32_t i, j, k ;

	// we will keep track, for each variable, which variables are its neighbors.
	// given i and j (i < j) will keep track of if i and j are neighbors.
	// adj info for i is kept as a list of adj variables.
	// head of the list is stored in '_MinDegreeOrdering_VarPos' -> if starts at m, then we store -(m+2); 
	// if _MinDegreeOrdering_VarPos[n] >= 0, then that var is already in the ordering.
//	int32_t adj_space_size = 30000 ;
//	int32_t adj_space_size = 16384 ;
#define ORD_COMP_ADJ_SPACE_SIZE 32768
	int32_t adj_space_size = ORD_COMP_ADJ_SPACE_SIZE ;
	// adj space should be counted in pairs :
	//		[i] is the value,
	//		[i+1] is a pointer to the next in the list.
	int32_t *degrees = new int32_t[_nVars] ;
	int32_t *adj_space = new int32_t[adj_space_size] ;
	if (NULL == degrees || NULL == adj_space) goto failed ;

	{

	_MinDegreeOrdering_VarList = new int32_t[_nVars] ;
	_MinDegreeOrdering_VarPos = new int32_t[_nVars] ;
	if (NULL == _MinDegreeOrdering_VarPos || NULL == _MinDegreeOrdering_VarList) goto failed ;

	// forget the induced width computed earlier
	_MinDegreeOrdering_InducedWidth = -1 ;

	// initialize the space
	int32_t empty_adj_space = 0 ;
	int32_t number_of_empty_cells = adj_space_size >> 1 ;
	for (i = 1 ; i < adj_space_size - 2 ; i += 2) 
		adj_space[i] = i + 1 ;
	// NIL at the end
	adj_space[adj_space_size - 1] = -1 ;

	// during each iteration, we will pick a variable that has the smallest degree

	// mark each variable as not in the ordering yet
	for (i = 0 ; i < _nVars ; i++) {
		_MinDegreeOrdering_VarPos[i] = -1 ;
		// use this array for counting the degree of unordered variables
		_MinDegreeOrdering_VarList[i] = 0 ;
		}

	// count degrees
	for (i = 0 ; i < _nVars ; i++) 
		degrees[i] = 0 ;

	// fill in initial graph
	for (j = 0 ; j < _nVars ; j++) {
		ARPGraphNode & node_j = _GraphAdjacencyMatrix[j] ;
		for (i = j+1 ; i < _nVars ; i++) {
			if (! node_j.IsAdjacent_QueryReleventFunctionsOnly(i)) continue ;
			// check space
			if (number_of_empty_cells < 2) {
				int32_t new_size = adj_space_size + ORD_COMP_ADJ_SPACE_SIZE ;
				int32_t *new_adj_space = new int32_t[new_size] ;
				if (NULL == new_adj_space) 
					goto failed ;
				memcpy(new_adj_space, adj_space, adj_space_size*sizeof(int32_t)) ;
				// generate structure for new space
				number_of_empty_cells += ORD_COMP_ADJ_SPACE_SIZE >> 1 ;
				for (k = adj_space_size + 1 ; k < new_size - 2 ; k += 2) {
					new_adj_space[k] = k + 1 ;
					}
				new_adj_space[new_size - 1] = empty_adj_space ;
				empty_adj_space = adj_space_size ;
				// let go of old space
				delete [] adj_space ;
				adj_space = new_adj_space ;
				adj_space_size = new_size ;
				}

			// add a node.
			// ***************************
			// get a cell for i.
			// ***************************
			int32_t cell = empty_adj_space ;
			empty_adj_space = adj_space[empty_adj_space + 1] ;
			--number_of_empty_cells ;
			adj_space[cell] = j ;
			// check if i has nothing
			if (-1 == _MinDegreeOrdering_VarPos[i]) 
				adj_space[cell + 1] = -1 ;
			else 
				adj_space[cell + 1] = -(_MinDegreeOrdering_VarPos[i] + 2) ;
			_MinDegreeOrdering_VarPos[i] = -(cell + 2) ;
			// count degree
			++degrees[i] ;
			// ***************************
			// get a cell for j.
			// ***************************
			cell = empty_adj_space ;
			empty_adj_space = adj_space[empty_adj_space + 1] ;
			--number_of_empty_cells ;
			adj_space[cell] = i ;
			// check if j has nothing
			if (-1 == _MinDegreeOrdering_VarPos[j]) 
				adj_space[cell + 1] = -1 ;
			else 
				adj_space[cell + 1] = -(_MinDegreeOrdering_VarPos[j] + 2) ;
			_MinDegreeOrdering_VarPos[j] = -(cell + 2) ;
			// count degree
			++degrees[j] ;
			}
		}

	// find the max degree
	int32_t max_degree = -INT_MAX ;

	// we will fill in the ordering from the end.
	// from the end, we will put min-degree nodes, except degree 0,
	// which go in the beginning.
	int32_t end = _nVars - 1 ;
	while (end >= 0) {
		// among variables that are not yet in the ordering, find one that has the smallest degree.
		int32_t min_degree_value = INT_MAX ;
		int32_t min_degree_var = -1 ;
		for (j = 0 ; j < _nVars ; j++) {
			if (_MinDegreeOrdering_VarPos[j] >= 0) continue ; // variable already in ordering
			if (degrees[j] < min_degree_value) {
				min_degree_value = degrees[j] ;
				min_degree_var = j ;
				}
			}
		// check if nothing left
		if (min_degree_var < 0) 
			break ; // no variable found
		if (max_degree < min_degree_value) 
			max_degree = min_degree_value ;

		// remove all references from others to this
		for (j = 0 ; j < _nVars ; j++) {
			if (min_degree_var == j) continue ;
			if (_MinDegreeOrdering_VarPos[j] >= 0) continue ; // variable already in ordering
			int32_t previous = -1 ;
			for (i = -_MinDegreeOrdering_VarPos[j] - 2 ; -1 != i ; i = adj_space[i+1]) {
				if (min_degree_var == adj_space[i]) {
					// remove min_degree_var from list
					if (-1 == previous)
						_MinDegreeOrdering_VarPos[j] = -(adj_space[i+1] + 2) ;
					else
						adj_space[previous+1] = adj_space[i+1] ;
					adj_space[i+1] = empty_adj_space ;
					empty_adj_space = i ;
					++number_of_empty_cells ;
					// reduce degree
					--degrees[j] ;
					break ;
					}
				previous = i ;
				}
			}

		// connect neighbors
		for (j = -_MinDegreeOrdering_VarPos[min_degree_var] - 2 ; -1 != j ; j = adj_space[j+1]) {
			for (i = adj_space[j+1] ; -1 != i ; i = adj_space[i+1]) {
				int32_t u = adj_space[i] ;
				int32_t v = adj_space[j] ;
				// check if u and v are already connected
				for (k = -_MinDegreeOrdering_VarPos[u] - 2 ; -1 != k ; k = adj_space[k+1]) {
					if (adj_space[k] == v) break ;
					}
				if (k >= 0) continue ; // already connected

				// check space
				if (number_of_empty_cells < 2) {
					int32_t new_size = adj_space_size + ORD_COMP_ADJ_SPACE_SIZE ;
					int32_t *new_adj_space = new int32_t[new_size] ;
					if (NULL == new_adj_space) goto failed ;
					memcpy(new_adj_space, adj_space, adj_space_size*sizeof(int32_t)) ;
					// generate structure for new space
					number_of_empty_cells += ORD_COMP_ADJ_SPACE_SIZE >> 1 ;
					for (k = adj_space_size + 1 ; k < new_size - 2 ; k += 2) {
						new_adj_space[k] = k + 1 ;
						}
					new_adj_space[new_size - 1] = empty_adj_space ;
					empty_adj_space = adj_space_size ;
					// let go of old space
					delete [] adj_space ;
					adj_space = new_adj_space ;
					adj_space_size = new_size ;
					}

				// ***************************
				// get a cell for u.
				// ***************************
				int32_t cell = empty_adj_space ;
				empty_adj_space = adj_space[empty_adj_space + 1] ;
				--number_of_empty_cells ;
				adj_space[cell] = v ;
				// check if u has nothing
				if (-1 == _MinDegreeOrdering_VarPos[u]) 
					adj_space[cell + 1] = -1 ;
				else 
					adj_space[cell + 1] = -_MinDegreeOrdering_VarPos[u] - 2 ;
				_MinDegreeOrdering_VarPos[u] = -(cell + 2) ;
				// count degree
				++degrees[u] ;
				// ***************************
				// get a cell for v.
				// ***************************
				cell = empty_adj_space ;
				empty_adj_space = adj_space[empty_adj_space + 1] ;
				--number_of_empty_cells ;
				adj_space[cell] = u ;
				// check if v has nothing
				if (-1 == _MinDegreeOrdering_VarPos[v]) 
					adj_space[cell + 1] = -1 ;
				else 
					adj_space[cell + 1] = -_MinDegreeOrdering_VarPos[v] - 2 ;
				_MinDegreeOrdering_VarPos[v] = -(cell + 2) ;
				// count degree
				++degrees[v] ;
				}
			}

		// release space from this var
		int32_t first = -_MinDegreeOrdering_VarPos[min_degree_var] - 2 ;
		int32_t last = -1 ;
		int32_t count = 0 ;
		for (j = first ; -1 != j ; j = adj_space[j+1]) {
			last = j ;
			++count ;
			}
		if (last >= 0) {
			adj_space[last+1] = empty_adj_space ;
			empty_adj_space = first ;
			number_of_empty_cells += count ;
			}

		// min_degree_var is the next variable in the ordering
		_MinDegreeOrdering_VarPos[min_degree_var] = end ;
		_MinDegreeOrdering_VarList[end] = min_degree_var ;
		--end ;
		}

	delete [] degrees ;
	delete [] adj_space ;

	ComputeInducedWidth(_MinDegreeOrdering_VarList, _MinDegreeOrdering_VarPos, _MinDegreeOrdering_InducedWidth) ;

#ifdef _DEBUG
	if (NULL != ARE::fpLOG) {
		char s[256] ;
		sprintf(s, "\nMin-Degree ordering width = %d varlist : ", (int32_t) _MinDegreeOrdering_InducedWidth) ;
		fwrite(s, 1, strlen(s), ARE::fpLOG) ;
		for (j = 0 ; j < _nVars ; j++) {
			sprintf(s, " %d", (int32_t) _MinDegreeOrdering_VarList[j]) ;
			fwrite(s, 1, strlen(s), ARE::fpLOG) ;
			}
		}
#endif // _DEBUG

	return 0 ;

	}

failed :
	if (NULL != degrees) 
		delete [] degrees ;
	if (NULL != adj_space) 
		delete [] adj_space ;
	DestroyMinDegreeOrdering() ;
	return 1 ;
}


int32_t ARE::ARP::ComputeMinFillOrdering(int32_t nRandomTries, int32_t BestKnownWidth)
{
	if (_nVars < 1) 
		return 0 ;
	if (nRandomTries < 2) {
		ComputeMinFillOrdering() ;
		goto done ;
		}

	// run n tries; save the best
	{
	int32_t *VarList = new int32_t[_nVars] ;
	if (NULL == VarList) 
		return 1 ;
	int32_t w = INT_MAX ;
	for (int32_t i = 0 ; i < nRandomTries ; i++) {
		if (0 != ComputeMinFillOrdering()) 
			continue ;
		if (_MinFillOrdering_InducedWidth < w) {
			w = _MinFillOrdering_InducedWidth ;
			for (int32_t j = 0 ; j < _nVars ; j++) VarList[j] = _MinFillOrdering_VarList[j] ;
			if (w <= BestKnownWidth) 
				break ;
			}
		}
	if (w < _MinFillOrdering_InducedWidth) {
		for (int32_t j = 0 ; j < _nVars ; j++) {
			_MinFillOrdering_VarList[j] = VarList[j] ;
			_MinFillOrdering_VarPos[_MinFillOrdering_VarList[j]] = j ;
			}
		_MinFillOrdering_InducedWidth = w ;
		}
	delete [] VarList ;
	}

done :
#ifdef _DEBUG
	if (NULL != ARE::fpLOG) {
		char s[256] ;
		sprintf(s, "\nMin-Fill   ordering width = %d varlist : ", (int32_t) _MinFillOrdering_InducedWidth) ;
		fwrite(s, 1, strlen(s), ARE::fpLOG) ;
		for (int32_t j = 0 ; j < _nVars ; j++) {
			sprintf(s, " %d", (int32_t) _MinFillOrdering_VarList[j]) ;
			fwrite(s, 1, strlen(s), ARE::fpLOG) ;
			}
		}
#endif // _DEBUG

	return 0 ;
}


int32_t ARE::ARP::ComputeMinFillOrdering(void)
{
	DestroyMinFillOrdering() ;
	if (_nVars < 1) 
		return 0 ;
	if (NULL == _GraphAdjacencyMatrix) 
		return 1 ;

	int32_t i, j, k, l ;

	// we will keep track, for each variable, which variables are its neighbors.
	// given i and j (i < j) will keep track of if i and j are neighbors.
	// adj info for i is kept as a list of adj variables.
	// head of the list is stored in '_MinFillOrdering_VarPos' -> if starts at m, then we store -(m+2); 
	// if _MinFillOrdering_VarPos[n] >= 0, then that var is already in the ordering.
//	int32_t adj_space_size = 30000 ;
//	int32_t adj_space_size = 16384 ;
#define ORD_COMP_ADJ_SPACE_SIZE 32768
	int32_t adj_space_size = ORD_COMP_ADJ_SPACE_SIZE ;
	// adj space should be counted in pairs :
	//		[i] is the value,
	//		[i+1] is a pointer to the next in the list.
	int32_t *degrees = new int32_t[_nVars] ;
	int32_t *adj_space = new int32_t[adj_space_size] ;
	if (NULL == degrees || NULL == adj_space) goto failed ;

	{

	_MinFillOrdering_VarList = new int32_t[_nVars] ;
	_MinFillOrdering_VarPos = new int32_t[_nVars] ;
	if (NULL == _MinFillOrdering_VarPos || NULL == _MinFillOrdering_VarList) goto failed ;

	// forget the induced width computed earlier
	_MinFillOrdering_InducedWidth = -1 ;

	// during each iteration, we will pick a variable that has the smallest number of fill-in edges

	// mark each variable as not in the ordering yet
	for (i = 0 ; i < _nVars ; i++) {
		_MinFillOrdering_VarPos[i] = -1 ;
		// use this array for counting the fill-in number of unordered variables
		_MinFillOrdering_VarList[i] = 0 ;
		}

	// initialize the space
	int32_t empty_adj_space = 0 ;
	int32_t number_of_empty_cells = adj_space_size >> 1 ;
	for (i = 1 ; i < adj_space_size - 2 ; i += 2) 
		adj_space[i] = i + 1 ;
	// NIL at the end
	adj_space[adj_space_size - 1] = -1 ;

	// count degrees
	for (i = 0 ; i < _nVars ; i++) 
		degrees[i] = 0 ;

	// fill in initial graph
	for (j = 0 ; j < _nVars ; j++) {
		ARPGraphNode & node_j = _GraphAdjacencyMatrix[j] ;
		for (i = j+1 ; i < _nVars ; i++) {
			if (! node_j.IsAdjacent_QueryReleventFunctionsOnly(i)) continue ;
			// check space
			if (number_of_empty_cells < 2) {
				int32_t new_size = adj_space_size + ORD_COMP_ADJ_SPACE_SIZE ;
				int32_t *new_adj_space = new int32_t[new_size] ;
				if (NULL == new_adj_space) 
					goto failed ;
				memcpy(new_adj_space, adj_space, adj_space_size*sizeof(int32_t)) ;
				// generate structure for new space
				number_of_empty_cells += ORD_COMP_ADJ_SPACE_SIZE >> 1 ;
				for (k = adj_space_size + 1 ; k < new_size - 2 ; k += 2) {
					new_adj_space[k] = k + 1 ;
					}
				new_adj_space[new_size - 1] = empty_adj_space ;
				empty_adj_space = adj_space_size ;
				// let go of old space
				delete [] adj_space ;
				adj_space = new_adj_space ;
				adj_space_size = new_size ;
				}

			// add a node.
			// ***************************
			// get a cell for i.
			// ***************************
			int32_t cell = empty_adj_space ;
			empty_adj_space = adj_space[empty_adj_space + 1] ;
			--number_of_empty_cells ;
			adj_space[cell] = j ;
			// check if i has nothing
			if (-1 == _MinFillOrdering_VarPos[i]) 
				adj_space[cell + 1] = -1 ;
			else 
				adj_space[cell + 1] = -(_MinFillOrdering_VarPos[i] + 2) ;
			_MinFillOrdering_VarPos[i] = -(cell + 2) ;
			// count degree
			++degrees[i] ;
			// ***************************
			// get a cell for j.
			// ***************************
			cell = empty_adj_space ;
			empty_adj_space = adj_space[empty_adj_space + 1] ;
			--number_of_empty_cells ;
			adj_space[cell] = i ;
			// check if j has nothing
			if (-1 == _MinFillOrdering_VarPos[j]) 
				adj_space[cell + 1] = -1 ;
			else 
				adj_space[cell + 1] = -(_MinFillOrdering_VarPos[j] + 2) ;
			_MinFillOrdering_VarPos[j] = -(cell + 2) ;
			// count degree
			++degrees[j] ;
			}
		}

	// we will fill in the ordering from the end.
	// from the end, we will put min-fillin nodes.
	int32_t end = _nVars - 1 ;
	while (end >= 0) {
		// among variables that are not yet in the ordering, find one that has the smallest fillin.
		int32_t min_fillin_value = INT_MAX ;
		int32_t min_fillin_var = -1 ;
		for (j = 0 ; j < _nVars ; j++) {
			if (_MinFillOrdering_VarPos[j] >= 0) continue ; // variable already in ordering
			// compute for all pairs of adjacent nodes, how many are not connected; this is the fillin value
			int32_t fillinvalue = 0 ;
			if (degrees[j] > 1) {
//				fillinvalue = (degrees[j]*(degrees[j]-1)) >> 1 ;
				for (i = -_MinFillOrdering_VarPos[j] - 2 ; -1 != i ; i = adj_space[i+1]) {
					for (k = adj_space[i+1] ; -1 != k ; k = adj_space[k+1]) {
						int32_t u = adj_space[i] ;
						int32_t v = adj_space[k] ;
						// u and v are two variables adjacent to j; check if u and v are already connected.
						for (l = -_MinFillOrdering_VarPos[u] - 2 ; -1 != l ; l = adj_space[l+1]) {
							if (adj_space[l] == v) break ;
							}
						if (l < 0) 
							// not connected
							fillinvalue++ ;
						}
					}
				}
			if (fillinvalue < min_fillin_value) {
				min_fillin_value = fillinvalue ;
				min_fillin_var = j ;
				}
			// if equal, break tries randomly
			else if (fillinvalue == min_fillin_value) {
				if (RNG.rand() > 0.5) {
					min_fillin_value = fillinvalue ;
					min_fillin_var = j ;
					}
				}
			}
		// check if nothing left
		if (min_fillin_var < 0) 
			break ; // no variable found

		// remove all references from others to this
		for (j = 0 ; j < _nVars ; j++) {
			if (min_fillin_var == j) continue ;
			if (_MinFillOrdering_VarPos[j] >= 0) continue ; // variable already in ordering
			int32_t previous = -1 ;
			for (i = -_MinFillOrdering_VarPos[j] - 2 ; -1 != i ; i = adj_space[i+1]) {
				if (min_fillin_var == adj_space[i]) {
					// remove min_fillin_var from list
					if (-1 == previous)
						_MinFillOrdering_VarPos[j] = -(adj_space[i+1] + 2) ;
					else
						adj_space[previous+1] = adj_space[i+1] ;
					adj_space[i+1] = empty_adj_space ;
					empty_adj_space = i ;
					++number_of_empty_cells ;
					// reduce degree
					--degrees[j] ;
					break ;
					}
				previous = i ;
				}
			}

		// connect neighbors
		for (j = -_MinFillOrdering_VarPos[min_fillin_var] - 2 ; -1 != j ; j = adj_space[j+1]) {
			for (i = adj_space[j+1] ; -1 != i ; i = adj_space[i+1]) {
				int32_t u = adj_space[i] ;
				int32_t v = adj_space[j] ;
				// check if u and v are already connected
				for (k = -_MinFillOrdering_VarPos[u] - 2 ; -1 != k ; k = adj_space[k+1]) {
					if (adj_space[k] == v) break ;
					}
				if (k >= 0) continue ; // already connected

				// check space
				if (number_of_empty_cells < 2) {
					int32_t new_size = adj_space_size + ORD_COMP_ADJ_SPACE_SIZE ;
					int32_t *new_adj_space = new int32_t[new_size] ;
					if (NULL == new_adj_space) goto failed ;
					memcpy(new_adj_space, adj_space, adj_space_size*sizeof(int32_t)) ;
					// generate structure for new space
					number_of_empty_cells += ORD_COMP_ADJ_SPACE_SIZE >> 1 ;
					for (k = adj_space_size + 1 ; k < new_size - 2 ; k += 2) {
						new_adj_space[k] = k + 1 ;
						}
					new_adj_space[new_size - 1] = empty_adj_space ;
					empty_adj_space = adj_space_size ;
					// let go of old space
					delete [] adj_space ;
					adj_space = new_adj_space ;
					adj_space_size = new_size ;
					}

				// ***************************
				// get a cell for u.
				// ***************************
				int32_t cell = empty_adj_space ;
				empty_adj_space = adj_space[empty_adj_space + 1] ;
				--number_of_empty_cells ;
				adj_space[cell] = v ;
				// check if u has nothing
				if (-1 == _MinFillOrdering_VarPos[u]) 
					adj_space[cell + 1] = -1 ;
				else 
					adj_space[cell + 1] = -_MinFillOrdering_VarPos[u] - 2 ;
				_MinFillOrdering_VarPos[u] = -(cell + 2) ;
				// count degree
				++degrees[u] ;
				// ***************************
				// get a cell for v.
				// ***************************
				cell = empty_adj_space ;
				empty_adj_space = adj_space[empty_adj_space + 1] ;
				--number_of_empty_cells ;
				adj_space[cell] = u ;
				// check if v has nothing
				if (-1 == _MinFillOrdering_VarPos[v]) 
					adj_space[cell + 1] = -1 ;
				else 
					adj_space[cell + 1] = -_MinFillOrdering_VarPos[v] - 2 ;
				_MinFillOrdering_VarPos[v] = -(cell + 2) ;
				// count degree
				++degrees[v] ;
				}
			}

		// release space from this var
		int32_t first = -_MinFillOrdering_VarPos[min_fillin_var] - 2 ;
		int32_t last = -1 ;
		int32_t count = 0 ;
		for (j = first ; -1 != j ; j = adj_space[j+1]) {
			last = j ;
			++count ;
			}
		if (last >= 0) {
			adj_space[last+1] = empty_adj_space ;
			empty_adj_space = first ;
			number_of_empty_cells += count ;
			}

		// min_fillin_var is the next variable in the ordering
		_MinFillOrdering_VarPos[min_fillin_var] = end ;
		_MinFillOrdering_VarList[end] = min_fillin_var ;
		--end ;
		}

	delete [] degrees ;
	delete [] adj_space ;

	ComputeInducedWidth(_MinFillOrdering_VarList, _MinFillOrdering_VarPos, _MinFillOrdering_InducedWidth) ;

	return 0 ;

	}

failed :
	DestroyMinFillOrdering() ;
	return 1 ;
}
*/

static int32_t mark_two_variables_as_adj_in_induced_width_computation(int32_t i, int32_t j, int32_t **adj_matrix, int32_t *adj_list_length, int32_t adj_list_initial_length) 
{
	int32_t k, l ;

	// check if adj list exists.
	if (NULL == adj_matrix[i]) {
		// adj list does not exist. allocate memory, mark j as adjacent to i.
		adj_matrix[i] = new int32_t[adj_list_initial_length] ;
		if (NULL == adj_matrix[i]) return 1 ;
		adj_list_length[i] = adj_list_initial_length ;
		(adj_matrix[i])[0] = j ;
		for (k = 1 ; k < adj_list_length[i] ; k++) (adj_matrix[i])[k] = -1 ;
		}
	else {
		// adj list does exist.
		// check if j is already there.
		for (k = 0 ; k < adj_list_length[i] ; k++) {
			if ((adj_matrix[i])[k] < 0) break ;
			if ((adj_matrix[i])[k] == j) return 0 ;
			}
		// j not there
		if (k >= adj_list_length[i]) {
			// out of space
			int32_t new_length = adj_list_length[i] + adj_list_initial_length ;
			int32_t *temp_adj = new int32_t[new_length] ;
			if (NULL == temp_adj) return 1 ;
			// copy existing list
			for (l = 0 ; l < adj_list_length[i] ; l++) temp_adj[l] = (adj_matrix[i])[l] ;
			for (; l < new_length ; l++) temp_adj[l] = -1 ;
			// free old memory
			delete [] adj_matrix[i] ;
			// set new length
			adj_list_length[i] = new_length ;
			// attach new memory
			adj_matrix[i] = temp_adj ;
			}
		// store adj info
		(adj_matrix[i])[k] = j ;
		}

	return 0 ;
}


int32_t ARE::ARP::ComputeInducedWidth(const int32_t *VarList, const int32_t *Var2PosMap, int32_t & InducedWidth)
{
	InducedWidth = -1 ;

	// check that the number of variables is positive
	if (_nVars < 1) 
		return 0 ;
	// check that an ordering is given
	if (NULL == VarList || NULL == Var2PosMap) 
		return 1 ;

	int32_t return_value = 1 ;
	int32_t i, j, k ;

	// we will use this variable to allocate memory to an adjacency list, and
	// if memory runs out, reallocate.
	const int32_t adj_list_initial_length = 32 ;
	// for every variable, we will keep track of the adjacent variables
	int32_t **adj_matrix = NULL ;
	// for every variable, we will keep track of the length of the adjacency list
	int32_t *adj_list_length = NULL ;
	// this is the maximum width we are computing
	int32_t max_width = 0 ;
	int32_t max_width_variable = -1 ;

	// allocate memory
	adj_matrix = new int32_t*[_nVars] ;
	if (NULL == adj_matrix) goto done ;
	for (i = 0 ; i < _nVars ; i++) adj_matrix[i] = NULL ;
	adj_list_length = new int32_t[_nVars] ;
	if (NULL == adj_list_length) goto done ;
	for (i = 0 ; i < _nVars ; i++) adj_list_length[i] = 0 ;

	// ***************************************************************************************************
	// for every variable, compute its initial adjacency list, based on the given constraint graph.
	// once this is done, we don't need the original constraint graph any more and we won't reference it.
	// ***************************************************************************************************
	for (i = 0 ; i < _nVars ; i++) {
		for (k = 0 ; k < nAdjFunctions(i) ; k++) {
			Function *f = this->AdjFunction(i, k) ;
			if (NULL == f) 
				continue ;
			if (f->IsQueryIrrelevant()) 
				continue ;
			for (j = 0 ; j < f->N() ; j++) {
				if (j == i) continue ;
				// now we need to add this bit of information that i and j are adjacent to the
				// adjcency lists. but first we need to check if there is enough memory.
				// do this for both variables
				// for variable i.
				if (mark_two_variables_as_adj_in_induced_width_computation(i, j, adj_matrix, adj_list_length, adj_list_initial_length)) goto done ;
				// for variable j.
				if (mark_two_variables_as_adj_in_induced_width_computation(j, i, adj_matrix, adj_list_length, adj_list_initial_length)) goto done ;
				}
			}
		}

	// ***************************************************************************************************
	// go backwards in the ordering, and for every variable, connect all of its adjacent variables
	// that come before it in the ordering.
	// at the same time, keep track of the width of every variable.
	// ***************************************************************************************************
	for (i = _nVars - 1 ; i >= 0 ; i--) {
		// get the name of the variable
		int32_t v = VarList[i] ;
		// check if this variable has an adj list
		if (NULL == adj_matrix[v]) continue ;

		// compute the width of this variable
		int32_t width = 0 ;
		for (j = 0 ; j < adj_list_length[v] ; j++) {
			int32_t vr = (adj_matrix[v])[j] ;
			if (vr < 0) break ;
			// check if this variable comes in the ordering before v
			if (Var2PosMap[vr] < i) width++ ;
			}
		// store the width, if maximum
		if (max_width < width) {
			max_width = width ;
			max_width_variable = v ;
			}

		// for any two parents of this variable, such that both are before it in the ordering, connect them.
		for (j = 0 ; j < adj_list_length[v] ; j++) {
			if ((adj_matrix[v])[j] < 0) break ;
			// check if this variable comes in the ordering before v
			if (Var2PosMap[(adj_matrix[v])[j]] >= i) continue ;
			for (int32_t k = 0 ; k < adj_list_length[v] ; k++) {
				if ((adj_matrix[v])[k] < 0) break ;
				if (k == j) continue ;
				// check if this variable comes in the ordering before v
				if (Var2PosMap[(adj_matrix[v])[k]] >= i) continue ;

				// connect these two variables.
				// NOTE : we only need to add adjacency information to the variable that comes
				// second in the ordering
				int32_t v1 = (adj_matrix[v])[j] ;
				int32_t v2 = (adj_matrix[v])[k] ;
				if (Var2PosMap[v1] < Var2PosMap[v2]) {
					// add for v2
					if (mark_two_variables_as_adj_in_induced_width_computation(v2, v1, adj_matrix, adj_list_length, adj_list_initial_length)) goto done ;
					}
				else {
					// add for v1
					if (mark_two_variables_as_adj_in_induced_width_computation(v1, v2, adj_matrix, adj_list_length, adj_list_initial_length)) goto done ;
					}
				}
			}

		// delete space owned by this variable. it is no longer needed.
		if (adj_matrix[v]) {
			delete [] adj_matrix[v] ;
			adj_matrix[v] = NULL ;
			}
		}

	// store max width
	InducedWidth = max_width ;

	// successfully done
	return_value = 0 ;

done :

	// free memory
	if (adj_matrix) {
		for (i = 0 ; i < _nVars ; i++) {
			if (adj_matrix[i]) delete [] adj_matrix[i] ;
			}
		delete [] adj_matrix ;
		}
	if (adj_list_length) delete [] adj_list_length ;

	return return_value ;
}


int32_t ARE::ARP::ComputeSingletonConsistency(int32_t & nNewSingletonDomainVariables)
{
	nNewSingletonDomainVariables = 0 ;

	if (NULL != ARE::fpLOG) {
		fprintf(ARE::fpLOG, "\nARP::ComputeSingletonConsistency() begin") ;
		}

	vector<vector<bool>> is_consistent ;
	SingletonConsistencyHelper(is_consistent) ;

	for (int32_t i = 0 ; i < _nVars ; i++) {
		int32_t old_k = _K[i] ;
		for (int32_t j = _K[i] - 1 ; j >= 0 ; j--) {
			if (is_consistent[i][j]) 
				continue ;
			// assignment variable[i] = value[j] is inconsistent with the rest of the problem; eliminate this value from the domain of the variable.
			for (int32_t k = 0 ; k < nAdjFunctions(i) ; k++) {
				ARE::Function *f = AdjFunction(i, k) ;
				if (NULL == f) continue ;
				f->RemoveVariableValue(i, j) ;
				if (NULL != ARE::fpLOG) {
					fprintf(ARE::fpLOG, "\nComputeSingletonConsistency(): removing var %d value %d f=%d", i, j, (int32_t) f->IDX()) ;
					fflush(ARE::fpLOG) ;
					}
				}
			// remove this value from the domain
			// check if domain size became 0
			if (--_K[i] < 1) 
				return -1 ;
			}
		if (old_k > 1 && 1 == _K[i]) 
			nNewSingletonDomainVariables++ ;
		}

	if (NULL != ARE::fpLOG) {
		fprintf(ARE::fpLOG, "\nARP::ComputeSingletonConsistency() end") ;
		}

	return 0 ;
}


void ARE::ARP::SingletonConsistencyHelper(vector<vector<bool> > & is_consistent)
{
	int32_t i, j, k ;

	// Constructing the CNF for singleton consistency
	// We have \sum_{i=1}^{N} #values(X_i) Boolean variables, namely a variable for each variable-value pair
	
	// minisat solver
	Solver S ;

	// The vector stores mapping from variable-value pairs to Boolean variables
	vector<vector<int32_t> > var2cnfvar(_nVars) ;
	int32_t count = 0 ;
	for (i = 0 ; i < _nVars ; i++) {
		var2cnfvar[i] = vector<int32_t>(_K[i]) ;
		for (j = 0 ; j < _K[i] ; j++) 
			var2cnfvar[i][j] = count++ ;
		}
	// Initialize the variables of minisat
	for (i = 0 ; i < count ; i++) 
		S.newVar() ;

	// Each solution of the CNF corresponds to a unique consistent assignment to the variables.
	// Add clauses to ensure that each variable takes exactly one value (i.e. at least and at most one)
	for (i = 0 ; i < _nVars ; i++) {
		vec<Lit> lits ;
		// For each variable, we have a clause which ensures that each variable takes at least one value
		for (j = 0 ; j < _K[i] ; j++) 
			lits.push(mkLit(var2cnfvar[i][j])) ;
		S.addClause(lits) ;
		// We have K[i]*(K[i]-1)/2 clauses which ensure that each variable takes at most one value
		for (j = 0 ; j < _K[i] ; j++) {
			for (k = j + 1 ; k < _K[i] ; k++) {
				lits.clear() ;
				lits.push(~mkLit(var2cnfvar[i][j])) ;
				lits.push(~mkLit(var2cnfvar[i][k])) ;
				S.addClause(lits) ;
				}
			}
		}

	// For each zero probability in each function, add a clause
	int32_t value_assignment[MAX_NUM_ARGUMENTS_PER_FUNCTION] ;
	for (i = 0 ; i < _nFunctions ; i++) {
		ARE::Function *f = _Functions[i] ;
		if (NULL == f) continue ;
		ARE_Function_TableType *data = f->TableData() ;
		vec<Lit> lits ;
		for (j = 0 ; j < f->TableSize() ; j++) {
			if (f->TableEntry(j) <= 0.000001) {
				lits.clear() ;
				// convert table address j to an assignment to the scope
				ComputeArgCombinationFromFnTableAdr(j, f->N(), f->Arguments(), value_assignment, _K) ;
				// Add a clause which is a negation of the current tuple assignment
				// E.g. If the tuple is (X=2,Y=3), then add a clause !X2 V !Y3 where X2 and Y3 are Boolean variables
				// corresponding to X=2 and Y=3 respectively
				// This will ensure that X=2 and Y=3 never occur together in a solution
				for (k = 0 ; k < f->N() ; k++) 
					lits.push(~mkLit(var2cnfvar[f->Arguments()[k]][value_assignment[k]])) ;
				S.addClause(lits) ;
				}
			}
		}

	// solve for each variable-value pair, by forcing the variable to have the given value
	is_consistent = vector<vector<bool> >(_nVars) ;
	vec<Lit> assumps ;
	for (i = 0 ; i < _nVars ; i++) {
		is_consistent[i] = vector<bool>(_K[i]) ;
		for (j = 0 ; j < _K[i] ; j++) {
			assumps.clear() ;
			assumps.push(mkLit(var2cnfvar[i][j])) ;
			if (S.solve(assumps))
				is_consistent[i][j] = true ;
			}
		}

	// done
}

