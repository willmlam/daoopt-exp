// daoopt_lookahead.cpp : Defines the entry point for the console application.
//
/*
daoopt_lookahead C:\UCI\linkage\pedigree30.uai C:\UCI\problems\0.evid 9
daoopt_lookahead_p1 C:\UCI\linkage\pedigree1.uai C:\UCI\problems\0.evid 5
*/

//#include "stdafx.h"
#include <fstream>


#include <Problem/Problem.hxx>
#include <CVO/VariableOrderComputation.hxx>

#include "Main.h"
#include "DaooptInterface.h"
#include "MiniBucketElimLH.h"

#ifdef LINUX
#include <unistd.h>
#endif

/*
daoopt_lookahead C:\UCI\linkage\pedigree23.uai C:\UCI\problems\0.evid C:\UCI\lookahead\experiments\pedigree23_var_elim_order.txt 6 3 0 logG1-6-3.txt
*/

INT64 nd1SpecialCalls = 0 ;
INT64 nd1GeneralCalls = 0 ;

// Hack to read a file into a string.
std::string getFileContents(const char* filename)
{
	std::ifstream in(filename, std::ios::in | std::ios::binary);
	std::string contents;
	in.seekg(0, std::ios::end);
	contents.resize(in.tellg());
	in.seekg(0, std::ios::beg);
	in.read(&contents[0], contents.size());
	in.close();
	return(contents);
}

int main(int argc, char** argv)
//int _tmain(int argc, _TCHAR* argv[])
{
/* double x, y, z = -1.0 ;
x = ELEM_ZERO ;
y = ELEM_ZERO ;
if (x < y) {
	int here = 1 ;
	}
if (x == y) {
	int here = 1 ;
	}
if (x < z) {
	int here = 1 ;
	}
double zz = z - x ;*/

	int nParams = argc ;
	if (10 != nParams) {
		printf("\nUSAGE:\n $ %s uaifile evidencefile varelimorderfile ibound lookaheadDepth BEtablesizetotal BEsingletablesize BuckerErrorIgnoreThreshold logfile", argv[0]) ;
		return 0 ;
		}
	std::string uaifile = (NULL != argv[1] ? argv[1] : "") ;
	std::string evidencefile = (NULL != argv[2] ? argv[2] : "") ;
	std::string varelimorderfile = (NULL != argv[3] ? argv[3] : "") ;
	std::string ibound = (NULL != argv[4] ? argv[4] : "") ;
	std::string sLookaheadDepth = (NULL != argv[5] ? argv[5] : "") ;
	std::string sBEtablesizetotal = (NULL != argv[6] ? argv[6] : "") ;
	std::string sBEsingletablesize = (NULL != argv[7] ? argv[7] : "") ;
	std::string sBuckerErrorIgnoreThreshold = (NULL != argv[8] ? argv[8] : "") ;
	std::string sLogFile = (NULL != argv[9] ? argv[9] : "") ; if (0 == sLogFile.length()) sLogFile = "daoopt_lookahead.txt" ;

	// Load problem specification and evidence into memory
	std::string problem = getFileContents(uaifile.c_str());
	std::string evidence = getFileContents(evidencefile.c_str());

	int maxNumProcessorThreads = 1 ;
#if defined _WINDOWS || WINDOWS
	{
		SYSTEM_INFO sysinfo ;
		GetSystemInfo(&sysinfo) ;
		maxNumProcessorThreads = sysinfo.dwNumberOfProcessors ;
	}
#endif

	daoopt::ProgramOptions options; // Set various options here, as documented in ProgramOptions.h
	options.ibound = ibound.length() > 0 ? atoi(ibound.c_str()) : 2 ;
  options.cbound = 1000;
	options.lookaheadDepth = sLookaheadDepth.length() > 0 ? atoi(sLookaheadDepth.c_str()) : 0 ;
	options.order_iterations = 0;
	options.seed = 2323;
	options.lds = -1;
	options.match = true;
  options.subprobOrder = 0;
  options.jglp = 5;
	options.lookahead_LE_AllTablesTotalLimit = sBEtablesizetotal.length() > 0 ? atof(sBEtablesizetotal.c_str()) : -DBL_MIN ;
	if (options.lookahead_LE_AllTablesTotalLimit < 0.0) options.lookahead_LE_AllTablesTotalLimit = -DBL_MIN ;
	else if (options.lookahead_LE_AllTablesTotalLimit > 20.0) options.lookahead_LE_AllTablesTotalLimit = 20.0 ;
    options.lookahead_LE_SingleTableLimit = sBEsingletablesize.length() > 0 ? atof(sBEsingletablesize.c_str()) : 0.0;
//	options.lookahead_LE_SingleTableLimit = options.lookahead_LE_AllTablesTotalLimit - 0.5 ;
	// even if tables are not generated, sample tables
	if (options.lookahead_LE_SingleTableLimit < 5.0) options.lookahead_LE_SingleTableLimit = 5.0 ;
options.lookahead_LE_SingleTableLimit = 6.0 ;
	options.lookahead_LE_IgnoreThreshold = sBuckerErrorIgnoreThreshold.length() > 0 ? atof(sBuckerErrorIgnoreThreshold.c_str()) : DBL_MIN ;
	if (options.lookahead_LE_IgnoreThreshold < DBL_MIN) options.lookahead_LE_IgnoreThreshold = DBL_MIN ;

	// store problem/evid pointers in options object.
	options.problemSpec = &problem[0];  // char* (non-const)
	options.problemSpec_len = problem.size();  // size_t
	options.evidSpec = &evidence[0];  // char*
	options.evidSpec_len = evidence.size();  // size_t
	options._fpLogFile = fopen(sLogFile.c_str(), "w") ;

	ARE::VarElimOrderComp::Order cvoBestOrder ; 
	if (varelimorderfile.length() > 0 ? 0 != stricmp("NULL", varelimorderfile.c_str()) : false) {
		FILE *fp = fopen(varelimorderfile.c_str(), "r") ;
		if (NULL != fp) 
			fclose(fp) ;
		else 
			varelimorderfile.erase() ;
		}
	if (0 == varelimorderfile.length() ? true : 0 == stricmp("NULL", varelimorderfile.c_str())) {
		if (NULL != options._fpLogFile && NULL != uaifile.c_str()) {
			fprintf(options._fpLogFile, "\nWILL COMPUTE VAR ELIM ORDER...") ;
			}
		printf("\nWILL COMPUTE VAR ELIM ORDER...") ;

		ARE::ARP *p = new ARE::ARP("DAOOPT_LA_TEST") ;
		int res_load = p->LoadFromBuffer("UAI", problem.c_str(), problem.length()) ;
		int res_pca = p->PerformPostConstructionAnalysis() ;
		if (evidence.length() > 0) {
			int nEvidenceVars = 0 ;
			int res_evid = p->LoadFromBuffer_Evidence("UAI", evidence.c_str(), evidence.length(), nEvidenceVars) ;
			int res_elim_evid = p->EliminateEvidence() ;
			}

		// initialize CVO stuff
		FILE *fpCVOLOG = fopen("cvo_log.txt", "w") ;
		ARE::VarElimOrderComp::CVOcontext cvoContext ;
		cvoContext._fpLOG = fpCVOLOG ;
		cvoContext._Problem = p ;
		cvoContext._BestOrder = &cvoBestOrder ;
		cvoContext._FindPracticalVariableOrder = true ;
		cvoContext._AlgCode = ARE::VarElimOrderComp::MinFill ;
		cvoContext._ObjCode = ARE::VarElimOrderComp::Width ;
		cvoContext._nRunsToDoMin = 1000 ;
		cvoContext._nRunsToDoMax = 1000000 ;
		cvoContext._TimeLimitInMilliSeconds = 3600000 ;
		cvoContext._nRandomPick = 8 ;
		cvoContext._eRandomPick = 0.5 ;
		cvoContext._EarlyTerminationOfBasic_W = true ;
		cvoContext._EarlyTerminationOfBasic_C = false ;
		cvoContext._nThreads = maxNumProcessorThreads > 1 ? maxNumProcessorThreads-1 : 1 ;

		// run CVO
		cvoContext.CreateCVOthread() ;
		if (0 == cvoContext._ThreadHandle) {
			return 1 ;
			}
		while (true) {
#if defined WINDOWS || _WINDOWS
			Sleep(50) ;
#else
      usleep(50000);
#endif
			if (0 == cvoContext._ThreadHandle) 
				break ;
			}
		cvoContext.StopCVOthread() ;
		if (cvoBestOrder._Width < 0 && cvoBestOrder._Width >= INT_MAX) {
			if (NULL != options._fpLogFile && NULL != uaifile.c_str()) {
				fprintf(options._fpLogFile, "\n\nFAILED TO COMPUTE VAR ELIM ORDER...") ;
				}
			printf("\n\nFAILED TO COMPUTE VAR ELIM ORDER...") ;
			return 2 ;
			}
		if (NULL != options._fpLogFile && NULL != uaifile.c_str()) {
			fprintf(options._fpLogFile, "\nVAR ELIM ORDER : w=%d", (int) cvoBestOrder._Width) ;
			}
		printf("\nVAR ELIM ORDER : w=%d", (int) cvoBestOrder._Width) ;

		int res_copyorder = p->SetVarElimOrdering(cvoBestOrder._VarListInElimOrder, cvoBestOrder._Width) ;

		// save order
		std::string vofn(uaifile) ;
		std::string::size_type ext_pos = vofn.rfind('.') ;
		vofn.erase(ext_pos) ;
		vofn += "_var_elim_order.txt" ;
		cvoBestOrder.SerializeAsElimOrder(vofn.c_str()) ;
		varelimorderfile = vofn ;

		if (NULL != options._fpLogFile && NULL != uaifile.c_str()) {
			fprintf(options._fpLogFile, "problem = %s", uaifile.c_str()) ;
			fprintf(options._fpLogFile, "\n N=%d", (int) p->N()) ;
			fprintf(options._fpLogFile, "\n variable elimination order (%s): \nint elim_order[%d]={", 
				varelimorderfile.c_str(), (int) cvoBestOrder._nVars) ;
			for (int i = 0 ; i < cvoBestOrder._nVars ; i++) {
				if (0 == i) 
					fprintf(options._fpLogFile, "%d", (int) cvoBestOrder._VarListInElimOrder[i]) ;
				else 
					fprintf(options._fpLogFile, ",%d", (int) cvoBestOrder._VarListInElimOrder[i]) ;
				}
			fprintf(options._fpLogFile, "} ;") ;
			}
		}

/* pedigree1
int elim_order[334]={22,307,20,305,306,114,304,113,112,111,70,69,68,67,49,23,50,21,0,302,300,298,296,294,292,290,288,286,284,275,273,271,269,267,265,263,261,259,250,248,246,244,242,240,238,236,234,232,230,201,200,199,139,137,135,133,131,129,127,125,123,121,119,118,117,115,116,79,76,75,32,65,59,57,55,53,51,40,38,36,34,63,30,28,26,24,18,16,14,12,10,8,6,4,2,33,71,72,25,77,78,46,80,81,82,62,74,73,54,41,45,29,64,61,58,37,42,19,186,301,66,303,122,108,52,1,3,283,285,27,120,17,150,151,153,152,154,155,189,190,162,149,148,260,258,157,156,262,192,191,256,194,193,84,227,197,198,281,145,231,196,195,110,48,43,233,182,188,229,164,205,86,83,287,158,159,160,39,147,308,309,11,95,47,203,293,9,321,322,166,92,5,280,100,245,130,170,132,243,138,124,310,311,324,323,268,274,255,295,228,88,276,217,31,176,237,251,209,219,178,141,207,252,326,325,102,226,104,142,56,291,235,312,313,146,94,222,87,91,60,266,173,254,282,44,327,328,289,278,257,253,224,277,7,330,329,168,211,316,317,136,126,315,314,107,103,264,239,96,35,128,99,241,213,270,299,165,225,143,184,15,333,318,332,331,320,319,106,215,297,175,90,161,169,171,185,177,181,101,105,144,140,109,85,89,93,98,13,174,134,97,223,247,272,221,249,279,208,202,179,206,220,218,216,214,212,210,187,183,180,172,167,163,204} ;
if (334 == cvoBestOrder._nVars) {
for (int ii = 0 ; ii < cvoBestOrder._nVars ; ii++) {
	cvoBestOrder._VarListInElimOrder[ii] = elim_order[ii] ;
	}}*/
/* pedigree18
int elim_order[1184]={14,1066,1065,704,691,683,679,675,663,659,651,587,576,573,569,577,570,579,566,563,561,554,550,542,551,536,528,517,529,504,500,496,493,489,497,490,501,486,505,479,476,472,474,468,464,469,457,465,446,433,425,421,418,408,422,400,426,396,428,392,430,384,434,385,370,358,350,346,342,330,326,318,273,272,271,270,269,266,256,252,249,253,247,246,237,236,234,230,225,223,222,208,207,203,202,204,199,189,183,180,176,178,174,170,165,163,139,138,124,118,108,107,105,101,67,54,46,42,38,26,22,0,206,461,483,558,584,1063,1061,1059,1057,1053,1051,1047,1045,1043,1041,1039,1037,1035,1033,1029,1027,1025,1023,1021,1019,1017,1015,1013,1011,1009,1007,1005,1003,1001,999,997,995,989,987,985,983,978,974,972,970,968,966,964,962,960,958,956,954,952,950,948,946,944,942,940,938,936,934,932,930,928,926,924,922,920,918,916,835,834,833,832,831,785,830,828,827,826,825,772,824,822,821,764,820,818,817,756,816,814,813,812,811,744,810,808,807,806,805,732,804,802,801,800,799,798,797,787,829,780,774,823,766,819,760,762,758,815,746,809,740,742,734,803,706,699,693,689,687,685,681,677,673,671,665,661,653,649,647,645,643,637,635,590,588,585,580,574,571,567,564,559,557,555,552,548,546,544,540,534,532,530,526,520,518,515,510,508,506,502,498,494,491,487,484,482,480,477,470,466,462,460,458,455,449,447,444,439,437,435,431,423,419,416,414,412,406,404,402,398,394,390,388,386,382,372,356,354,352,348,344,332,320,296,267,263,261,259,257,254,250,243,241,239,232,228,226,219,217,215,213,211,209,205,200,196,194,192,190,187,185,181,172,168,166,161,159,157,155,153,151,149,147,141,136,131,127,125,122,120,116,114,112,110,103,99,97,93,91,89,87,85,83,81,79,73,71,69,62,60,58,56,52,50,48,44,40,36,34,32,30,28,24,20,18,16,12,10,6,4,2,47,15,8,9,27,102,43,55,96,119,128,106,109,23,68,39,164,175,95,171,184,130,140,198,193,177,224,173,104,221,231,129,245,238,260,265,280,281,287,286,295,294,233,235,179,319,284,285,248,371,327,401,351,347,393,331,409,410,429,397,411,427,473,359,537,343,475,543,538,562,603,602,607,606,539,616,617,622,623,608,609,610,611,621,620,626,627,692,680,660,733,578,676,664,705,1055,757,761,684,765,741,773,794,745,717,786,652,1032,80,1150,1149,1050,115,78,142,1134,1135,976,1049,1137,1136,1031,977,82,1056,1,144,162,998,3,994,169,229,1124,1125,197,264,992,220,210,244,996,268,383,148,201,299,300,195,212,262,298,297,1096,1095,1098,1097,415,381,292,293,150,301,389,452,509,456,547,438,1086,1085,527,521,450,523,917,593,636,638,591,919,1062,715,915,113,1154,1153,135,913,1160,1159,92,275,274,277,276,133,191,1064,152,214,282,283,76,289,288,279,278,84,182,291,290,88,156,160,218,379,98,227,167,303,117,242,395,258,186,420,305,407,1166,1165,123,251,121,188,255,126,137,132,143,1167,1168,1170,1169,545,436,499,1179,1180,478,134,837,596,597,599,598,1171,1172,568,790,601,600,459,604,605,403,549,839,467,614,615,613,612,618,619,100,531,387,391,463,535,541,399,471,553,481,796,753,624,625,777,413,556,417,485,488,560,572,424,492,495,565,729,630,631,713,633,632,894,1100,1099,575,1139,1138,432,503,719,240,111,1152,1151,29,31,533,507,1105,1144,1143,1142,11,1104,1103,49,690,51,682,345,37,57,1141,1140,1028,1012,1022,721,90,70,1067,1068,1070,1069,17,1030,674,1006,662,1145,1146,1127,1126,969,353,899,686,1147,1148,767,873,644,1130,1131,1128,1129,1038,45,984,935,1000,921,1072,1071,1132,1133,1107,1106,1121,1120,1077,1078,957,1109,1108,1079,1080,897,1117,1116,1073,1074,313,1112,1113,688,1014,349,357,646,949,933,743,955,445,845,841,1082,1081,589,749,843,788,925,927,1173,1174,943,355,1092,1091,516,759,1046,986,53,654,146,694,852,864,669,882,77,595,362,583,875,582,451,514,443,522,907,1004,454,701,1054,1155,1156,361,891,975,967,1044,855,1176,1175,867,321,525,668,7,154,883,1123,1122,727,1020,1157,1158,629,628,333,857,751,941,64,889,311,512,441,1088,1087,1089,1090,782,703,1102,1101,881,335,380,373,909,779,1002,1181,1182,878,923,1111,1110,847,586,367,865,711,519,511,725,1118,1119,1084,1083,988,1058,979,634,440,448,25,905,581,592,666,442,513,216,892,1162,1161,1164,1163,338,902,848,735,329,965,872,858,910,914,322,1114,1115,648,368,65,896,314,642,309,863,888,906,376,453,524,594,860,377,947,886,334,1026,655,991,840,41,678,838,86,1076,1075,324,364,341,5,959,328,850,980,316,853,836,842,1040,312,304,308,846,723,844,35,639,887,75,769,63,366,1042,868,339,961,712,158,667,898,72,378,13,317,650,900,981,990,792,911,775,912,904,993,908,903,784,1060,982,901,310,716,700,708,709,74,374,1183,315,306,707,1177,145,1178,375,714,702,710,365,61,695,1093,885,1094,874,851,931,1010,307,302,854,731,1008,929,849,657,866,1016,937,736,862,861,739,856,1018,939,859,945,747,1024,870,672,369,66,698,697,360,724,728,720,21,658,641,640,325,94,789,776,793,781,340,363,59,336,696,768,752,748,337,19,737,670,656,405,323,33,871,755,951,869,754,718,750,738,730,726,722,953,770,1048,971,895,1034,1052,1036,973,963,893,890,884,880,879,795,791,783,778,771,763,876,877} ;
if (1184 == cvoBestOrder._nVars) {
for (int ii = 0 ; ii < cvoBestOrder._nVars ; ii++) {
	cvoBestOrder._VarListInElimOrder[ii] = elim_order[ii] ;
	}}
cvoBestOrder.SerializeAsElimOrder("pedigree18_var_elim_order.txt") ;*/
/* pedigree30
int elim_order[1289]={101,1171,1170,888,875,867,863,859,847,843,835,771,760,757,753,761,754,763,750,747,745,738,734,726,735,720,712,701,713,688,684,680,677,673,681,674,685,670,689,663,660,656,658,652,648,653,641,649,630,617,609,605,602,592,606,584,610,580,612,576,614,568,618,569,554,541,533,529,525,513,509,501,448,433,430,426,434,427,423,419,410,420,412,406,398,407,388,384,370,385,358,389,350,391,346,342,330,326,318,273,272,271,270,269,266,256,252,249,253,247,246,237,236,234,230,225,223,222,208,207,203,202,204,199,189,183,180,176,178,174,170,165,163,139,138,124,118,108,107,105,0,206,416,645,667,742,768,1168,1166,1164,1162,1158,1156,1152,1150,1148,1146,1144,1142,1140,1138,1134,1132,1130,1128,1126,1124,1122,1120,1118,1116,1114,1112,1110,1108,1106,1104,1102,1100,1019,1018,1017,1016,1015,969,1014,1012,1011,1010,1009,956,1008,1006,1005,948,1004,1002,1001,940,1000,998,997,996,995,928,994,992,991,990,989,916,988,986,985,984,983,982,981,971,1013,964,958,1007,950,1003,944,946,942,999,930,993,924,926,918,987,890,883,877,873,871,869,865,861,857,855,849,845,837,833,831,829,827,821,819,774,772,769,764,758,755,751,748,743,741,739,736,732,730,728,724,718,716,714,710,704,702,699,694,692,690,686,682,678,675,671,668,666,664,661,654,650,646,644,642,639,633,631,628,623,621,619,615,607,603,600,598,596,590,588,586,582,578,574,572,570,566,560,558,556,549,543,539,537,535,531,527,519,517,515,511,503,493,491,489,487,486,485,484,446,483,481,480,435,479,477,476,428,475,473,472,421,471,469,468,413,467,465,464,404,463,461,460,396,459,457,456,386,455,482,439,437,478,431,474,424,470,417,415,466,408,462,400,402,458,392,394,454,382,383,387,390,393,397,401,405,409,411,414,418,422,425,429,432,436,372,356,354,352,348,344,332,320,296,267,263,261,259,257,254,250,243,241,239,232,228,226,219,217,215,213,211,209,205,200,196,194,192,190,187,185,181,172,168,166,161,159,157,155,153,151,149,147,141,136,131,127,125,122,120,116,114,112,110,103,99,97,93,91,89,87,85,83,81,79,73,71,69,67,62,58,56,54,52,50,48,46,44,42,40,38,36,34,32,30,28,26,24,22,20,18,16,14,12,10,8,6,4,2,104,102,128,60,129,95,109,164,173,175,61,171,184,179,106,221,198,140,130,235,193,119,238,224,231,245,265,233,260,280,281,285,284,287,286,177,319,96,248,327,347,359,295,294,444,447,449,441,451,440,450,331,343,510,445,526,371,502,351,514,495,581,496,534,593,577,611,594,595,1161,613,399,585,659,1160,555,1154,727,657,542,723,762,746,1137,721,790,791,722,787,786,794,795,1136,792,793,806,807,811,810,1155,804,805,848,800,801,844,868,864,876,901,929,836,941,925,945,970,889,530,1229,1230,957,978,80,111,82,142,148,917,169,229,860,162,197,264,949,3,78,1,210,1254,1255,76,150,244,268,1200,1201,144,381,1202,1203,403,1242,1241,297,298,201,299,300,195,1239,1240,1190,1191,565,262,293,292,301,1167,212,488,490,634,567,599,693,1169,707,711,438,640,775,573,717,220,777,636,705,622,822,899,133,1097,191,135,182,1099,966,1101,258,274,275,277,276,291,290,214,820,1103,100,152,160,84,278,279,113,156,88,937,218,303,123,282,283,92,98,227,379,1259,1258,288,289,921,167,442,1265,1264,521,117,242,453,452,909,587,604,903,579,563,1271,1270,251,186,121,188,255,126,137,132,143,1021,1272,1273,1274,1275,134,679,1276,1277,305,1023,733,980,897,780,781,783,782,571,715,662,643,784,785,647,1284,1285,575,796,797,789,788,756,798,799,729,719,651,583,725,655,591,737,740,669,597,665,815,814,601,816,817,803,802,691,752,809,808,1204,1205,1249,683,1210,395,1208,1209,240,115,1243,1244,1076,616,1248,1247,15,676,744,672,749,608,687,846,731,974,17,620,27,41,33,49,47,39,68,55,1232,1231,1139,830,9,759,1079,1198,1199,90,858,70,1113,1225,1226,927,1214,1213,1147,25,1221,1222,832,1129,947,1215,1216,1212,1211,1119,1227,1228,838,1059,1081,700,773,629,905,766,866,333,874,345,1057,953,1183,1182,856,1195,1194,7,536,1145,1279,1278,1027,1235,1236,1107,355,1251,1250,373,357,1173,1172,1175,1174,520,518,1253,1252,538,544,1127,1153,528,919,972,1105,512,1082,353,1062,767,698,706,313,853,709,627,1070,638,146,870,635,321,779,349,532,1048,547,364,885,540,1087,813,812,911,625,696,939,1025,1176,1177,1245,1246,1089,887,878,57,1073,1067,5,770,361,1045,1029,1178,1179,1180,1181,1256,1257,1157,1206,1207,557,1109,1121,1031,1111,1217,1218,1219,1220,498,551,955,1043,11,1135,1286,1287,1095,494,311,154,1037,335,1065,1189,1188,74,1051,1197,1196,1192,1193,380,1035,1115,19,703,695,818,624,632,1224,1223,1165,1185,1184,765,776,1186,1187,626,697,505,564,75,1269,1268,963,504,367,322,336,216,1262,1263,895,1084,1053,850,1260,1261,1092,309,826,368,314,1233,1234,824,913,872,1088,1237,1238,1040,35,545,362,1044,1055,1078,37,497,851,840,1039,708,637,778,86,1061,862,1036,329,376,516,1091,1071,562,1032,66,550,377,1072,1066,559,1080,935,334,961,1028,1131,31,63,1096,316,308,304,312,1024,341,823,158,896,1022,378,892,893,884,1143,45,1149,51,374,366,552,561,1288,1281,1280,879,891,338,375,145,1283,1282,894,886,898,1090,553,65,77,1098,1086,522,443,369,881,1093,72,1094,976,307,64,324,1133,492,339,834,317,500,1085,968,1163,1083,908,13,915,1033,1117,328,977,960,973,965,907,1020,1042,1038,1034,1030,1026,21,923,1041,839,1060,1074,951,1054,1058,1075,59,1159,1077,959,1151,53,1069,1068,1064,1052,43,1141,1063,1056,1050,943,931,23,1125,1123,1049,1047,1046,29,882,360,548,365,880,310,302,306,315,828,825,546,523,363,340,952,936,900,499,912,904,94,325,920,842,841,508,507,979,938,589,975,1267,1266,934,933,932,922,914,910,906,902,854,852,524,506,337,323,967,954,962} ;
if (1289 == cvoBestOrder._nVars) {
for (int ii = 0 ; ii < cvoBestOrder._nVars ; ii++) {
	cvoBestOrder._VarListInElimOrder[ii] = elim_order[ii] ;
	}}*/

	if (NULL != options._fpLogFile && NULL != uaifile.c_str()) {
		fprintf(options._fpLogFile, "\n\nproblem = %s", uaifile.c_str()) ;
		fprintf(options._fpLogFile, "\n i-bound=%d", (int) options.ibound) ;
		fprintf(options._fpLogFile, "\n") ;
		}

	// fetch var elim order
	options.in_orderingFile = varelimorderfile ;
/*	static std::vector<int> varElimOrder ;
	varElimOrder.clear() ;
	int induced_width = -1 ;
	p->GetVarElimOrdering(varElimOrder, induced_width) ;
	options.varOrder = &varElimOrder ;*/

	printf("\nWILL INITIALIZE ...") ;

	// Example for using DaooptInterface
	daoopt::DaooptInterface daoopt ;
	daoopt.initialize() ;
	// Pass in options object created above.
	printf("\nWILL PREPROCESS ...") ;
	daoopt.preprocess(options) ;

  /*
	double Estimate_NumSearchSpaceNodesExpanded = daoopt.estimate() ;
	if (0.0 == Estimate_NumSearchSpaceNodesExpanded) { // solved by preprocessing
		double Query_Answer = daoopt.getSolution(NULL) ;
		printf("\nSOLVED by preprocessing... answer=%g", Query_Answer) ;
		}
    */

/*	if (options.lookaheadDepth > 0) {
		if (NULL != options._fpLogFile && NULL != uaifile.c_str()) {
			fprintf(options._fpLogFile, "\n\nWILL COMPUTE LOCAL ERROR for each bucket ... error=(exact value - MB value)\n") ;
			fflush(options._fpLogFile) ;
			}
		daoopt::MiniBucketElimLH *mbEx = dynamic_cast<daoopt::MiniBucketElimLH*>(daoopt.Heuristic()) ;
		if (NULL != mbEx) {
			double avgError = DBL_MAX, avgExact = DBL_MAX, stdDev = -DBL_MAX ;
			INT64 memory_limit = 1000000 ;
			mbEx->computeLocalErrorTable(memory_limit, true, avgError, avgExact, stdDev) ;
			}
		}
	fflush(options._fpLogFile) ;*/
/*	double avgError_H = DBL_MAX ;
	if (NULL != mbEx) {
		mbEx->localError_H(avgError_H, true) ;
		}
	fflush(options._fpLogFile) ;*/

	try {
		static bool runSearch = true ;
		if (runSearch) {
			printf("\nWILL SOLVE ...") ;
			size_t nodesdone = 0 ;
			size_t nodeLimit = 10000 ;
			time_t ttNow, ttLastLog = 0 ;
			while (! daoopt.solve(nodeLimit)) {
				double sol = daoopt.getSolution(NULL) ;
				time(&ttNow) ;
				if ((ttNow - ttLastLog) >= 5 && NULL != options._fpLogFile && NULL != uaifile.c_str()) {
					ttLastLog = ttNow ;
					nodesdone += nodeLimit ;
					time_t ttNow ; time(&ttNow) ; char strDT[32] ;
					struct tm *pTime = localtime(&ttNow) ;
					char *sDT = asctime(pTime) ; // from manual : The string result produced by asctime contains exactly 26 characters and has the form Wed Jan 02 02:03:55 1980\n\0.
					memcpy(strDT, sDT, 26) ;
					strDT[24] = 0 ;
					fprintf(options._fpLogFile, "\n Computing solution, %lld steps (done %lld), current solution is %g, ttNow=%s", (__int64) nodeLimit, (__int64) nodesdone, sol, strDT) ;
					fflush(options._fpLogFile) ;
					}
				}
			}
		}
	catch(...){
		printf("\n\nmain exception...\n") ;
		}

	daoopt.outputStatistics() ;

	if (NULL != options._fpLogFile && NULL != uaifile.c_str()) {
		fprintf(options._fpLogFile, "\n\nDONE : nG=%lld nS=%lld\n", nd1GeneralCalls, nd1SpecialCalls) ;
		fflush(options._fpLogFile) ;
		}

	return 0;
}

