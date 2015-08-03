#include "Factor.h"

#ifdef SPARSE_FACTORS
#include "FactorSparse.cpp"
#else
#include "FactorDense.cpp"
#endif

