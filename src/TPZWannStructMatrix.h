#pragma once

#if defined(PZ_USING_MKL)
#include "TPZSSpStructMatrix.h"
using TPZWannStructMatrix = TPZSSpStructMatrix<STATE>;
#elif defined(PZ_USING_MUMPS) || defined(USING_MUMPS)
#include "TPZSSpStructMatrixMumps.h"
using TPZWannStructMatrix = TPZSSpStructMatrixMumps<STATE>;
#else
#include "pzskylstrmatrix.h"
using TPZWannStructMatrix = TPZSkylineStructMatrix<STATE>;
#endif
