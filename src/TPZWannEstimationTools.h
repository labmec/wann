#pragma once

#include <TPZMultiphysicsCompMesh.h>
#include "ProblemData.h"

class TPZWannEstimationTools {

protected:
    static constexpr REAL estimator_tol = 1e-3;

public:
  static void EstimateAndRefine(TPZMultiphysicsCompMesh* cmeshHdiv, TPZCompMesh* cmeshH1, ProblemData* SimData, int nthreads = 1);
  static TPZVec<int64_t> ErrorEstimation(TPZMultiphysicsCompMesh* cmeshHdiv, TPZCompMesh* cmeshH1, ProblemData* SimData, int nthreads = 1);
  static TPZVec<int64_t> ErrorEstimationOld(TPZMultiphysicsCompMesh* cmeshHdiv, TPZCompMesh* cmeshH1, ProblemData* SimData);
  static void CheckRef(TPZCompMesh* cmesh);
  static void FakeRefine2(TPZGeoMesh* gmesh, ProblemData* SimData);
  // Maybe something between? MarkElements?
  
private:
  static REAL ForcingFunctionWellbore(TPZCompEl* celMixed, const TPZManVector<REAL,3>& pt);
  static void hRefinement(TPZGeoMesh* gmesh, TPZVec<int64_t>& needRefinement, ProblemData* SimData);
  static void FakeRefine(TPZGeoMesh* gmesh, TPZVec<int64_t>& needRefinement, ProblemData* SimData);
  static REAL ElementDiameter(TPZGeoEl* gel);
};