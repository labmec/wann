#pragma once

#include <TPZMultiphysicsCompMesh.h>
#include "ProblemData.h"

class TPZWannEstimationTools {

protected:
    static constexpr REAL estimator_tol = 1e-3;
    static constexpr REAL relative_estimator_tol = 1e-3;

public:
  static void EstimateAndRefine(TPZMultiphysicsCompMesh* cmeshHdiv, TPZCompMesh* cmeshH1, ProblemData* SimData, int nthreads = 1);
  static TPZVec<int> ErrorEstimation(TPZMultiphysicsCompMesh* cmeshHdiv, TPZCompMesh* cmeshH1, ProblemData* SimData, int nthreads = 0);
  static TPZVec<int> ErrorEstimationOld(TPZMultiphysicsCompMesh* cmeshHdiv, TPZCompMesh* cmeshH1, ProblemData* SimData, int nthreads = 0);
  static void CheckRef(TPZCompMesh* cmesh);
  static void FakeRefine(TPZGeoMesh* gmesh, ProblemData* SimData);
  
private:
  static REAL ForcingFunctionWellbore(TPZCompEl* celMixed, const TPZManVector<REAL,3>& pt);
  static void hRefinement(TPZGeoMesh* gmesh, TPZVec<int>& RefinementIndicator, ProblemData* SimData);
  static REAL ElementDiameter(TPZGeoEl* gel);
  static void meshSmoothing(TPZGeoMesh* gmesh, TPZVec<int>& RefinementIndicator);
};