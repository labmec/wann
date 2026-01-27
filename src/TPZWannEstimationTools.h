#pragma once

#include <TPZMultiphysicsCompMesh.h>
#include "ProblemData.h"

class TPZWannEstimationTools {

protected:
    static constexpr REAL estimator_tol = 1e-3;
    static constexpr REAL relative_estimator_tol = 1e-3;

public:
  static void EstimateAndRefine(TPZMultiphysicsCompMesh* cmeshHdiv, TPZCompMesh* cmeshH1, ProblemData* SimData, int nthreads = 1);
  static REAL PragerSynge(TPZMultiphysicsCompMesh* cmeshHdiv, TPZCompMesh* cmeshH1, ProblemData* SimData, TPZVec<REAL>& elementErrors, int nthreads = 0);
  static REAL ErrorEstimation(TPZMultiphysicsCompMesh* cmeshHdiv, TPZCompMesh* cmeshH1, ProblemData* SimData, TPZVec<REAL>& elementErrors, int nthreads = 0);
  static void MarkElementsForRefinement(const TPZVec<REAL>& elementErrors, TPZVec<int>& refinementIndicator, REAL tol);
  static void MeshSmoothing(TPZGeoMesh* gmesh, TPZVec<int>& RefinementIndicator);
  static void MeshWellCompatibility(TPZGeoMesh* gmesh, TPZVec<int>& RefinementIndicator, ProblemData* SimData);
  static void hRefinement(TPZGeoMesh* gmesh, TPZVec<int>& RefinementIndicator);
  static void RefineFromFile(TPZGeoMesh* og_gmesh, const std::string& filename);

  // Old or testing stuff
  static TPZVec<int> ErrorEstimationOld(TPZMultiphysicsCompMesh* cmeshHdiv, TPZCompMesh* cmeshH1, ProblemData* SimData, int nthreads = 0);
  static void CheckRef(TPZCompMesh* cmesh);
  static void FakeRefine(TPZGeoMesh* gmesh, ProblemData* SimData);

private:
  static REAL ForcingFunctionWellbore(TPZCompEl* celMixed, const TPZManVector<REAL,3>& pt);
  static REAL ElementDiameter(TPZGeoEl* gel);
};