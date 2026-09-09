#pragma once

#include <TPZMultiphysicsCompMesh.h>
#include "ProblemData.h"

class TPZWannAdaptivityTools {

protected:
    static constexpr REAL estimator_tol = 1e-3;
    static constexpr REAL relative_estimator_tol = 1e-3;

public:
  enum class MarkingStrategy {
    ESimple,
    EBulk
  };

  static REAL PragerSynge(TPZMultiphysicsCompMesh* cmeshHdiv, TPZCompMesh* cmeshH1, ProblemData* SimData, TPZVec<REAL>& elementErrors, int nthreads = 0); 
  static REAL GoalOriented(TPZMultiphysicsCompMesh* cmeshHdiv, TPZMultiphysicsCompMesh* cmeshDual, ProblemData* SimData, TPZVec<REAL>& elementErrors, int nthreads); 
  static void AdaptivityProcess(TPZGeoMesh* gmesh, ProblemData* SimData, TPZVec<REAL>& elementErrors, TPZVec<int64_t>& refinedElements, REAL theta, MarkingStrategy strategy = MarkingStrategy::EBulk);

  // --- Old or testing stuff ---
  static REAL ErrorEstimation(TPZMultiphysicsCompMesh* cmeshHdiv, TPZCompMesh* cmeshH1, ProblemData* SimData, TPZVec<REAL>& elementErrors, int nthreads = 0);

private:
  static REAL ForcingFunctionWellbore(TPZCompEl* celMixed, const TPZManVector<REAL,3>& pt);
  static REAL ElementDiameter(TPZGeoEl* gel);
  static REAL GoalContribution3D(TPZCompEl* celMixed, TPZCompEl* celDual, ProblemData* SimData); 
  static REAL GoalContributionInterface(TPZCompEl* celMixed, TPZCompEl* celDual, ProblemData* SimData);   
  static REAL GoalContribution1D(TPZCompEl* celMixed, TPZCompEl* celDual, ProblemData* SimData);
  static void MarkElementsForRefinement(const TPZVec<REAL>& elementErrors, TPZVec<int>& refinementIndicator, REAL theta, MarkingStrategy strategy = MarkingStrategy::EBulk);
  static void SimpleMarking(const TPZVec<REAL>& elementErrors, TPZVec<int>& refinementIndicator, REAL theta);
  static void BulkMarking(const TPZVec<REAL>& elementErrors, TPZVec<int>& refinementIndicator, REAL theta);
  static void MeshSmoothing(TPZGeoMesh* gmesh, TPZVec<int>& refinementIndicator);
  static void MeshWellCompatibility(TPZGeoMesh* gmesh, TPZVec<int>& refinementIndicator, ProblemData* SimData);
};