/**
 * @file
 * @brief Class and method contracts for TPZWannAdaptivityTools.
 */

#pragma once

#include <TPZMultiphysicsCompMesh.h>
#include "ProblemData.h"

/**
 * @brief Legacy error estimation and refinement experiments awaiting API migration.
 * Located in needRefactor and excluded from WANN_lib by the current source glob.
 * Several methods refer to material-ID fields absent from current ProblemData;
 * these declarations describe historical implementation, not a supported API.
 */
class TPZWannAdaptivityTools {

protected:
    static constexpr REAL estimator_tol = 1e-3;
    static constexpr REAL relative_estimator_tol = 1e-3;

public:
  /**
   * @brief Relative-maximum threshold (ESimple) or cumulative-sum marking (EBulk).
   */
  enum class MarkingStrategy {
    ESimple,
    EBulk
  };

  /**
   * @brief Legacy reservoir flux-discrepancy/balance estimator comparing mixed and H1
   * solutions on shared geometry. Fills geometrically indexed indicators and returns
   * their Euclidean norm; uses nthreads+1 workers and zero source.
   */
  static REAL PragerSynge(TPZMultiphysicsCompMesh* cmeshHdiv, TPZCompMesh* cmeshH1, ProblemData* SimData, TPZVec<REAL>& elementErrors, int nthreads = 0);
  /**
   * @brief Legacy primal/dual signed goal estimator with geometric-indexed absolute
   * local indicators. Uses nthreads+1 workers; well contributions are disabled.
   */
  static REAL GoalOriented(TPZMultiphysicsCompMesh* cmeshHdiv, TPZMultiphysicsCompMesh* cmeshDual, ProblemData* SimData, TPZVec<REAL>& elementErrors, int nthreads);
  /**
   * @brief Mark geometric indices, enforce well compatibility, smooth, recheck
   * compatibility, populate refinedElements, and divide the marked elements.
   * theta is the relative threshold or bulk fraction according to strategy.
   */
  static void AdaptivityProcess(TPZGeoMesh* gmesh, ProblemData* SimData, TPZVec<REAL>& elementErrors, TPZVec<int64_t>& refinedElements, REAL theta, MarkingStrategy strategy = MarkingStrategy::EBulk);

  // --- Old or testing stuff ---
  /**
   * @brief Historical mixed/H1 estimator with scalar permeability and well special case.
   * Fills geometric indicators and writes an error log; requires the old data API.
   */
  static REAL ErrorEstimation(TPZMultiphysicsCompMesh* cmeshHdiv, TPZCompMesh* cmeshH1, ProblemData* SimData, TPZVec<REAL>& elementErrors, int nthreads = 0);

private:
  /**
   * @brief Placeholder returning zero for every element and point.
   */
  static REAL ForcingFunctionWellbore(TPZCompEl* celMixed, const TPZManVector<REAL,3>& pt);
  /**
   * @brief Return the maximum nodal distance using only the first gel->Dimension()
   * coordinate components; not a general embedded-geometry diameter.
   */
  static REAL ElementDiameter(TPZGeoEl* gel);
  /**
   * @brief Integrate legacy 3D primal/dual residual and boundary contributions.
   */
  static REAL GoalContribution3D(TPZCompEl* celMixed, TPZCompEl* celDual, ProblemData* SimData);
  /**
   * @brief Reserved declaration with no current implementation.
   */
  static REAL GoalContributionInterface(TPZCompEl* celMixed, TPZCompEl* celDual, ProblemData* SimData);
  /**
   * @brief Integrate legacy well primal/dual residual terms; not called by GoalOriented().
   */
  static REAL GoalContribution1D(TPZCompEl* celMixed, TPZCompEl* celDual, ProblemData* SimData);
  /**
   * @brief Dispatch the selected marking scheme; unknown values fall back to simple marking.
   */
  static void MarkElementsForRefinement(const TPZVec<REAL>& elementErrors, TPZVec<int>& refinementIndicator, REAL theta, MarkingStrategy strategy = MarkingStrategy::EBulk);
  /**
   * @brief Mark indicators >= theta times the maximum indicator (including zero ties).
   */
  static void SimpleMarking(const TPZVec<REAL>& elementErrors, TPZVec<int>& refinementIndicator, REAL theta);
  /**
   * @brief Mark descending indicators until their sum reaches theta times the total.
   * Uses indicator values directly, not their squares.
   */
  static void BulkMarking(const TPZVec<REAL>& elementErrors, TPZVec<int>& refinementIndicator, REAL theta);
  /**
   * @brief Extend refinement marks according to geometric neighbor/refinement patterns.
   */
  static void MeshSmoothing(TPZGeoMesh* gmesh, TPZVec<int>& refinementIndicator);
  /**
   * @brief Extend refinement marks to maintain legacy well/surface/reservoir compatibility.
   */
  static void MeshWellCompatibility(TPZGeoMesh* gmesh, TPZVec<int>& refinementIndicator, ProblemData* SimData);
};