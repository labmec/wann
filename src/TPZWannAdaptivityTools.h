/**
 * @file
 * @brief Class and method contracts for TPZWannAdaptivityTools.
 */

#pragma once

#include <TPZMultiphysicsCompMesh.h>
#include <string>
#include "ProblemData.h"

/**
 * @brief Reservoir error estimation and refinement compatible with well coupling.
 */
class TPZWannAdaptivityTools {

public:
  /**
   * @brief Reservoir flux-discrepancy/balance estimator comparing mixed and H1
   * solutions on shared geometry. Fills geometrically indexed indicators and returns
   * their Euclidean norm; uses nthreads+1 workers and zero source.
   * Writes the error field to estimatorFile when post-processing verbosity is enabled.
   */
  static REAL PragerSynge(TPZMultiphysicsCompMesh* cmeshHdiv, TPZCompMesh* cmeshH1,
                        ProblemData* SimData, TPZVec<REAL>& elementErrors, int nthreads = 0,
                        const std::string& estimatorFile = "PragerSyngeEst.vtk");
  /**
   * @brief Mark geometric indices, enforce well compatibility, smooth, recheck
   * compatibility, populate refinedElements, and divide the marked elements.
   * theta is the marking threshold relative to the largest indicator.
   */
  static void AdaptivityProcess(TPZGeoMesh* gmesh, ProblemData* SimData, TPZVec<REAL>& elementErrors, TPZVec<int64_t>& refinedElements, REAL theta);

private:
  /**
   * @brief Return the maximum nodal distance using only the first gel->Dimension()
   * coordinate components; not a general embedded-geometry diameter.
   */
  static REAL ElementDiameter(TPZGeoEl* gel);
  /**
   * @brief Mark indicators >= theta times the maximum indicator for positive indicators.
   */
  static void SimpleMarking(const TPZVec<REAL>& elementErrors, TPZVec<int>& refinementIndicator, REAL theta);
  /**
   * @brief Extend refinement marks according to geometric neighbor/refinement patterns.
   */
  static void MeshSmoothing(TPZGeoMesh* gmesh, TPZVec<int>& refinementIndicator);
  /**
   * @brief Extend refinement marks to maintain well/surface/reservoir compatibility.
   */
  static void MeshWellCompatibility(TPZGeoMesh* gmesh, TPZVec<int>& refinementIndicator, ProblemData* SimData);
};
