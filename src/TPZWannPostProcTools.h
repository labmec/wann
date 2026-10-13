/**
 * @file
 * @brief Class and method contracts for TPZWannPostProcTools.
 */

#pragma once

#include "ProblemData.h"
#include <TPZGeoMeshTools.h>
#include <TPZVTKGenerator.h>
#include <pzvec_extras.h>
#include <pzcompel.h>
#include <pzcmesh.h>
#include <TPZMultiphysicsCompMesh.h>
#include <Material/DarcyFlow/TPZMixedDarcyFlow.h>

/**
 * @brief VTK export, well inflow integration, and productivity calculations.
 * Uses the current computational solution and reloads geometric references where
 * needed. GenerateTrainingData() is a legacy first-well, x-aligned mixed-only
 * exporter; WriteVTKs() avoids that legacy path.
 */
class TPZWannPostProcTools {

public:
  /**
   * @brief Append x, inferred div(Q)/(farfield-pressure), and radius to TRAININGDIR.
   * @note Legacy mixed-only first-well exporter; assumes x alignment, a specific
   * circumferential sampling line, and a well BC named surface_farfield.
   */
  static void GenerateTrainingData(TPZGeoMesh* gmesh, ProblemData* SimData);
  /**
   * @brief Export pressure, axial flow, divergence, and pressure gradient for each well.
   * Names include Hdiv/H1 and a zero-based well index.
   */
  static void WriteWellboreVTK(TPZCompMesh* cmesh, ProblemData* SimData);
  /**
   * @brief Export reservoir Pressure/Flux for dimension 3 using configured output
   * name, resolution, threads, and ToCylindrical safety-check flag.
   */
  static void WriteReservoirVTK(TPZCompMesh* cmesh, ProblemData* SimData);
  /**
   * @brief Write reservoir and per-well VTK files without invoking legacy training export.
   */
  static void WriteVTKs(TPZCompMesh* cmesh, ProblemData* SimData);
  /**
   * @brief Reload computational references, write VTKs, and invoke GenerateTrainingData().
   * Inherits the legacy exporter restrictions; use WriteVTKs() for VTK-only output.
   */
  static void PostProcessAllData(TPZCompMesh* cmesh, TPZGeoMesh* gmesh, ProblemData* SimData);
  /**
   * @brief Integrate signed reservoir normal flux into equal-length heel-to-toe bins.
   * @param[in] cmesh Solved mixed or H1 computational mesh used to evaluate reservoir flux.
   * @param[in] SimData Simulation configuration defining the well, reservoir IDs, and integration order.
   * @param[in] wellIndex Zero-based index of the selected well.
   * @param[in] nSegments Positive number of bins.
   * @return One integrated flux per bin, oriented by the geometric face normal.
   * Entire faces belong to the bin containing their axial centroid; faces are not
   * split at bin boundaries. Uses meshvec[0] for mixed flux or H1 variable 7.
   */
  static TPZVec<REAL> ComputeWellFluxes(TPZCompMesh* cmesh, ProblemData* SimData, int wellIndex, int nSegments);
  /**
   * @brief Compute endpoint flow divided by signed endpoint-minus-farfield pressure
   * difference, unscale pressure, and convert to STB/(day psi). Requires named
   * point_heel/point_toe well BCs and reservoir surface_farfield BC.
   * @note The current loop accepts heel or toe endpoints and overwrites values;
   * the result can depend on element order. No zero-pressure-difference check exists.
   */
  static REAL ProductivityIndex(TPZCompMesh* cmesh, ProblemData* SimData, int wellIndex);
};