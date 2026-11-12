/**
 * @file
 * @brief Class and method contracts for TPZWannPostProcTools.
 */

#pragma once

#include "ProblemData.h"
#include <TPZVTKGenerator.h>
#include <pzvec_extras.h>
#include <pzcompel.h>
#include <pzcmesh.h>
#include <TPZMultiphysicsCompMesh.h>
#include <Material/DarcyFlow/TPZMixedDarcyFlow.h>

/**
 * @brief VTK export, well inflow integration, and productivity calculations.
 * Uses the current computational solution and reloads geometric references where
 * needed.
 */
class TPZWannPostProcTools {

public:
  /**
   * @brief Export pressure and axial flow for each well, plus divergence for mixed wells.
   * Names include Hdiv/H1 and a zero-based well index.
   */
  static void WriteWellboreVTK(TPZCompMesh* cmesh, ProblemData* SimData);
  /**
   * @brief Write paper-format x,Curve1 CSVs from the solved mixed well elements.
   * Exports flow magnitude in STB/day, inflow normalized by the heel rate, and
   * far-field pressure drop in psi against normalized heel-to-toe distance.
   * Interior endpoints use the toe-side trace. A zero heel rate leaves the
   * dimensionless inflow column empty.
   */
  static bool WriteWellProfilesCSV(TPZMultiphysicsCompMesh* cmesh, const ProblemData& data,
                                   int subdivisionsPerElement = 4);
  /**
   * @brief Export reservoir Pressure/Flux for dimension 3 using configured output
   * name, resolution, threads, and ToCylindrical safety-check flag.
   */
  static void WriteReservoirVTK(TPZCompMesh* cmesh, ProblemData* SimData);
  /**
   * @brief Write reservoir and per-well VTK files for the solved mesh.
   */
  static void WriteVTKs(TPZCompMesh* cmesh, ProblemData* SimData);
  /**
   * @brief Integrate signed reservoir normal flux into equal-length heel-to-toe bins.
   * @param[in] cmesh Solved mixed or H1 computational mesh used to evaluate reservoir flux.
   * @param[in] SimData Simulation configuration defining the well, reservoir IDs, and integration order.
   * @param[in] wellIndex Zero-based index of the selected well.
   * @param[in] nSegments Positive number of bins.
   * @return One integrated flux per bin, positive from the reservoir into the well.
   * Entire faces belong to the bin containing their axial centroid; faces are not
   * split at bin boundaries.
   */
  static TPZVec<REAL> ComputeWellFluxes(TPZCompMesh* cmesh, ProblemData* SimData, int wellIndex, int nSegments);
  /**
   * @brief Compute endpoint flow divided by signed endpoint-minus-farfield pressure
   * difference, unscale pressure, and convert to STB/(day psi). Requires named
   * point_heel/point_toe well BCs and reservoir surface_farfield BC.
   * @note The well flow and pressure are evaluated at the heel.
   */
  static REAL ProductivityIndex(TPZCompMesh* cmesh, ProblemData* SimData, int wellIndex);
};
