/**
 * @file
 * @brief Class and method contracts for TPZWannApproxTools.
 */

#pragma once

#include "TPZWannGeometryTools.h"
#include <TPZMultiphysicsCompMesh.h>
#include <TPZHDivApproxCreator.h>
#include <Material/DarcyFlow/TPZMixedDarcyFlow.h>
#include <Material/DarcyFlow/TPZDarcyFlow.h>
#include <Material/TPZNullMaterial.h>
#include <Material/TPZNullMaterialCS.h>
#include <Material/TPZLagrangeMultiplierCS.h>
#include <Material/TPZLagrangeMultiplier.h>
#include <Mesh/pzintel.h>
#include <TPZAnalyticSolution.h>

/**
 * @brief Builds computational spaces and reservoir/well coupling connects.
 * The mixed builder supports multiple wells; the H1 builder currently uses
 * only the first well and x-coordinate-based pressure constraints. Both builders
 * use the first fluid and diagonal reservoir permeability divided by viscosity.
 * Mesh factories return newly allocated objects for caller-managed lifetime.
 */
class TPZWannApproxTools {

public:
  /**
   * @brief Build reservoir flux/pressure atomic meshes, pressure skins, 1D wells,
   * shared pressure connects, and multiphysics Lagrange interfaces.
   * Returns a newly allocated mesh; condensation is disabled in this builder.
   */
  static TPZMultiphysicsCompMesh* CreateMixedCompMesh(TPZGeoMesh* gmesh, ProblemData* SimData);
  /**
   * @brief Build the pressure-only coupled mesh for the first well with shared
   * x-coordinate pressure dependencies. An optional TLaplaceExample1 exact
   * solution/forcing is attached to the well material. Returns a new mesh.
   */
  static TPZCompMesh* CreateH1CompMesh(TPZGeoMesh *gmesh, ProblemData *SimData, TPZAnalyticSolution *exact);

private:
  /**
   * @brief Reserved private helper declaration; no implementation or calls currently exist.
   */
  static void ComputeTotalMobility(TPZCompMesh* cmesh, ProblemData* SimData);
  /**
   * @brief Build disconnected quadrilateral skins, retain axial pressure modes,
   * and share connects circumferentially within each segment of each well.
   * Uses endpoint-projected coordinates, requires wellPorder >= 1 and unique ownership.
   */
  static void AddPressureSkinElements(TPZCompMesh* cmesh, ProblemData* SimData, const int laglevel);
  /**
   * @brief Collapse circumferential skin pressure modes through identity connect
   * dependencies grouped by global x; remove circumferential/face interior modes
   * and rebuild blocks. Current implementation assumes a single x-aligned well.
   */
  static void EqualizeH1PressureConnects(TPZCompMesh *cmesh, ProblemData *SimData);
  /**
   * @brief Add disconnected 1D pressure elements to meshvec[1] and H(div) flow
   * elements plus heel/toe endpoints to meshvec[0], for every configured well.
   */
  static void AddWellboreElements(TPZVec<TPZCompMesh*>& meshvec, ProblemData* SimData, const int laglevel);
  /**
   * @brief Make each 1D well pressure element reuse its adjacent pressure-skin
   * edge and endpoint connects, accounting for reversed node order.
   */
  static void EqualizeL2PressureConnects(TPZCompMesh* cmesh, ProblemData* SimData);
  /**
   * @brief Build reservoir normal-flux trace elements on EHDivBoundInterface geometry.
   */
  static void AddHDivBoundInterfaceElements(TPZCompMesh* cmesh, ProblemData* SimData);
  /**
   * @brief Create multiphysics interface elements between pressure-skin and
   * H(div)-boundary traces. The laglevel parameter is currently unused.
   */
  static void AddInterfaceElements(TPZMultiphysicsCompMesh* cmesh, ProblemData* SimData, const int laglevel);
};