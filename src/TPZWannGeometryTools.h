/**
 * @file
 * @brief Class and method contracts for TPZWannGeometryTools.
 */

#pragma once

#include "ProblemData.h"
#include <TPZGmshReader.h>
#include <TPZCylinderMap.h>
#include <tpzgeoelrefpattern.h>
#include <TPZVTKGeoMesh.h>
#include <tpzchangeel.h>
#include <pzvec_extras.h>
#include <TPZGeoMeshTools.h>
#include <TPZRefPatternDataBase.h>
#include <TPZRefPatternTools.h>
#include <pzcheckgeom.h>
#include <pzgeoel.h>
#include <pzgeoelbc.h>

/**
 * @brief Static utilities for importing and preparing coupled reservoir/well geometry.
 * Creates auxiliary surface elements and orders node IDs along each well axis.
 * The coupled mesh verifier assumes a 3D reservoir with 2D boundary faces.
 * Mesh factories return newly allocated objects for caller-managed lifetime.
 */
class TPZWannGeometryTools {

public:
  /**
   * @brief Import, verify, add coupling surfaces, and order well node IDs.
   * @param[in] simData Populated simulation configuration.
   * @return Newly allocated geometry mesh; caller manages its lifetime.
   * Cylindrical conversion and refinement are not invoked by this factory.
   */
  static TPZGeoMesh* CreateGeoMesh(ProblemData* simData);
  /**
   * @brief Read m_Mesh.file relative to INPUTDIR; rebuild codimension-one boundary
   * elements with TPZGeoElBC for consistent orientation. Returns a newly allocated mesh.
   */
  static TPZGeoMesh* ReadMeshFromGmsh(ProblemData* simData);
  /**
   * @brief Project selected surface nodes to an x-aligned cylinder and change their
   * elements to cylindrical maps; convert adjoining elements to GeoBlend.
   * Requires unrefined selected elements; modifies geometry/connectivity in place.
   */
  static void ModifyGeometricMeshToCylWell(TPZGeoMesh *gmesh, int matid, REAL cylradius, TPZManVector<REAL,3> &cylcenter);
  /**
   * @brief Divide leaf elements listed by geometric index; skip already divided
   * elements and trigger DebugStop() for a missing element.
   */
  static void hRefinement(TPZGeoMesh* gmesh, TPZVec<int64_t>& toRefine);
  /**
   * @brief Replay refinement levels from filename. Each line contains a count
   * followed by that many geometric element indices; indices are used after each level.
   */
  static void RefineFromFile(TPZGeoMesh* og_gmesh, const std::string& filename);
  /**
   * @brief Insert a scalar coordinate only when no neighboring set entry is within tol.
   * Despite the name, also used for projected axial coordinates.
   */
  static void InsertXCoorInSet(const REAL x, std::set<REAL>& nodeCoordsX, const REAL tol);
  /**
   * @brief Return a neighboring set representative within tol; an empty set or
   * no acceptable representative triggers DebugStop().
   */
  static REAL FindClosestX(const REAL x, const std::set<REAL>& nodeCoordsX, const REAL tol);
  /**
   * @brief Return whether a neighboring set entry is within tol of x.
   */
  static bool CheckXInSet(const REAL x, const std::set<REAL>& nodeCoordsX, const REAL tol);
  /**
   * @brief Locate each well heel/toe point using named BC IDs and order its surface
   * node IDs along the heel-to-toe axis.
   */
  static void OrderWellIds(TPZGeoMesh *gmesh, ProblemData *SimData);
  /**
   * @brief Assign increasing unique node IDs to axial groups of leaf surface elements
   * with material wellid. Coordinates and node-vector indices are unchanged.
   */
  static void OrderIdsSingleWell(TPZGeoMesh *gmesh, int wellid, const TPZManVector<REAL,3>& axisPoint, const TPZManVector<REAL,3>& axis);
  /**
   * @brief Return dot(point-axisPoint, axis/Norm(axis)).
   * Requires three-component vectors and a finite nonzero axis.
   */
  static REAL ComputeAxialCoordinate(const TPZManVector<REAL,3>& point, const TPZManVector<REAL,3>& axisPoint, const TPZManVector<REAL,3>& axis);

private:
  /**
   * @brief Add pressure-skin, pressure-interface, and H(div)-boundary geometry
   * on every configured well cylindrical surface.
   */
  static void CreateCouplingEls(TPZGeoMesh *gmesh, ProblemData *SimData);
  /**
   * @brief Check 3D reservoir IDs/face neighbors and 2D boundary neighbors.
   * @return True when errors were found (not a success flag).
   * Marks faulty elements with diagnostic material IDs 1001 through 1005.
   */
  static bool VerifyMesh(TPZGeoMesh *gmesh, ProblemData *SimData);
};
