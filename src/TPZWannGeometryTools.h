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

class TPZWannGeometryTools {

public:
  static TPZGeoMesh* CreateGeoMesh(ProblemData* simData);
  static TPZGeoMesh* ReadMeshFromGmsh(ProblemData* simData);
  static void ModifyGeometricMeshToCylWell(TPZGeoMesh *gmesh, int matid, REAL cylradius, TPZManVector<REAL,3> &cylcenter);
  static void hRefinement(TPZGeoMesh* gmesh, TPZVec<int64_t>& toRefine);
  static void RefineFromFile(TPZGeoMesh* og_gmesh, const std::string& filename);
  static void InsertXCoorInSet(const REAL x, std::set<REAL>& nodeCoordsX, const REAL tol);
  static REAL FindClosestX(const REAL x, const std::set<REAL>& nodeCoordsX, const REAL tol);
  static bool CheckXInSet(const REAL x, const std::set<REAL>& nodeCoordsX, const REAL tol);
  static void OrderIds(TPZGeoMesh *gmesh, ProblemData *SimData);
  static void OrderIdsSingleWell(TPZGeoMesh *gmesh, int wellid, const TPZManVector<REAL,3>& axisPoint, const TPZManVector<REAL,3>& axis);
  static REAL ComputeAxialCoordinate(const TPZManVector<REAL,3>& point, const TPZManVector<REAL,3>& axisPoint, const TPZManVector<REAL,3>& axis);

private:
  static void CreateCouplingEls(TPZGeoMesh *gmesh, ProblemData *SimData);
  static bool VerifyMesh(TPZGeoMesh *gmesh, ProblemData *SimData);
};
