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
  static void hRefinement(TPZGeoMesh* gmesh, TPZVec<int>& RefinementIndicator);
  static void RefineFromFile(TPZGeoMesh* og_gmesh, const std::string& filename);
  static void InsertXCoorInSet(const REAL x, std::set<REAL>& nodeCoordsX, const REAL tol);
  static REAL FindClosestX(const REAL x, const std::set<REAL>& nodeCoordsX, const REAL tol);
  static bool CheckXInSet(const REAL x, const std::set<REAL>& nodeCoordsX, const REAL tol);
  static void OrderIds(TPZGeoMesh *gmesh, ProblemData *SimData);
  static void DividePyramids(TPZGeoMesh *gmesh);

private:
  static bool SetBC(ProblemData* simData, const std::string& bcName, int matid);
  static void CreatePressure2DEls(TPZGeoMesh *gmesh, ProblemData *SimData);
  static bool VerifyMesh(TPZGeoMesh *gmesh, ProblemData *SimData);
};