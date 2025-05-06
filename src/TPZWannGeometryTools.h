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
  static void InsertXCoorInSet(const REAL x, std::set<REAL>& nodeCoordsX, const REAL tol);
  static REAL FindClosestX(const REAL x, const std::set<REAL>& nodeCoordsX, const REAL tol);

private:
  static TPZGeoMesh* ReadMeshFromGmsh(ProblemData* simData);
  static void ModifyGeometricMeshToCylWell(TPZGeoMesh *gmesh, ProblemData *SimData);
  static void CreatePressure2DElsAndOrderIds(TPZGeoMesh *gmesh, ProblemData *SimData);
};