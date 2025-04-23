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
#include <TPZMultiphysicsCompMesh.h>
#include <TPZHDivApproxCreator.h>
#include <Material/DarcyFlow/TPZMixedDarcyFlow.h>
#include <Material/TPZNullMaterial.h>
#include <Material/TPZNullMaterialCS.h>
#include <Material/TPZLagrangeMultiplierCS.h>
#include <Material/TPZLagrangeMultiplier.h>
#include <Mesh/pzintel.h>

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

class TPZWannApproxTools {

public:
  static TPZMultiphysicsCompMesh* CreateMultiphysicsCompMesh(TPZGeoMesh* gmesh, ProblemData* SimData);
  
private:
  static void AddPressureSkinElements(TPZCompMesh* cmesh, ProblemData* SimData, const int laglevel);
  static void AddWellboreElements(TPZVec<TPZCompMesh*>& meshvec, ProblemData* SimData, const int laglevel);
  static void EqualizePressureConnects(TPZCompMesh* cmesh, ProblemData* SimData);
  static void AddHDivBoundInterfaceElements(TPZCompMesh* cmesh, ProblemData* SimData);
  static void AddInterfaceElements(TPZMultiphysicsCompMesh* cmesh, ProblemData* SimData, const int laglevel);
};