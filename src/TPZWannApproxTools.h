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

class TPZWannApproxTools {

public:
  static TPZMultiphysicsCompMesh* CreateMultiphysicsCompMesh(TPZGeoMesh* gmesh, ProblemData* SimData);
  static TPZCompMesh* CreateH1CompMesh(TPZGeoMesh *gmesh, ProblemData *SimData) ;

private:
  static void AddPressureSkinElements(TPZCompMesh* cmesh, ProblemData* SimData, const int laglevel);
  static void EqualizePressureSkinConnects(TPZCompMesh* cmesh, ProblemData* SimData);
  static void EqualizeH1Connects(TPZCompMesh *cmesh, ProblemData *SimData);
  static void AddWellboreElements(TPZVec<TPZCompMesh*>& meshvec, ProblemData* SimData, const int laglevel);
  static void EqualizePressureConnects(TPZCompMesh* cmesh, ProblemData* SimData);
  static void AddHDivBoundInterfaceElements(TPZCompMesh* cmesh, ProblemData* SimData);
  static void AddInterfaceElements(TPZMultiphysicsCompMesh* cmesh, ProblemData* SimData, const int laglevel);
};