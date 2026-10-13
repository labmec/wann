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

class TPZWannApproxTools {

public:
  static TPZMultiphysicsCompMesh* CreateMixedCompMesh(TPZGeoMesh* gmesh, ProblemData* SimData);
  static TPZCompMesh* CreateH1CompMesh(TPZGeoMesh *gmesh, ProblemData *SimData, TPZAnalyticSolution *exact);

private:
  static void ComputeTotalMobility(TPZCompMesh* cmesh, ProblemData* SimData);
  static void AddPressureSkinElements(TPZCompMesh* cmesh, ProblemData* SimData, const int laglevel);
  static void EqualizeH1PressureConnects(TPZCompMesh *cmesh, ProblemData *SimData);
  static void AddWellboreElements(TPZVec<TPZCompMesh*>& meshvec, ProblemData* SimData, const int laglevel);
  static void EqualizeL2PressureConnects(TPZCompMesh* cmesh, ProblemData* SimData);
  static void AddHDivBoundInterfaceElements(TPZCompMesh* cmesh, ProblemData* SimData);
  static void AddInterfaceElements(TPZMultiphysicsCompMesh* cmesh, ProblemData* SimData, const int laglevel);
};