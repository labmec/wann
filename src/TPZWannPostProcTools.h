#pragma once

#include "ProblemData.h"
#include <TPZGeoMeshTools.h>
#include <TPZVTKGenerator.h>
#include <pzvec_extras.h>
#include <pzcompel.h>
#include <pzcmesh.h>
#include <TPZMultiphysicsCompMesh.h>
#include <Material/DarcyFlow/TPZMixedDarcyFlow.h>

class TPZWannPostProcTools {

public:
  static void GenerateTrainingData(TPZGeoMesh* gmesh, ProblemData* SimData);
  static void WriteWellboreVTK(TPZCompMesh* cmesh, ProblemData* SimData);
  static void WriteReservoirVTK(TPZCompMesh* cmesh, ProblemData* SimData);
  static void WriteVTKs(TPZCompMesh* cmesh, ProblemData* SimData);
  static void PostProcessAllData(TPZCompMesh* cmesh, TPZGeoMesh* gmesh, ProblemData* SimData);
  static TPZVec<REAL> ComputeWellFluxes(TPZCompMesh* cmesh, ProblemData* SimData, TPZVec<REAL> segmentPoints);
};