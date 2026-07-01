#ifdef HAVE_CONFIG_H
#include <pz_config.h>
#endif

#include <iostream>
#include "TPZWannGeometryTools.h"
#include "TPZWannApproxTools.h"
#include "TPZWannPostProcTools.h"
#include "TPZWannAdaptivityTools.h"
#include <TPZLinearAnalysis.h>
#include <TPZSSpStructMatrix.h>
#include <pzskylstrmatrix.h>
#include <pzstepsolver.h>
#include <TPZAnalyticSolution.h>
#include "TPZWannAnalysis.h"

using namespace std;

// ===================
// Auxiliary functions
// ===================

ProblemData DualizeProblemData(const ProblemData &SimData);

// ===================

int main(int argc, char *argv[]) {
  
  // === Setup ===

  TLaplaceExample1 exact; // Global variable to be used in the material objects
  exact.fDimension = 3;
  exact.fExact = TLaplaceExample1::ENone;

  REAL thetaRef = 0.3; // Relative tolerance for goal-oriented refinement
  REAL absoluteRefTol = 1e-4; // Absolute tolerance for goal-oriented refinement
  int maxRefSteps = 3; // Maximum number of refinement steps

  std::string jsonfile = "easyMesh.json";

  if (argc > 2) {
    std::cout << argv[0] << " being called with too many arguments." << std::endl;
    DebugStop();
  } else if (argc == 2) {
    jsonfile = argv[1];
  }

  if (jsonfile.find(".json") == std::string::npos) {
    jsonfile += ".json";
  }

  std::cout << "Using json file: " << jsonfile << std::endl;
  std::cout << "\n--------- Starting simulation ---------" << std::endl;

#ifdef PZ_LOG
  TPZLogger::InitializePZLOG();
#endif

  // Problem data
  ProblemData SimData;  
  SimData.ReadJson(jsonfile);

  // For some reason we have to change the sign of boundary condition 
  // when using cylindrical map. TODO: fix it
  if (SimData.m_Mesh.ToCylindrical) {
    SimData.m_Wellbore.BCs["point_heel"].value *= -1;
  }

  // Read original geometric mesh. Can perform some predefined custom refinement
  // if specified in the json file.
  // Some fileds of ProblemData are only filled in CreateGeoMesh (bc.second.matid)
  // TODO: Can this be avoided?
  TPZGeoMesh* gmesh = TPZWannGeometryTools::CreateGeoMesh(&SimData);

  // Generate SimData for the dual problem.
  // Right now this is hardcoded for the chosen functional. 
  // TODO: make it more general.
  ProblemData SimDataDual = DualizeProblemData(SimData);

  // === Refinement loop ===

  int refStep = 0;
  REAL goalError = 2*absoluteRefTol;
  TPZVec<REAL> goalErrors;

  while (refStep < maxRefSteps && std::abs(goalError) > absoluteRefTol) {
    std::cout << "\n====== Refinement step " << refStep + 1 << " ======" << std::endl;

    // H(div) Simulation ---

    TPZMultiphysicsCompMesh* cmeshMixed = TPZWannApproxTools::CreateMultiphysicsCompMesh(gmesh, &SimData, &exact);
    TPZWannAnalysis anMixed(cmeshMixed, RenumType::EMetis);
    anMixed.SetProblemData(&SimData);
    anMixed.Initialize();
    anMixed.NewtonIteration();

    // Dual problem simulation ---

    TPZMultiphysicsCompMesh* cmeshMixedDual = TPZWannApproxTools::CreateMultiphysicsCompMesh(gmesh, &SimDataDual, &exact, true);
    TPZWannAnalysis anDual(cmeshMixedDual, RenumType::EMetis);
    anDual.SetProblemData(&SimDataDual);
    anDual.Initialize();
    anDual.NewtonIteration();

    std::cout << "\n--------- Simulation finished ---------" << std::endl;
    std::cout << "\n--------- Starting post-processing ---------" << std::endl;

    TPZWannPostProcTools::WriteVTKs(cmeshMixed, &SimData);
    TPZWannPostProcTools::WriteVTKs(cmeshMixedDual, &SimDataDual);

    std::cout << "\n--------- Post-processing finished ---------" << std::endl;
    std::cout << "\n--------- Starting goal-oriented refinement ---------" << std::endl;

    // Goal oriented error estimation on the 1D well
    TPZVec<REAL> elementErrors(gmesh->NElements(), 0.0);

    REAL totalError = TPZWannAdaptivityTools::GoalOriented(cmeshMixed, cmeshMixedDual, &SimData, elementErrors, SimData.m_Numerics.nthreads);
    goalErrors.push_back(totalError);
    std::cout << "Total goal-oriented error: " << totalError << std::endl;

    // Wrapper for the refinement process
    TPZWannAdaptivityTools::AdaptivityProcess(gmesh, &SimData, elementErrors, thetaRef);

    // Print refined mesh
    if (SimData.m_PostProc.verbosityLevel > 0) {
      std::ofstream out("gmeshGoalRef_" + std::to_string(refStep+1) + ".vtk");
      TPZVTKGeoMesh::PrintGMeshVTK(gmesh, out);
    }

    std::cout << "\n--------- Goal-oriented refinement finished ---------" << std::endl;

    refStep++;
    delete cmeshMixed;
    delete cmeshMixedDual;
  }

  std::cout << "\n\nEstimated error at each refinement step: " << std::endl;
  for (size_t i = 0; i < goalErrors.size(); ++i) {
    std::cout << "  Refinement step " << i+1 << ": " << goalErrors[i] << std::endl;
  }
  delete gmesh;
}

ProblemData DualizeProblemData(const ProblemData &SimData) {
    ProblemData dual = SimData;
    dual.m_PostProc.reservoir_vtk = SimData.m_PostProc.reservoir_vtk + "_dual";
    dual.m_PostProc.wellbore_vtk = SimData.m_PostProc.wellbore_vtk + "_dual";
    dual.m_Reservoir.pOrder = SimData.m_Reservoir.pOrder+1;
    dual.m_Wellbore.pOrder = SimData.m_Wellbore.pOrder+1;

    // Essential boundary conditions need to be homogeneous.
    // Here we assume that we are using only mixed formulations,
    // so the essential boundary conditions are the Neumann ones (type 1).

    // For the chosen functional, the natural boundary conditions
    // (Dirichlet ones, type 0) are also set to zero.
    // For other functionals, different Dirichlet boundary conditions may be used.

    for (auto &bcpair : dual.m_Wellbore.BCs) {
        if (bcpair.second.type == 1) {
            // bcpair.second.type = 0;
            bcpair.second.value = 0.0;
        } else if (bcpair.second.type == 0) {
            bcpair.second.value = 0.0;
        }
    }

    for (auto &bcpair : dual.m_Reservoir.BCs) {
        if (bcpair.second.type == 1) {
            bcpair.second.value = 0.0;
        } else if (bcpair.second.type == 0) {
            bcpair.second.value = 0.0;
        }
    }
    return dual;
}