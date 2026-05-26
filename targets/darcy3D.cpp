#ifdef HAVE_CONFIG_H
#include <pz_config.h>
#endif

#include <iostream>

#include <TPZLinearAnalysis.h>
#include <TPZSSpStructMatrix.h>
#include <pzskylstrmatrix.h>
#include <pzstepsolver.h>
#include <TPZAnalyticSolution.h>
#include <TPZH1ApproxCreator.h>
#include <TPZHDivApproxCreator.h>
#include <TPZVTKGenerator.h>
#include <Material/DarcyFlow/TPZMixedDarcyFlow.h>
#include <Material/DarcyFlow/TPZDarcyFlow.h>
#include <TPZFileStream.h>

#include "TPZWannGeometryTools.h"
#include "TPZWannPostProcTools.h"
#include "ProblemData.h"

// ================
// Global variables
// ================

const int globalNthreads = 0;
const bool shouldPlot = true;

// Refinement parameters
const int maxIterations = 2;
const REAL errorTolerance = 1e-2;
const REAL relativeRefTol = 0.1; 

// ===================
// Function prototypes
// ===================

// Computational mesh for mixed Darcy formulation
TPZMultiphysicsCompMesh* MixedDarcyCompMesh(TPZGeoMesh* gmesh, ProblemData* SimData, TLaplaceExample1* exact, bool isCondensed = false);

void VerifyMesh(TPZGeoMesh* gmesh, ProblemData* SimData);

// ============
// Main program
// ============

int main(int argc, char *argv[]) {
  
  TLaplaceExample1 exact; // Global variable to be used in the material objects
  exact.fDimension = 3;
  exact.fExact = TLaplaceExample1::ENone;

  std::string jsonfile = "penmatcha1999.json";

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

  // Read original geometric mesh and perform the refinement process 
  // described in refinementProcess.txt file
  TPZGeoMesh* gmesh = TPZWannGeometryTools::CreateGeoMesh(&SimData);

  // Plot geometric mesh
  {
    std::ofstream out("gmeshVerified.vtk");
    TPZVTKGeoMesh::PrintGMeshVTK(gmesh, out);
  }

  // H(div) Simulation

  // For some reason we have to change the sign of boundary condition 
  // when using cylindrical map. TODO: fix it
  if (SimData.m_Mesh.ToCylindrical) {
    SimData.m_Wellbore.BCs["point_heel"].value *= -1;
  }

  // Mixed Darcy simulation
  TPZMultiphysicsCompMesh* cmeshMixed = MixedDarcyCompMesh(gmesh, &SimData, &exact);
  TPZLinearAnalysis anMixed(cmeshMixed);
#ifdef PZ_USING_MKL
  TPZSSpStructMatrix<STATE> matMixed(cmeshMixed);
#else
  TPZSkylStrMatrix matMixed(cmeshMixed);
#endif
  matMixed.SetNumThreads(globalNthreads);
  anMixed.SetStructuralMatrix(matMixed);
  TPZStepSolver<STATE> stepMixed;
  stepMixed.SetDirect(ELDLt);
  anMixed.SetSolver(stepMixed);
  anMixed.Run();

  // Plotting
  TPZWannPostProcTools::WriteReservoirVTK(cmeshMixed, &SimData);

  // Clean up
  delete cmeshMixed;
  delete gmesh;
}

// ========================
// Functions implementation
// ========================

TPZMultiphysicsCompMesh *MixedDarcyCompMesh(TPZGeoMesh *gmesh, ProblemData *SimData, TLaplaceExample1 *exact, bool isCondensed) {
  TPZHDivApproxCreator hdivCreator(gmesh);
  hdivCreator.ProbType() = ProblemType::EDarcy;
  hdivCreator.HdivFamily() = HDivFamily::EHDivStandard;
  hdivCreator.SetDefaultOrder(SimData->m_Reservoir.pOrder);
  hdivCreator.SetShouldCondense(isCondensed);

  TLaplaceExample1* exactsol = dynamic_cast<TLaplaceExample1 *>(exact);
  bool hasAnalyticSol = (exactsol != nullptr && exactsol->fExact != TLaplaceExample1::ENone);

  // Reservoir data
  auto &ReservoirData = SimData->m_Reservoir;

  TPZMixedDarcyFlow *reservoirMat = new TPZMixedDarcyFlow(SimData->EDomain, gmesh->Dimension());
  if (hasAnalyticSol) {
    reservoirMat->SetExactSol(exact->ExactSolution(), 3);
    reservoirMat->SetForcingFunction(exact->ForceFunc(), 3);
    reservoirMat->SetConstantPermeability(1.0);
  } else {
    reservoirMat->SetConstantPermeability(ReservoirData.perm/SimData->m_Fluid.viscosity);
  }
  hdivCreator.InsertMaterialObject(reservoirMat);

  // Boundary conditions ---
  TPZFMatrix<STATE> val1(1, 1, 0.);
  TPZManVector<STATE> val2(1, 0);
  TPZBndCondT<STATE> *BCond;
  for (auto &bcpair : ReservoirData.BCs) {
    auto &bc = bcpair.second;
    val2[0] = bc.value;
    BCond = reservoirMat->CreateBC(reservoirMat, bc.matid, bc.type, val1, val2);
    if (hasAnalyticSol)
      BCond->SetForcingFunctionBC(exact->ExactSolution(), 1);
    hdivCreator.InsertMaterialObject(BCond);
  }

  // Add a boundary condition on the well surface
  // We choose to add the the same boundary condition imposed in the toe point
  int type = SimData->m_Wellbore.BCs["point_toe"].type;
  val2[0] = SimData->m_Wellbore.BCs["point_toe"].value;
  BCond = reservoirMat->CreateBC(reservoirMat, SimData->ESurfWellCyl, type, val1, val2);
  hdivCreator.InsertMaterialObject(BCond); 

  TPZMultiphysicsCompMesh *cmesh = hdivCreator.CreateApproximationSpace();
  return cmesh;
}

// void VerifyMesh(TPZGeoMesh *gmesh, ProblemData *SimData) {
//     bool hasErrors = false;

//     // Check if every 3D element has a neighbour on each face.
//     for (int iel = 0; iel < gmesh->NElements(); iel++) {
//     TPZGeoEl *gel = gmesh->ElementVec()[iel];

//     if (gel->HasSubElement()) continue;
//     if (gel->MaterialId() != SimData->EDomain) continue;
    
//     int firstFace = gel->FirstSide(gel->Dimension() - 1);
//     int lastFace = gel->FirstSide(gel->Dimension()) - 1;

//     for (int side = firstFace; side <= lastFace; side++) {
//       TPZGeoElSide gelside(gel, side);
//       TPZGeoElSide neighbour = gelside.Neighbour();
//       if (!neighbour.Exists()) {
//         std::cout << "Element " << gel->Index() << " side " << side << " has no neighbour!" << std::endl;
//         gel->SetMaterialId(2000);
//         hasErrors = true;
//       } else { 
//         if (neighbour.Element()->Dimension() == 2) {
//           int matid = neighbour.Element()->MaterialId();
//           if (matid != SimData->ESurfWellCyl && matid != SimData->ESurfHeel && matid != SimData->ESurfToe && matid != SimData->EFarField && matid != SimData->ECapRock) {
//             std::cout << "Element " << gel->Index() << " side " << side << " has a neighbour with wrong material id: " << matid << std::endl;
//             gel->SetMaterialId(3000);
//             hasErrors = true;
//           }
//         }
//       }
//     }
//   }

//   // Plot mesh if it has errors
//   if (hasErrors) {
//     std::ofstream out("gmeshWithErrors.vtk");
//     TPZVTKGeoMesh::PrintGMeshVTK(gmesh, out);
//     std::cout << "Mesh has errors! Plotted elements with issues in "
//                  "gmeshWithErrors.vtk"
//               << std::endl;
//   } else {
//     std::cout << "Mesh verification passed! No issues found." << std::endl;
//   }
// }