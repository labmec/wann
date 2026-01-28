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
#include "TPZWannAdaptivityTools.h"
#include "ProblemData.h"

// ================
// Global variables
// ================

const int globalNthreads = 8;
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

// Computational mesh for H1 Darcy formulation
TPZCompMesh* H1DarcyCompMesh(TPZGeoMesh* gmesh, ProblemData* SimData, TLaplaceExample1* exact, bool isCondensed = false);

// ============
// Main program
// ============

int main(int argc, char *argv[]) {
  
  TLaplaceExample1 exact; // Exact solution (if known)
  exact.fDimension = 3;
  exact.fExact = TLaplaceExample1::ENone;

  std::string jsonfile = "wann3d_test.json";

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

#ifdef PZ_LOG
  TPZLogger::InitializePZLOG();
#endif

  // Problem data
  ProblemData SimData;  
  SimData.ReadJson(jsonfile);

  // Initial geometric mesh
  TPZGeoMesh* gmesh = TPZWannGeometryTools::ReadMeshFromGmsh(&SimData);

  // Convert to cylindrical coordinates if needed
  if (SimData.m_Mesh.ToCylindrical) {
    TPZWannGeometryTools::ModifyGeometricMeshToCylWell(gmesh, &SimData);
  }

  // --- Adaptive refinement loop ---

  int refIt = 0;
  REAL estimatedError = errorTolerance + 1.0;

  // Open file to store refinement process
  std::string file = SimData.m_Mesh.file;
  std::string baseName = file.substr(0, file.find_last_of('.'));
  std::string path = std::string(INPUTDIR) + "/" + baseName + "_refProcess.txt";
  std::ofstream refinementLog(path);
  if (!refinementLog) {
    std::cerr << "Error: Could not open refinement_log.txt for writing." << std::endl;
    return 1;
  }

  std::cout << "\n=== Starting adaptive refinement loop ===" << std::endl;

  while (refIt < maxIterations && estimatedError > errorTolerance) {
    // --- Computing Hdiv and H1 approximations ---
    std::cout << "\n--- Computing Hdiv and H1 approximations at iteration " << refIt << " ---" << std::endl;

    // Create computational meshes
    TPZMultiphysicsCompMesh* cmeshMixed = MixedDarcyCompMesh(gmesh, &SimData, &exact, false);
    TPZCompMesh* cmeshH1 = H1DarcyCompMesh(gmesh, &SimData, &exact);

    // H(div) analysis
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

    // H1 analysis
    TPZLinearAnalysis anH1(cmeshH1);
#ifdef PZ_USING_MKL
    TPZSSpStructMatrix<STATE> matH1(cmeshH1);
#else
    TPZSkylStrMatrix matH1(cmeshH1);
#endif
    matH1.SetNumThreads(globalNthreads);
    anH1.SetStructuralMatrix(matH1);

    TPZStepSolver<STATE> stepH1;
    stepH1.SetDirect(ECholesky);
    anH1.SetSolver(stepH1);
    anH1.Run();

    // Plot solutions (for internal checking only)
    if (shouldPlot) {
      // GeoMesh
      {
        std::ofstream plotfile("geomesh_" + std::to_string(refIt) + ".vtk");
        TPZVTKGeoMesh::PrintGMeshVTK(gmesh, plotfile);
      }

      // Mixed Darcy
      {
        const std::string plotfile = "mixedDarcy_plot_" + std::to_string(refIt);
        constexpr int vtkRes{0};
        TPZManVector<std::string, 2> fields = {"Flux", "Pressure"};
        auto vtk = TPZVTKGenerator(cmeshMixed, fields, plotfile, vtkRes);
        vtk.Do();
      }

      // H1 Darcy
      {
        const std::string plotfile = "h1Darcy_plot_" + std::to_string(refIt);
        constexpr int vtkRes{0};
        TPZManVector<std::string, 2> fields = {"Flux", "Pressure"};
        auto vtk = TPZVTKGenerator(cmeshH1, fields, plotfile, vtkRes);
        vtk.Do();
      }
    }

    // --- Error estimation and adaptive refinement ---
    std::cout << "\n--- Error estimation and adaptive refinement at iteration " << refIt << " ---" << std::endl;

    int64_t ngel = cmeshMixed->Reference()->NElements();
    TPZVec<REAL> elementErrors(ngel, 0.0);
    TPZVec<int> refinementIndicator(ngel, 0);

    REAL estimatedError = TPZWannAdaptivityTools::PragerSynge(cmeshMixed, cmeshH1, &SimData, elementErrors, globalNthreads);
    std::cout << "Estimated error: " << estimatedError << std::endl;

    // Mark elements for refinement
    TPZWannAdaptivityTools::MarkElementsForRefinement(elementErrors, refinementIndicator, relativeRefTol);
    REAL sumRef = std::accumulate(refinementIndicator.begin(), refinementIndicator.end(), 0);
    std::cout << "Initial number of to-refine elements: " << sumRef << std::endl;

    // Smooth the refinement map
    TPZWannAdaptivityTools::MeshSmoothing(cmeshMixed->Reference(), refinementIndicator);
    sumRef = std::accumulate(refinementIndicator.begin(), refinementIndicator.end(), 0);
    std::cout << "Number of to-refine elements after mesh smoothing: " << sumRef << std::endl;

    // Compatibilize refinement to well geometry
    TPZWannAdaptivityTools::MeshWellCompatibility(cmeshMixed->Reference(), refinementIndicator, &SimData);
    sumRef = std::accumulate(refinementIndicator.begin(), refinementIndicator.end(), 0);
    std::cout << "Final number of to-refine elements: " << sumRef << std::endl;

    // Finally perform the refinement
    TPZWannGeometryTools::hRefinement(cmeshMixed->Reference(), refinementIndicator);

    // Export refinementIndicators vector to file "inputs/mesh_name_refProcess.txt"
    // Each line corresponds to one refinement iteration
    // On each line: <number of elements> <refinementIndicator[0]> ... <refinementIndicator[n-1]>
    refinementLog << refinementIndicator.size() << " ";
    for (size_t i = 0; i < refinementIndicator.size(); ++i) {
      refinementLog << refinementIndicator[i];
      if (i != refinementIndicator.size() - 1) refinementLog << " ";
    }
    refinementLog << "\n";

    // --- Clean up ---

    // Remove dependencies before deleting H1 mesh
    int ncon = cmeshH1->NConnects();
    for (int i = 0; i < ncon; ++i) {
      cmeshH1->ConnectVec()[i].RemoveDepend();
    }

    delete cmeshMixed;
    delete cmeshH1;
    refIt++;
  }

  refinementLog.close();
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
    reservoirMat->SetConstantPermeability(ReservoirData.perm);
  }
  hdivCreator.InsertMaterialObject(reservoirMat);

  // Boundary conditions --- 
  // TODO: add version with analytic solution

  TPZFMatrix<STATE> val1(1, 1, 0.);
  TPZManVector<STATE> val2(1, 0);

  // No flux on cylinder bases
  val2[0] = 0.;
  TPZBndCondT<STATE> *BCond = reservoirMat->CreateBC(reservoirMat, SimData->ESurfHeel, 1, val1, val2);
  hdivCreator.InsertMaterialObject(BCond);

  BCond = reservoirMat->CreateBC(reservoirMat, SimData->ESurfToe, 1, val1, val2);
  hdivCreator.InsertMaterialObject(BCond);

  // Zero pressure on farfield
  val2[0] = 0.;
  BCond = reservoirMat->CreateBC(reservoirMat, SimData->EFarField, 0, val1, val2);
  hdivCreator.InsertMaterialObject(BCond);

  // Prescribed pressure on well surface
  val2[0] = 1.;
  BCond = reservoirMat->CreateBC(reservoirMat, SimData->ESurfWellCyl, 0, val1, val2);
  hdivCreator.InsertMaterialObject(BCond); 

  TPZMultiphysicsCompMesh *cmesh = hdivCreator.CreateApproximationSpace();
  return cmesh;
}

TPZCompMesh *H1DarcyCompMesh(TPZGeoMesh *gmesh, ProblemData *SimData, TLaplaceExample1 *exact, bool isCondensed) {
  TPZH1ApproxCreator h1Creator(gmesh);
  h1Creator.SetDefaultOrder(SimData->m_Reservoir.pOrder+2);
  h1Creator.ProbType() = ProblemType::EDarcy;
  h1Creator.SetShouldCondense(isCondensed);

  TLaplaceExample1* exactsol = dynamic_cast<TLaplaceExample1 *>(exact);
  bool hasAnalyticSol = (exactsol != nullptr && exactsol->fExact != TLaplaceExample1::ENone);

  // Insert material
  TPZDarcyFlow *reservoirMat = new TPZDarcyFlow(SimData->EDomain, gmesh->Dimension());

  if (hasAnalyticSol) {
    reservoirMat->SetExactSol(exact->ExactSolution(), 3);
    reservoirMat->SetForcingFunction(exact->ForceFunc(), 3);
    reservoirMat->SetConstantPermeability(1.0);
  } else {
    reservoirMat->SetConstantPermeability(SimData->m_Reservoir.perm);
  }
  h1Creator.InsertMaterialObject(reservoirMat);

  // Boundary conditions --- 
  // TODO: add version with analytic solution

  TPZFMatrix<STATE> val1(1, 1, 0.);
  TPZManVector<STATE> val2(1, 0);

  // No flux on cylinder bases
  val2[0] = 0.;
  TPZBndCondT<STATE> *BCond = reservoirMat->CreateBC(reservoirMat, SimData->ESurfHeel, 1, val1, val2);
  h1Creator.InsertMaterialObject(BCond);

  BCond = reservoirMat->CreateBC(reservoirMat, SimData->ESurfToe, 1, val1, val2);
  h1Creator.InsertMaterialObject(BCond);

  // Zero pressure on farfield
  val2[0] = 0.;
  BCond = reservoirMat->CreateBC(reservoirMat, SimData->EFarField, 0, val1, val2);
  h1Creator.InsertMaterialObject(BCond);

  // Prescribed pressure on well surface
  val2[0] = 1.;
  BCond = reservoirMat->CreateBC(reservoirMat, SimData->ESurfWellCyl, 0, val1, val2);
  h1Creator.InsertMaterialObject(BCond); 

  // Create the H1 computational mesh
  TPZCompMesh *cmesh = h1Creator.CreateClassicH1ApproximationSpace();
  return cmesh;
}