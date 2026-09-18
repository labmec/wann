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
#include "TPZMixedDarcyAnisotropic.h"
#include "TPZDarcyAnisotropic.h"
#include "ProblemData.h"

// ================
// Global variables
// ================

const bool shouldPlot = true;

// Refinement parameters
const int maxIterations = 3;
const REAL errorTolerance = 1e-3;
const REAL relativeRefTol = 0.5; 

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

  std::string jsonfile = "ozkan1999adaptive.json";

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

  // Convert to cylindrical coordinates if set in json
  if (SimData.m_Mesh.ToCylindrical) {
    TPZManVector<REAL,3> cylcenter = {0.,0.,0.};
    REAL hr = SimData.m_Reservoir.height;
    REAL lr = SimData.m_Wellbore.height;
    cylcenter[2] = lr - hr/2.;
    TPZWannGeometryTools::ModifyGeometricMeshToCylWell(gmesh, SimData.ESurfWellCyl, SimData.m_Wellbore.radius, cylcenter);
  }

  // --- Adaptive refinement loop ---

  int refIt = 0;
  REAL estimatedError = errorTolerance + 1.0;
  TPZVec estimatedErrorVec(maxIterations+1, 0.0);

  // Open file to store refinement process
  std::string file = SimData.m_Mesh.file;
  std::string baseName = file.substr(0, file.find_last_of('.'));
  std::string path = std::string(INPUTDIR) + "/" + baseName + "_refProcess.txt";
  std::ofstream refinementLog(path);
  if (!refinementLog) {
    std::cerr << "Error: Could not open refinement_log.txt for writing." << std::endl;
    DebugStop();
  }

  std::cout << "\n=== Starting adaptive refinement loop ===" << std::endl;

  while (refIt <= maxIterations && estimatedError > errorTolerance) {
    std::cout << "\n--- Computing Hdiv and H1 approximations at iteration " << refIt << " ---" << std::endl;

    // Create computational meshes
    TPZMultiphysicsCompMesh* cmeshMixed = MixedDarcyCompMesh(gmesh, &SimData, &exact, false);
    TPZCompMesh* cmeshH1 = H1DarcyCompMesh(gmesh, &SimData, &exact);

    if (cmeshMixed->NEquations() >= 2000000) {
      std::cerr << "Error: Number of equations exceeds 2 million. Stopping refinement." << std::endl;
      break;
    }

    // H(div) analysis
    TPZLinearAnalysis anMixed(cmeshMixed);
#ifdef PZ_USING_MKL
    TPZSSpStructMatrix<STATE> matMixed(cmeshMixed);
#else
    TPZSkylStrMatrix matMixed(cmeshMixed);
#endif
    matMixed.SetNumThreads(SimData.m_Numerics.nthreads);
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
    matH1.SetNumThreads(SimData.m_Numerics.nthreads);
    anH1.SetStructuralMatrix(matH1);

    TPZStepSolver<STATE> stepH1;
    stepH1.SetDirect(ECholesky);
    anH1.SetSolver(stepH1);
    anH1.Run();

    // Plot solutions (for internal control only)
    if (shouldPlot) {
      // GeoMesh
      {
        std::ofstream plotfile("geomesh_" + std::to_string(refIt) + ".vtk");
        TPZVTKGeoMesh::PrintGMeshVTK(gmesh, plotfile);
      }

      // H1 Darcy
      {
        const std::string plotfile = "h1Darcy_plot_" + std::to_string(refIt);
        constexpr int vtkRes{0};
        TPZManVector<std::string, 2> fields = {"Flux", "Pressure"};
        bool safe_check = SimData.m_Mesh.ToCylindrical? true : false;
        auto vtk = TPZVTKGenerator(cmeshH1, fields, plotfile, vtkRes, 3, safe_check);
        vtk.Do();
      }

      // Mixed Darcy
      {
        const std::string plotfile = "mixedDarcy_plot_" + std::to_string(refIt);
        constexpr int vtkRes{0};
        TPZManVector<std::string, 2> fields = {"Flux", "Pressure"};
        bool safe_check = SimData.m_Mesh.ToCylindrical? true : false;
        auto vtk = TPZVTKGenerator(cmeshMixed, fields, plotfile, vtkRes, 3, safe_check);
        vtk.Do();
      }
    }

    std::cout << "\n--- Error estimation and adaptive refinement at iteration " << refIt << " ---" << std::endl;

    int64_t ngel = cmeshMixed->Reference()->NElements();
    TPZVec<REAL> elementErrors(ngel, 0.0);
    TPZVec<int> refinementIndicator(ngel, 0);

    // Prager-Synge error estimation
    estimatedError = TPZWannAdaptivityTools::PragerSynge(cmeshMixed, cmeshH1, &SimData, elementErrors, SimData.m_Numerics.nthreads);
    estimatedErrorVec[refIt] = estimatedError;
    std::cout << "Estimated error: " << estimatedErrorVec[refIt] << std::endl;

    // Mesh adaptive refinement (in the last iteration we only compute the error, no refinement)
    if (refIt < maxIterations) {
      TPZVec<int64_t> refinedElements;
      TPZWannAdaptivityTools::AdaptivityProcess(gmesh, &SimData, elementErrors, refinedElements, relativeRefTol, TPZWannAdaptivityTools::MarkingStrategy::ESimple);

      // Export refinement information to file
      refinementLog << refinedElements.size() << " ";
      for (size_t i = 0; i < refinedElements.size(); ++i) {
        refinementLog << refinedElements[i];
        if (i != refinedElements.size() - 1) refinementLog << " ";
      }
      refinementLog << "\n";
    }

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

  std::cout << "\n=== Adaptive refinement loop completed ===" << std::endl;
  std::cout << "Error history: ";
  for (int i = 0; i < refIt; ++i) {
    std::cout << estimatedErrorVec[i];
    if (i != refIt - 1) std::cout << ", ";
  }
  std::cout << std::endl;

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

  TPZMixedDarcyAnisotropic *reservoirMat = new TPZMixedDarcyAnisotropic(SimData->EDomain, gmesh->Dimension());
  TPZFNMatrix<9, STATE> perm(3, 3, 0.);
  for (int i = 0; i < 3; i++) {
    perm(i, i) = ReservoirData.perm[i];
  }
  reservoirMat->SetConstantPermeability(perm);
  hdivCreator.InsertMaterialObject(reservoirMat);

  // Boundary conditions --- 

  TPZFMatrix<STATE> val1(1, 1, 0.);
  TPZManVector<STATE> val2(1, 0);

  // No flux on cylinder bases and caprock
  val2[0] = 0.;
  TPZBndCondT<STATE> *BCond = reservoirMat->CreateBC(reservoirMat, SimData->ESurfHeel, 1, val1, val2);
  hdivCreator.InsertMaterialObject(BCond);
  BCond = reservoirMat->CreateBC(reservoirMat, SimData->ESurfToe, 1, val1, val2);
  hdivCreator.InsertMaterialObject(BCond);
  BCond = reservoirMat->CreateBC(reservoirMat, SimData->ECapRock, 1, val1, val2);
  hdivCreator.InsertMaterialObject(BCond);

  // Zero pressure on farfield
  val2[0] = 0.;
  BCond = reservoirMat->CreateBC(reservoirMat, SimData->EFarField, 0, val1, val2);
  hdivCreator.InsertMaterialObject(BCond);

  // Prescribed flux on well surface
  REAL Q = SimData->m_Wellbore.BCs["point_heel"].value;
  REAL Area = 2.0 * M_PI * SimData->m_Wellbore.radius * SimData->m_Wellbore.height;
  val2[0] = -Q / Area;
  BCond = reservoirMat->CreateBC(reservoirMat, SimData->ESurfWellCyl, 1, val1, val2);
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
  TPZDarcyAnisotropic *reservoirMat = new TPZDarcyAnisotropic(SimData->EDomain, gmesh->Dimension());
  TPZFNMatrix<9, STATE> perm(3, 3, 0.);
  for (int i = 0; i < 3; i++) {
    perm(i, i) = SimData->m_Reservoir.perm[i];
  }
  reservoirMat->SetConstantPermeability(perm);
  h1Creator.InsertMaterialObject(reservoirMat);

  // Bondary conditions ---

  TPZFMatrix<STATE> val1(1, 1, 0.);
  TPZManVector<STATE> val2(1, 0);

  // No flux on cylinder bases and caprock
  val2[0] = 0.;
  TPZBndCondT<STATE> *BCond = reservoirMat->CreateBC(reservoirMat, SimData->ESurfHeel, 1, val1, val2);
  h1Creator.InsertMaterialObject(BCond);
  BCond = reservoirMat->CreateBC(reservoirMat, SimData->ESurfToe, 1, val1, val2);
  h1Creator.InsertMaterialObject(BCond);
  BCond = reservoirMat->CreateBC(reservoirMat, SimData->ECapRock, 1, val1, val2);
  h1Creator.InsertMaterialObject(BCond);

  // Zero pressure on farfield
  val2[0] = 0.;
  BCond = reservoirMat->CreateBC(reservoirMat, SimData->EFarField, 0, val1, val2);
  h1Creator.InsertMaterialObject(BCond);

  // Prescribed flux on well surface
  REAL Q = SimData->m_Wellbore.BCs["point_heel"].value;
  REAL Area = 2.0 * M_PI * SimData->m_Wellbore.radius * SimData->m_Wellbore.height;
  val2[0] = Q / Area;
  BCond = reservoirMat->CreateBC(reservoirMat, SimData->ESurfWellCyl, 1, val1, val2);
  h1Creator.InsertMaterialObject(BCond); 

  // Create the H1 computational mesh
  TPZCompMesh *cmesh = h1Creator.CreateClassicH1ApproximationSpace();
  return cmesh;
}