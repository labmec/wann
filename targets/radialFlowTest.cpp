#ifdef HAVE_CONFIG_H
#include <pz_config.h>
#endif

// NeoPZ includes
#include <iostream>
#include <TPZLinearAnalysis.h>
#include <TPZSSpStructMatrix.h>
#include <pzskylstrmatrix.h>
#include <pzstepsolver.h>
#include <TPZAnalyticSolution.h>
#include <TPZH1ApproxCreator.h>
#include <TPZHDivApproxCreator.h>

// WannLib includes
#include "TPZWannAnalysis.h"
#include "TPZWannGeometryTools.h"
#include "TPZWannApproxTools.h"
#include "TPZWannPostProcTools.h"

// ====================
// Functions prototypes
// ====================

TPZMultiphysicsCompMesh *MixedDarcyCompMesh(TPZGeoMesh *gmesh, ProblemData *SimData, TLaplaceExample1 *exact, bool isCondensed = false);
TPZCompMesh *H1DarcyCompMesh(TPZGeoMesh *gmesh, ProblemData *SimData, TLaplaceExample1 *exact, bool isCondensed = false);
TPZGeoMesh *CreateRadialMesh(ProblemData *SimData);

// =============
// Main function
// =============

int main(int argc, char *argv[]) {
  
  TLaplaceExample1 exact; // Global variable to be used in the material objects
  exact.fDimension = 3;
  exact.fExact = TLaplaceExample1::ENone;

  std::string jsonfile = "radialTest.json";

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
  TPZGeoMesh* gmesh = CreateRadialMesh(&SimData);

  // plot gmesh
  {
    std::ofstream plotfile("RadialMesh.vtk");
    TPZVTKGeoMesh::PrintGMeshVTK(gmesh, plotfile);
  }

  // Create computational meshes
  TPZMultiphysicsCompMesh* cmeshMixed = MixedDarcyCompMesh(gmesh, &SimData, &exact, false);
  TPZCompMesh* cmeshH1 = H1DarcyCompMesh(gmesh, &SimData, &exact, false);

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

  std::cout << "\n--------- Simulation finished ---------" << std::endl;
  std::cout << "\n--------- Starting post-processing ---------" << std::endl;

  // Plot mesh
  {
    std::ofstream plotfile("RadialMesh.vtk");
    TPZVTKGeoMesh::PrintGMeshVTK(gmesh, plotfile);
  }

  // Plotting solution
  TPZWannPostProcTools::WriteReservoirVTK(cmeshMixed, &SimData);
  TPZWannPostProcTools::WriteReservoirVTK(cmeshH1, &SimData);

  // Expected flow entering the wellbore
  // For this test, the reservoir height represent the external radius
  REAL pr = SimData.m_Reservoir.BCs.at("surface_farfield").value;
  REAL pw = SimData.m_Reservoir.BCs.at("surface_wellbore_cylinder").value;
  REAL Q = 2.0*M_PI*SimData.m_Reservoir.perm*SimData.m_Reservoir.length / (SimData.m_Fluid.viscosity * log(SimData.m_Reservoir.height/SimData.m_Wellbore.radius));
  Q = Q * (pr - pw);
  
  // Computed flow entering the wellbore
  TPZVec<REAL> segmentPoints = {0.0, SimData.m_Wellbore.length};
  TPZVec<REAL> QHdiv = TPZWannPostProcTools::ComputeWellFluxes(cmeshMixed, &SimData, segmentPoints);
  TPZVec<REAL> QH1 = TPZWannPostProcTools::ComputeWellFluxes(cmeshH1, &SimData, segmentPoints);
  
  std::cout << "\nExpected flow entering the wellbore: " << Q << std::endl;
  std::cout << "Computed flow entering the wellbore (H(div) mesh): " << QHdiv
            << " (relative error: " << std::abs(QHdiv[0] - Q)/std::abs(Q) * 100 << " %)" << std::endl;
  std::cout << "Computed flow entering the wellbore (H1 mesh): " << QH1
            << " (relative error: " << std::abs(QH1[0] - Q)/std::abs(Q) * 100 << " %)" << std::endl;

   std::cout << "\n--------- Post-processing finished ---------" << std::endl;

  delete cmeshMixed;
  delete cmeshH1;
  delete gmesh;
}

// ========================
// Functions implementation
// ========================

// Basically a copy of TPZWannGeometryTools::CreateGeoMesh, but implementing
// the cylindrical mapping on the far field as well.
TPZGeoMesh* CreateRadialMesh(ProblemData* SimData) {

  TPZGeoMesh* gmesh = TPZWannGeometryTools::ReadMeshFromGmsh(SimData);

  if (1) {
    std::ofstream out("gmeshorig.vtk");
    TPZVTKGeoMesh::PrintGMeshVTK(gmesh, out);
  }

  if (SimData->m_Mesh.ToCylindrical) {
    // Internal and external cylindrical surfaces
    TPZWannGeometryTools::ModifyGeometricMeshToCylWell(gmesh, SimData->ESurfWellCyl, SimData->m_Wellbore.radius);
    TPZWannGeometryTools::ModifyGeometricMeshToCylWell(gmesh, SimData->EFarField, SimData->m_Reservoir.height);
  }

  if (SimData->m_Mesh.customRefinement != 0) {
    std::string file = SimData->m_Mesh.file;
    std::string baseName = file.substr(0, file.find_last_of('.'));
    std::string refProcessFile = baseName + "_refProcess.txt";
    std::string path(std::string(INPUTDIR) + "/" + refProcessFile);
    TPZWannGeometryTools::RefineFromFile(gmesh, path);
    if (SimData->m_PostProc.verbosityLevel) {
      std::ofstream out("gmesh_customref.vtk");
      TPZVTKGeoMesh::PrintGMeshVTK(gmesh, out);
    }
  
  // Only perform uniform and directional refinements if no custom refinement is specified
  } else {
    if (SimData->m_Mesh.NumUniformRef) {
      TPZCheckGeom checkgeom(gmesh);
      checkgeom.UniformRefine(SimData->m_Mesh.NumUniformRef);

      if (SimData->m_PostProc.verbosityLevel) {
        std::ofstream out("gmeshnonlin.vtk");
        TPZVTKGeoMesh::PrintGMeshVTK(gmesh, out);
      }
    }

    if (SimData->m_Mesh.NumDirRef) {
      gRefDBase.InitializeRefPatterns(gmesh->Dimension());
      for (int i = 0; i < SimData->m_Mesh.NumDirRef; i++) {
        TPZRefPatternTools::RefineDirectional(gmesh, {SimData->ECurveHeel, SimData->ECurveToe});
      }
      if (SimData->m_PostProc.verbosityLevel) {
        std::ofstream out("gmeshnonlin_ref.vtk");
        TPZVTKGeoMesh::PrintGMeshVTK(gmesh, out);
      }
    }
  }

  return gmesh;
}

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
  TPZNullMaterial<STATE> *mat = new TPZNullMaterial<>(SimData->EHDivBoundInterface, 2);

  // Boundary conditions --- 

  for (auto &bcpair : ReservoirData.BCs) {
    auto &bc = bcpair.second;
    TPZFMatrix<STATE> val1(1, 1, 0.);
    TPZManVector<STATE> val2(1, 0);
    val2[0] = bc.value;
    TPZBndCondT<STATE> *BCond = reservoirMat->CreateBC(reservoirMat, bc.matid, bc.type, val1, val2);
    if (hasAnalyticSol) BCond->SetForcingFunctionBC(exact->ExactSolution(), 3);
    hdivCreator.InsertMaterialObject(BCond);
  }

  TPZMultiphysicsCompMesh *cmesh = hdivCreator.CreateApproximationSpace();
  return cmesh;
}

TPZCompMesh *H1DarcyCompMesh(TPZGeoMesh *gmesh, ProblemData *SimData, TLaplaceExample1 *exact, bool isCondensed) {
  TPZH1ApproxCreator h1Creator(gmesh);
  h1Creator.SetDefaultOrder(SimData->m_Reservoir.pOrder);
  h1Creator.ProbType() = ProblemType::EDarcy;
  h1Creator.SetShouldCondense(isCondensed);

  TLaplaceExample1* exactsol = dynamic_cast<TLaplaceExample1 *>(exact);
  bool hasAnalyticSol = (exactsol != nullptr && exactsol->fExact != TLaplaceExample1::ENone);

  // Reservoir data
  auto &ReservoirData = SimData->m_Reservoir;

  // Insert material
  TPZDarcyFlow *reservoirMat = new TPZDarcyFlow(SimData->EDomain, gmesh->Dimension());

  if (hasAnalyticSol) {
    reservoirMat->SetExactSol(exact->ExactSolution(), 3);
    reservoirMat->SetForcingFunction(exact->ForceFunc(), 3);
    reservoirMat->SetConstantPermeability(1.0);
  } else {
    reservoirMat->SetConstantPermeability(ReservoirData.perm/SimData->m_Fluid.viscosity);
  }
  h1Creator.InsertMaterialObject(reservoirMat);

  // Boundary conditions ---

  for (auto &bcpair : ReservoirData.BCs) {
    auto &bc = bcpair.second;
    TPZFMatrix<STATE> val1(1, 1, 0.);
    TPZManVector<STATE> val2(1, 0);
    val2[0] = bc.value;
    TPZBndCondT<STATE> *BCond = reservoirMat->CreateBC(reservoirMat, bc.matid, bc.type, val1, val2);
    if (hasAnalyticSol) BCond->SetForcingFunctionBC(exact->ExactSolution(), 3);
    h1Creator.InsertMaterialObject(BCond);
  }

  // Create the H1 computational mesh
  TPZCompMesh *cmesh = h1Creator.CreateClassicH1ApproximationSpace();
  return cmesh;
}
