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
REAL TotalSurfaceFlow(TPZCompMesh *cmesh, ProblemData *SimData);

// =============
// Main function
// =============

int main(int argc, char *argv[]) {
  
  TLaplaceExample1 exact; // Global variable to be used in the material objects
  exact.fDimension = 3;
  exact.fExact = TLaplaceExample1::ENone;

  std::string jsonfile = "wann3D_Cyl.json";

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

  // Aply cylindrical transformation in the exterior part of the reservoir as well!
  // if (SimData.m_Mesh.ToCylindrical) {
  //   CylindricalFarField(gmesh, &SimData);
  // }

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
  REAL QHdiv = TotalSurfaceFlow(cmeshMixed, &SimData);
  REAL QH1 = TotalSurfaceFlow(cmeshH1, &SimData);
  
  std::cout << "\nExpected flow entering the wellbore: " << Q << std::endl;
  std::cout << "Computed flow entering the wellbore (H(div) mesh): " << QHdiv
            << " (relative error: " << std::abs(QHdiv - Q)/std::abs(Q) * 100 << " %)" << std::endl;
  std::cout << "Computed flow entering the wellbore (H1 mesh): " << QH1
            << " (relative error: " << std::abs(QH1 - Q)/std::abs(Q) * 100 << " %)" << std::endl;

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

REAL TotalSurfaceFlow(TPZCompMesh *cmesh, ProblemData *SimData) {
  // Check if cmesh is Hdiv or H1
  // We are assuming that only the Hdiv mesh is multiphyiscs
  bool isHdiv;
  TPZMultiphysicsCompMesh *cmeshMult = dynamic_cast<TPZMultiphysicsCompMesh *>(cmesh);
  if (cmeshMult) {
    isHdiv = true;
  } else {
    isHdiv = false;
  }

  int matid = SimData->EDomain; 
  REAL totalFlow = 0.;

  // Ensure references point to current mesh
  cmesh->Reference()->ResetReference();
  cmesh->LoadReferences();
  TPZGeoMesh *gmesh = cmesh->Reference();

  int64_t ngel = cmesh->Reference()->NElements();

  for (auto gel : gmesh->ElementVec()) {
    if (gel->HasSubElement()) continue; // Skip non-leaf elements
    if (gel->MaterialId() != matid) continue;

    // Loop over sides to see if we are in a boundary element
    bool isBoundaryElement = false;
    int side = -1;
    TPZGeoEl *faceGel = nullptr;
    int firstFace = gel->FirstSide(gel->Dimension() - 1); 
    int lastFace = gel->FirstSide(gel->Dimension());
    for (int iside = firstFace; iside < lastFace; iside++) {
      TPZGeoElSide gelside(gel, iside);
      TPZGeoElSide neighside = gelside.HasNeighbour(SimData->ESurfWellCyl);
      if (neighside) {
        isBoundaryElement = true;
        side = iside;
        faceGel = neighside.Element();
        break;
      }
    }

    if (!isBoundaryElement) continue; // Skip elements that are not on the cylindrical surface

    // Center of the element
    TPZManVector<REAL, 3> qsi(faceGel->Dimension());
    TPZManVector<REAL, 3> xCenter(3, 0.);
    faceGel->CenterPoint(faceGel->NSides() - 1, qsi); // center of the element interior
    faceGel->X(qsi, xCenter);

    // Compute contributions
    const TPZIntPoints *intrule = nullptr;
    intrule = gel->CreateSideIntegrationRule(side, SimData->m_Reservoir.pOrder + 2);
    for (int ip = 0; ip < intrule->NPoints(); ip++) {
      int dimF = faceGel->Dimension();
      TPZManVector<REAL, 3> ptOnSide(faceGel->Dimension());
      TPZFNMatrix<9, REAL> jacobian, axes, jacinv;
      REAL weight, detjac;
      faceGel->Jacobian(ptOnSide, jacobian, axes, detjac, jacinv);

      // Compute normal
      // Must done in the integration point to account for curved sides
      TPZManVector<REAL, 3> v1(3), v2(3), normal(3);
      v1[0] = axes(0, 0);
      v1[1] = axes(0, 1);
      v1[2] = axes(0, 2);
      v2[0] = axes(1, 0);
      v2[1] = axes(1, 1);
      v2[2] = axes(1, 2);

      normal[0] = v1[1] * v2[2] - v1[2] * v2[1];
      normal[1] = v1[2] * v2[0] - v1[0] * v2[2];
      normal[2] = v1[0] * v2[1] - v1[1] * v2[0];

      REAL norm = sqrt(normal[0] * normal[0] + normal[1] * normal[1] +
                       normal[2] * normal[2]);
      if (norm > 1e-12) {
        normal[0] /= norm;
        normal[1] /= norm;
        normal[2] /= norm;
      }

      intrule->Point(ip, ptOnSide, weight);
      TPZManVector<REAL, 3> ptInElement(gel->Dimension());
      TPZTransform<> trans = gel->SideToSideTransform(side, gel->NSides() - 1);
      trans.Apply(ptOnSide, ptInElement);
      weight *= fabs(detjac);
      TPZManVector<REAL, 3> sigh(gel->Dimension(), 0.0);
      TPZCompEl *cel = gel->Reference();
      if (isHdiv) {
        auto mfcel = dynamic_cast<TPZMultiphysicsElement *>(cel);
        if (!mfcel)
          DebugStop();
        TPZCompEl *celHdiv = mfcel->Element(
            0); // We are assuming that the first mesh is the H(div) mesh
        celHdiv->Solution(ptInElement, 1, sigh);
      } else {
        cel->Solution(ptInElement, 7, sigh);
      }

      // get normal flux
      REAL normalFLux =
          sigh[0] * normal[0] + sigh[1] * normal[1] + sigh[2] * normal[2];
      totalFlow += normalFLux * weight;
    }
  }
  // std::cout << "Number of boundary elements: " << count << std::endl;
  return totalFlow;
}