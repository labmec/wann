#ifdef HAVE_CONFIG_H
#include <pz_config.h>
#endif

// NeoPZ includes
#include <charconv>
#include <cmath>
#include <iostream>
#include <sstream>
#include <string_view>
#include <Material/DarcyFlow/TPZDarcyFlow.h>
#include <TPZLinearAnalysis.h>
#include "TPZWannStructMatrix.h"
#include <pzstepsolver.h>
#include <TPZH1ApproxCreator.h>
#include <TPZHDivApproxCreator.h>

// WannLib includes
#include "TPZWannGeometryTools.h"
#include "TPZWannPostProcTools.h"

// =====================
// Global exact solution
// =====================

ProblemData SimData;

REAL RadialMobility(const ProblemData& data) {
  // ReadJson has scaled viscosity so this mobility acts on scaled pressure.
  return data.m_Reservoir.perm[0] / data.m_Fluid[0].viscosity;
}

void ExactSolution(const TPZVec<REAL> &loc, TPZVec<STATE> &result, TPZFMatrix<STATE> &deriv) {
  result.Resize(1);
  deriv.Redim(3, 1);

  // loc[0], loc[1], loc[2] are x, y, z coordinates
  // Return pressure and its gradient in the same scaled units as the solve.

  REAL y = loc[1];
  REAL z = loc[2];

  // Example: radial solution for cylindrical coordinates
  REAL r = std::sqrt(y * y + z * z);
  REAL r0 = SimData.m_Wellbore[0].radius;  // wellbore radius
  REAL r_ext = SimData.m_Reservoir.height;  // external radius (reservoir height)
  REAL p_ext = SimData.m_Reservoir.BCs.at("surface_farfield").value;  // pressure at the far field
  REAL p0 = SimData.m_Reservoir.BCs.at("surface_wellbore_cylinder").value;  // pressure at inner boundary (wellbore)

  // Pressure: logarithmic profile p = p0 + (p_ext - p0) * ln(r/r0) / ln(r_ext/r0)
  REAL ln_ratio = std::log(r / r0) / std::log(r_ext / r0);
  result[0] = p0 + (p_ext - p0) * ln_ratio;

  // Gradient
  REAL dPdr = (p_ext - p0) / (r * std::log(r_ext / r0));

  // Chain rule: dp/dx = dP/dr * dr/dx = dP/dr * x/r
  deriv(0, 0) = 0.0;  // dP/dx
  deriv(1, 0) = dPdr * y / r;  // dP/dy
  deriv(2, 0) = dPdr * z / r;  // dP/dz
}

void ZeroSource(const TPZVec<REAL> &, TPZVec<STATE> &result) {
  result.Resize(1);
  result[0] = 0.0;  // No source term for the radial flow problem
}

// ====================
// Functions prototypes
// ====================

TPZMultiphysicsCompMesh *MixedDarcyCompMesh(TPZGeoMesh *gmesh, ProblemData *SimData, bool isCondensed = false);
TPZCompMesh *H1DarcyCompMesh(TPZGeoMesh *gmesh, ProblemData *SimData, bool isCondensed = false);
TPZGeoMesh *CreateRadialMesh(ProblemData *SimData);
void ModifyToCylMesh(TPZGeoMesh *gmesh, int matid, REAL radius, const TPZManVector<REAL, 3> &center);

// =============
// Main function
// =============

int main(int argc, char *argv[]) {
  
  int maxRef = 4;
  std::string jsonfile = argc >= 2 ? argv[1] : "radialTest.json";
  if (argc > 3) {
    std::cerr << "Usage: " << argv[0] << " [input.json [max-refinements]]\n";
    return 1;
  }
  if (argc == 3) {
    const std::string_view argument(argv[2]);
    const auto [end, error] = std::from_chars(argument.data(), argument.data() + argument.size(), maxRef);
    if (error != std::errc{} || end != argument.data() + argument.size() || maxRef < 0) {
      std::cerr << "The maximum number of refinements must be a nonnegative integer.\n";
      return 1;
    }
  }

  if (jsonfile.find(".json") == std::string::npos) {
    jsonfile += ".json";
  }

  std::cout << "Using json file: " << jsonfile << std::endl;
  std::cout << "\n======= Starting simulation =======\n" << std::endl;

#ifdef PZ_LOG
  TPZLogger::InitializePZLOG();
#endif

  // Problem data
  SimData.ReadJson(jsonfile);

  if (SimData.m_Wellbore[0].radius <= 0.0 || SimData.m_Reservoir.height <= SimData.m_Wellbore[0].radius ||
      SimData.m_Reservoir.BCs.count("surface_farfield") == 0 ||
      SimData.m_Reservoir.BCs.count("surface_wellbore_cylinder") == 0) {
    std::cerr << "Invalid simulation data for exact solution." << std::endl;
    DebugStop();
  }

  // The radial input uses reservoir height to store the outer radius.
  const REAL outerRadius = SimData.m_Reservoir.height;
  const REAL wellRadius = SimData.m_Wellbore[0].radius;
  const REAL logRadiusRatio = std::log(outerRadius / wellRadius);
  const REAL pressureDrop = SimData.m_Reservoir.BCs.at("surface_farfield").value -
                            SimData.m_Reservoir.BCs.at("surface_wellbore_cylinder").value;
  const REAL mobility = RadialMobility(SimData);
  const REAL length = SimData.m_Reservoir.length;
  const REAL Q = 2.0 * M_PI * mobility * length * pressureDrop / logRadiusRatio;

  // Integrate the exact reference fields over the annulus with volume measure 2*pi*L*r dr.
  // Pressure is relative to the inner boundary. NeoPZ weights the flux norm by 1/mobility.
  const REAL pressureIntegral = outerRadius * outerRadius *
      (logRadiusRatio * (logRadiusRatio - 1.) + 0.5) - 0.5 * wellRadius * wellRadius;
  const REAL exactPressureNorm = std::abs(pressureDrop) / logRadiusRatio *
      std::sqrt(M_PI * length * pressureIntegral);
  const REAL exactWeightedFluxNorm = std::abs(pressureDrop) *
      std::sqrt(2.0 * M_PI * length * mobility / logRadiusRatio);

  struct Formulation {
    const char* name;
    bool mixed;
    int order;
  };
  const int order = SimData.m_Numerics.reservoirPorder;
  const Formulation formulations[] = {
    {"H(div) mixed", true, order},
    {"H1", false, order},
    {"H1", false, order + 2}
  };

  std::ofstream results("radialTestResults.txt");
  if (!results) {
    std::cerr << "Could not open radialTestResults.txt for writing.\n";
    return 1;
  }
  TPZGeoMesh* gmesh = CreateRadialMesh(&SimData);
  TPZCheckGeom checkgeom(gmesh);

  for (int ref = 0; ref <= maxRef; ++ref) {
    for (const auto& formulation : formulations) {
      ProblemData caseData = SimData;
      caseData.m_Numerics.reservoirPorder = formulation.order;
      gmesh->ResetReference();

      // Build the approximation for this formulation.
      TPZCompMesh* cmesh = nullptr;
      TPZVec<TPZCompMesh*> atomicMeshes;
      if (formulation.mixed) {
        TPZMultiphysicsCompMesh* mixed = MixedDarcyCompMesh(gmesh, &caseData, true);
        cmesh = mixed;
        atomicMeshes = mixed->MeshVector();
      } else {
        cmesh = H1DarcyCompMesh(gmesh, &caseData, true);
      }

      std::ostringstream heading;
      heading << "\n===== " << formulation.name << " degree " << formulation.order
              << ", refinement level " << ref << " =====\n"
              << "Number of equations: " << cmesh->NEquations() << '\n';
      std::cout << heading.str() << std::flush;
      results << heading.str() << std::flush;

      // Solve and compute the errors against the exact radial solution.
      REAL computedFlow;
      TPZVec<REAL> errors(formulation.mixed ? 5 : 3, 0.);
      {
        TPZLinearAnalysis an(cmesh);
        TPZWannStructMatrix mat(cmesh);
        mat.SetNumThreads(caseData.m_Numerics.nthreads);
        an.SetStructuralMatrix(mat);
        TPZStepSolver<STATE> step;
        step.SetDirect(formulation.mixed ? ELDLt : ECholesky);
        an.SetSolver(step);
        an.Run();

        computedFlow = TPZWannPostProcTools::ComputeWellFluxes(cmesh, &caseData, 0, 1)[0];
        std::fstream errorfile("errors.txt", std::ios::app);
        an.SetThreadsForError(caseData.m_Numerics.nthreads);
        an.PostProcessError(errors, false, errorfile);
      }

      const REAL pressureError = errors[formulation.mixed ? 0 : 1];
      const REAL pressureErrorPa = pressureError / caseData.m_Numerics.pressureScale;
      const REAL weightedFluxError = errors[formulation.mixed ? 1 : 2];
      std::ostringstream report;
      report << "Expected flow entering the wellbore: " << Q << '\n'
             << "Computed flow entering the wellbore: " << computedFlow << '\n'
             << "L2 norm of pressure error (using pressure in Pa): " << pressureErrorPa << '\n'
             << "Mobility-weighted L2 norm of flux error: " << weightedFluxError << '\n';
      if (formulation.mixed) {
        report << "H(div) norm of flux error: " << errors[4] << '\n'
               << "L2 norm of divergence error: " << errors[2] << '\n';
      }
      if (exactPressureNorm > 0. && exactWeightedFluxNorm > 0.) {
        report << "Relative pressure error (E_p): " << pressureError / exactPressureNorm << '\n'
               << "Relative Darcy flux error (E_u): " << weightedFluxError / exactWeightedFluxNorm << '\n';
      } else {
        report << "Relative pressure and Darcy flux errors are undefined for zero reference norms.\n";
      }
      if (Q != 0.) {
        report << "Relative total flow error (E_Q): " << std::abs((computedFlow - Q) / Q) << '\n';
      } else {
        report << "Relative total flow error is undefined for zero exact flow.\n";
      }
      std::cout << report.str();
      results << report.str();

      if (ref == maxRef) {
        caseData.m_PostProc.reservoir_vtk += "_p" + std::to_string(formulation.order);
        TPZWannPostProcTools::WriteReservoirVTK(cmesh, &caseData);
      }

      delete cmesh;
      for (auto* mesh : atomicMeshes) delete mesh;
    }

    if (ref < maxRef) checkgeom.UniformRefine(1);
  }
  delete gmesh;
  return 0;
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
    TPZManVector<REAL,3> cylcenter = {0., 0., 0.};
    int matid_ext = SimData->m_Reservoir.BCs.at("surface_farfield").matid;
    int matid_int = SimData->m_Reservoir.BCs.at("surface_wellbore_cylinder").matid;
    ModifyToCylMesh(gmesh, matid_int, SimData->m_Wellbore[0].radius, cylcenter);
    ModifyToCylMesh(gmesh, matid_ext, SimData->m_Reservoir.height, cylcenter);
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
  }

  return gmesh;
}

TPZMultiphysicsCompMesh *MixedDarcyCompMesh(TPZGeoMesh *gmesh, ProblemData *SimData, bool isCondensed) {
  TPZHDivApproxCreator hdivCreator(gmesh);
  hdivCreator.SetProbType(ProblemType::EDarcy);
  hdivCreator.SetHdivFamily(HDivFamily::EHDivStandard);
  hdivCreator.SetDefaultOrder(SimData->m_Numerics.reservoirPorder);
  hdivCreator.SetShouldCondense(isCondensed);

  // Reservoir data
  auto &ReservoirData = SimData->m_Reservoir;

  TPZMixedDarcyFlow *reservoirMat = new TPZMixedDarcyFlow(ReservoirData.matid, gmesh->Dimension());
  reservoirMat->SetExactSol(ExactSolution, 20);
  reservoirMat->SetConstantPermeability(RadialMobility(*SimData));
  reservoirMat->SetForcingFunction(ZeroSource, 20);
  hdivCreator.InsertMaterialObject(reservoirMat);

  // Boundary conditions --- 

  for (auto &bcpair : ReservoirData.BCs) {
    auto &bc = bcpair.second;
    TPZFMatrix<STATE> val1(1, 1, 0.);
    TPZManVector<STATE> val2(1, 0);
    val2[0] = bc.value;
    TPZBndCondT<STATE> *BCond = reservoirMat->CreateBC(reservoirMat, bc.matid, bc.type, val1, val2);
    hdivCreator.InsertMaterialObject(BCond);
  }

  TPZMultiphysicsCompMesh *cmesh = hdivCreator.CreateApproximationSpace();
  return cmesh;
}

TPZCompMesh *H1DarcyCompMesh(TPZGeoMesh *gmesh, ProblemData *SimData, bool isCondensed) {
  TPZH1ApproxCreator h1Creator(gmesh);
  h1Creator.SetDefaultOrder(SimData->m_Numerics.reservoirPorder);
  h1Creator.SetProbType(ProblemType::EDarcy);
  h1Creator.SetShouldCondense(isCondensed);

  // Reservoir data
  auto &ReservoirData = SimData->m_Reservoir;

  // Insert material
  TPZDarcyFlow *reservoirMat = new TPZDarcyFlow(ReservoirData.matid, gmesh->Dimension());

  reservoirMat->SetExactSol(ExactSolution, 20);
  reservoirMat->SetConstantPermeability(RadialMobility(*SimData));
  reservoirMat->SetForcingFunction(ZeroSource, 20);
  h1Creator.InsertMaterialObject(reservoirMat);

  // Boundary conditions ---

  for (auto &bcpair : ReservoirData.BCs) {
    auto &bc = bcpair.second;
    TPZFMatrix<STATE> val1(1, 1, 0.);
    TPZManVector<STATE> val2(1, 0);
    val2[0] = bc.value;
    TPZBndCondT<STATE> *BCond = reservoirMat->CreateBC(reservoirMat, bc.matid, bc.type, val1, val2);
    h1Creator.InsertMaterialObject(BCond);
  }

  // Create the H1 computational mesh
  TPZCompMesh *cmesh = h1Creator.CreateClassicH1ApproximationSpace();
  return cmesh;
}

void ModifyToCylMesh(TPZGeoMesh *gmesh, int matid, REAL radius, const TPZManVector<REAL, 3> &center) {
  const TPZManVector<REAL,3> cylcenter = center;
  const TPZManVector<REAL,3> cylaxis = {1.,0.,0.};
  int64_t nel = gmesh->NElements();
  for(int64_t iel = 0; iel < nel ; iel++) {
    TPZGeoEl* gel = gmesh->Element(iel);
    if(!gel) continue;
    if(gel->HasSubElement()) DebugStop();
    if(gel->MaterialId() != matid) continue;
    
    TPZManVector<int64_t, 4> nodeindices;
    gel->GetNodeIndices(nodeindices);
    const int nnodes = gel->NCornerNodes();
    //Moving the nodes to the cylinder surface
    TPZManVector<REAL,3> xnode(3,0);
    for (int in = 0; in < nnodes; in++)
    {
      gel->NodePtr(in)->GetCoordinates(xnode);
      // component of xnode that is orthogonal to cyl axis
      TPZManVector<REAL, 3> x_orth = xnode - cylcenter;
      const REAL dax = Dot(x_orth, cylaxis);
      for (int ix = 0; ix < 3; ix++)
      {
        x_orth[ix] -= dax * cylaxis[ix];
      }
      const auto normdiff = fabs(Norm(x_orth) - radius);
      if (normdiff > 1e-10)
      {
        // Moving the node to the cylinder shell
        REAL computed_radius = Norm(x_orth);
        for (int ix = 0; ix < 3; ix++)
        {
          xnode[ix] = cylcenter[ix] + (x_orth[ix] / computed_radius) * radius + dax * cylaxis[ix];
        }
        x_orth = xnode - cylcenter;
        for (int ix = 0; ix < 3; ix++)
        {
          x_orth[ix] -= dax * cylaxis[ix];
        }
        REAL new_radius = Norm(x_orth);
        gel->NodePtr(in)->SetCoord(xnode);

        PZError << __PRETTY_FUNCTION__
                << "\nNode not on cylinder shell: " << gel->NodePtr(in)->Id()
                << "\nComputed radius: " << computed_radius
                << "\nGiven radius: " << radius
                << "\nElement index: " << iel
                << "\nNew coordinates: " << xnode
                << "\nNew radius: " << new_radius
                << std::endl;
      }
    }

    TPZChangeEl::ChangeToCylinder(gmesh, iel, cylcenter, cylaxis, radius);
  }

  gmesh->BuildConnectivity();

  nel = gmesh->NElements();
  for(int64_t iel = 0; iel < nel ; iel++) {
    TPZGeoEl* gel = gmesh->Element(iel);
    if(!gel) continue;
    if(gel->HasSubElement()) DebugStop();
    // if(gel->MaterialId() != ESurfWellCyl && gel->MaterialId() != ESurfHeel && gel->MaterialId() != ESurfToe) continue;
    if(gel->MaterialId() != matid) continue;
    int nsides = gel->NSides();
    // for (int iside = gel->NCornerNodes(); iside < nsides; iside++) {
    for (int iside = gel->FirstSide(1); iside < nsides; iside++) {
      TPZGeoElSide gelSide(gel,iside);
      TPZStack<TPZGeoElSide> allneigh;
      for(auto neigh = gelSide.Neighbour(); neigh != gelSide ; neigh++) {      
        if(neigh.Element()->IsGeoBlendEl()) {
          continue;
        }
        if (neigh.Element()->MaterialId() == matid) {
          continue;
        }
        // if(neigh.Element()->Dimension() < 2) continue;
        // if(neigh.Element()->MaterialId() != EDomain) DebugStop();
        allneigh.Push(neigh);
      }
    //   std::cout << "Element " << iel << " side " << iside << " has " << allneigh.size() << " neighbours" << std::endl;
      for(auto it : allneigh){
        TPZChangeEl::ChangeToGeoBlend(gmesh, it.Element()->Index());
      }
    }
  }
}
