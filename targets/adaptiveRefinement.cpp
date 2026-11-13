#ifdef HAVE_CONFIG_H
#include <pz_config.h>
#endif

#include <charconv>
#include <cmath>
#include <iomanip>
#include <iostream>
#include <string_view>

#include <TPZLinearAnalysis.h>
#include "TPZWannStructMatrix.h"
#include <pzstepsolver.h>
#include <TPZH1ApproxCreator.h>
#include <TPZHDivApproxCreator.h>
#include <TPZVTKGenerator.h>
#include <Material/DarcyFlow/TPZMixedDarcyFlow.h>
#include <Material/DarcyFlow/TPZDarcyFlow.h>

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
const REAL errorTolerance = 1e-4;
const REAL relativeRefTol = 0.5; 

// ===================
// Function prototypes
// ===================

TPZMultiphysicsCompMesh* MixedDarcyCompMesh(TPZGeoMesh* gmesh, const ProblemData& data,
                                           const std::map<int, REAL>& wallFlux);
TPZCompMesh* H1DarcyCompMesh(TPZGeoMesh* gmesh, const ProblemData& data,
                            const std::map<int, REAL>& wallFlux);

void WriteDarcyVTK(TPZCompMesh* mesh, const std::string& filename, bool cylindrical);
void DeleteMixedMesh(TPZMultiphysicsCompMesh* mesh);

// ============
// Main program
// ============

int main(int argc, char *argv[]) {
  
  std::string jsonfile = argc >= 2 ? argv[1] : "ozkan1999adaptive.json";
  int maxRefinements = 3;
  if (argc > 3) {
    std::cerr << "Usage: " << argv[0] << " [input.json [max-refinements]]\n";
    return 1;
  }
  if (argc == 3) {
    const std::string_view argument(argv[2]);
    const auto [end, error] = std::from_chars(argument.data(), argument.data() + argument.size(),
                                               maxRefinements);
    if (error != std::errc{} || end != argument.data() + argument.size() || maxRefinements < 0) {
      std::cerr << "The maximum number of refinements must be a nonnegative integer.\n";
      return 1;
    }
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

  const REAL pressureScale = SimData.m_Numerics.pressureScale;
  if (!std::isfinite(pressureScale) || pressureScale <= 0.) {
    std::cerr << "pressureScale must be finite and positive.\n";
    return 1;
  }
  // The Prager–Synge energy norm scales with the square root of the pressure scale.
  const REAL errorToPhysical = 1. / std::sqrt(pressureScale);
  std::cout << "Reference pressure (Pa): " << SimData.m_Numerics.refPressure << '\n'
            << "Pressure scale: " << pressureScale << '\n';

  TPZGeoMesh* gmesh = TPZWannGeometryTools::CreateGeoMesh(&SimData);
  std::map<int, REAL> wallFlux;
  for (const auto& well : SimData.m_Wellbore) {
    const auto& heel = well.BCs.at("point_heel");
    const auto& toe = well.BCs.at("point_toe");
    if (heel.type != 1 || toe.type != 1) {
      std::cerr << "Adaptive Darcy benchmark requires prescribed rates at both endpoints of "
                << well.name << '\n';
      delete gmesh;
      return 1;
    }
    const REAL area = 2. * M_PI * well.radius * well.length;
    if (!(area > 0.)) DebugStop();
    wallFlux.emplace(well.matidSurf, (toe.value - heel.value) / area);
  }

  // --- Adaptive refinement loop ---

  int refIt = 0;
  REAL estimatedError = errorTolerance + 1.0;
  TPZVec<REAL> estimatedErrorVec;
  TPZVec<REAL> maxErrorVec;

  // Open file to store refinement process
  std::string file = SimData.m_Mesh.file;
  std::string baseName = file.substr(0, file.find_last_of('.'));
  std::string path = baseName + "_refProcess.txt";
  std::ofstream refinementLog(path);
  if (!refinementLog) {
    std::cerr << "Could not open the refinement log for writing." << std::endl;
    delete gmesh;
    return 1;
  }

  std::cout << "\n=== Starting adaptive refinement loop ===" << std::endl;

  while (refIt <= maxRefinements && estimatedError > errorTolerance) {
    std::cout << "\n--- Computing Hdiv and H1 approximations at iteration " << refIt << " ---" << std::endl;

    TPZVec<REAL> elementErrors(gmesh->NElements(), 0.);
    TPZMultiphysicsCompMesh* cmeshMixed = MixedDarcyCompMesh(gmesh, SimData, wallFlux);
    if (cmeshMixed->NEquations() >= 2000000) {
      std::cerr << "Error: Number of equations exceeds 2 million. Stopping refinement." << std::endl;
      DeleteMixedMesh(cmeshMixed);
      break;
    }

    // H(div) analysis
    {
      TPZLinearAnalysis anMixed(cmeshMixed);
      TPZWannStructMatrix matMixed(cmeshMixed);
      matMixed.SetNumThreads(SimData.m_Numerics.nthreads);
      anMixed.SetStructuralMatrix(matMixed);
      TPZStepSolver<STATE> stepMixed;
      stepMixed.SetDirect(ELDLt);
      anMixed.SetSolver(stepMixed);
      anMixed.Run();
    }

    TPZCompMesh* cmeshH1 = H1DarcyCompMesh(gmesh, SimData, wallFlux);
    // H1 analysis
    {
      TPZLinearAnalysis anH1(cmeshH1);
      TPZWannStructMatrix matH1(cmeshH1);
      matH1.SetNumThreads(SimData.m_Numerics.nthreads);
      anH1.SetStructuralMatrix(matH1);

      TPZStepSolver<STATE> stepH1;
      stepH1.SetDirect(ECholesky);
      anH1.SetSolver(stepH1);
      anH1.Run();
    }

    // Export the geometry and both reservoir solutions for this iteration.
    const std::string iteration = std::to_string(refIt);
    if (shouldPlot) {
      std::ofstream plotfile("geomesh_" + iteration + ".vtk");
      TPZVTKGeoMesh::PrintGMeshVTK(gmesh, plotfile);
      WriteDarcyVTK(cmeshH1, "h1Darcy_plot_" + iteration, SimData.m_Mesh.ToCylindrical);
      WriteDarcyVTK(cmeshMixed, "mixedDarcy_plot_" + iteration, SimData.m_Mesh.ToCylindrical);
    }

    std::cout << "\n--- Error estimation and adaptive refinement at iteration " << refIt << " ---" << std::endl;
    estimatedError = TPZWannAdaptivityTools::PragerSynge(
        cmeshMixed, cmeshH1, &SimData, elementErrors, SimData.m_Numerics.nthreads,
        "PragerSyngeEst_" + iteration + ".vtk");

    for (auto& connect : cmeshH1->ConnectVec()) connect.RemoveDepend();
    delete cmeshH1;
    DeleteMixedMesh(cmeshMixed);

    estimatedErrorVec.push_back(estimatedError);
    const REAL maxError = *std::max_element(elementErrors.begin(), elementErrors.end());
    maxErrorVec.push_back(maxError);
    std::cout << "Estimated error: " << estimatedError << " (scaled), "
              << estimatedError * errorToPhysical << " (physical)\n"
              << "Max element error: " << maxError << " (scaled), "
              << maxError * errorToPhysical << " (physical)" << std::endl;

    // Mesh adaptive refinement (in the last iteration we only compute the error, no refinement)
    if (refIt < maxRefinements && estimatedError > errorTolerance) {
      TPZVec<int64_t> refinedElements;
      TPZWannAdaptivityTools::AdaptivityProcess(gmesh, &SimData, elementErrors, refinedElements, relativeRefTol);

      // Export refinement information to file
      refinementLog << refinedElements.size() << " ";
      for (size_t i = 0; i < refinedElements.size(); ++i) {
        refinementLog << refinedElements[i];
        if (i != refinedElements.size() - 1) refinementLog << " ";
      }
      refinementLog << std::endl;
    }

    refIt++;
  }

  std::cout << "\n=== Adaptive refinement loop completed ===" << std::endl;
  std::cout << "\nError history\n"
            << std::setw(10) << "Iteration"
            << std::setw(20) << "Total (scaled)"
            << std::setw(20) << "Total (physical)"
            << std::setw(20) << "Max elem. (scaled)"
            << std::setw(22) << "Max elem. (physical)" << '\n';
  for (int i = 0; i < refIt; ++i) {
    std::cout << std::setw(10) << i
              << std::setw(20) << estimatedErrorVec[i]
              << std::setw(20) << estimatedErrorVec[i] * errorToPhysical
              << std::setw(20) << maxErrorVec[i]
              << std::setw(22) << maxErrorVec[i] * errorToPhysical << '\n';
  }

  refinementLog.close();
  delete gmesh;
  return 0;
}

void WriteDarcyVTK(TPZCompMesh* mesh, const std::string& filename, bool cylindrical) {
  const int vtkResolution = 0;
  TPZManVector<std::string, 2> fields = {"Flux", "Pressure"};
  TPZVTKGenerator vtk(mesh, fields, filename, vtkResolution, 3, cylindrical);
  vtk.Do();
}

void DeleteMixedMesh(TPZMultiphysicsCompMesh* mesh) {
  const auto atomicMeshes = mesh->MeshVector();
  for (auto& connect : mesh->ConnectVec()) connect.RemoveDepend();
  delete mesh;
  for (auto* atomicMesh : atomicMeshes) {
    for (auto& connect : atomicMesh->ConnectVec()) connect.RemoveDepend();
    delete atomicMesh;
  }
}

namespace {
template<class Creator, class Material>
void SetDarcyData(Creator& creator, Material* material, const ProblemData& data,
                  const std::map<int, REAL>& wallFlux, REAL sign) {
  TPZFNMatrix<9, STATE> mobility(3, 3, 0.);
  for (int i = 0; i < 3; ++i)
    mobility(i, i) = data.m_Reservoir.perm[i] / data.m_Fluid[0].viscosity;
  material->SetConstantPermeability(mobility);
  creator.InsertMaterialObject(material);

  TPZFMatrix<STATE> val1(1, 1, 0.);
  TPZManVector<STATE, 1> val2(1);
  for (const auto& [name, bc] : data.m_Reservoir.BCs) {
    val2[0] = bc.value;
    creator.InsertMaterialObject(material->CreateBC(material, bc.matid, bc.type, val1, val2));
  }
  for (const auto& [id, flux] : wallFlux) {
    val2[0] = sign * flux;
    creator.InsertMaterialObject(material->CreateBC(material, id, 1, val1, val2));
  }
}
}

TPZMultiphysicsCompMesh* MixedDarcyCompMesh(TPZGeoMesh* gmesh, const ProblemData& data,
                                           const std::map<int, REAL>& wallFlux) {
  TPZHDivApproxCreator creator(gmesh);
  creator.SetProbType(ProblemType::EDarcy);
  creator.SetHdivFamily(HDivFamily::EHDivStandard);
  creator.SetDefaultOrder(data.m_Numerics.reservoirPorder);
  creator.SetShouldCondense(true);
  auto* material = new TPZMixedDarcyAnisotropic(data.m_Reservoir.matid, gmesh->Dimension());
  SetDarcyData(creator, material, data, wallFlux, 1.);
  return creator.CreateApproximationSpace();
}

TPZCompMesh* H1DarcyCompMesh(TPZGeoMesh* gmesh, const ProblemData& data,
                            const std::map<int, REAL>& wallFlux) {
  TPZH1ApproxCreator creator(gmesh);
  creator.SetProbType(ProblemType::EDarcy);
  creator.SetDefaultOrder(data.m_Numerics.reservoirPorder + 2);
  creator.SetShouldCondense(true);
  auto* material = new TPZDarcyAnisotropic(data.m_Reservoir.matid, gmesh->Dimension());
  // Keep the pressure penalty relative to the Darcy stiffness at the input pressure scale.
  const auto& permeability = data.m_Reservoir.perm;
  const REAL maxMobility = *std::max_element(permeability.begin(), permeability.end()) /
                           data.m_Fluid[0].viscosity;
  material->SetBigNumber(material->BigNumber() * maxMobility);
  SetDarcyData(creator, material, data, wallFlux, -1.);
  return creator.CreateClassicH1ApproximationSpace();
}
