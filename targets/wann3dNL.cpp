#ifdef HAVE_CONFIG_H
#include <pz_config.h>
#endif

#include <iostream>
#include "TPZWannGeometryTools.h"
#include "TPZWannApproxTools.h"
#include "TPZWannPostProcTools.h"
#include "TPZWannAnalysis.h"

int main(int argc, char *argv[]) {
  std::string jsonfile = "ozkan1999.json";

  if (argc > 2) {
    std::cout << argv[0] << " being called with too many arguments." << std::endl;
    DebugStop();
  } else if (argc == 2) {
    jsonfile = argv[1];
  }

  if (jsonfile.find(".json") == std::string::npos) {
    jsonfile += ".json";
  }

  std::cout << "\n--------- Starting simulation ---------" << std::endl;
  std::cout << "Using json file: " << jsonfile << std::endl;

#ifdef PZ_LOG
  TPZLogger::InitializePZLOG();
#endif

  // Problem data
  ProblemData SimData;  
  SimData.ReadJson(jsonfile);

  TPZGeoMesh* gmesh = nullptr;
  TPZMultiphysicsCompMesh* cmeshMixed = nullptr;
  TPZCompMesh* cmeshH1 = nullptr;

  // Read the geometric mesh
  gmesh = TPZWannGeometryTools::CreateGeoMesh(&SimData);

  // Mixed simulation
  cmeshMixed = TPZWannApproxTools::CreateMixedCompMesh(gmesh, &SimData);
  TPZWannAnalysis anMixed(cmeshMixed, RenumType::EMetis);
  anMixed.SetProblemData(&SimData);
  anMixed.Initialize();
  anMixed.NewtonIteration();

  // H1 simulation for comparison
  cmeshH1 = TPZWannApproxTools::CreateH1CompMesh(gmesh, &SimData);
  TPZWannAnalysis anH1(cmeshH1, RenumType::EMetis);
  anH1.SetProblemData(&SimData);
  anH1.Initialize();
  anH1.NewtonIteration();

  // --- Post-processing ---

  // Export VTK files for both simulations
  TPZWannPostProcTools::WriteVTKs(cmeshMixed, &SimData);
  TPZWannPostProcTools::WriteVTKs(cmeshH1, &SimData);

  // For mixed results, export CSV files with well profiles (flow, inflow, pressure drop)
  TPZWannPostProcTools::WriteWellProfilesCSV(cmeshMixed, SimData, 4);

  // Compute and display well fluxes and productivity indices
  for (int iw = 0; iw < SimData.m_Wellbore.size(); ++iw) {
    const auto &well = SimData.m_Wellbore[iw];
    std::cout << "\nWell " << well.name << '\n';
    const auto mixedFlux = TPZWannPostProcTools::ComputeWellFluxes(cmeshMixed, &SimData, iw, 1);
    const auto h1Flux = TPZWannPostProcTools::ComputeWellFluxes(cmeshH1, &SimData, iw, 1);
    std::cout << "Integrated well flux H(div): " << mixedFlux[0] << '\n'
              << "Integrated well flux H1: " << h1Flux[0] << '\n';

    const auto &heel = well.BCs.at("point_heel");
    const auto &toe = well.BCs.at("point_toe");
    if (heel.type == 1 && toe.type == 1) {
      const REAL expected = toe.value - heel.value;
      std::cout << "Fluid loss H(div): " << mixedFlux[0] - expected << '\n'
                << "Fluid loss H1: " << h1Flux[0] - expected << '\n';
    }

    std::cout << "Productivity Index H(div): "
              << TPZWannPostProcTools::ProductivityIndex(cmeshMixed, &SimData, iw) << '\n'
              << "Productivity Index H1: "
              << TPZWannPostProcTools::ProductivityIndex(cmeshH1, &SimData, iw) << '\n';
  }
  
  // --- Clean up ---

  // Delete H(div) meshes. The multiphysics mesh owns the atomic meshes, so we delete them first.
  const auto atomicMeshes = cmeshMixed->MeshVector();
  for (auto* mesh : atomicMeshes) delete mesh;
  delete cmeshMixed;

  // Delete H1 mesh. Delete the connect dependencies first to avoid warnings
  for (auto& connect : cmeshH1->ConnectVec()) connect.RemoveDepend();
  delete cmeshH1;

  return 0;
}
