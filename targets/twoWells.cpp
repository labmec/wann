#ifdef HAVE_CONFIG_H
#include <pz_config.h>
#endif

#include <iostream>
#include "TPZWannGeometryTools.h"
#include "TPZWannApproxTools.h"
#include "TPZWannPostProcTools.h"
#include <TPZAnalyticSolution.h>
#include "TPZWannAnalysis.h"

using namespace std;

int main(int argc, char *argv[]) {
  
  TLaplaceExample1 exact; // Global variable to be used in the material objects
  exact.fDimension = 3;
  exact.fExact = TLaplaceExample1::ENone;

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

  // Read original geometric mesh and perform the refinement process 
  // described in refinementProcess.txt file
  TPZGeoMesh* gmesh = TPZWannGeometryTools::CreateGeoMesh(&SimData);
  std::ofstream out("gmeshWannFinal.vtk");
  TPZVTKGeoMesh::PrintGMeshVTK(gmesh, out);

  // --- H(div) Simulation ---

  // Mixed and H1 comp meshes
  TPZMultiphysicsCompMesh* cmeshMixed = TPZWannApproxTools::CreateMixedCompMesh(gmesh, &SimData);
  TPZCompMesh* cmeshH1 = TPZWannApproxTools::CreateH1CompMesh(gmesh, &SimData);

  // Mixed analysis
  TPZWannAnalysis anMixed(cmeshMixed, RenumType::EMetis);
  anMixed.SetProblemData(&SimData);
  anMixed.Initialize();
  anMixed.NewtonIteration();

  // H1 analysis
  TPZWannAnalysis anH1(cmeshH1, RenumType::EMetis);
  anH1.SetProblemData(&SimData);
  anH1.Initialize();
  anH1.NewtonIteration();

  // --- Post-processing ----
  TPZWannPostProcTools::WriteVTKs(cmeshMixed, &SimData);
  TPZWannPostProcTools::WriteVTKs(cmeshH1, &SimData);

  // ---- Clean up ----

  // To avoid pesky messages, delete all dependencies in the H1 mesh before deleting it.
  for (int64_t ic = 0; ic < cmeshH1->NConnects(); ic++) {
    cmeshH1->ConnectVec()[ic].RemoveDepend();
  }

  delete cmeshMixed;
  delete cmeshH1;
  delete gmesh;
}