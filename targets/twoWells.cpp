#ifdef HAVE_CONFIG_H
#include <pz_config.h>
#endif

#include <iostream>
#include "TPZWannGeometryTools.h"
#include "TPZWannApproxTools.h"
#include "TPZWannPostProcTools.h"
#include <TPZLinearAnalysis.h>
#include <TPZSSpStructMatrix.h>
#include <pzskylstrmatrix.h>
#include <pzstepsolver.h>
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
  std::ofstream out("gmeshWannFinal.vtk");
  TPZVTKGeoMesh::PrintGMeshVTK(gmesh, out);

  // --- H(div) Simulation ---

  // For some reason we have to change the sign of boundary condition when using cylindrical map. 
  // TODO: fix it
  // TODO: If unable to fix it, adapt the workaround to multiple wellbores
  for (auto &Wellbore : SimData.m_Wellbore) {
    for (auto &BC : Wellbore.BCs) {
      BC.second.value *= -1.0;
    }
  }

  TPZMultiphysicsCompMesh* cmeshMixed = TPZWannApproxTools::CreateMultiphysicsCompMesh(gmesh, &SimData, &exact);
  TPZWannAnalysis anMixed(cmeshMixed, RenumType::EMetis);
  anMixed.SetProblemData(&SimData);
  anMixed.SetRescaling(false); // Enable matrix rescaling to improve conditioning
  anMixed.Initialize();
  anMixed.NewtonIteration();

  // Revert sign changes
  for (auto &Wellbore : SimData.m_Wellbore) {
    for (auto &BC : Wellbore.BCs) {
      BC.second.value *= -1.0;
    }
  }

  // --- Post-processing ----
  TPZWannPostProcTools::WriteVTKs(cmeshMixed, &SimData);

  delete cmeshMixed;
  delete gmesh;
}