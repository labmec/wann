#ifdef HAVE_CONFIG_H
#include <pz_config.h>
#endif

#include <iostream>
#include <TPZAnalyticSolution.h>
#include "TPZWannGeometryTools.h"
#include "TPZWannApproxTools.h"
#include "TPZWannPostProcTools.h"
#include "TPZWannAdaptivityTools.h"
#include "TPZWannAnalysis.h"

using namespace std;

const int global_nthread = 12;

int main(int argc, char *argv[]) {
  
  TLaplaceExample1 exact; // Global variable to be used in the material objects
  exact.fDimension = 3;
  exact.fExact = TLaplaceExample1::ENone;

  std::string jsonfile = "wann3D_parameterStudy.json";

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

  // Pressure scale factor
  // REAL pFactor = 1e-5;
  // SimData.m_Fluid.viscosity *= pFactor;
  // SimData.m_Fluid.density *= pFactor;

  // Read original geometric mesh and perform the refinement process 
  // described in refinementProcess.txt file
  TPZGeoMesh* gmesh = TPZWannGeometryTools::CreateGeoMesh(&SimData);

  // H(div) Simulation

  // For some reason we have to change the sign of boundary condition 
  // when using cylindrical map. TODO: fix it
  if (SimData.m_Mesh.ToCylindrical) {
    SimData.m_Wellbore.BCs["point_heel"].value *= -1;
  }

  TPZMultiphysicsCompMesh* cmeshMixed = TPZWannApproxTools::CreateMultiphysicsCompMesh(gmesh, &SimData, &exact);
  TPZWannAnalysis anMixed(cmeshMixed, RenumType::EMetis);
  anMixed.SetProblemData(&SimData);
  anMixed.SetRescaling(false); // Enable matrix rescaling to improve conditioning
  anMixed.Initialize();
  anMixed.NewtonIteration();

  // Reverse sign changes
  if (SimData.m_Mesh.ToCylindrical) {
    SimData.m_Wellbore.BCs["point_heel"].value *= -1;
  }

  std::cout << "\n--------- Simulation finished ---------" << std::endl;
  std::cout << "\n--------- Starting post-processing ---------" << std::endl;

  TPZWannPostProcTools::WriteVTKs(cmeshMixed, &SimData);
  REAL ProductivityIndex = TPZWannPostProcTools::ProductivityIndex(cmeshMixed, &SimData);
  std::cout << "Computed productivity index: " << ProductivityIndex << std::endl;

  std::cout << "\n--------- Post-processing finished ---------" << std::endl;

  // --- Clean up ---

  delete cmeshMixed;
  delete gmesh;
}