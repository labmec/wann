#ifdef HAVE_CONFIG_H
#include <pz_config.h>
#endif

#include <iostream>
#include "TPZWannGeometryTools.h"
#include "TPZWannApproxTools.h"
#include "TPZWannPostProcTools.h"
#include "TPZWannAdaptivityTools.h"
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

  std::string jsonfile = "ozkan1999Box.json";

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

  // H1 Simulation
  TPZCompMesh* cmeshH1 = TPZWannApproxTools::CreateH1CompMesh(gmesh, &SimData, &exact);
  TPZWannAnalysis anH1(cmeshH1, RenumType::EMetis);
  anH1.SetProblemData(&SimData);
  anH1.SetRescaling(false); // Enable matrix rescaling to improve conditioning
  anH1.Initialize();
  anH1.NewtonIteration();

  std::cout << "\n--------- Simulation finished ---------" << std::endl;
  std::cout << "\n--------- Starting post-processing ---------" << std::endl;

  TPZWannPostProcTools::WriteVTKs(cmeshMixed, &SimData);
  TPZWannPostProcTools::WriteVTKs(cmeshH1, &SimData);

  // Integrated flux along segments of the well
  TPZVec<REAL> segmentPoints = {0.0, SimData.m_Wellbore.length};

  TPZVec<REAL> fluxes = TPZWannPostProcTools::ComputeWellFluxes(cmeshMixed, &SimData, segmentPoints);
  TPZVec<REAL> fluxesH1 = TPZWannPostProcTools::ComputeWellFluxes(cmeshH1, &SimData, segmentPoints);

  std::cout << "Integrated well fluxes H(div): ";
  REAL inflow = 0.0;
  for (int iseg = 0; iseg < fluxes.size(); iseg++) {
    std::cout << fluxes[iseg] << " ";
    inflow += fluxes[iseg];
  }
  std::cout << std::endl;
  std::cout << "Fluid loss H(div): " << std::abs(inflow) - std::abs(SimData.m_Wellbore.BCs["point_heel"].value) << std::endl;

  std::cout << "Integrated well fluxes H1: ";
  REAL inflowH1 = 0.0;
  for (int iseg = 0; iseg < fluxesH1.size(); iseg++) {
    std::cout << fluxesH1[iseg] << " ";
    inflowH1 += fluxesH1[iseg];
  }
  std::cout << std::endl;
  std::cout << "Fluid loss H1: " << std::abs(inflowH1) - std::abs(SimData.m_Wellbore.BCs["point_heel"].value) << std::endl;

  // Computing productivity index
  REAL PI = TPZWannPostProcTools::ProductivityIndex(cmeshMixed, &SimData);
  REAL PI_H1 = TPZWannPostProcTools::ProductivityIndex(cmeshH1, &SimData);

  std::cout << "Productivity Index H(div): " << PI << std::endl;
  std::cout << "Productivity Index H1: " << PI_H1 << std::endl;

  std::cout << "\n--------- Post-processing finished ---------" << std::endl;

  // --- Clean up ---

  // Remove dependencies before deleting H1 mesh
  int ncon = cmeshH1->NConnects();
  for (int i = 0; i < ncon; ++i) {
    cmeshH1->ConnectVec()[i].RemoveDepend();
  }

  delete cmeshMixed;
  delete cmeshH1;
  delete gmesh;
}