#ifdef HAVE_CONFIG_H
#include <pz_config.h>
#endif

#include <iostream>
#include "TPZWannGeometryTools.h"
#include "TPZWannApproxTools.h"
#include "TPZWannPostProcTools.h"
#include "TPZWannEstimationTools.h"
#include <TPZLinearAnalysis.h>
#include <TPZSSpStructMatrix.h>  //symmetric sparse matrix storage>
#include <pzskylstrmatrix.h>
#include <pzstepsolver.h>
#include <TPZAnalyticSolution.h>
#include "TPZWannAnalysis.h"

using namespace std;

const int global_nthread = 10;

int main(int argc, char *argv[]) {
  
  TLaplaceExample1 exact; // Global variable to be used in the material objects
  exact.fDimension = 3;
  exact.fExact = TLaplaceExample1::ENone;

  std::string jsonfile = "wann3d_test.json";

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

  // Create computational meshes (H(div) mixed and H1)
  TPZGeoMesh* gmesh = TPZWannGeometryTools::CreateGeoMesh(&SimData);
  TPZMultiphysicsCompMesh* cmeshMixed = TPZWannApproxTools::CreateMultiphysicsCompMesh(gmesh, &SimData, &exact);
  TPZCompMesh* cmeshH1 = TPZWannApproxTools::CreateH1CompMesh(gmesh, &SimData);

  // Non-linear H(div) analysis
  TPZWannAnalysis anMixed(cmeshMixed, RenumType::EMetis);
  anMixed.SetProblemData(&SimData);
  anMixed.Initialize();
  anMixed.NewtonIteration();

  // Non-linear H1 analysis
  // TODO
  // Equation filter must be applied in the Dirichelet BCs

  std::cout << "\n--------- Simulation finished ---------" << std::endl;
  std::cout << "\n--------- Starting post-processing ---------" << std::endl;

  TPZWannPostProcTools::PostProcessAllData(cmeshMixed, gmesh, &SimData);

  // Integrated flux along segments of the well
  TPZVec<REAL> segmentPoints = {0.0, SimData.m_Wellbore.length / 2.0,
                                SimData.m_Wellbore.length};
  TPZVec<REAL> fluxes = TPZWannPostProcTools::ComputeWellFluxes(
      cmeshMixed, &SimData, segmentPoints);

  std::cout << "Integrated well fluxes H(div): ";
  for (int iseg = 0; iseg < fluxes.size(); iseg++) {
    std::cout << fluxes[iseg] << " ";
  }
  std::cout << std::endl;

  std::cout << "\n--------- Post-processing finished ---------" << std::endl;

  delete cmeshMixed;
  delete cmeshH1;
  delete gmesh;
}