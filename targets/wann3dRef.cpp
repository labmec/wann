#ifdef HAVE_CONFIG_H
#include <pz_config.h>
#endif

#include <iostream>
#include "TPZWannGeometryTools.h"
#include "TPZWannApproxTools.h"
#include "TPZWannPostProcTools.h"
#include "TPZWannEstimationTools.h"
#include <TPZLinearAnalysis.h>
#include <TPZSSpStructMatrix.h>  //symmetric sparse matrix storage
#include <pzstepsolver.h>

using namespace std;

const int global_nthread = 1;

int main(int argc, char *argv[]) {

  std::string jsonfile = "case_1.json";
  jsonfile = "wann3d.json";
  // jsonfile = "wann3d_test.json";

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
  std::cout << "--------- Starting simulation ---------" << std::endl;
#ifdef PZ_LOG
  TPZLogger::InitializePZLOG();
#endif

  // Problem data
  ProblemData SimData;  
  SimData.ReadJson(jsonfile);
  
  TPZGeoMesh* gmesh = TPZWannGeometryTools::CreateGeoMesh(&SimData);
  
  // Refinement loop (experimental)
  for (int refIt = 0; refIt < 2; refIt++) {
    TPZMultiphysicsCompMesh* cmesh = TPZWannApproxTools::CreateMultiphysicsCompMesh(gmesh, &SimData);
    TPZCompMesh* cmeshH1 = TPZWannApproxTools::CreateH1CompMesh(gmesh, &SimData);

    // ---- Hdiv analysis ----
    TPZLinearAnalysis analysis(cmesh);
    TPZSSpStructMatrix<STATE> skylstr(cmesh);
    skylstr.SetNumThreads(global_nthread);
    analysis.SetStructuralMatrix(skylstr);

    TPZStepSolver<STATE> step;
    step.SetDirect(ELDLt);
    analysis.SetSolver(step);
    analysis.Run();

    // ---- H1 analysis ----
    TPZLinearAnalysis analysisH1(cmeshH1);
    TPZSSpStructMatrix<STATE> skylstrH1(cmeshH1);
    skylstrH1.SetNumThreads(global_nthread);
    analysisH1.SetStructuralMatrix(skylstrH1);

    TPZStepSolver<STATE> stepH1;
    stepH1.SetDirect(ECholesky);
    analysisH1.SetSolver(stepH1);
    analysisH1.Run();

    std::cout << "--------- Simulation finished ---------" << std::endl;

    // ----
  
    TPZWannPostProcTools::PostProcessAllData(cmesh, gmesh, &SimData);
    //TPZWannPostProcTools::WriteVTKs(cmeshH1, &SimData);

    std::cout << "\n--------- Starting estimate and refine ---------" << std::endl;

    TPZWannEstimationTools::EstimateAndRefine(cmesh, cmeshH1, &SimData, global_nthread);

    std::cout << "--------- Estimate and refine finished ---------" << std::endl;

    delete cmesh;
    delete cmeshH1;
  }
  delete gmesh;
}