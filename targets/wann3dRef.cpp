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

using namespace std;

const int global_nthread = 0;


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
  std::cout << "--------- Starting simulation ---------" << std::endl;
#ifdef PZ_LOG
  TPZLogger::InitializePZLOG();
#endif

  // Problem data
  ProblemData SimData;  
  SimData.ReadJson(jsonfile);
  
  TPZGeoMesh* gmesh = TPZWannGeometryTools::CreateGeoMesh(&SimData);
  
  // Refinement loop (experimental)
  for (int refIt = 0; refIt < 1; refIt++) {
    TPZMultiphysicsCompMesh* cmesh = TPZWannApproxTools::CreateMultiphysicsCompMesh(gmesh, &SimData, &exact);
    TPZCompMesh* cmeshH1 = TPZWannApproxTools::CreateH1CompMesh(gmesh, &SimData);

    // {
    //   std::ofstream out("cmesh.txt");
    //   cmesh->Print(out);
    // }

    if (refIt >= 0) {
    // ---- Hdiv analysis ----
    TPZLinearAnalysis analysis(cmesh);
    #ifdef PZ_USING_MKL
    TPZSSpStructMatrix<STATE> skylstr(cmesh);
    #else
    TPZSkylStrMatrix skylstr(cmesh);
    #endif
    // TPZSkylStrMatrix skylstr(cmesh);
    // TPZSkylineStructMatrix<STATE> skylstr(cmesh);
    skylstr.SetNumThreads(global_nthread);
    analysis.SetStructuralMatrix(skylstr);

    TPZStepSolver<STATE> step;
    step.SetDirect(ELDLt);
    analysis.SetSolver(step);
    analysis.Run();

    {
      std::ofstream out("solution.txt");
      TPZFMatrix<REAL>& sol = analysis.Solution();
      sol.Print("Solution", out, EMathematicaInput);
    }
    }

    // if (refIt == 1) {

    //   int64_t ncon = cmesh->NConnects();
    //   for (int64_t icon = 0; icon < ncon; icon++)
    //   {
    //     TPZConnect &c = cmesh->ConnectVec()[icon];
    //     if (c.HasDependency())
    //       continue;
    //     if (!c.LagrangeMultiplier()) // Flux
    //     {
    //       int nshape = c.NShape();
    //       int blocksize = c.NShape() * c.NState();
    //       int64_t seqnum = c.SequenceNumber();
    //       if (seqnum < 0)
    //         continue;
    //       int64_t pos = cmesh->Block().Position(seqnum);
    //       for (int i = 0; i < blocksize; i++)
    //       {
    //         TPZFMatrix<REAL> &sol = analysis.Solution();
    //         sol(pos + i, 0) = 0.0;
    //       }
    //     }
    //     else // Pressure
    //     {
    //       int nshape = c.NShape();
    //       int blocksize = c.NShape() * c.NState();
    //       int64_t seqnum = c.SequenceNumber();
    //       if (seqnum < 0)
    //         continue;
    //       int64_t pos = cmesh->Block().Position(seqnum);
    //       if (c.Order() == 1) // vertex functions
    //       {
    //         for (int i = 0; i < blocksize; i++)
    //         {
    //           TPZFMatrix<REAL> &sol = analysis.Solution();
    //           sol(pos + i, 0) = 1.0;
    //         }
    //       }
    //       else // internal functions
    //       {
    //         for (int i = 0; i < blocksize; i++)
    //         {
    //           TPZFMatrix<REAL> &sol = analysis.Solution();
    //           sol(pos + i, 0) = 0.0;
    //         }
    //       }
    //     }
    //   }

    //   TPZFMatrix<REAL>& sol = analysis.Solution();
    //   {
    //     std::ofstream out("solution2.txt");
    //     sol.Print("Solution", out, EMathematicaInput);
    //   }
    //   TPZFMatrix<REAL>& rhs = analysis.Rhs();
    //   auto mat = analysis.MatrixSolver<REAL>().Matrix();

    //   TPZFMatrix<REAL> residual = rhs;
    //   mat->MultAdd(sol, rhs, residual, 1., -1.);
    //   analysis.LoadSolution(residual);
    //   cmesh->TransferMultiphysicsSolution();
    // }


    // ---- H1 analysis ----
    TPZLinearAnalysis analysisH1(cmeshH1);
    #ifdef PZ_USING_MKL
    TPZSSpStructMatrix<STATE> skylstrH1(cmeshH1);
    #else
    TPZSkylStrMatrix skylstrH1(cmeshH1);
    #endif
    skylstrH1.SetNumThreads(global_nthread);
    analysisH1.SetStructuralMatrix(skylstrH1);

    TPZStepSolver<STATE> stepH1;
    stepH1.SetDirect(ECholesky);
    analysisH1.SetSolver(stepH1);
    analysisH1.Run();

    std::cout << "--------- Simulation finished ---------" << std::endl;
    std::cout << "\n--------- Starting post-processing ---------" << std::endl;

    TPZWannPostProcTools::PostProcessAllData(cmesh, gmesh, &SimData);
    TPZWannPostProcTools::WriteVTKs(cmeshH1, &SimData);

    // Integrated flux along segments of the well
    TPZVec<REAL> segmentPoints = {0.0, SimData.m_Wellbore.length/2.0, SimData.m_Wellbore.length};
    TPZVec<REAL> fluxes = TPZWannPostProcTools::ComputeWellFluxes(cmesh, &SimData, segmentPoints);
    TPZVec<REAL> fluxesH1 = TPZWannPostProcTools::ComputeWellFluxes(cmeshH1, &SimData, segmentPoints);
    
    std::cout << "Integrated well fluxes H(div): ";
    for (int iseg = 0; iseg < fluxes.size(); iseg++) {
      std::cout << fluxes[iseg] << " ";
    }
    std::cout << std::endl;
    std::cout << "Integrated well fluxes H1:    ";
    for (int iseg = 0; iseg < fluxesH1.size(); iseg++) {
      std::cout << fluxesH1[iseg] << " ";
    }
    std::cout << std::endl;

    std::cout << "\n--------- Post-processing finished ---------" << std::endl;
    std::cout << "\n--------- Starting estimate and refine ---------" << std::endl;

    TPZWannEstimationTools::EstimateAndRefine(cmesh, cmeshH1, &SimData, global_nthread);

    std::cout << "\n--------- Estimate and refine finished ---------" << std::endl;

    delete cmesh;
    delete cmeshH1;
  }
  delete gmesh;
}