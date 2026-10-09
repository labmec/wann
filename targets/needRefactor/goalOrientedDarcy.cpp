#ifdef HAVE_CONFIG_H
#include <pz_config.h>
#endif

#include <iostream>

#include <TPZLinearAnalysis.h>
#include <TPZSSpStructMatrix.h>
#include <pzskylstrmatrix.h>
#include <pzstepsolver.h>
#include <TPZAnalyticSolution.h>
#include <TPZH1ApproxCreator.h>
#include <TPZHDivApproxCreator.h>
#include <TPZVTKGenerator.h>
#include <Material/DarcyFlow/TPZMixedDarcyFlow.h>
#include <Material/DarcyFlow/TPZDarcyFlow.h>
#include <TPZFileStream.h>
#include <pzcondensedcompel.h>

#include "TPZWannGeometryTools.h"
#include "TPZWannPostProcTools.h"
#include "TPZWannAdaptivityTools.h"
#include "ProblemData.h"

// ================
// Global variables
// ================

REAL gLenght = -1.0; // To save the well lenght in a global variable (gambiarra)

// Dirichlet condition for the dual problem
auto DirichletFunctionDual = [](const TPZVec<REAL> &pt, TPZVec<STATE> &result, TPZFMatrix<STATE> &deriv) {
  result.Resize(1); // Ensure proper size
  deriv.Resize(3, 1); // 2D gradient but needs 3 components

  REAL x = pt[0];
  REAL y = pt[1];
  REAL z = pt[2];

  // Test 1
//   REAL a = gLenght / 3.0;
//   REAL b = 2.0 * gLenght / 3.0;
//   REAL k = 0.2;
//   result[0] = (1.0 / (1.0 + std::exp(-k * (x - a)))) *
//               (1.0 / (1.0 + std::exp(k * (x - b))));

  result[0] = 1.0; // Test 2

  // Dummy argument, not used in this application
  deriv(0, 0) = 0.0; // du/dx
  deriv(1, 0) = 0.0; // du/dy
  deriv(2, 0) = 0.0; // empty
};

// Forcing function for the dual problem
auto ForcingFunctionDual = [](const TPZVec<REAL> &pt, TPZVec<STATE> &result) {
  result.Resize(1); // Ensure proper size
  REAL x = pt[0];
  REAL y = pt[1];
  REAL z = pt[2];

  result[0] = 0.; // Test 1

  // Test 2
//   if (x >= gLenght/2.0 && x <= 3.0 * gLenght/4.0 && z >= 8.0 && z <= 12.0 && y >=-10.0 && y <= 10.0) {
//     result[0] = 1.; // divide by area of the square
//   }
};

// ===================
// Function prototypes
// ===================

// Computational mesh for mixed Darcy formulation
TPZMultiphysicsCompMesh* MixedDarcyCompMesh(TPZGeoMesh* gmesh, ProblemData* SimData, TLaplaceExample1* exact, bool isCondensed = false);

// Computational mesh for the Dual problem
TPZMultiphysicsCompMesh* MixedDarcyDualCompMesh(TPZGeoMesh* gmesh, ProblemData* SimData, bool isCondensed = false);

// Goal-oriented error estimation for 3D Darcy flow
REAL GoalOrientedDarcy(TPZMultiphysicsCompMesh* cmeshHdiv, TPZMultiphysicsCompMesh* cmeshDual, ProblemData* SimData, TPZVec<REAL>& elementErrors, int nthreads);

// Generate Problem data for the dual problem
ProblemData DualizeProblemData(const ProblemData &SimData);

// ============
// Main program
// ============

int main(int argc, char *argv[]) {
  
  // --- Set up ---

  REAL errorTolerance = 1e-3; // Error tolerance for adaptive refinement
  const int maxIterations = 5; // Maximum number of adaptive refinement iterations
  REAL theta = 0.2; // Marking parameter for adaptive refinement

  TLaplaceExample1 exact; // Exact solution (if known)
  exact.fDimension = 3;
  exact.fExact = TLaplaceExample1::ENone;

  std::string jsonfile = "simplestGoalDarcy.json";

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

#ifdef PZ_LOG
  TPZLogger::InitializePZLOG();
#endif

  // Problem data
  ProblemData SimData;  
  SimData.ReadJson(jsonfile);

  // Save well length in global variable (gambiarra)
  gLenght = SimData.m_Wellbore.length;

  int nThreads = SimData.m_Numerics.nthreads; // Number of threads for parallel computations

  // Initial geometric mesh
  TPZGeoMesh* gmesh = TPZWannGeometryTools::ReadMeshFromGmsh(&SimData);

  // Convert to cylindrical coordinates if set in json
  if (SimData.m_Mesh.ToCylindrical) {
    TPZManVector<REAL,3> cylcenter = {0.,0.,0.};
    REAL hr = SimData.m_Reservoir.height;
    REAL lr = SimData.m_Wellbore.height;
    cylcenter[2] = lr - hr/2.;
    TPZWannGeometryTools::ModifyGeometricMeshToCylWell(gmesh, SimData.ESurfWellCyl, SimData.m_Wellbore.radius, cylcenter);
  }

  // Problem data for the dual problem
  // Need to be done after reading the mesh (some fields are set in the mesh reading process)
  ProblemData SimDataDual = DualizeProblemData(SimData);

  // --- Adaptive refinement loop ---

  int refIt = 0;
  REAL estimatedError = errorTolerance + 1.0;
  TPZManVector<REAL, maxIterations> errorHistory(maxIterations, 0.0);
  TPZManVector<int64_t, maxIterations> nEquations(maxIterations, 0);

  std::cout << "\n=== Starting adaptive refinement loop ===" << std::endl;

  while (refIt < maxIterations && estimatedError > errorTolerance) {
    std::cout << "\n--- Computing Hdiv and Dual approximations at iteration " << refIt << " ---" << std::endl;

    // Create computational meshes
    TPZMultiphysicsCompMesh* cmeshMixed = MixedDarcyCompMesh(gmesh, &SimData, &exact, true);
    TPZMultiphysicsCompMesh* cmeshDual = MixedDarcyDualCompMesh(gmesh, &SimDataDual, true);

    // H(div) analysis
    TPZLinearAnalysis anMixed(cmeshMixed);
#ifdef PZ_USING_MKL
    TPZSSpStructMatrix<STATE> matMixed(cmeshMixed);
#else
    TPZSkylStrMatrix matMixed(cmeshMixed);
#endif
    matMixed.SetNumThreads(nThreads);
    anMixed.SetStructuralMatrix(matMixed);
    TPZStepSolver<STATE> stepMixed;
    stepMixed.SetDirect(ELDLt);
    anMixed.SetSolver(stepMixed);
    anMixed.Run();

    // Dual analysis
    TPZLinearAnalysis anDual(cmeshDual);
#ifdef PZ_USING_MKL
    TPZSSpStructMatrix<STATE> matDual(cmeshDual);
#else
    TPZSkylStrMatrix matDual(cmeshDual);
#endif
    matDual.SetNumThreads(nThreads);
    anDual.SetStructuralMatrix(matDual);
    TPZStepSolver<STATE> stepDual;
    stepDual.SetDirect(ELDLt);
    anDual.SetSolver(stepDual);
    anDual.Run();

    // Plot solutions (for internal control only)
    if (SimData.m_PostProc.verbosityLevel) {
      TPZWannPostProcTools::WriteReservoirVTK(cmeshMixed, &SimData);
      TPZWannPostProcTools::WriteReservoirVTK(cmeshDual, &SimDataDual);
    }

    std::cout << "\n--- Error estimation and adaptive refinement at iteration " << refIt << " ---" << std::endl;

    int64_t ngel = cmeshMixed->Reference()->NElements();
    TPZVec<REAL> elementErrors(ngel, 0.0);
    TPZVec<int> refinementIndicator(ngel, 0);

    REAL estimatedError = GoalOrientedDarcy(cmeshMixed, cmeshDual, &SimData, elementErrors, nThreads);
    std::cout << "Estimated error: " << estimatedError << std::endl;
    errorHistory[refIt] = estimatedError;
    nEquations[refIt] = cmeshMixed->NEquations();

    // Wrapper for the refinement process
    TPZWannAdaptivityTools::AdaptivityProcess(gmesh, &SimData, elementErrors, theta);

    // Uniform refinement for testing
    // TPZCheckGeom checkgeom(gmesh);
    // checkgeom.UniformRefine(1);

    // --- Clean up ---

    delete cmeshMixed;
    delete cmeshDual;
    refIt++;
  }

  std::cout << "Number of equations and Goal errors per iteration: " << std::endl;
  for (int i = 0; i < maxIterations; ++i) {
    std::cout << "    Iteration " << i << ": " << nEquations[i] << " equations, error: " << errorHistory[i] << std::endl;
  }
}

// ========================
// Functions implementation
// ========================

TPZMultiphysicsCompMesh *MixedDarcyCompMesh(TPZGeoMesh *gmesh, ProblemData *SimData, TLaplaceExample1 *exact, bool isCondensed) {
  TPZHDivApproxCreator hdivCreator(gmesh);
  hdivCreator.ProbType() = ProblemType::EDarcy;
  hdivCreator.HdivFamily() = HDivFamily::EHDivStandard;
  hdivCreator.SetDefaultOrder(SimData->m_Reservoir.pOrder);
  hdivCreator.SetShouldCondense(isCondensed);

  TLaplaceExample1* exactsol = dynamic_cast<TLaplaceExample1 *>(exact);
  bool hasAnalyticSol = (exactsol != nullptr && exactsol->fExact != TLaplaceExample1::ENone);

  TPZMixedDarcyFlow *reservoirMat = new TPZMixedDarcyFlow(SimData->EDomain, gmesh->Dimension());
  reservoirMat->SetConstantPermeability(SimData->m_Reservoir.perm[0]/SimData->m_Fluid.viscosity);
  if (hasAnalyticSol) {
    reservoirMat->SetExactSol(exact->ExactSolution(), 3);
    reservoirMat->SetForcingFunction(exact->ForceFunc(), 3);
    reservoirMat->SetConstantPermeability(1.0);
  }
  
  hdivCreator.InsertMaterialObject(reservoirMat);

  // Boundary conditions --- 

  // Hardcoded. We are not using the values from the .json file
  TPZFMatrix<STATE> val1(1, 1, 0.);
  TPZManVector<STATE> val2(1, 0);

  for (auto &bcpair : SimData->m_Reservoir.BCs) {
    auto &bc = bcpair.second;
    TPZFMatrix<STATE> val1(1, 1, 0.);
    TPZManVector<STATE> val2(1, 0);
    val2[0] = bc.value;
    TPZBndCondT<STATE> *BCond = reservoirMat->CreateBC(reservoirMat, bc.matid, bc.type, val1, val2);
    if (hasAnalyticSol) BCond->SetForcingFunctionBC(exact->ExactSolution(), 3);
    hdivCreator.InsertMaterialObject(BCond);
  }

  TPZMultiphysicsCompMesh *cmesh = hdivCreator.CreateApproximationSpace();
  return cmesh;
}

TPZMultiphysicsCompMesh *MixedDarcyDualCompMesh(TPZGeoMesh *gmesh, ProblemData *SimDataDual, bool isCondensed) {
  TPZHDivApproxCreator hdivCreator(gmesh);
  hdivCreator.ProbType() = ProblemType::EDarcy;
  hdivCreator.HdivFamily() = HDivFamily::EHDivStandard;
  hdivCreator.SetDefaultOrder(SimDataDual->m_Reservoir.pOrder);
  hdivCreator.SetShouldCondense(isCondensed);

  TPZMixedDarcyFlow *reservoirMat = new TPZMixedDarcyFlow(SimDataDual->EDomain, gmesh->Dimension());
  reservoirMat->SetConstantPermeability(SimDataDual->m_Reservoir.perm[0]/SimDataDual->m_Fluid.viscosity);
  reservoirMat->SetForcingFunction(ForcingFunctionDual, 3);
  hdivCreator.InsertMaterialObject(reservoirMat);

  TPZFMatrix<STATE> val1(1, 1, 0.);
  TPZManVector<STATE> val2(1, 0);

  for (auto &bcpair : SimDataDual->m_Reservoir.BCs)
    {
      auto &bc = bcpair.second;
      TPZFMatrix<STATE> val1(1, 1, 0.);
      TPZManVector<STATE> val2(1, 0);
      val2[0] = bc.value;
      TPZBndCondT<STATE> *BCond = reservoirMat->CreateBC(reservoirMat, bc.matid, bc.type, val1, val2);

      // We only put a dirichlet dual condition on the well surface
      if (bc.type == 0 && bc.matid == SimDataDual->ESurfWellCyl) {
        BCond->SetForcingFunctionBC(DirichletFunctionDual, 3);
      }
      hdivCreator.InsertMaterialObject(BCond);
    }

  TPZMultiphysicsCompMesh *cmesh = hdivCreator.CreateApproximationSpace();
  return cmesh;
}

REAL GoalOrientedDarcy(TPZMultiphysicsCompMesh* cmeshHdiv, TPZMultiphysicsCompMesh* cmeshDual, ProblemData* SimData, TPZVec<REAL>& elementErrors, int nthreads) {    
  nthreads++; // If nthreads = 0, we use 1 thread
  int dim = cmeshHdiv->Dimension();

  // Ensure references point to dual cmesh
  cmeshDual->Reference()->ResetReference();
  cmeshDual->LoadReferences();

  // Resizes elementContributions and fill with zeros
  int64_t ngel = cmeshHdiv->Reference()->NElements();
  elementErrors.Resize(ngel);
  elementErrors.Fill(0);

  int64_t ncel = cmeshHdiv->NElements();
  TPZAdmChunkVector<TPZCompEl *> &elementvec_m = cmeshHdiv->ElementVec();
  TPZVec<REAL> elementContributionsAux(ncel, 0.0); // Auxiliary vector for plotting

  // Parallelization setup
  std::vector<std::thread> threads(nthreads);
  TPZManVector<REAL> partialErrors(nthreads, 0.0);

  auto worker = [&](int tid, int64_t start, int64_t end) {
    REAL localTotalError = 0.0;
    for (int64_t icel = start; icel < end; ++icel) {
      TPZCompEl *celMixed = elementvec_m[icel];

      // Check if element is condensed
      TPZCondensedCompEl *condEl = dynamic_cast<TPZCondensedCompEl *>(celMixed);
      if (condEl) {
        // If compel is condensed, load solution on the unconsensed compel
        condEl->LoadSolution();
        celMixed = condEl->ReferenceCompEl();
      }

      int matid = celMixed->Material()->Id();
      if (matid != SimData->EDomain) continue;

      TPZGeoEl *gel = celMixed->Reference();
      if (gel->HasSubElement()) continue;
      
      TPZCompEl *celDual = gel->Reference();

      // Check if element is condensed
      TPZCondensedCompEl *condElDual = dynamic_cast<TPZCondensedCompEl *>(celDual);
      if (condElDual) {
        // If compel is condensed, load solution on the unconsensed compel
        condElDual->LoadSolution();
        celDual = condElDual->ReferenceCompEl();
      }

      if (!celMixed || !celDual) continue;
      if (celDual->Material()->Id() != matid) DebugStop();

      TPZManVector<REAL, 3> perm = SimData->m_Reservoir.perm;
      for (int i = 0; i < 3; ++i) {
        perm[i] = perm[i]/SimData->m_Fluid.viscosity;
      }

      REAL goalError = 0.0;

      // Set integration rule
      const TPZIntPoints* intrule = gel->CreateSideIntegrationRule(gel->NSides() - 1, 5);

      for (int ip = 0; ip < intrule->NPoints(); ++ip) {
        TPZManVector<REAL,3> ptInElement(gel->Dimension());
        REAL weight, detjac;
        intrule->Point(ip, ptInElement, weight);
        TPZFNMatrix<9, REAL> jacobian, axes, jacinv;
        gel->Jacobian(ptInElement, jacobian, axes, detjac, jacinv);
        weight *= detjac; // fabs(detjac);

        TPZManVector<REAL, 3> x(3, 0.0);
        gel->X(ptInElement, x); // Real coordinates for x

        // No source function on reservoir
        TPZManVector<REAL,3> force(1,0.0);

        // Compute mixed solution ph and sigh
        TPZManVector<REAL,3> ph(gel->Dimension(),0.0);
        TPZManVector<REAL,3> sigh(gel->Dimension(),0.0);
        TPZManVector<REAL,1> divsigh(1,0.0);
        celMixed->Solution(ptInElement, 2, ph);
        celMixed->Solution(ptInElement, 1, sigh);
        celMixed->Solution(ptInElement, 5, divsigh);

        // Compute dual solution zh and psih
        TPZManVector<REAL,3> zh(gel->Dimension(),0.0);
        TPZManVector<REAL,3> psih(gel->Dimension(),0.0);
        TPZManVector<REAL,1> divpsih(1,0.0);
        celDual->Solution(ptInElement, 2, zh);
        celDual->Solution(ptInElement, 1, psih);
        celDual->Solution(ptInElement, 5, divpsih);

        // F contribution
        REAL FTerm = -force[0]*zh[0];

        // Flux contribution
        REAL FluxTerm = 0.0;
        for (int d = 0; d < dim; ++d) {
          FluxTerm += (1./perm[d])*sigh[d]*psih[d];
        }

        // div contributions
        REAL divTermA = -divsigh[0]*zh[0];
        REAL divTermB = -divpsih[0]*ph[0];

        REAL ATerm = FluxTerm + divTermA + divTermB;

        // Minus sign due to the definition of the dual problem
        goalError += (FTerm - ATerm) * weight;
      }

      // Boundary contributions TODO: Improve this
      std::set<int> bcIds;
      std::map<int, REAL> bcValues;

      // Gather the Dirichlet boundary conditions of direct problem
       for (auto& bc : SimData->m_Reservoir.BCs) {
        if (bc.second.type == 0) { // Dirichlet
          bcIds.insert(bc.second.matid);
          bcValues[bc.second.matid] = bc.second.value;
        }
      }

      int firstSide = gel->FirstSide(dim-1);
      int lastSide = gel->FirstSide(dim);

      for (int side = firstSide; side < lastSide; ++side) {
        TPZGeoElSide gelSide(gel, side);
        TPZGeoElSide neigh = gelSide.HasNeighbour(bcIds);
        if (neigh) {
          TPZGeoEl *gelB = neigh.Element();
          if (gelB->HasSubElement()) DebugStop();
          TPZCompEl *celDualB = gelB->Reference();
          if (!celDualB) DebugStop();
          if (gelB->HasSubElement()) DebugStop();

          int neighMatid = gelB->MaterialId();
          REAL pe_cte = bcValues[neighMatid];

          // Check if element is condensed
          TPZCondensedCompEl *condEl = dynamic_cast<TPZCondensedCompEl *>(celDualB);
          if (condEl) {
            // If compel is condensed, load solution on the unconsensed compel
            condEl->LoadSolution();
            celDualB = condEl->ReferenceCompEl();
          }

          const TPZIntPoints* intruleSide = gelB->CreateSideIntegrationRule(gelB->NSides() - 1, 5);
          // const TPZIntPoints* intruleSide = &celDualB->GetIntegrationRule();

          for (int ip = 0; ip < intruleSide->NPoints(); ++ip) {
            TPZManVector<REAL,3> ptInElement(gelB->Dimension());
            REAL weight, detjac;
            intruleSide->Point(ip, ptInElement, weight);
            TPZFNMatrix<9, REAL> jacobian, axes, jacinv;
            gelB->Jacobian(ptInElement, jacobian, axes, detjac, jacinv);
            weight *= detjac;

            TPZManVector<REAL, 3> x(3, 0.0);
            gelB->X(ptInElement, x); // Real coordinates for x

            // Compute dual solution psih
            TPZManVector<REAL, 3> psih(1, 0.0);
            celDualB->Solution(ptInElement, 18, psih);

            goalError += -pe_cte*psih[0]*weight;
          }
        }
      }

      int64_t igeo = cmeshHdiv->Element(icel)->Reference()->Index();
      elementErrors[igeo] = std::abs(goalError);
      elementContributionsAux[icel] = std::abs(goalError); // For plotting
      localTotalError += goalError;
    }
    partialErrors[tid] = localTotalError;
  };

  int64_t chunk = ncel/nthreads;
  for (int t = 0; t < nthreads; ++t) {
    int64_t start = t * chunk;
    int64_t end = (t == nthreads - 1) ? ncel : (t + 1) * chunk;
    threads[t] = std::thread(worker, t, start, end);
  }
  for (auto& th : threads) th.join();

  REAL totalError = 0.0;
  for (auto val : partialErrors) totalError += val;

  // VTK output
  if (SimData->m_PostProc.verbosityLevel) {
    std::ofstream out_estimator("goalEstimation.vtk");
    TPZVTKGeoMesh::PrintCMeshVTK(cmeshHdiv, out_estimator, elementContributionsAux,
                                 "EstimatedError");
  }

  return totalError;
}

ProblemData DualizeProblemData(const ProblemData &SimData) {
    ProblemData dual = SimData;
    dual.m_PostProc.reservoir_vtk = SimData.m_PostProc.reservoir_vtk + "_dual";
    dual.m_PostProc.wellbore_vtk = SimData.m_PostProc.wellbore_vtk + "_dual";
    dual.m_Reservoir.pOrder = SimData.m_Reservoir.pOrder;
    dual.m_Wellbore.pOrder = SimData.m_Wellbore.pOrder;

    // Essential boundary conditions need to be homogeneous.
    // Here we assume that we are using only mixed formulations,
    // so the essential boundary conditions are the Neumann ones (type 1).

    // For the chosen functional, the natural boundary conditions
    // (Dirichlet ones, type 0) are also set to zero.
    // For other functionals, different Dirichlet boundary conditions may be used.

    for (auto &bcpair : dual.m_Wellbore.BCs) {
        if (bcpair.second.type == 1) {
            // bcpair.second.type = 0;
            bcpair.second.value = 0.0;
        } else if (bcpair.second.type == 0) {
            bcpair.second.value = 0.0;
        }
    }

    for (auto &bcpair : dual.m_Reservoir.BCs) {
        if (bcpair.second.type == 1) {
            bcpair.second.value = 0.0;
        } else if (bcpair.second.type == 0) {
            bcpair.second.value = 0.0;
        }
    }
    return dual;
}