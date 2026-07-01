#include <iostream>
#include <random>
#include <thread>
#include <numeric>
#include "TPZRefPatternDataBase.h"
#include "TPZVTKGeoMesh.h"
#include "TPZVTKGenerator.h"
#include "DarcyFlow/TPZDarcyFlow.h"
#include "DarcyFlow/TPZMixedDarcyFlow.h"
#include "TPZMultiphysicsCompMesh.h"
#include "pzmultiphysicselement.h"
#include "pzcondensedcompel.h"
#include "pzlog.h"
#include "pzintel.h"
#include "TPZVTKGeoMesh.h"

// Wann imports
#include "TPZWannAdaptivityTools.h"
#include "TPZWannGeometryTools.h"

REAL TPZWannAdaptivityTools::ErrorEstimation(TPZMultiphysicsCompMesh* cmeshMixed, TPZCompMesh* cmesh, ProblemData* SimData, TPZVec<REAL>& elementErrors, int nthreads) {

  {
    std::ofstream clearlog("error_estimation.txt", std::ios_base::trunc);
  }

  nthreads++; // If nthreads = 0, we use 1 thread

  // Ensure references point to H1 cmesh
  cmesh->Reference()->ResetReference();
  cmesh->LoadReferences();

  int64_t ncel = cmeshMixed->NElements();
  int64_t ngel = cmeshMixed->Reference()->NElements();
  TPZAdmChunkVector<TPZCompEl *> &elementvec_m = cmeshMixed->ElementVec();
  elementErrors.Resize(ngel); // Ensure proper size
  elementErrors.Fill(0.0); 

  TPZVec<REAL> elementErrorsAux(ncel, 0.0); // Aux vector for plotting

  // Parallelization setup
  std::vector<std::thread> threads(nthreads);
  TPZManVector<REAL> partialErrors(nthreads, 0.0);

  auto worker = [&](int tid, int64_t start, int64_t end) {
    REAL localTotalError = 0.0;
    for (int64_t icel = start; icel < end; ++icel) {
      TPZCompEl *celMixed = elementvec_m[icel];
      int matid = celMixed->Material()->Id();
      if (matid != SimData->EDomain && matid != SimData->ECurveWell) continue;

      TPZGeoEl *gel = celMixed->Reference();
      TPZCompEl *celH1 = gel->Reference();
      if (!celMixed || !celH1) continue;
      if (gel->HasSubElement()) continue;
      if (celH1->Material()->Id() != matid) DebugStop();

      REAL hk = ElementDiameter(gel);
      REAL perm = SimData->m_Reservoir.perm[0];
      if (matid == SimData->ECurveWell) {
        perm = SimData->m_Wellbore.perm[0];
      }
      REAL sqrtPerm = sqrt(perm);

      REAL fluxError = 0.0;
      REAL balanceError = 0.0;

      // Set integration rule
      const TPZIntPoints* intrule = nullptr;
      const TPZIntPoints &intruleMixed = celMixed->GetIntegrationRule();
      const TPZIntPoints &intruleH1 = celH1->GetIntegrationRule();
      if (intruleMixed.NPoints() < intruleH1.NPoints()) {
        intrule = &intruleH1;
      } else {
        intrule = &intruleMixed;
      }

      for (int ip = 0; ip < intrule->NPoints(); ++ip) {
        TPZManVector<REAL,3> ptInElement(gel->Dimension());
        REAL weight, detjac;
        intrule->Point(ip, ptInElement, weight);
        TPZFNMatrix<9, REAL> jacobian, axes, jacinv;
        gel->Jacobian(ptInElement, jacobian, axes, detjac, jacinv);
        weight *= fabs(detjac);

        TPZManVector<REAL,3> force(1,0.0);

        // Compute H1 term K^(1/2) grad(u)
        TPZManVector<REAL,3> termH1(gel->Dimension(),0.0);
        celH1->Solution(ptInElement, 2, termH1);
        for (int d = 0; d < termH1.size(); ++d) {
          termH1[d] = sqrtPerm * termH1[d];
        }

        // Compute Hdiv term -K^(-1/2) sig and div(sig)
        TPZManVector<REAL,3> termHdiv(3,0.0);
        TPZManVector<REAL,1> divFluxMixed(1,0.0);
        TPZMultiphysicsElement *celMulti = dynamic_cast<TPZMultiphysicsElement*>(celMixed);
        if (celMulti) {
          celMulti->Solution(ptInElement, 1, termHdiv);
          celMulti->Solution(ptInElement, 5, divFluxMixed);
        }
        for (int d = 0; d < termHdiv.size(); ++d) {
          termHdiv[d] = (-1./sqrtPerm) * termHdiv[d];
        }

        // Flux contribution
        REAL diffFlux = 0.0;
        for (int d = 0; d < termH1.size(); ++d) {
          REAL diff = termH1[d] - termHdiv[d];
          diffFlux += diff * diff;
        }

        // Balance contribution
        REAL diffBalance = (divFluxMixed[0] - force[0]) * (divFluxMixed[0] - force[0]);

        fluxError += diffFlux * weight;
        balanceError += diffBalance * weight;
      }

      if (SimData->m_PostProc.verbosityLevel) {
        std::ofstream errorlog("error_estimation.txt", std::ios_base::app);
        if (errorlog.is_open()) {
          errorlog << "Element " << icel << " - MatId " << matid
            << ": Flux error = " << fluxError
            << ", Balance error = " << balanceError << std::endl;
          errorlog.close();
        }
      }

      // Balance error is zero on wellbore (f belong to divergence space)
      if (matid == SimData->ECurveWell) {
        balanceError = 0.0;
      }

      REAL contribution = sqrt(fluxError) + (hk/(M_PI*sqrt(perm)))*sqrt(balanceError);
      int64_t igeo = cmeshMixed->Element(icel)->Reference()->Index();
      elementErrors[igeo] = (contribution * contribution);

      localTotalError += elementErrors[igeo];
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
  
  totalError = sqrt(totalError);

  // VTK output
  if (SimData->m_PostProc.verbosityLevel) {
    std::ofstream out_estimator("EstimatedError.vtk");
    TPZVTKGeoMesh::PrintCMeshVTK(cmeshMixed, out_estimator, elementErrors, "EstimatedError");
  }

  return totalError;
}

REAL TPZWannAdaptivityTools::PragerSynge(TPZMultiphysicsCompMesh* cmeshMixed, TPZCompMesh* cmeshH1, ProblemData* SimData, TPZVec<REAL>& elementErrors, int nthreads) {
  // If nthreads = 0, we use 1 thread
  nthreads++;

  // Ensure references point to H1 cmesh
  cmeshH1->Reference()->ResetReference();
  cmeshH1->LoadReferences();

  REAL dim = cmeshMixed->Dimension();
  int64_t ncel = cmeshMixed->NElements();
  int64_t ngel = cmeshMixed->Reference()->NElements();
  TPZAdmChunkVector<TPZCompEl *> &elementvec_m = cmeshMixed->ElementVec();
  elementErrors.Resize(ngel); // Ensure proper size
  elementErrors.Fill(0.0); 

  // Auxiliary vector for plotting
  TPZVec<REAL> elementErrorsAux(ncel, 0.0);

  // Parallelization setup
  std::vector<std::thread> threads(nthreads);
  TPZManVector<REAL> threadErrors(nthreads, 0.0);

  auto worker = [&](int tid, int64_t start, int64_t end) {
    REAL totalError = 0.0;
    for (int64_t icel = start; icel < end; ++icel) {
      TPZCompEl *celMixed = elementvec_m[icel];

      // Check if mixed element is condensed
      TPZCondensedCompEl *condEl = dynamic_cast<TPZCondensedCompEl *>(celMixed);
      if (condEl) {
        // If compel is condensed, load solution on the unconsensed compel
        condEl->LoadSolution();
        celMixed = condEl->ReferenceCompEl();
      }

      int matid = celMixed->Material()->Id();
      if (matid != SimData->EDomain) continue;

      TPZGeoEl *gel = celMixed->Reference();
      TPZCompEl *celH1 = gel->Reference();

      // Check if H1 element is condensed
      condEl = dynamic_cast<TPZCondensedCompEl *>(celH1);
      if (condEl) {
        // If compel is condensed, load solution on the unconsensed compel
        condEl->LoadSolution();
        celH1 = condEl->ReferenceCompEl();
      }

      if (!celMixed || !celH1) continue;
      if (gel->HasSubElement()) continue;
      if (celH1->Material()->Id() != matid) DebugStop();

      REAL hk = ElementDiameter(gel);
      REAL perm = SimData->m_Reservoir.perm[0];
      REAL sqrtPerm = sqrt(perm);

      REAL fluxError = 0.0;
      REAL balanceError = 0.0;

      // Set integration rule
      const TPZIntPoints* intrule = nullptr;
      const TPZIntPoints &intruleMixed = celMixed->GetIntegrationRule();
      const TPZIntPoints &intruleH1 = celH1->GetIntegrationRule();
      if (intruleMixed.NPoints() < intruleH1.NPoints()) {
        intrule = &intruleH1;
      } else {
        intrule = &intruleMixed;
      }

      for (int ip = 0; ip < intrule->NPoints(); ++ip) {
        TPZManVector<REAL,3> ptInElement(gel->Dimension());
        REAL weight, detjac;
        intrule->Point(ip, ptInElement, weight);
        TPZFNMatrix<9, REAL> jacobian, axes, jacinv;
        gel->Jacobian(ptInElement, jacobian, axes, detjac, jacinv);
        weight *= fabs(detjac);

        TPZManVector<REAL,3> force(1, 0.0);
        // TODO: Change for analytical force if available

        // Compute H1 term K^(1/2) grad(u)
        TPZManVector<REAL,3> termH1(dim, 0.0);
        celH1->Solution(ptInElement, 2, termH1);
        for (int d = 0; d < dim; ++d) {
          termH1[d] = sqrtPerm * termH1[d];
        }

        // Compute Hdiv term -K^(-1/2) sig and div(sig)
        TPZManVector<REAL,3> termHdiv(dim, 0.0);
        TPZManVector<REAL,1> divFluxMixed(1, 0.0);
        celMixed->Solution(ptInElement, 1, termHdiv);
        celMixed->Solution(ptInElement, 5, divFluxMixed);
        for (int d = 0; d < dim; ++d) {
          termHdiv[d] = (-1./sqrtPerm) * termHdiv[d];
        }

        // Flux contribution
        REAL diffFlux = 0.0;
        for (int d = 0; d < dim; ++d) {
          REAL diff = termH1[d] - termHdiv[d];
          diffFlux += diff * diff;
        }

        // Balance contribution
        REAL diffBalance = (divFluxMixed[0] - force[0]) * (divFluxMixed[0] - force[0]);

        fluxError += diffFlux * weight;
        balanceError += diffBalance * weight;
      }

      int64_t igeo = cmeshMixed->Element(icel)->Reference()->Index();
      elementErrors[igeo] = sqrt(fluxError) + (hk/(M_PI*sqrt(perm)))*sqrt(balanceError);
      elementErrorsAux[icel] = elementErrors[igeo];

      totalError += elementErrors[igeo] * elementErrors[igeo];
    }
    threadErrors[tid] = totalError;
  };

  int64_t chunk = ncel/nthreads;
  for (int t = 0; t < nthreads; ++t) {
    int64_t start = t * chunk;
    int64_t end = (t == nthreads - 1) ? ncel : (t + 1) * chunk;
    threads[t] = std::thread(worker, t, start, end);
  }
  for (auto& th : threads) th.join();

  REAL finalErrorSquared = 0.0;
  for (auto val : threadErrors) finalErrorSquared += val;

  // VTK output
  if (SimData->m_PostProc.verbosityLevel) {
    std::ofstream out_estimator("PragerSyngeEst.vtk");
    TPZVTKGeoMesh::PrintCMeshVTK(cmeshMixed, out_estimator, elementErrorsAux, "EstimatedError");
  }

  return std::sqrt(finalErrorSquared);
}

REAL TPZWannAdaptivityTools::GoalOriented(TPZMultiphysicsCompMesh* cmeshHdiv, TPZMultiphysicsCompMesh* cmeshDual, ProblemData* SimData, TPZVec<REAL>& elementErrors, int nthreads) {    
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

      // Only consider elements in the reservoir and wellbore domains
      int matid = celMixed->Material()->Id();
      if (matid != SimData->EDomain && matid != SimData->ECurveWell) continue;

      TPZGeoEl *gel = celMixed->Reference();
      if (gel->HasSubElement()) DebugStop();
      
      TPZCompEl *celDual = gel->Reference();

      // Check if element is condensed
      TPZCondensedCompEl *condElDual = dynamic_cast<TPZCondensedCompEl *>(celDual);
      if (condElDual) {
        // If compel is condensed, load solution on the unconsensed compel
        condElDual->LoadSolution();
        celDual = condElDual->ReferenceCompEl();
      }

      if (!celMixed || !celDual) DebugStop();
      if (celDual->Material()->Id() != matid) DebugStop();

      REAL goalError = 0.0;

      if (matid == SimData->EDomain) {
        goalError = GoalContribution3D(celMixed, celDual, SimData);
      } else if (matid == SimData->ECurveWell) {
        goalError = 0.0; // GoalContribution1D(celMixed, celDual, SimData);
      } else {
        DebugStop();
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

REAL TPZWannAdaptivityTools::GoalContribution1D(TPZCompEl* celMixed, TPZCompEl* celDual, ProblemData* SimData) {
  // Some useful constants
  REAL diameter = SimData->m_Wellbore.radius * 2.0;
  REAL rho = SimData->m_Fluid.density;
  REAL mu = SimData->m_Fluid.viscosity;
  REAL c = (2.252610888 * pow(diameter, 19. / 7.)) / (pow(mu, 1. / 7.) * pow(rho, 3. / 7.));
  REAL CNL = pow(c, -7. / 4.);
  REAL CLin = 128. * mu / (M_PI * pow(diameter, 4));
  REAL goalError = 0.0;

  // Set integration rule
  TPZGeoEl *gel = celMixed->Reference();
  const TPZIntPoints *intrule = nullptr;
  intrule = gel->CreateSideIntegrationRule(gel->NSides() - 1, 5);

  for (int ip = 0; ip < intrule->NPoints(); ++ip) {
    TPZManVector<REAL, 3> ptInElement(gel->Dimension());
    REAL weight, detjac;
    intrule->Point(ip, ptInElement, weight);
    TPZFNMatrix<9, REAL> jacobian, axes, jacinv;
    gel->Jacobian(ptInElement, jacobian, axes, detjac, jacinv);
    weight *= detjac; // fabs(detjac);

    TPZManVector<REAL, 3> x(3, 0.0);
    gel->X(ptInElement, x); // Real coordinates for x

    // No source function on reservoir
    TPZManVector<REAL, 3> force(1, 0.0);

    // Compute mixed solution ph and qh
    TPZManVector<REAL, 3> ph(gel->Dimension(), 0.0);
    TPZManVector<REAL, 3> qh(gel->Dimension(), 0.0);
    TPZManVector<REAL, 1> divqh(1, 0.0);
    celMixed->Solution(ptInElement, 2, ph);
    celMixed->Solution(ptInElement, 1, qh);
    celMixed->Solution(ptInElement, 3, divqh);

    // Compute dual solution zh and psih
    TPZManVector<REAL, 3> zh(gel->Dimension(), 0.0);
    TPZManVector<REAL, 3> psih(gel->Dimension(), 0.0);
    TPZManVector<REAL, 1> divpsih(1, 0.0);
    celDual->Solution(ptInElement, 2, zh);
    celDual->Solution(ptInElement, 1, psih);
    celDual->Solution(ptInElement, 3, divpsih);

    // From qh compute G^(-1) qh
    REAL velocity = 4.0 * qh[0] / (M_PI * diameter * diameter);
    REAL reynolds = (rho * std::abs(velocity) * diameter) / mu;
    bool turbulent = reynolds > 1187.38;

    turbulent = false;

    REAL gh = 0.0;
    if (turbulent) {
      gh = -CNL * qh[0] * (pow(std::abs(qh[0]), 3. / 4.));
    } else {
      gh = -CLin * qh[0];
    }

    // Flux contribution
    REAL FluxTerm = gh * psih[0];

    // div contributions
    REAL divTermA = -divqh[0] * zh[0];
    REAL divTermB = -divpsih[0] * ph[0];

    // Minus sign due to the definition of the dual problem
    goalError += (-FluxTerm - divTermA - divTermB) * weight;
  }
  return goalError;
}

REAL TPZWannAdaptivityTools::GoalContribution3D(TPZCompEl* celMixed, TPZCompEl* celDual, ProblemData* SimData) {
  TPZManVector<REAL, 3> perm = SimData->m_Reservoir.perm;
  TPZGeoEl *gel = celMixed->Reference();
  int dim = gel->Dimension();
  if (dim != 3) DebugStop();

  REAL goalError = 0.0;

  // Set integration rule
  const TPZIntPoints *intrule = nullptr;
  intrule = gel->CreateSideIntegrationRule(gel->NSides() - 1, 5);

  for (int ip = 0; ip < intrule->NPoints(); ++ip) {
    TPZManVector<REAL, 3> ptInElement(gel->Dimension());
    REAL weight, detjac;
    intrule->Point(ip, ptInElement, weight);
    TPZFNMatrix<9, REAL> jacobian, axes, jacinv;
    gel->Jacobian(ptInElement, jacobian, axes, detjac, jacinv);
    weight *= detjac;

    TPZManVector<REAL, 3> x(3, 0.0);
    gel->X(ptInElement, x); // Real coordinates for x

    // No source function on reservoir
    TPZManVector<REAL, 3> force(1, 0.0);

    // Compute mixed solution ph and sigh
    TPZManVector<REAL, 3> ph(gel->Dimension(), 0.0);
    TPZManVector<REAL, 3> sigh(gel->Dimension(), 0.0);
    TPZManVector<REAL, 1> divsigh(1, 0.0);
    celMixed->Solution(ptInElement, 2, ph);
    celMixed->Solution(ptInElement, 1, sigh);
    celMixed->Solution(ptInElement, 5, divsigh);

    // Compute dual solution zh and psih
    TPZManVector<REAL, 3> zh(gel->Dimension(), 0.0);
    TPZManVector<REAL, 3> psih(gel->Dimension(), 0.0);
    TPZManVector<REAL, 1> divpsih(1, 0.0);
    celDual->Solution(ptInElement, 2, zh);
    celDual->Solution(ptInElement, 1, psih);
    celDual->Solution(ptInElement, 5, divpsih);

    // F contribution
    REAL FTerm = -force[0] * zh[0];

    // Flux contribution
    REAL FluxTerm = 0.0;
    for (int d = 0; d < dim; ++d) {
      FluxTerm += (SimData->m_Fluid.density / perm[d]) * sigh[d] * psih[d];
    }

    // div contributions
    REAL divTermA = -divsigh[0] * zh[0];
    REAL divTermB = -divpsih[0] * ph[0];

    REAL ATerm = FluxTerm + divTermA + divTermB;

    goalError += (FTerm - ATerm) * weight;
  }

  // --- Boundary contributions ---
  std::set<int> bcIds;
  std::map<int, REAL> bcValues;

  // Gather the Dirichlet boundary conditions
  for (auto &bc : SimData->m_Reservoir.BCs) {
    if (bc.second.type == 0) { // Dirichlet
      bcIds.insert(bc.second.matid);
      bcValues[bc.second.matid] = bc.second.value;
    }
  }

  int firstSide = gel->FirstSide(dim - 1);
  int lastSide = gel->FirstSide(dim);

  for (int side = firstSide; side < lastSide; ++side) {
    TPZGeoElSide gelSide(gel, side);
    TPZGeoElSide neigh = gelSide.HasNeighbour(bcIds);
    if (neigh) {
      TPZGeoEl *gelB = neigh.Element();
      if (gelB->HasSubElement()) DebugStop();
      TPZCompEl *celDualB = gelB->Reference();
      if (!celDualB) DebugStop();

      int neighMatid = gelB->MaterialId();
      REAL pe_cte = bcValues[neighMatid];

      // Check if element is condensed
      TPZCondensedCompEl *condEl = dynamic_cast<TPZCondensedCompEl *>(celDualB);
      if (condEl) {
        // If compel is condensed, load solution on the unconsensed compel
        condEl->LoadSolution();
        celDualB = condEl->ReferenceCompEl();
      }

      TPZIntPoints *intruleSide = gelB->CreateSideIntegrationRule(gelB->NSides() - 1, 5);
      for (int ip = 0; ip < intruleSide->NPoints(); ++ip) {
        TPZManVector<REAL, 3> ptInElement(gelB->Dimension());
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

        goalError += -pe_cte * psih[0] * weight;
      }
    }
  }

  // Interface contributions ---

  std::set<int> interfaceIds;
  interfaceIds.insert(SimData->EPressureInterface);

  for (int side = firstSide; side < lastSide; ++side) {
    TPZGeoElSide gelSide(gel, side);
    TPZGeoElSide neigh = gelSide.HasNeighbour(interfaceIds);
    if (neigh) {
      TPZGeoEl *faceGel = neigh.Element();
      REAL dimF = faceGel->Dimension();
      if (dimF != 2) DebugStop();
      intrule = gel->CreateSideIntegrationRule(side, 5);
      for (int ip = 0; ip < intrule->NPoints(); ip++) {
        TPZManVector<REAL, 3> ptOnSide(faceGel->Dimension());
        TPZFNMatrix<9, REAL> jacobian, axes, jacinv;
        REAL weight, detjac;
        faceGel->Jacobian(ptOnSide, jacobian, axes, detjac, jacinv);

        // Compute normal
        // Must be done in the integration point to account for curved sides
        TPZManVector<REAL, 3> v1(3), v2(3), normal(3);
        v1[0] = axes(0, 0);
        v1[1] = axes(0, 1);
        v1[2] = axes(0, 2);
        v2[0] = axes(1, 0);
        v2[1] = axes(1, 1);
        v2[2] = axes(1, 2);

        normal[0] = v1[1] * v2[2] - v1[2] * v2[1];
        normal[1] = v1[2] * v2[0] - v1[0] * v2[2];
        normal[2] = v1[0] * v2[1] - v1[1] * v2[0];

        REAL norm = sqrt(normal[0] * normal[0] + normal[1] * normal[1] +
                         normal[2] * normal[2]);
        if (norm > 1e-12) {
          normal[0] /= norm;
          normal[1] /= norm;
          normal[2] /= norm;
        }

        intrule->Point(ip, ptOnSide, weight);
        TPZManVector<REAL, 3> ptInElement(gel->Dimension());
        TPZTransform<> trans = gel->SideToSideTransform(side, gel->NSides() - 1);
        trans.Apply(ptOnSide, ptInElement);
        weight *= detjac;

        TPZManVector<REAL, 3> sigh(gel->Dimension(), 0.0);
        TPZManVector<REAL, 3> ph(gel->Dimension(), 0.0);
        TPZManVector<REAL, 3> psih(gel->Dimension(), 0.0);
        TPZManVector<REAL, 1> zh(gel->Dimension(), 0.0);

        celMixed->Solution(ptInElement, 2, ph);
        celMixed->Solution(ptInElement, 1, sigh);
        celDual->Solution(ptInElement, 2, zh);
        celDual->Solution(ptInElement, 1, psih);

        REAL fluxMixed = 0.0;
        REAL fluxDual = 0.0;
        for (int d = 0; d < dim; ++d) {
          fluxMixed += sigh[d] * normal[d];
          fluxDual +=  psih[d] * normal[d];
        }

        goalError += (ph[0]*fluxDual) * weight;
      }   
    }
  }
  return goalError;
}

void TPZWannAdaptivityTools::MarkElementsForRefinement(const TPZVec<REAL>& elementErrors, TPZVec<int>& refinementIndicator, REAL theta) {
  int64_t nels = elementErrors.size();
  refinementIndicator.Resize(nels);
  refinementIndicator.Fill(0);
  REAL totalError = 0.0;
  for (int64_t iel = 0; iel < nels; ++iel) {
      totalError += elementErrors[iel];
  }
  totalError = totalError * theta;

  // Sort the element errors
  std::vector<std::pair<REAL, int64_t>> errorIndexPairs;
  for (int64_t iel = 0; iel < nels; ++iel) {
      errorIndexPairs.emplace_back(elementErrors[iel], iel);
  }
  std::sort(errorIndexPairs.begin(), errorIndexPairs.end());

  REAL accumulatedError = 0.0;
  while (accumulatedError < totalError && !errorIndexPairs.empty()) {
      auto [error, index] = errorIndexPairs.back();
      errorIndexPairs.pop_back();
      accumulatedError += error;
      refinementIndicator[index] = 1;
  }
}

void TPZWannAdaptivityTools::MeshSmoothing(TPZGeoMesh* gmesh, TPZVec<int>& refinementIndicator) {
  // If an element has most of its neighbors refined, then refine it too
  for (int64_t iel = 0; iel < refinementIndicator.size(); ++iel) {
    if (refinementIndicator[iel] != 1) continue;
    TPZGeoEl* gel = gmesh->Element(iel);
    if (!gel) DebugStop();
    int firstside = gel->FirstSide(gel->Dimension()-1);
    int lastside = gel->FirstSide(gel->Dimension());
    for (int side = firstside; side < lastside; ++side) {
      TPZGeoElSide gelSide(gel, side);
      TPZGeoElSide neigh = gelSide.Neighbour();
      TPZGeoEl* neighGel = neigh.Element();
      if (!neighGel) DebugStop();
      int threshold = neighGel->Dimension() == 1 ? 2 : 5;
      int neighIndex = neighGel->Index();
      if (refinementIndicator[neighIndex] == 1) continue;
      int firstsideNeigh = neighGel->FirstSide(neighGel->Dimension()-1);
      int lastsideNeigh = neighGel->FirstSide(neighGel->Dimension());
      int countRefNeigh = 0;
      for (int sideNeigh = firstsideNeigh; sideNeigh < lastsideNeigh; ++sideNeigh) {
        TPZGeoElSide neighSide(neighGel, sideNeigh);
        TPZGeoElSide neigh2 = neighSide.Neighbour();
        if (!neigh2) DebugStop();
        int64_t neigh2Index = neigh2.Element()->Index();
        if (refinementIndicator[neigh2Index] == 1) countRefNeigh++;
      }
      if (countRefNeigh >= threshold) {
        refinementIndicator[neighGel->Index()] = 1;
      }
    }  
  }

  // If an element is refined twice, ensure its neighbors are refined at least once
  for (int64_t iel = 0; iel < refinementIndicator.size(); ++iel) {
    if (refinementIndicator[iel] != 1) continue;
    TPZGeoEl* gel = gmesh->Element(iel);

    // Get corner indexes of gel
    TPZManVector<int64_t> cornerIndexes(gel->NCornerNodes());
    for (int i = 0; i < gel->NCornerNodes(); ++i) {
      cornerIndexes[i] = gel->NodeIndex(i);
    }
    TPZGeoEl* fatherGel = gel->Father();
    if (!fatherGel) continue;
    int firstside = fatherGel->FirstSide(fatherGel->Dimension()-1);
    int lastside = fatherGel->FirstSide(fatherGel->Dimension());
    for (int side = firstside; side < lastside; ++side) {
      TPZGeoElSide gelSide(fatherGel, side);
      TPZGeoElSide neigh = gelSide.Neighbour();
      TPZGeoEl* neighGel = neigh.Element();
      if (!neighGel) DebugStop();
      if (neighGel->Dimension() != gel->Dimension()) continue;
      if (neighGel->HasSubElement()) continue;

      // Get corner idexes of neighGel
      TPZManVector<int64_t> neighCornerIndexes(neighGel->NCornerNodes());
      for (int i = 0; i < neighGel->NCornerNodes(); ++i) {
        neighCornerIndexes[i] = neighGel->NodeIndex(i);
      }

      // Verify if there are common corners
      int commonCorners = 0;
      for (int i = 0; i < gel->NCornerNodes(); ++i) {
        for (int j = 0; j < neighGel->NCornerNodes(); ++j) {
          if (cornerIndexes[i] == neighCornerIndexes[j]) {
            commonCorners++;
            break;
          }
        }
      }
      
      if (commonCorners > 0) refinementIndicator[neighGel->Index()] = 1;
    }
  }
}

void TPZWannAdaptivityTools::MeshWellCompatibility(TPZGeoMesh* gmesh, TPZVec<int>& refinementIndicator, ProblemData* SimData) {
  REAL dim = gmesh->Dimension();
  REAL tol = 1e-6;
  std::set<REAL> xcoords2D;
  TPZVec<int64_t> needRefinement2D;

  std::set<int> matIds = {
      SimData->ESurfWellCyl, SimData->EPressure2DSkin, 
      SimData->EPressureInterface, SimData->EHDivBoundInterface
    };

  // From the initial list of "to-refine" elements, get the the centroids
  // of the faces that intersect the wellbore surface (fill xcoords2D)
  for (int64_t i = 0; i < refinementIndicator.size(); ++i) {
    if (refinementIndicator[i] != 1) continue;
    TPZGeoEl* gel = gmesh->Element(i);
    if (!gel) DebugStop();
    if (gel->Dimension() == 1) {
      TPZManVector<REAL,3> center(gel->Dimension(), 0.0);
      gel->CenterPoint(gel->NSides()-1, center); 
      TPZManVector<REAL,3> rcenter(3, 0.0);
      gel->X(center, rcenter);
      TPZWannGeometryTools::InsertXCoorInSet(rcenter[0], xcoords2D, tol);
    } else if (gel->Dimension() == 3) {
      int firstside = gel->FirstSide(dim-1);
      int lastside = gel->FirstSide(dim);
      for (int side = firstside; side < lastside; ++side) {
        TPZGeoElSide gelSide(gel, side);
        TPZGeoElSide neigh = gelSide.HasNeighbour(matIds);
        if (neigh) {
          TPZManVector<REAL,3> center(neigh.Element()->Dimension(), 0.0);
          neigh.Element()->CenterPoint(neigh.Side(), center);
          TPZManVector<REAL,3> rcenter(3, 0.0);
          neigh.Element()->X(center, rcenter);
          TPZWannGeometryTools::InsertXCoorInSet(rcenter[0], xcoords2D, tol);
        }
      }
    }
  }

  // Manually set BC elements to be refined
  for (int64_t i = 0; i < refinementIndicator.size(); ++i) {
    if (refinementIndicator[i] == 0) continue;
    TPZGeoEl* gel = gmesh->Element(i);
    if (!gel) DebugStop();
    if (gel->Dimension() != 3) continue; 
    int firstside = gel->FirstSide(dim-1);
    int lastside = gel->FirstSide(dim);
    for (int side = firstside; side < lastside; ++side) {
      TPZGeoElSide gelSide(gel, side);
      std::set<int> bcIds = {SimData->EFarField, SimData->ESurfHeel, SimData->ESurfToe};
      TPZGeoElSide neigh = gelSide.HasNeighbour(bcIds);
      if (neigh) {
        refinementIndicator[neigh.Element()->Index()] = 1;
      }
    }
  }

  // From xcoords2D, get the 2D elements that need refinement (fill needRefinement2D)
  for (int64_t iel = 0; iel < gmesh->NElements(); ++iel) {
    TPZGeoEl* gel = gmesh->Element(iel);
    if (matIds.find(gel->MaterialId()) == matIds.end()) continue;
    if (gel->HasSubElement()) continue;
    TPZManVector<REAL,3> center(gel->Dimension(), 0.0);
    gel->CenterPoint(gel->NSides()-1, center);
    TPZManVector<REAL,3> rcenter(3, 0.0);
    gel->X(center, rcenter);
    if (TPZWannGeometryTools::CheckXInSet(rcenter[0], xcoords2D, tol)) {
      needRefinement2D.push_back(gel->Index());
      refinementIndicator[gel->Index()] = 1;
    }
  }

  // From 2D elements, guarantee that 1D and 3D neighbors will be refined
  for (int64_t i = 0; i < needRefinement2D.size(); ++i) {
    TPZGeoEl* gel = gmesh->Element(needRefinement2D[i]);
    if (gel->Dimension() != 2) DebugStop();
    // Pick the 3D neighbor
    TPZGeoElSide gelSide(gel);
    TPZGeoElSide neigh3D = gelSide.HasNeighbour(SimData->EDomain);
    if (neigh3D) {
      refinementIndicator[neigh3D.Element()->Index()] = 1;
    }
    // Pick the 1D neighbor
    int firstside = gel->FirstSide(1);
    int lastside = gel->FirstSide(2);
    for (int side = firstside; side < lastside; ++side) {
      TPZGeoElSide gelSide(gel, side);
      TPZGeoElSide neigh = gelSide.HasNeighbour(SimData->ECurveWell);
      if (neigh) {
        refinementIndicator[neigh.Element()->Index()] = 1;
      }
    }
  }
}

void TPZWannAdaptivityTools::AdaptivityProcess(TPZGeoMesh* gmesh, ProblemData* SimData, TPZVec<REAL>& elementErrors, REAL tol) {
  if (elementErrors.size() != gmesh->NElements()) {
    std::cerr << "Error: elementErrors size does not match number of geometric elements." << std::endl;
    DebugStop();
  }

  TPZVec<int> refinementIndicator(gmesh->NElements(), 0);
    
  // Mark elements for refinement
  TPZWannAdaptivityTools::MarkElementsForRefinement(elementErrors, refinementIndicator, tol);
  REAL sumRef = std::accumulate(refinementIndicator.begin(), refinementIndicator.end(), 0);
  std::cout << "Initial number of to-refine elements: " << sumRef << std::endl;

  // Compatibilize refinement to well geometry
  TPZWannAdaptivityTools::MeshWellCompatibility(gmesh, refinementIndicator, SimData);
  sumRef = std::accumulate(refinementIndicator.begin(), refinementIndicator.end(), 0);
  std::cout << "Number of to-refine elements after well compatibility: " << sumRef << std::endl;

  // Smooth the refinement map
  TPZWannAdaptivityTools::MeshSmoothing(gmesh, refinementIndicator);
  sumRef = std::accumulate(refinementIndicator.begin(), refinementIndicator.end(), 0);
  std::cout << "Number of to-refine elements after mesh smoothing: " << sumRef << std::endl;

  // Compatibilize refinement to well geometry
  TPZWannAdaptivityTools::MeshWellCompatibility(gmesh, refinementIndicator, SimData);
  sumRef = std::accumulate(refinementIndicator.begin(), refinementIndicator.end(), 0);
  std::cout << "Final number of to-refine elements: " << sumRef << std::endl;

  // From refinementIndicator, get the list of elements to refine
  TPZVec<int64_t> toRefine;
  for (int64_t iel = 0; iel < refinementIndicator.size(); ++iel) {
    if (refinementIndicator[iel] == 1) {
      toRefine.push_back(iel);
    }
  }

  // Finally perform the refinement.
  TPZWannGeometryTools::hRefinement(gmesh, toRefine);
}


REAL TPZWannAdaptivityTools::ForcingFunctionWellbore(TPZCompEl* celMixed, const TPZManVector<REAL,3>& pt) {
  // This function should compute the forcing function value at the given point
  // For now, we return a constant value as a placeholder
  return 0.0;
}

REAL TPZWannAdaptivityTools::ElementDiameter(TPZGeoEl* gel) {
  REAL maxdist = 0.;
  int nnodes = gel->NNodes();
  for (int i = 0; i < nnodes; ++i) {
    TPZManVector<REAL,3> xi(3,0.), xj(3,0.);
    gel->Node(i).GetCoordinates(xi);
    for (int j = i+1; j < nnodes; ++j) {
      gel->Node(j).GetCoordinates(xj);
      REAL dist = 0.;
      for (int d = 0; d < gel->Dimension(); ++d) {
        dist += (xi[d] - xj[d]) * (xi[d] - xj[d]);
      }
      dist = sqrt(dist);
      if (dist > maxdist) maxdist = dist;
    }
  }
  return maxdist;
}