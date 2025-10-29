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
#include "pzlog.h"
#include "pzintel.h"
#include "TPZVTKGeoMesh.h"

// Wann imports
#include "TPZWannEstimationTools.h"
#include "TPZWannGeometryTools.h"

void TPZWannEstimationTools::EstimateAndRefine(TPZMultiphysicsCompMesh* cmeshHdiv, TPZCompMesh* cmeshH1, ProblemData* SimData, int nthreads) {
  TPZVec<int> RefinementIndicator = ErrorEstimation(cmeshHdiv, cmeshH1, SimData, nthreads);
  REAL sumRef = std::accumulate(RefinementIndicator.begin(), RefinementIndicator.end(), 0);
  std::cout << "Initial number of to-refine elements: " << sumRef << std::endl;

  meshSmoothing(cmeshHdiv->Reference(), RefinementIndicator);
  sumRef = std::accumulate(RefinementIndicator.begin(), RefinementIndicator.end(), 0);
  std::cout << "Number of to-refine elements after mesh smoothing: " << sumRef << std::endl;

  hRefinement(cmeshHdiv->Reference(), RefinementIndicator, SimData);
  sumRef = std::accumulate(RefinementIndicator.begin(), RefinementIndicator.end(), 0);
  std::cout << "Final number of to-refine elements: " << sumRef << std::endl;

  // After refining we have to reorder the well IDs again
  TPZWannGeometryTools::OrderIds(cmeshHdiv->Reference(), SimData);
}

TPZVec<int> TPZWannEstimationTools::ErrorEstimation(TPZMultiphysicsCompMesh* cmeshMixed, TPZCompMesh* cmesh, ProblemData* SimData, int nthreads) {

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
  TPZVec<REAL> elementErrors(ncel, 0.0);
  TPZVec<int> RefinementIndicator(ngel, 0);

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
      REAL perm = SimData->m_Reservoir.perm;
      if (matid == SimData->ECurveWell) {
        perm = SimData->m_Wellbore.perm;
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

      if (SimData->m_VerbosityLevel) {
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
      elementErrors[icel] = (contribution * contribution);

      localTotalError += elementErrors[icel];
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
  std::ofstream out_estimator("EstimatedError.vtk");
  TPZVTKGeoMesh::PrintCMeshVTK(cmeshMixed, out_estimator, elementErrors, "EstimatedError");

  std::cout << "\nTotal estimated error: " << totalError << std::endl;

  // Indices of geoel elements that need refinement 
  for (int64_t i = 0; i < ncel; ++i) {
    if (elementErrors[i] > estimator_tol) {
      int64_t igeo = cmeshMixed->Element(i)->Reference()->Index();
      RefinementIndicator[igeo] = 1;
    }
  }

  return RefinementIndicator;
}

void TPZWannEstimationTools::hRefinement(TPZGeoMesh* gmesh, TPZVec<int>& RefinementIndicator, ProblemData* SimData) {
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
  for (int64_t i = 0; i < RefinementIndicator.size(); ++i) {
    if (RefinementIndicator[i] == 0) continue;
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
  for (int64_t i = 0; i < RefinementIndicator.size(); ++i) {
    if (RefinementIndicator[i] == 0) continue;
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
        RefinementIndicator[neigh.Element()->Index()] = 1;
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
      RefinementIndicator[gel->Index()] = 1;
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
      RefinementIndicator[neigh3D.Element()->Index()] = 1;
    }
    // Pick the 1D neighbor
    int firstside = gel->FirstSide(1);
    int lastside = gel->FirstSide(2);
    for (int side = firstside; side < lastside; ++side) {
      TPZGeoElSide gelSide(gel, side);
      TPZGeoElSide neigh = gelSide.HasNeighbour(SimData->ECurveWell);
      if (neigh) {
        RefinementIndicator[neigh.Element()->Index()] = 1;
      }
    }
  }

  // Perform h-refinement on the specified elements
  for (int64_t i = 0; i < RefinementIndicator.size(); ++i) {
    if (RefinementIndicator[i] == 0) continue;
    TPZVec<TPZGeoEl *> pv;
    TPZGeoEl* gel = gmesh->Element(i);
    if (!gel) DebugStop();
    if (gel->HasSubElement()) continue; // We have to do this because of duplicated entries
    gel->Divide(pv);
  }
}

void TPZWannEstimationTools::meshSmoothing(TPZGeoMesh* gmesh, TPZVec<int>& RefinementIndicator) {
  // If an element has most of its neighbors refined, then refine it too
  for (int64_t iel = 0; iel < RefinementIndicator.size(); ++iel) {
    if (RefinementIndicator[iel] != 1) continue;
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
      if (RefinementIndicator[neighIndex] == 1) continue;
      int firstsideNeigh = neighGel->FirstSide(neighGel->Dimension()-1);
      int lastsideNeigh = neighGel->FirstSide(neighGel->Dimension());
      int countRefNeigh = 0;
      for (int sideNeigh = firstsideNeigh; sideNeigh < lastsideNeigh; ++sideNeigh) {
        TPZGeoElSide neighSide(neighGel, sideNeigh);
        TPZGeoElSide neigh2 = neighSide.Neighbour();
        if (!neigh2) DebugStop();
        int64_t neigh2Index = neigh2.Element()->Index();
        if (RefinementIndicator[neigh2Index] == 1) countRefNeigh++;
      }
      if (countRefNeigh >= threshold) {
        RefinementIndicator[neighGel->Index()] = 1;
      }
    }  
  }

  // If an element has a neighbor with refinement two levels higher, then refine it too
  for (int64_t iel = 0; iel < RefinementIndicator.size(); ++iel) {
    if (RefinementIndicator[iel] != 1) continue;
    TPZGeoEl* gel = gmesh->Element(iel);
    if (!gel) DebugStop();
    int matid = gel->MaterialId();
    int firstside = gel->FirstSide(gel->Dimension()-1);
    int lastside = gel->FirstSide(gel->Dimension());
    for (int side = firstside; side < lastside; ++side) {
      TPZGeoElSide gelSide(gel, side);
      TPZGeoElSide lowerSide = gelSide.HasLowerLevelNeighbour(matid);
      if (lowerSide) {
        TPZGeoEl* neighGel = lowerSide.Element();
        if (!neighGel->HasSubElement()) RefinementIndicator[neighGel->Index()] = 1;
      }
    }
  }
  
  // Old version
  // for (int64_t iel = 0; iel < RefinementIndicator.size(); ++iel) {
  //   if (RefinementIndicator[iel] != 1) continue;
  //   TPZGeoEl* gel = gmesh->Element(iel);
  //   TPZGeoEl* fatherGel = gel->Father();
  //   if (!fatherGel) continue;
  //   int firstside = fatherGel->FirstSide(fatherGel->Dimension()-1);
  //   int lastside = fatherGel->FirstSide(fatherGel->Dimension());
  //   for (int side = firstside; side < lastside; ++side) {
  //     TPZGeoElSide gelSide(fatherGel, side);
  //     TPZGeoElSide neigh = gelSide.Neighbour();
  //     TPZGeoEl* neighGel = neigh.Element();
  //     if (!neighGel) DebugStop();
  //     if (neighGel->Dimension() != gel->Dimension()) continue;
  //     if (!neighGel->HasSubElement()) {
  //       RefinementIndicator[neighGel->Index()] = 1;
  //     }
  //   }
  // }
}

REAL TPZWannEstimationTools::ForcingFunctionWellbore(TPZCompEl* celMixed, const TPZManVector<REAL,3>& pt) {
  // This function should compute the forcing function value at the given point
  // For now, we return a constant value as a placeholder
  return 0.0;
}

REAL TPZWannEstimationTools::ElementDiameter(TPZGeoEl* gel) {
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

// ==============================
// Old or experimental stuff
// ==============================

TPZVec<int> TPZWannEstimationTools::ErrorEstimationOld(TPZMultiphysicsCompMesh* cmeshMixed, TPZCompMesh* cmesh, ProblemData* SimData, int nthreads) {
  {
    std::ofstream clearlog("error_estimation.txt", std::ios_base::trunc);
  }

  nthreads++;
  // Ensure references point to H1 cmesh
  cmesh->Reference()->ResetReference();
  cmesh->LoadReferences();

  int64_t ncel = cmeshMixed->NElements();
  int64_t ngel = cmeshMixed->Reference()->NElements();
  TPZAdmChunkVector<TPZCompEl *> &elementvec_m = cmeshMixed->ElementVec();
  TPZVec<REAL> elementErrors(ncel, 0.0);
  TPZVec<int> RefinementIndicator(ngel, 0);

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
      REAL perm = SimData->m_Reservoir.perm;
      if (matid == SimData->ECurveWell) {
        perm = SimData->m_Wellbore.perm;
      }
      REAL sqrtPerm = sqrt(perm);

      REAL fluxError = 0.0;
      REAL balanceError = 0.0;
      REAL conformityError = 0.0;

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

        // Compute L2 term K^(1/2) grad(u)
        TPZManVector<REAL,3> termL2(3,0.0);
        celMixed->Solution(ptInElement, 9, termL2);
        for (int d = 0; d < termL2.size(); ++d) {
          termL2[d] = sqrtPerm * termL2[d];
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
        for (int d = 0; d < termL2.size(); ++d) {
          REAL diff = termL2[d] - termHdiv[d];
          diffFlux += diff * diff;
        }

        // Balance contribution
        REAL diffBalance = (divFluxMixed[0] - force[0]) * (divFluxMixed[0] - force[0]);

        // Conformity contribution
        REAL diffConformity = 0.0;
        for (int d = 0; d < termH1.size(); ++d) {
          REAL diff = termH1[d] - termL2[d];
          diffConformity += diff * diff;
        }

        fluxError += diffFlux * weight;
        balanceError += diffBalance * weight;
        conformityError += diffConformity * weight;
      }

      if (SimData->m_VerbosityLevel) {
        std::ofstream errorlog("error_estimation.txt", std::ios_base::app);
        if (errorlog.is_open()) {
          errorlog << "Element " << icel << " - MatId " << matid
            << ": Flux error = " << fluxError
            << ", Balance error = " << balanceError 
            << ", Conformity error = " << conformityError << std::endl;
          errorlog.close();
        }
      }

      // Balance error is zero on wellbore (f belong to divergence space)
      if (matid == SimData->ECurveWell) {
        balanceError = 0.0;
      }

      REAL contribution = sqrt(fluxError) + (hk/(M_PI*sqrt(perm)))*sqrt(balanceError);
      elementErrors[icel] = (contribution * contribution) + conformityError;

      localTotalError += elementErrors[icel];
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

  if (SimData->m_VerbosityLevel > 0) {
    std::ofstream errorlog_total("error_estimation.txt", std::ios_base::app);
    if (errorlog_total.is_open()) {
      errorlog_total << "\nError estimation in energy norm: " << totalError << std::endl;
      errorlog_total.close();
    }
  }

  // VTK output
  std::ofstream out_estimator("EstimatedError.vtk");
  TPZVTKGeoMesh::PrintCMeshVTK(cmeshMixed, out_estimator, elementErrors, "EstimatedError");

  // Indices of geoel elements that need refinement 
  for (int64_t i = 0; i < ncel; ++i) {
    if (elementErrors[i] > estimator_tol) {
      int64_t igeo = cmeshMixed->Element(i)->Reference()->Index();
      RefinementIndicator[igeo] = 1;
    }
  }

  return RefinementIndicator;
}

void TPZWannEstimationTools::FakeRefine(TPZGeoMesh* gmesh, ProblemData* SimData) {
  REAL dim = gmesh->Dimension();
  REAL tol = 1e-6;
  int n_refined = 0;

  std::set<int> matIds = {
      SimData->EFarField
    };

  int64_t nel = gmesh->NElements();

  // Perform h-refinement on the specified elements
  for (int64_t i = 0; i < nel; ++i) {
    TPZVec<TPZGeoEl *> pv;
    TPZGeoEl* gel = gmesh->Element(i);
    if (gel->Dimension() != 3) continue;
    bool willRefine = false;
    int firstside = gel->FirstSide(2);
    int lastside = gel->FirstSide(3);
    for (int side = firstside; side < lastside; ++side){
      TPZGeoElSide gelSide(gel,side);
      TPZGeoElSide neigh = gelSide.HasNeighbour(matIds);
      if (neigh) {
        willRefine = true;
      }
    }
    if (!gel) DebugStop();
    if (gel->HasSubElement()) continue; // We have to do this because of duplicated entries
    if (willRefine) {
      n_refined++;
      if (n_refined >= 2) {
        gel->Divide(pv);
        break;
      }
    }
  }
  std::cout << n_refined << " elements were refined in the fake refinement." << std::endl;
}

// void TPZWannEstimationTools::CheckRef(TPZCompMesh* cmesh) {
//   std::map<int, int> matid_to_count;
//   for (int64_t i = 0; i < cmesh->NElements(); ++i) {
//     TPZCompEl* cel = cmesh->Element(i);
//     if (!cel) continue;
//     TPZGeoEl* gel = cel->Reference();
//     if (gel->Father()) {
//       matid_to_count[gel->MaterialId()]++;
//       std::cout << "Search for index " << i << std::endl;
//     }
//   }
//   std::cout << "Count of subelements by material id:" << std::endl;
//   for (const auto& pair : matid_to_count) {
//     std::cout << "MatId " << pair.first << ": " << pair.second << std::endl;
//   }
// }