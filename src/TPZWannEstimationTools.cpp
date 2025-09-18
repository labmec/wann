#include <iostream>
#include <random>
#include <thread>
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
  TPZVec<int64_t> needRefinement = ErrorEstimation(cmeshHdiv, cmeshH1, SimData, nthreads);
  std::cout << "Number of elements to refine: " << needRefinement.size() << std::endl;
  hRefinement(cmeshHdiv->Reference(), needRefinement, SimData);

  // After refining we have to reorder the well IDs again
  TPZWannGeometryTools::OrderIds(cmeshHdiv->Reference(), SimData);
}

TPZVec<int64_t> TPZWannEstimationTools::ErrorEstimation(TPZMultiphysicsCompMesh* cmeshMixed, TPZCompMesh* cmesh, ProblemData* SimData, int nthreads) {
  {
    std::ofstream clearlog("error_estimation.txt", std::ios_base::trunc);
  }

  // Ensure references point to H1 cmesh
  cmesh->Reference()->ResetReference();
  cmesh->LoadReferences();

  int64_t ncel = cmeshMixed->NElements();
  TPZAdmChunkVector<TPZCompEl *> &elementvec_m = cmeshMixed->ElementVec();
  TPZVec<REAL> elementErrors(ncel, 0.0);

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

      // // Gambiarra
      // bool flag = true;
      // int firstside = gel->FirstSide(1);
      // int lastside = gel->FirstSide(2);
      // for (int side = firstside; side < lastside; ++side) {
      //   if (gel->Dimension() != 3)
      //     continue;
      //   TPZGeoElSide gelSide(gel, side);
      //   std::set<int> matIds = {
      //                           SimData->EFarField,
      //                         };
      //   TPZGeoElSide neigh = gelSide.HasNeighbour(matIds);
      //   if (neigh)
      //     flag = false;
      // }
      // if (!flag) elementErrors[icel] = 0.0; continue;
      // // End gambiarra

      localTotalError += elementErrors[icel];
    }
    partialErrors[tid] = localTotalError;
  };

  int64_t chunk = ncel / nthreads;
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

  std::cout << "Estimated error in energy norm: " << totalError << std::endl;

  // VTK output
  std::ofstream out_estimator("EstimatedError.vtk");
  TPZVTKGeoMesh::PrintCMeshVTK(cmeshMixed, out_estimator, elementErrors, "EstimatedError");

  // Indices of geoel elements that need refinement 
  TPZVec<int64_t> needRefinement;
  for (int64_t i = 0; i < ncel; ++i) {
    if (elementErrors[i] > estimator_tol) {
      int64_t igeo = cmeshMixed->Element(i)->Reference()->Index();
      needRefinement.push_back(igeo);
    }
  }

  return needRefinement;
}

void TPZWannEstimationTools::hRefinement(TPZGeoMesh* gmesh, TPZVec<int64_t>& needRefinement, ProblemData* SimData) {
  REAL dim = gmesh->Dimension();
  REAL tol = 1e-6;
  std::set<REAL> xcoords2D;
  TPZVec<int64_t> needRefinement1D;
  TPZVec<int64_t> needRefinement2D;
  TPZVec<int64_t> needRefinement3D;

  std::set<int> matIds = {
      SimData->ESurfWellCyl, SimData->EPressure2DSkin, 
      SimData->EPressureInterface, SimData->EHDivBoundInterface
    };

  // From the initial list of "to-refine" elements, get the the centroids
  // of the faces that intersect the wellbore surface (fill xcoords2D)
  for (int64_t i = 0; i < needRefinement.size(); ++i) {
    TPZGeoEl* gel = gmesh->Element(needRefinement[i]);
    if (!gel) DebugStop();
    if (gel->Dimension() == 1) {
      needRefinement1D.push_back(gel->Index());
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

  // From initial list of "to-refine" elements, add the boundary elements that 
  // need refinement
  for (int64_t i = 0; i < needRefinement.size(); ++i) {
    TPZGeoEl* gel = gmesh->Element(needRefinement[i]);
    if (!gel) DebugStop();
    if (gel->Dimension() != 3) continue; 
    int firstside = gel->FirstSide(dim-1);
    int lastside = gel->FirstSide(dim);
    for (int side = firstside; side < lastside; ++side) {
      TPZGeoElSide gelSide(gel, side);
      TPZGeoElSide neigh = gelSide.HasNeighbour(SimData->EFarField);
      if (neigh) {
        needRefinement.push_back(neigh.Element()->Index());
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
      needRefinement.push_back(gel->Index());
    }
  }

  // From 2D elements, guarantee the 1D and 3D neighbors will be refined
  for (int64_t i = 0; i < needRefinement2D.size(); ++i) {
    TPZGeoEl* gel = gmesh->Element(needRefinement2D[i]);
    if (gel->Dimension() != 2) DebugStop();
    // Pick the 3D neighbor
    TPZGeoElSide gelSide(gel);
    TPZGeoElSide neigh3D = gelSide.HasNeighbour(SimData->EDomain);
    if (neigh3D) {
      needRefinement.push_back(neigh3D.Element()->Index()); // Will have duplicated elements
    }
    // Pick the 1D neighbor
    int firstside = gel->FirstSide(1);
    int lastside = gel->FirstSide(2);
    for (int side = firstside; side < lastside; ++side) {
      TPZGeoElSide gelSide(gel, side);
      TPZGeoElSide neigh = gelSide.HasNeighbour(SimData->ECurveWell);
      if (neigh) {
        needRefinement.push_back(neigh.Element()->Index()); // Will have duplicated elements!
      }
    }
  }

  // // Print elements to be refined per dimension
  // std::cout << needRefinement.size() << " elements to be refined initially." << std::endl;

  // std::cout << needRefinement1D.size() << " 1D elements to be refined." << std::endl;
  // for (int64_t i = 0; i < needRefinement1D.size(); ++i) {
  //   std::cout << needRefinement1D[i] << " ";
  // }
  // std::cout << std::endl;

  // std::cout << needRefinement2D.size() << " 2D elements to be refined." << std::endl;
  // for (int64_t i = 0; i < needRefinement2D.size(); ++i) {
  //   std::cout << needRefinement2D[i] << " ";
  // }
  // std::cout << std::endl;

  // std::cout << needRefinement3D.size() << " 3D elements to be refined." << std::endl;
  // for (int64_t i = 0; i < needRefinement3D.size(); ++i) {
  //   std::cout << needRefinement3D[i] << " ";
  // }
  // std::cout << std::endl;

  // Perform h-refinement on the specified elements
  for (int64_t i = 0; i < needRefinement.size(); ++i) {
    TPZVec<TPZGeoEl *> pv;
    TPZGeoEl* gel = gmesh->Element(needRefinement[i]);
    if (!gel) DebugStop();
    if (gel->HasSubElement()) continue; // We have to do this because of duplicated entries
    gel->Divide(pv);
    int banana = 0;
  }
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

void TPZWannEstimationTools::FakeRefine(TPZGeoMesh* gmesh, TPZVec<int64_t>& needRefinement, ProblemData* SimData) {
  REAL dim = gmesh->Dimension();
  REAL tol = 1e-6;
  int n_refined = 0;
  std::set<REAL> xcoords2D;
  TPZVec<int64_t> needRefinement1D;
  TPZVec<int64_t> needRefinement2D;
  TPZVec<int64_t> needRefinement3D;

  std::set<int> matIds = {
      SimData->ESurfWellCyl, SimData->EPressure2DSkin, 
      SimData->EPressureInterface, SimData->EHDivBoundInterface
    };

  // Perform h-refinement on the specified elements
  for (int64_t i = 0; i < needRefinement.size(); ++i) {
    TPZVec<TPZGeoEl *> pv;
    TPZGeoEl* gel = gmesh->Element(needRefinement[i]);
    if (gel->Dimension() != 3) continue;
    bool willRefine = true;
    int firstside = gel->FirstSide(2);
    int lastside = gel->FirstSide(3);
    for (int side = firstside; side < lastside; ++side){
      TPZGeoElSide gelSide(gel,side);
      TPZGeoElSide neigh = gelSide.HasNeighbour(matIds);
      if (neigh) {
        willRefine = false;
      }
    }
    if (!gel) DebugStop();
    if (gel->HasSubElement()) continue; // We have to do this because of duplicated entries
    if (!willRefine) continue;
    gel->Divide(pv);
    n_refined++;
  }
  std::cout << n_refined << " elements were refined in the fake refinement." << std::endl;
}

void TPZWannEstimationTools::FakeRefine2(TPZGeoMesh* gmesh, ProblemData* SimData) {
  std::set<int> matIds = {
    SimData->ESurfWellCyl, SimData->EPressure2DSkin, 
    SimData->EPressureInterface, SimData->EHDivBoundInterface,
    SimData->EFarField, SimData->ESurfHeel, SimData->ESurfToe, SimData->ECapRock
  };

  bool flag = true;

  while (flag == true) {
    std::random_device rd;
    std::mt19937 gen(rd());
    std::uniform_int_distribution<> dis(1, gmesh->NElements());
    int iel = dis(gen);

    TPZGeoEl* gel = gmesh->Element(iel);
    if (gel->Dimension() != 3) continue;
    bool willRefine = true;
    int firstside = gel->FirstSide(2);
    int lastside = gel->FirstSide(3);
    for (int side = firstside; side < lastside; ++side){
      TPZGeoElSide gelSide(gel,side);
      TPZGeoElSide neigh = gelSide.HasNeighbour(matIds);
      if (neigh) {
        willRefine = false;
      }
    }
    if (gel->HasSubElement()) continue;
    if (!willRefine) continue;
    TPZVec<TPZGeoEl *> pv;
    gel->Divide(pv);
    std::cout << "Refining element " << gel->Index() << std::endl;
    flag = false;
  }
}

void TPZWannEstimationTools::CheckRef(TPZCompMesh* cmesh) {
  std::map<int, int> matid_to_count;
  for (int64_t i = 0; i < cmesh->NElements(); ++i) {
    TPZCompEl* cel = cmesh->Element(i);
    if (!cel) continue;
    TPZGeoEl* gel = cel->Reference();
    if (gel->Father()) {
      matid_to_count[gel->MaterialId()]++;
      std::cout << "Search for index " << i << std::endl;
    }
  }
  std::cout << "Count of subelements by material id:" << std::endl;
  for (const auto& pair : matid_to_count) {
    std::cout << "MatId " << pair.first << ": " << pair.second << std::endl;
  }
}

TPZVec<int64_t> TPZWannEstimationTools::ErrorEstimationOld(TPZMultiphysicsCompMesh* cmeshMixed, TPZCompMesh* cmesh, ProblemData* SimData) {
  {
    std::ofstream clearlog("error_estimation.txt", std::ios_base::trunc);
  }

  // Ensure references point to H1 cmesh
  cmesh->Reference()->ResetReference();
  cmesh->LoadReferences();

  int64_t ncel = cmeshMixed->NElements();
  TPZAdmChunkVector<TPZCompEl *> &elementvec_m = cmeshMixed->ElementVec();
  TPZVec<REAL> elementErrors(ncel, 0.0);
  REAL totalError = 0.0;
  
  const int numIntPts = 5; // Number of integration points. Can be adjusted.

  // Loop on elements
  for (int64_t icel = 0; icel < ncel; ++icel) {

    TPZCompEl *celMixed = elementvec_m[icel];
    int matid = celMixed->Material()->Id();
    if (matid != SimData->EDomain && matid != SimData->ECurveWell) continue;

    TPZGeoEl *gel = celMixed->Reference();
    TPZCompEl *celH1 = gel->Reference();

    if (!celMixed || !celH1) continue;
    if (gel->HasSubElement()) continue;
    if (celH1->Material()->Id() != matid) DebugStop();

    REAL hk = ElementDiameter(gel); // Element size

    // Permeability
    REAL perm = SimData->m_Reservoir.perm;
    if (matid == SimData->ECurveWell) {
      perm = SimData->m_Wellbore.perm;
    }

    REAL fluxError = 0.0;
    REAL balanceError = 0.0;
    REAL conformityError = 0.0;

    TPZIntPoints *intrule = gel->CreateSideIntegrationRule(gel->NSides()-1, numIntPts);
    if (!intrule) continue;
  
    int npts = intrule->NPoints();
    for (int ip = 0; ip < npts; ++ip) {
      TPZManVector<REAL,3> ptInElement(gel->Dimension());
      REAL weight;
      intrule->Point(ip, ptInElement, weight);
    
      // Get Jacobian for integration
      TPZFMatrix<REAL> jacobian, axes, jacinv;
      REAL detjac;
      gel->Jacobian(ptInElement, jacobian, axes, detjac, jacinv);
      weight *= fabs(detjac);
    
      TPZManVector<REAL,3> x(3,0.0);
      gel->X(ptInElement, x); // Real coordinates for x

      TPZManVector<REAL,3> force(1,0.0); // f = 0 in this case

      // Compute H1 term
      TPZVec<REAL> termH1(3,0.0);
      celH1->Solution(ptInElement, 7, termH1);

      // Compute L2 term
      TPZVec<REAL> termL2(3,0.0);
      celMixed->Solution(ptInElement, 9, termL2);
      for (int d = 0; d < 3; ++d) {
        termL2[d] = -perm * termL2[d];
      }

      // Compute Hdiv term and Hdiv flux divergence
      TPZVec<REAL> termHdiv(3,0.0);
      TPZVec<REAL> divFluxMixed(1,0.0);
      TPZMultiphysicsElement *celMulti = dynamic_cast<TPZMultiphysicsElement*>(celMixed);
      if (celMulti) {
        celMulti->Solution(ptInElement, 1, termHdiv);
        celMulti->Solution(ptInElement, 5, divFluxMixed);
      }

      // Flux contribution
      REAL diffFlux = 0.0;
      for (int d = 0; d < gel->Dimension(); ++d) {
        REAL diff = termL2[d] - termHdiv[d];
        diffFlux += diff * diff;
      }

      // Balance contribution
      REAL diffBalance = (divFluxMixed[0] - force[0]) * (divFluxMixed[0] - force[0]);

      // Conformity contribution
      REAL diffConformity = 0.0;
      for (int d = 0; d < 3; ++d) {
        REAL diff = termH1[d] - termL2[d];
        diffConformity += diff * diff;
      }

      // Add weighted contribution to element error
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

    // In the wellbore, forcing function is exactly represented by the divergence space
    if (matid == SimData->ECurveWell) {
      balanceError = 0.0;
    }

    elementErrors[icel] = pow((sqrt(fluxError) + (hk/M_PI)*sqrt(balanceError)), 2)
      + conformityError;
    totalError += elementErrors[icel];

    delete intrule;
  }
  
  totalError = sqrt(totalError);

  if (SimData->m_VerbosityLevel > 0) {
    std::ofstream errorlog_total("error_estimation.txt", std::ios_base::app);
    if (errorlog_total.is_open()) {
      errorlog_total << "\nError estimation in energy norm: " << totalError << std::endl;
      errorlog_total.close();
    }
  }

  // VTK output
  std::ofstream out_estimator("EstimatedErrorReservoir.vtk");
  TPZVTKGeoMesh::PrintCMeshVTK(cmeshMixed, out_estimator, elementErrors, "EstimatedErrorReservoir");

  // Indices of geoel elements that need refinement 
  TPZVec<int64_t> needRefinement;
  for (int64_t i = 0; i < ncel; ++i) {
    if (elementErrors[i] > estimator_tol) {
      int64_t igeo = cmeshMixed->Element(i)->Reference()->Index();
      needRefinement.push_back(igeo);
    }
  }

  return needRefinement;
}