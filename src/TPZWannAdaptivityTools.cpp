/**
 * @file
 * @brief Error estimation and well-compatible refinement.
 * See the matching header for function contracts and docs/ for the architecture.
 */

#include "TPZWannAdaptivityTools.h"
#include "TPZWannGeometryTools.h"
#include "TPZMaterial.h"
#include "TPZVTKGeoMesh.h"
#include "pzcondensedcompel.h"

#include <algorithm>
#include <cmath>
#include <fstream>
#include <iostream>
#include <map>
#include <numeric>
#include <set>
#include <thread>
#include <vector>

REAL TPZWannAdaptivityTools::PragerSynge(TPZMultiphysicsCompMesh* cmeshMixed, TPZCompMesh* cmeshH1, ProblemData* SimData, TPZVec<REAL>& elementErrors, int nthreads) {
  // If nthreads = 0, we use 1 thread
  nthreads++;

  // Ensure references point to H1 cmesh
  cmeshH1->Reference()->ResetReference();
  cmeshH1->LoadReferences();

  const int dim = cmeshMixed->Dimension();
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

      if (!celMixed || !celMixed->Material()) continue;
      int matid = celMixed->Material()->Id();
      if (matid != SimData->m_Reservoir.matid) continue;

      TPZGeoEl *gel = celMixed->Reference();
      TPZCompEl *celH1 = gel->Reference();

      // Check if H1 element is condensed
      condEl = dynamic_cast<TPZCondensedCompEl *>(celH1);
      if (condEl) {
        // If compel is condensed, load solution on the unconsensed compel
        condEl->LoadSolution();
        celH1 = condEl->ReferenceCompEl();
      }

      if (!celH1 || gel->HasSubElement()) DebugStop();
      if (celH1->Material()->Id() != matid) DebugStop();

      REAL hk = ElementDiameter(gel);
      TPZManVector<REAL,3> perm = SimData->m_Reservoir.perm;
      for (auto& value : perm) value /= SimData->m_Fluid[0].viscosity;
      TPZManVector<REAL,3> sqrtPerm(3,0.0);
      for (int d = 0; d < dim; ++d) {
        sqrtPerm[d] = sqrt(perm[d]);
      }
      REAL minPerm = *std::min_element(perm.begin(), perm.end());

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
          termH1[d] = sqrtPerm[d] * termH1[d];
        }

        // Compute Hdiv term -K^(-1/2) sig and div(sig)
        TPZManVector<REAL,3> termHdiv(dim, 0.0);
        TPZManVector<REAL,1> divFluxMixed(1, 0.0);
        celMixed->Solution(ptInElement, 1, termHdiv);
        celMixed->Solution(ptInElement, 5, divFluxMixed);
        for (int d = 0; d < dim; ++d) {
          termHdiv[d] = (-1./sqrtPerm[d]) * termHdiv[d];
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
      elementErrors[igeo] = sqrt(fluxError) + (hk/(M_PI*sqrt(minPerm)))*sqrt(balanceError);
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

void TPZWannAdaptivityTools::SimpleMarking(const TPZVec<REAL>& elementErrors, TPZVec<int>& refinementIndicator, REAL theta) {
  const int64_t nels = elementErrors.size();
  refinementIndicator.Resize(nels);
  refinementIndicator.Fill(0);

  if (nels == 0) {
    return;
  }

  const REAL maxError = *std::max_element(elementErrors.begin(), elementErrors.end());
  const REAL threshold = theta * maxError;

  for (int64_t iel = 0; iel < nels; ++iel) {
    if (elementErrors[iel] > 0. && elementErrors[iel] >= threshold) {
      refinementIndicator[iel] = 1;
    }
  }
}

void TPZWannAdaptivityTools::MeshSmoothing(TPZGeoMesh* gmesh, TPZVec<int>& refinementIndicator) {
  // If an element has most of its neighbors refined, then refine it too
  for (int64_t iel = 0; iel < refinementIndicator.size(); ++iel) {
    if (refinementIndicator[iel] != 1) continue;
    TPZGeoEl* gel = gmesh->Element(iel);
    if (!gel) DebugStop();
    if (gel->Dimension() == 0 || gel->HasSubElement()) continue;
    int firstside = gel->FirstSide(gel->Dimension()-1);
    int lastside = gel->FirstSide(gel->Dimension());
    for (int side = firstside; side < lastside; ++side) {
      TPZGeoElSide gelSide(gel, side);
      TPZGeoElSide neigh = gelSide.Neighbour();
      TPZGeoEl* neighGel = neigh.Element();
      if (!neighGel) DebugStop();
      if (neighGel->Dimension() == 0 || neighGel->HasSubElement()) continue;
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
  const int reservoirId = SimData->m_Reservoir.matid;
  const std::set<int> couplingIds = {
      SimData->EPressure2DSkin, SimData->EPressureInterface, SimData->EHDivBoundInterface
  };
  std::vector<std::vector<int64_t>> groups;

  for (const auto& well : SimData->m_Wellbore) {
    const int heelId = well.BCs.at("point_heel").matid;
    const int toeId = well.BCs.at("point_toe").matid;
    TPZGeoEl *heel = nullptr, *toe = nullptr;
    std::vector<TPZGeoEl*> lines, walls;
    for (auto* gel : gmesh->ElementVec()) {
      if (!gel || gel->HasSubElement()) continue;
      const int id = gel->MaterialId();
      if (gel->Dimension() == 0 && id == heelId) heel = gel;
      if (gel->Dimension() == 0 && id == toeId) toe = gel;
      if (gel->Dimension() == 1 && id == well.matid) lines.push_back(gel);
      if (gel->Dimension() == 2 && id == well.matidSurf) walls.push_back(gel);
    }
    if (!heel || !toe || lines.empty() || walls.empty()) DebugStop();

    TPZManVector<REAL,3> origin(3), endpoint(3);
    heel->NodePtr(0)->GetCoordinates(origin);
    toe->NodePtr(0)->GetCoordinates(endpoint);
    const TPZManVector<REAL,3> axis = endpoint - origin;
    auto axialCenter = [&](TPZGeoEl* gel) {
      TPZManVector<REAL,3> center(3);
      TPZGeoElSide(gel).CenterX(center);
      return TPZWannGeometryTools::ComputeAxialCoordinate(center, origin, axis);
    };

    constexpr REAL tolerance = 1.e-6;
    std::set<REAL> centers;
    std::map<REAL, size_t> groupByCenter;
    for (auto* line : lines) {
      const REAL axial = axialCenter(line);
      TPZWannGeometryTools::InsertXCoorInSet(axial, centers, tolerance);
      if (!groupByCenter.emplace(axial, groups.size()).second) DebugStop();
      groups.push_back({line->Index()});
    }
    for (auto* wall : walls) {
      const REAL axial = axialCenter(wall);
      const REAL center = TPZWannGeometryTools::FindClosestX(axial, centers, tolerance);
      auto& group = groups.at(groupByCenter.at(center));
      group.push_back(wall->Index());
      TPZGeoElSide face(wall);
      int reservoirNeighbors = 0;
      for (auto neighbor = face.Neighbour(); neighbor != face; ++neighbor) {
        auto* gel = neighbor.Element();
        if (gel->HasSubElement()) continue;
        const int id = gel->MaterialId();
        if (id == reservoirId) ++reservoirNeighbors;
        if (id == reservoirId || couplingIds.count(id)) group.push_back(gel->Index());
      }
      if (reservoirNeighbors != 1) DebugStop();
    }
  }

  // A reservoir element can touch more than one group. Close the marks over
  // all such groups, including groups on different wells if they share a cell.
  bool changed;
  do {
    changed = false;
    for (const auto& group : groups) {
      if (!std::any_of(group.begin(), group.end(), [&](int64_t index) {
            return refinementIndicator[index] != 0;
          })) continue;
      for (const auto index : group) {
        if (refinementIndicator[index]) continue;
        refinementIndicator[index] = 1;
        changed = true;
      }
    }
  } while (changed);

  std::set<int> boundaryIds;
  for (const auto& [name, bc] : SimData->m_Reservoir.BCs)
    boundaryIds.insert(bc.matid);
  for (int64_t i = 0; i < refinementIndicator.size(); ++i) {
    if (!refinementIndicator[i]) continue;
    auto* gel = gmesh->Element(i);
    if (gel->MaterialId() != reservoirId) continue;
    for (int side = gel->FirstSide(gel->Dimension()-1); side < gel->FirstSide(gel->Dimension()); ++side) {
      const auto boundary = TPZGeoElSide(gel, side).HasNeighbour(boundaryIds);
      if (boundary && !boundary.Element()->HasSubElement())
        refinementIndicator[boundary.Element()->Index()] = 1;
    }
  }
}

void TPZWannAdaptivityTools::AdaptivityProcess(TPZGeoMesh* gmesh, ProblemData* SimData, TPZVec<REAL>& elementErrors, TPZVec<int64_t>& refinedElements, REAL tol) {
  if (elementErrors.size() != gmesh->NElements()) {
    std::cerr << "Error: elementErrors size does not match number of geometric elements." << std::endl;
    DebugStop();
  }

  TPZVec<int> refinementIndicator(gmesh->NElements(), 0);
    
  // Mark elements for refinement
  SimpleMarking(elementErrors, refinementIndicator, tol);
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

  // List of elements to refine
  refinedElements.resize(0);
  for (int64_t iel = 0; iel < refinementIndicator.size(); ++iel) {
    if (refinementIndicator[iel] == 1) {
      refinedElements.push_back(iel);
    }
  }

  // Finally perform the refinement.
  TPZWannGeometryTools::hRefinement(gmesh, refinedElements);
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
