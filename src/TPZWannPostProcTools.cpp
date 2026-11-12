/**
 * @file
 * @brief Solution export, surface-flux integration, and productivity calculations.
 * See the matching header for function contracts and docs/ for the architecture.
 */

#include "TPZWannPostProcTools.h"
#include "TPZWannGeometryTools.h"

#include <algorithm>
#include <cmath>
#include <fstream>
#include <iomanip>
#include <limits>
#include <memory>
#include <vector>
#include "pzcondensedcompel.h"
#include "TPZWannMixedWell.h"
#include "TPZWannH1Well.h"

#include <set>

namespace {
constexpr REAL barrelsPerCubicMeter = 6.28981;
constexpr REAL secondsPerDay = 86400.;
constexpr REAL pascalsPerPsi = 6894.76;

bool FindWellEndpoints(TPZGeoMesh* mesh, int heelId, int toeId,
                       TPZVec<REAL>& heelPoint, TPZVec<REAL>& toePoint) {
  bool foundHeel = false, foundToe = false;
  for (auto* element : mesh->ElementVec()) {
    if (!element || element->Dimension() != 0) continue;
    if (element->MaterialId() == heelId) {
      element->NodePtr(0)->GetCoordinates(heelPoint);
      foundHeel = true;
    } else if (element->MaterialId() == toeId) {
      element->NodePtr(0)->GetCoordinates(toePoint);
      foundToe = true;
    }
  }
  return foundHeel && foundToe;
}
}

void TPZWannPostProcTools::WriteWellboreVTK(TPZCompMesh* cmesh, ProblemData* SimData) {
  auto& PostProcData = SimData->m_PostProc;

  // Check whether the simulation is H1 or Hdiv
  bool isMultiphysics = true;
  TPZMultiphysicsCompMesh *cmeshMP = dynamic_cast<TPZMultiphysicsCompMesh *>(cmesh);
  if (!cmeshMP) isMultiphysics = false;

  std::string simType = isMultiphysics ? "Hdiv" : "H1";

  TPZStack<std::string> fieldnames; 
  fieldnames.Push("Pressure");
  fieldnames.Push("Flux");
  if (isMultiphysics) fieldnames.Push("Divergence");

  for (size_t wellIndex = 0; wellIndex < SimData->m_Wellbore.size(); wellIndex++) {
    const auto& well = SimData->m_Wellbore[wellIndex];
    const std::string filename = PostProcData.wellbore_vtk + "_" + simType
                               + "_" + std::to_string(wellIndex);
    TPZVTKGenerator vtk(cmesh, std::set<int>{well.matid}, fieldnames,
                       filename, PostProcData.vtk_resolution);
    vtk.SetNThreads(PostProcData.nthreads);
    vtk.Do();
  }
}

bool TPZWannPostProcTools::WriteWellProfilesCSV(TPZMultiphysicsCompMesh* cmesh, const ProblemData& data,
                                                int subdivisionsPerElement) {
  const auto farfield = data.m_Reservoir.BCs.find("surface_farfield");
  const REAL pressureScale = data.m_Numerics.pressureScale;
  if (farfield == data.m_Reservoir.BCs.end() || farfield->second.type != 0 ||
      !std::isfinite(pressureScale) || pressureScale <= 0. || subdivisionsPerElement < 1) {
    std::cerr << "Well CSV export requires a far-field pressure, positive pressure scale, "
                 "and at least one subdivision per well element.\n";
    return false;
  }

  struct WellSegment {
    TPZCompEl* element;
    REAL start;
    REAL end;
    REAL direction;
  };

  for (int iw = 0; iw < data.m_Wellbore.size(); ++iw) {
    const auto& well = data.m_Wellbore[iw];
    const auto& heel = well.BCs.at("point_heel");
    const auto& toe = well.BCs.at("point_toe");
    TPZManVector<REAL, 3> heelPoint(3), toePoint(3);
    if (!FindWellEndpoints(cmesh->Reference(), heel.matid, toe.matid, heelPoint, toePoint)) {
      std::cerr << "Heel/toe points not found for well " << well.name << '\n';
      return false;
    }
    const TPZManVector<REAL, 3> axis = toePoint - heelPoint;
    const REAL length = Norm(axis);
    if (!(length > 0.) || !std::isfinite(length)) return false;

    std::vector<WellSegment> segments;
    for (auto* element : cmesh->ElementVec()) {
      if (!element) continue;
      auto* geo = element->Reference();
      if (!geo || geo->Dimension() != 1 || geo->MaterialId() != well.matid ||
          geo->HasSubElement()) continue;
      TPZManVector<REAL, 3> first(3), second(3);
      geo->NodePtr(0)->GetCoordinates(first);
      geo->NodePtr(1)->GetCoordinates(second);
      const REAL a = TPZWannGeometryTools::ComputeAxialCoordinate(first, heelPoint, axis) / length;
      const REAL b = TPZWannGeometryTools::ComputeAxialCoordinate(second, heelPoint, axis) / length;
      segments.push_back({element, std::min(a, b), std::max(a, b), a < b ? 1. : -1.});
    }
    if (segments.empty()) {
      std::cerr << "No mixed well elements found for " << well.name << '\n';
      return false;
    }
    std::sort(segments.begin(), segments.end(), [](const auto& a, const auto& b) {
      return a.start < b.start;
    });

    const auto* material = segments.front().element->Material();
    const int fluxVar = material->VariableIndex("Flux");
    const int divergenceVar = material->VariableIndex("Divergence");
    const int pressureVar = material->VariableIndex("Pressure");
    TPZManVector<REAL, 1> qsi(1, -segments.front().direction);
    TPZManVector<STATE, 1> value(1);
    segments.front().element->Solution(qsi, fluxVar, value);
    const REAL heelRate = std::abs(heel.type == 1 ? heel.value : value[0]);

    std::string suffix = data.m_PostProc.wellbore_vtk;
    if (data.m_Wellbore.size() > 1) suffix += "_" + std::to_string(iw);
    const std::string fluxPath = "flux" + suffix + ".csv";
    const std::string inflowPath = "divergence" + suffix + ".csv";
    const std::string pressurePath = "pressure" + suffix + ".csv";
    std::ofstream flux(fluxPath), inflow(inflowPath), pressure(pressurePath);
    if (!flux || !inflow || !pressure) {
      std::cerr << "Could not open well profile CSV files for " << well.name << '\n';
      return false;
    }
    for (auto* stream : {&flux, &inflow, &pressure}) {
      *stream << "x,Curve1\n" << std::setprecision(std::numeric_limits<REAL>::max_digits10);
    }

    for (size_t i = 0; i < segments.size(); ++i) {
      const auto& segment = segments[i];
      // Keep a single sample at interfaces, using the element on the toe side.
      const int last = i + 1 == segments.size() ? subdivisionsPerElement : subdivisionsPerElement - 1;
      for (int point = 0; point <= last; ++point) {
        const REAL fraction = REAL(point) / subdivisionsPerElement;
        const REAL position = segment.start + fraction * (segment.end - segment.start);
        const REAL x = point == subdivisionsPerElement ? 1. : std::clamp(position, REAL(0), REAL(1));
        qsi[0] = segment.direction * (2. * fraction - 1.);
        segment.element->Solution(qsi, fluxVar, value);
        flux << x << ',' << std::abs(value[0]) * barrelsPerCubicMeter * secondsPerDay << '\n';
        segment.element->Solution(qsi, divergenceVar, value);
        inflow << x << ',';
        if (heelRate > 0.) inflow << std::abs(value[0]) * length / heelRate;
        inflow << '\n';
        segment.element->Solution(qsi, pressureVar, value);
        pressure << x << ',' << (farfield->second.value - value[0]) / (pressureScale * pascalsPerPsi) << '\n';
      }
    }
    flux.close();
    inflow.close();
    pressure.close();
    if (!flux || !inflow || !pressure) {
      std::cerr << "Could not write well profile CSV files for " << well.name << '\n';
      return false;
    }
  }
  return true;
}

void TPZWannPostProcTools::WriteReservoirVTK(TPZCompMesh* cmesh, ProblemData* SimData) {
  auto& PostProcData = SimData->m_PostProc;

  // Check whether the simulation is H1 or Hdiv
  bool isMultiphysics = true;
  TPZMultiphysicsCompMesh *cmeshMP = dynamic_cast<TPZMultiphysicsCompMesh *>(cmesh);
  if (!cmeshMP) isMultiphysics = false;

  std::string simType = isMultiphysics ? "Hdiv" : "H1";
  std::string filename = PostProcData.reservoir_vtk + "_" + simType;

  TPZStack<std::string> fieldnames; 
  fieldnames.Push("Pressure");
  fieldnames.Push("Flux");

  bool safe_check = SimData->m_Mesh.ToCylindrical? true : false;

  TPZVTKGenerator vtk(cmesh, fieldnames, filename, PostProcData.vtk_resolution, 3, safe_check);
  vtk.SetNThreads(PostProcData.nthreads);
  vtk.Do();
}

void TPZWannPostProcTools::WriteVTKs(TPZCompMesh* cmesh, ProblemData* SimData) {
  WriteReservoirVTK(cmesh, SimData);
  WriteWellboreVTK(cmesh, SimData);
}

TPZVec<REAL> TPZWannPostProcTools::ComputeWellFluxes(TPZCompMesh *cmesh, ProblemData *SimData, int wellIndex, int nSegments) {
  if (wellIndex < 0 || wellIndex >= SimData->m_Wellbore.size()) {
    std::cerr << "Invalid well index: " << wellIndex << std::endl;
    DebugStop();
  }

  if (nSegments <= 0) {
    std::cerr << "Invalid number of segments: " << nSegments << std::endl;
    DebugStop();
  }

  TPZVec<REAL> fluxes(nSegments, 0.);

  int matid = SimData->m_Reservoir.matid;
  const auto &well = SimData->m_Wellbore[wellIndex];
  const int wellMatid = well.matidSurf;

  // Ensure references point to current mesh
  cmesh->Reference()->ResetReference();
  cmesh->LoadReferences();
  TPZGeoMesh *gmesh = cmesh->Reference();

  // Locate this well's endpoints to compute axil coordinates and length
  const auto heelBC = well.BCs.find("point_heel");
  const auto toeBC = well.BCs.find("point_toe");
  if (heelBC == well.BCs.end() || toeBC == well.BCs.end() ||
      heelBC->second.matid < 0 || toeBC->second.matid < 0 ||
      heelBC->second.matid == toeBC->second.matid) {
    std::cerr << "Invalid heel/toe material IDs for well " << wellIndex << std::endl;
    DebugStop();
  }
  TPZManVector<REAL, 3> heelPoint(3, 0.), toePoint(3, 0.);
  if (!FindWellEndpoints(gmesh, heelBC->second.matid, toeBC->second.matid, heelPoint, toePoint)) {
    std::cerr << "Heel/toe points not found for well " << wellIndex << std::endl;
    DebugStop();
  }
  const TPZManVector<REAL, 3> axis = toePoint - heelPoint;
  const REAL length = Norm(axis);
  if (!(length > 0.) || !std::isfinite(length)) DebugStop();
  const REAL axialTolerance = 1.e-8 * length;

  for (auto* wall : gmesh->ElementVec()) {
    if (!wall || wall->HasSubElement() || wall->MaterialId() != wellMatid) continue;
    TPZGeoElSide face(wall);
    const auto reservoirSide = face.HasNeighbour(matid);
    if (!reservoirSide) DebugStop();
    auto* reservoir = reservoirSide.Element();
    auto* element = reservoir->Reference();
    if (auto* condensed = dynamic_cast<TPZCondensedCompEl*>(element))
      element = condensed->ReferenceCompEl();
    if (!element) DebugStop();
    const int fluxVar = element->Material()->VariableIndex("Flux");

    TPZManVector<REAL, 3> center(3);
    face.CenterX(center);
    const REAL axialCoord = TPZWannGeometryTools::ComputeAxialCoordinate(center, heelPoint, axis);
    if (!std::isfinite(axialCoord) || axialCoord < -axialTolerance || axialCoord > length + axialTolerance) {
      std::cerr << "Face centroid outside well " << wellIndex << std::endl;
      DebugStop();
    }
    const REAL boundedCoord = std::clamp(axialCoord, REAL(0), length);
    const int segment = std::min(nSegments - 1, static_cast<int>(boundedCoord / length * nSegments));

    TPZTransform<> faceToSide(wall->Dimension());
    face.SideTransform3(reservoirSide, faceToSide);
    auto sideToVolume = reservoir->SideToSideTransform(reservoirSide.Side(), reservoir->NSides() - 1);
    auto faceToVolume = sideToVolume.Multiply(faceToSide);
    std::unique_ptr<TPZIntPoints> rule(wall->CreateSideIntegrationRule(
        face.Side(), 2 * SimData->m_Numerics.reservoirPorder + 4));
    for (int ip = 0; ip < rule->NPoints(); ++ip) {
      TPZManVector<REAL,3> qsi(wall->Dimension()), volumePoint(reservoir->Dimension()), normal(3);
      TPZFNMatrix<9, REAL> jac, axes, inverse;
      REAL weight, determinant;
      rule->Point(ip, qsi, weight);
      wall->Jacobian(qsi, jac, axes, determinant, inverse);
      face.Normal(qsi, reservoir, wall, normal);
      faceToVolume.Apply(qsi, volumePoint);
      TPZManVector<STATE,3> flux(3, 0.);
      element->Solution(volumePoint, fluxVar, flux);
      fluxes[segment] += Dot(flux, normal) * weight * std::abs(determinant);
    }
  }
  return fluxes;
}

REAL TPZWannPostProcTools::ProductivityIndex(TPZCompMesh *cmesh, ProblemData *SimData, int wellIndex) {
  int matidWell = SimData->m_Wellbore[wellIndex].matid;
  int matidPointHeel = SimData->m_Wellbore[wellIndex].BCs["point_heel"].matid;
  int matidPointToe = SimData->m_Wellbore[wellIndex].BCs["point_toe"].matid;
  
  // Check if toe and heel matids are valid
  if (matidPointHeel < 0 || matidPointToe < 0) {
    std::cerr << "Invalid heel/toe material IDs for well " << wellIndex << std::endl;
    DebugStop();
  }

  REAL pff = SimData->m_Reservoir.BCs["surface_farfield"].value;

  cmesh->Reference()->ResetReference();
  cmesh->LoadReferences();

  // Check whether the simulation is H1 or Hdiv
  bool isMultiphysics = true;
  TPZMultiphysicsCompMesh *cmeshMP = dynamic_cast<TPZMultiphysicsCompMesh *>(cmesh);
  if (!cmeshMP) isMultiphysics = false;

  TPZGeoMesh* gmesh = cmesh->Reference();
  REAL totalFlux = 0.;
  REAL pheel = 0.;

  for (TPZCompEl* cel : cmesh->ElementVec()) {
    if (!cel) continue;
    TPZGeoEl* gel = cel->Reference();
    if (!gel) DebugStop();
    if (gel->HasSubElement()) DebugStop();
    if (gel->MaterialId() != matidWell) continue;

    TPZManVector<REAL, 3> qsi(gel->Dimension(), 0.);
    TPZGeoElSide gelsideA(gel, 0);
    TPZGeoElSide gelsideB(gel, 1);
    TPZGeoElSide neighsideA = gelsideA.HasNeighbour(matidPointHeel);
    TPZGeoElSide neighsideB = gelsideB.HasNeighbour(matidPointHeel);

    if (neighsideA) {
      qsi[0] = -1.0;
    } else if (neighsideB) {
      qsi[0] = 1.0;
    } else {
      continue; // We want an element that is connected to the heel point
    }

    TPZMaterial* mat = cel->Material();
    int pind;
    int qind;
    if (isMultiphysics) {
      TPZWannMixedWell* wellmat = dynamic_cast<TPZWannMixedWell*>(mat);
      if (!wellmat) DebugStop();
      pind = wellmat->VariableIndex("Pressure");
      qind = wellmat->VariableIndex("Flux");
    } else {
      TPZWannH1Well* wellmat = dynamic_cast<TPZWannH1Well*>(mat);
      if (!wellmat) DebugStop();
      pind = wellmat->VariableIndex("Pressure");
      qind = wellmat->VariableIndex("Flux");
    }

    TPZManVector<STATE, 3> output(1);
    cel->Solution(qsi,pind,output);
    pheel = output[0];
    cel->Solution(qsi, qind, output);
    totalFlux = output[0];
  }

  const REAL pressureScale = SimData->m_Numerics.pressureScale;

  // p_scaled = (p_Pa - refPressure) * pressureScale, so the reference
  // cancels in the difference and deltaP_Pa = deltaP_scaled / pressureScale.
  // Flux remains in m^3/s. Preserve the existing signed pressure difference.
  const REAL pressureDifferencePa = (pheel - pff) / pressureScale;
  const REAL PI = totalFlux / pressureDifferencePa;

  // Convert m^3/(s Pa) to STB/(day psi), using the existing volume convention.
  return PI * barrelsPerCubicMeter * secondsPerDay * pascalsPerPsi;
}
