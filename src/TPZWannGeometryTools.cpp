/**
 * @file
 * @brief Geometry import, coupling surfaces, axial ordering, and refinement helpers.
 * See the matching header for function contracts and docs/ for the architecture.
 */

#include "TPZWannGeometryTools.h"
#include "TPZRefPatternTools.h"
#include <cmath>

TPZGeoMesh* TPZWannGeometryTools::CreateGeoMesh(ProblemData* simData) {

  // Import mesh from gmsh file
  TPZGeoMesh* gmesh = ReadMeshFromGmsh(simData);

  // Verify consistency of initial mesh
  bool hasErrors = TPZWannGeometryTools::VerifyMesh(gmesh, simData);

  // Plot mesh if it has errors
  if (hasErrors) {
    std::ofstream out("gmeshWithErrors.vtk");
    TPZVTKGeoMesh::PrintGMeshVTK(gmesh, out);
    std::cout << "Mesh has errors! Plotted elements with issues in "
                 "gmeshWithErrors.vtk"
              << std::endl;
    DebugStop();
  }

  // This flag set all wellbore surfaces to be cylindrical
  // Should we define a well-wise flag?
  if (simData->m_Mesh.ToCylindrical) {
    for (int iw = 0; iw < simData->m_Wellbore.size(); iw++) {
      ModifyGeometricMeshToCylWell(gmesh, *simData, iw);
    }
  }

  // Create geo elements for wellbore-reservoir coupling
  // Such auxiliary elements live in the wellbore cylindrical surface
  CreateCouplingEls(gmesh, simData);

  // Guarantee that nodes in the wellbore surfaces are in a increasing order along the wellbore axis
  OrderWellIds(gmesh, simData);

  if (simData->m_PostProc.verbosityLevel > 0) {
    std::ofstream out("gmeshWannFinal.vtk");
    TPZVTKGeoMesh::PrintGMeshVTK(gmesh, out);
  }

  return gmesh;
}

TPZGeoMesh* TPZWannGeometryTools::ReadMeshFromGmsh(ProblemData* simData){

  std::string file = simData->m_Mesh.file;
  std::string path(std::string(INPUTDIR) + "/" + file);
  TPZGeoMesh* gmesh = new TPZGeoMesh();
  {
    TPZGmshReader reader;
    reader.GeometricGmshMesh(path, gmesh);

    // Remove gmsh boundary elements and create GeoElBC so normals are consistent
    int64_t nel = gmesh->NElements();
    for(int64_t el = 0; el < nel; el++){
        TPZGeoEl *gel = gmesh->Element(el);
        if(!gel || gel->Dimension() != gmesh->Dimension()-1) continue;
        TPZGeoElSide gelside(gel);
        TPZGeoElSide neigh = gelside.Neighbour();
        gel->RemoveConnectivities();
        int matid = gel->MaterialId();
        delete gel;
        TPZGeoElBC gbc(neigh, matid);
    }

    if (simData->m_PostProc.verbosityLevel > 0) {
      std::ofstream out("gmeshWannOrig.vtk");
      TPZVTKGeoMesh::PrintGMeshVTK(gmesh, out);
    }
  }
  return gmesh;
}

void TPZWannGeometryTools::ModifyGeometricMeshToCylWell(TPZGeoMesh *gmesh, const ProblemData &simData, int wellIndex) {
  if (wellIndex < 0 || wellIndex >= simData.m_Wellbore.size()) DebugStop();
  const auto &well = simData.m_Wellbore[wellIndex];
  TPZManVector<REAL,3> heelCoord(3, 0.), toeCoord(3, 0.);
  const int heelMatid = well.BCs.at("point_heel").matid;
  const int toeMatid = well.BCs.at("point_toe").matid;
  int heelCount = 0, toeCount = 0;
  for (auto gel : gmesh->ElementVec()) {
    if (!gel || gel->HasSubElement() || gel->Dimension() != 0) continue;
    if (gel->MaterialId() == heelMatid) {
      gel->NodePtr(0)->GetCoordinates(heelCoord);
      heelCount++;
    }
    if (gel->MaterialId() == toeMatid) {
      gel->NodePtr(0)->GetCoordinates(toeCoord);
      toeCount++;
    }
  }
  if (heelCount != 1 || toeCount != 1) DebugStop();

  TPZManVector<REAL,3> cylaxis = toeCoord - heelCoord;
  const REAL axisNorm = Norm(cylaxis);
  if (!(axisNorm > 0.) || !std::isfinite(axisNorm)) DebugStop();
  for (int d = 0; d < 3; d++) cylaxis[d] /= axisNorm;
  const int matid = well.matidSurf;
  const REAL cylradius = well.radius;
  if (!(cylradius > 0.) || !std::isfinite(cylradius)) DebugStop();

  // The tagged well curve lies on the wall. Find the actual heel circle center
  // from nodes shared by the heel cap and cylindrical surface.
  std::set<int64_t> capNodes, cylinderNodes;
  for (auto gel : gmesh->ElementVec()) {
    if (!gel || gel->HasSubElement() || gel->Dimension() != 2) continue;
    if (gel->MaterialId() != well.matidHeelSurf && gel->MaterialId() != matid) continue;
    auto &nodes = gel->MaterialId() == matid ? cylinderNodes : capNodes;
    for (int in = 0; in < gel->NCornerNodes(); in++) nodes.insert(gel->NodeIndex(in));
  }
  std::vector<TPZManVector<REAL,3>> rim;
  for (auto node : capNodes) {
    if (!cylinderNodes.count(node)) continue;
    TPZManVector<REAL,3> x(3, 0.);
    gmesh->NodeVec()[node].GetCoordinates(x);
    rim.push_back(x);
  }
  if (rim.size() < 3) DebugStop();

  const auto &origin = rim[0];
  TPZManVector<REAL,3> chord(3, 0.), transverse(3, 0.);
  REAL chordLength = 0.;
  for (const auto &x : rim) {
    TPZManVector<REAL,3> delta = x - origin;
    if (Norm(delta) > chordLength) {
      chord = delta;
      chordLength = Norm(delta);
    }
  }
  if (!(chordLength > 0.)) DebugStop();
  TPZManVector<REAL,3> chordDirection = chord;
  for (int d = 0; d < 3; d++) chordDirection[d] /= chordLength;
  REAL transverseLength = 0., centerOffset = 0.;
  for (const auto &x : rim) {
    TPZManVector<REAL,3> delta = x - origin;
    const REAL along = Dot(delta, chordDirection);
    TPZManVector<REAL,3> perpendicular = delta;
    for (int d = 0; d < 3; d++) perpendicular[d] -= along * chordDirection[d];
    const REAL distance = Norm(perpendicular);
    if (distance > transverseLength) {
      transverse = perpendicular;
      transverseLength = distance;
      centerOffset = (Dot(delta, delta) - chordLength * along) / (2. * distance);
    }
  }
  if (!(transverseLength > 1.e-10 * chordLength)) DebugStop();
  TPZManVector<REAL,3> cylCenter = origin;
  for (int d = 0; d < 3; d++) {
    cylCenter[d] += 0.5 * chord[d] + centerOffset * transverse[d] / transverseLength;
  }

  // Validate before changing any nodes: the rim must be a circle normal to the axis.
  const REAL tolerance = 1.e-10 + 1.e-8 * cylradius;
  for (const auto &x : rim) {
    TPZManVector<REAL,3> radial = x - cylCenter;
    const REAL radius = Norm(radial);
    if (!std::isfinite(radius) || std::abs(radius - cylradius) > tolerance ||
        std::abs(Dot(radial, cylaxis)) > tolerance) {
      PZError << "Inconsistent heel circle for well " << wellIndex
              << ": computed radius " << radius << ", expected " << cylradius << std::endl;
      DebugStop();
    }
  }

  int64_t nel = gmesh->NElements();
  for(int64_t iel = 0; iel < nel ; iel++) {
    TPZGeoEl* gel = gmesh->Element(iel);
    if(!gel) continue;
    if(gel->HasSubElement()) DebugStop();
    if(gel->MaterialId() != matid) continue;
    
    TPZManVector<int64_t, 4> nodeindices;
    gel->GetNodeIndices(nodeindices);
    const int nnodes = gel->NCornerNodes();
    //Moving the nodes to the cylinder surface
    TPZManVector<REAL,3> xnode(3,0);
    for (int in = 0; in < nnodes; in++)
    {
      gel->NodePtr(in)->GetCoordinates(xnode);
      // component of xnode that is orthogonal to cyl axis
      TPZManVector<REAL, 3> x_orth = xnode - cylCenter;
      const REAL dax = Dot(x_orth, cylaxis);
      for (int ix = 0; ix < 3; ix++)
      {
        x_orth[ix] -= dax * cylaxis[ix];
      }
      const auto normdiff = fabs(Norm(x_orth) - cylradius);
      if (normdiff > 1e-10)
      {
        // Moving the node to the cylinder shell
        REAL computed_radius = Norm(x_orth);
        if (!(computed_radius > 0.) || !std::isfinite(computed_radius)) DebugStop();
        for (int ix = 0; ix < 3; ix++)
        {
          xnode[ix] = cylCenter[ix] + (x_orth[ix] / computed_radius) * cylradius + dax * cylaxis[ix];
        }
        x_orth = xnode - cylCenter;
        for (int ix = 0; ix < 3; ix++)
        {
          x_orth[ix] -= dax * cylaxis[ix];
        }
        REAL new_radius = Norm(x_orth);
        gel->NodePtr(in)->SetCoord(xnode);

        PZError << __PRETTY_FUNCTION__
                << "\nNode not on cylinder shell: " << gel->NodePtr(in)->Id()
                << "\nComputed radius: " << computed_radius
                << "\nGiven radius: " << cylradius
                << "\nElement index: " << iel
                << "\nNew coordinates: " << xnode
                << "\nNew radius: " << new_radius
                << std::endl;
      }
    }

    TPZChangeEl::ChangeToCylinder(gmesh, iel, cylCenter, cylaxis, cylradius);
  }

  gmesh->BuildConnectivity();

  nel = gmesh->NElements();
  for(int64_t iel = 0; iel < nel ; iel++) {
    TPZGeoEl* gel = gmesh->Element(iel);
    if(!gel) continue;
    if(gel->HasSubElement()) DebugStop();
    // if(gel->MaterialId() != ESurfWellCyl && gel->MaterialId() != ESurfHeel && gel->MaterialId() != ESurfToe) continue;
    if(gel->MaterialId() != matid) continue;
    int nsides = gel->NSides();
    // for (int iside = gel->NCornerNodes(); iside < nsides; iside++) {
    for (int iside = gel->FirstSide(1); iside < nsides; iside++) {
      TPZGeoElSide gelSide(gel,iside);
      TPZStack<TPZGeoElSide> allneigh;
      for(auto neigh = gelSide.Neighbour(); neigh != gelSide ; neigh++) {      
        if(neigh.Element()->IsGeoBlendEl()) {
          continue;
        }
        if (neigh.Element()->MaterialId() == matid) {
          continue;
        }
        // if(neigh.Element()->Dimension() < 2) continue;
        // if(neigh.Element()->MaterialId() != EDomain) DebugStop();
        allneigh.Push(neigh);
      }
    //   std::cout << "Element " << iel << " side " << iside << " has " << allneigh.size() << " neighbours" << std::endl;
      for(auto it : allneigh){
        TPZChangeEl::ChangeToGeoBlend(gmesh, it.Element()->Index());
      }
    }
  }  
}

void TPZWannGeometryTools::CreateCouplingEls(TPZGeoMesh* gmesh, ProblemData* SimData) {

  // Gather all ids related to wellbore cylindrical surfaces
  std::set<int64_t> wellboreCylindricalSurfaces;
  for (int i = 0; i < SimData->m_Wellbore.size(); i++) {
    auto &WellboreData = SimData->m_Wellbore[i];
    wellboreCylindricalSurfaces.insert(WellboreData.matidSurf);
  }

  const int nel = gmesh->NElements();
  for (int64_t iel = 0; iel < nel; iel++) {
    TPZGeoEl* gel = gmesh->Element(iel);
    if (wellboreCylindricalSurfaces.find(gel->MaterialId()) == wellboreCylindricalSurfaces.end()) continue;
    TPZGeoElBC bc(gel,gel->NSides()-1,SimData->EPressure2DSkin); 
    TPZGeoElBC bc2(gel,gel->NSides()-1,SimData->EPressureInterface);
    TPZGeoElBC bc3(gel,gel->NSides()-1,SimData->EHDivBoundInterface);
  }
}

void TPZWannGeometryTools::OrderWellIds(TPZGeoMesh* gmesh, ProblemData* SimData) {
  for (int i = 0; i < SimData->m_Wellbore.size(); i++) {
    const int matid = SimData->m_Wellbore[i].matidSurf;
    // TODO: we are still picking bc matids by names. Should we change that?
    const int matidPointHeel = SimData->m_Wellbore[i].BCs["point_heel"].matid;
    const int matidPointToe = SimData->m_Wellbore[i].BCs["point_toe"].matid;
    if (matidPointHeel == -1 || matidPointToe == -1) {
      std::cout << "Error: matidPointHeel or matidPointToe not found for wellbore " << i << std::endl;
      DebugStop();
    }
    TPZManVector<REAL,3> heelPoint = {0., 0., 0.};
    TPZManVector<REAL,3> toePoint = {0., 0., 0.};

    for (int64_t iel = 0; iel < gmesh->NElements(); iel++) {
      TPZGeoEl* gel = gmesh->Element(iel);
      if (!gel) continue;
      if (gel->Dimension() != 0) continue;
      if (gel->MaterialId() == matidPointHeel) {
        gel->NodePtr(0)->GetCoordinates(heelPoint);
      } else if (gel->MaterialId() == matidPointToe) {
        gel->NodePtr(0)->GetCoordinates(toePoint);
      }
    }

    const TPZManVector<REAL,3> axis = toePoint - heelPoint;
    OrderIdsSingleWell(gmesh, matid, heelPoint, axis);
  }
}

REAL TPZWannGeometryTools::ComputeAxialCoordinate(const TPZManVector<REAL,3>& point, const TPZManVector<REAL,3>& axisPoint, const TPZManVector<REAL,3>& axis) {
  if (point.size() != 3 || axisPoint.size() != 3 || axis.size() != 3) {
    DebugStop();
  }
  const REAL axisNorm = Norm(axis);
  if (!(axisNorm > 0.) || !std::isfinite(axisNorm)) {
    DebugStop();
  }
  REAL axialCoord = 0.;
  for (int d = 0; d < 3; d++) {
    axialCoord += (point[d] - axisPoint[d]) * (axis[d] / axisNorm);
  }
  return axialCoord;
}

void TPZWannGeometryTools::OrderIdsSingleWell(TPZGeoMesh* gmesh, int wellid, const TPZManVector<REAL,3>& axisPoint, const TPZManVector<REAL,3>& axis) {
  const int nel = gmesh->NElements();
  const REAL tol = 1.e-4;
  std::set<int64_t> pressure2Dels;
  std::set<REAL> nodeCoordsAxial;
  for (int64_t iel = 0; iel < nel; iel++) {
    TPZGeoEl* gel = gmesh->Element(iel);
    if (!gel) continue;
    if (gel->MaterialId() != wellid) continue;  
    if (gel->HasSubElement()) continue;
    pressure2Dels.insert(iel);
    for (int i = 0; i < gel->NCornerNodes(); i++) {
      TPZManVector<REAL,3> coor(3);
      gel->NodePtr(i)->GetCoordinates(coor);
      REAL axialCoord = ComputeAxialCoordinate(coor, axisPoint, axis);
      InsertXCoorInSet(axialCoord, nodeCoordsAxial, tol);
    }
  }

  std::map<REAL,std::set<int64_t>> axialToNodes;
  for (auto iel : pressure2Dels) {
    TPZGeoEl* gel = gmesh->Element(iel);
    if (gel->MaterialId() != wellid) DebugStop();
    for (int i = 0; i < gel->NCornerNodes(); i++) {
      TPZManVector<REAL,3> coor(3);
      gel->NodePtr(i)->GetCoordinates(coor);
      REAL axialCoord = ComputeAxialCoordinate(coor, axisPoint, axis);
      REAL closestPoint = FindClosestX(axialCoord, nodeCoordsAxial, tol);
      axialToNodes[closestPoint].insert(gel->NodeIndex(i));
    }
  }

  for (auto& it : axialToNodes) {
    const auto& nodes = it.second;
    const int64_t nnodes = nodes.size();
    if (nnodes < 2) DebugStop();          
    for (auto& node : nodes) {
      int64_t maxid = gmesh->CreateUniqueNodeId();
      gmesh->NodeVec()[node].SetNodeId(maxid);
    }
  }
}

void TPZWannGeometryTools::InsertXCoorInSet(const REAL x, std::set<REAL>& nodeCoordsX, const REAL tol) {
  if (nodeCoordsX.size() == 0) {
    nodeCoordsX.insert(x);
    return;
  }
  auto it = std::lower_bound(nodeCoordsX.begin(), nodeCoordsX.end(), x);
  if (it == nodeCoordsX.end()) {
    REAL ref = *(--it);
    if (fabs(x - ref) > tol) {
      nodeCoordsX.insert(x);
    }
  } else {
    REAL ref1 = *it;
    if (fabs(x - ref1) <= tol) {
      return;
    }
    if (it != nodeCoordsX.begin()) {
      REAL ref2 = *(--it);
      if (fabs(x - ref2) <= tol) {
        return;
      }
    }
    nodeCoordsX.insert(x);
  }
}

REAL TPZWannGeometryTools::FindClosestX(const REAL x, const std::set<REAL>& nodeCoordsX, const REAL tol) {
  if (nodeCoordsX.empty()) DebugStop();
  REAL closestX = -1000;
  auto it = std::lower_bound(nodeCoordsX.begin(), nodeCoordsX.end(), x);
  if (it == nodeCoordsX.end()) {
    closestX = *(--it);
    if(fabs(x - closestX) >= tol) DebugStop();
    return closestX;
  }
  REAL ref = *it;
  if (fabs(x - ref) <= tol) {
    closestX = ref;
    return closestX;
  }
  if (it != nodeCoordsX.begin()) {
    ref = *(--it);
    if (fabs(x - ref) <= tol) {
      closestX = ref;
      return closestX;
    }
  }
  DebugStop();
  return -1;
}

bool TPZWannGeometryTools::CheckXInSet(const REAL x, const std::set<REAL>& nodeCoordsX, const REAL tol) {
  auto it = nodeCoordsX.lower_bound(x);
  if (it != nodeCoordsX.end() && fabs(*it - x) <= tol) {
    return true;
  }
  if (it != nodeCoordsX.begin()) {
    --it;
    if (fabs(*it - x) <= tol) {
      return true;
    }
  }
  return false;
}

void TPZWannGeometryTools::hRefinement(TPZGeoMesh* gmesh, TPZVec<int64_t>& toRefine) {
  for (int64_t i = 0; i < toRefine.size(); ++i) {
    TPZVec<TPZGeoEl *> pv;
    TPZGeoEl* gel = gmesh->Element(toRefine[i]);
    if (!gel) DebugStop();
    if (gel->HasSubElement()) continue;
    gel->Divide(pv);
  }
}

void TPZWannGeometryTools::RefineFromFile(TPZGeoMesh* og_gmesh, const std::string& filename) {
  // Open the file
  std::ifstream infile(filename);
  if (!infile) {
    std::cerr << "Error: Could not open file '" << filename << "' for reading." << std::endl;
    DebugStop();
  }

  std::string line;
  while (std::getline(infile, line)) {
    std::istringstream iss(line);
    int vecSize;
    if (!(iss >> vecSize)) {
      std::cerr << "Error: Could not read vector size from line: '" << line << "'\n";
      DebugStop();
    }
    TPZVec<int64_t> toRefine(vecSize, 0);
    for (int i = 0; i < vecSize; ++i) {
      if (!(iss >> toRefine[i])) {
        std::cerr << "Error: Not enough entries for vector of size " << vecSize << " in line: '" << line << "'\n";
        DebugStop();
      }
    }

    hRefinement(og_gmesh, toRefine);
  }
}

bool TPZWannGeometryTools::VerifyMesh(TPZGeoMesh *gmesh, ProblemData *SimData) {
  bool hasErrors = false;

  // Check neighbors of 3D elements: All 3D elements should have a neighbour on each face.
  // If the neighbour is 2D, it should have a material id corresponding to a surface of the simulaiton
  for (int iel = 0; iel < gmesh->NElements(); iel++) {
    TPZGeoEl *gel = gmesh->ElementVec()[iel];

    if (!gel) DebugStop();
    if (gel->HasSubElement()) continue;
    if (gel->Dimension() != 3) continue; // only check 3D elements first

    // All 3D elements should have the material id of the domain
    if (gel->MaterialId() != SimData->m_Reservoir.matid) {
      std::cout << "Element " << gel->Index() << " has wrong material id: " << gel->MaterialId() << std::endl;
      gel->SetMaterialId(1001);
      hasErrors = true;
    }

    // Check face neighborhood
    // A 3D element should have exacly one neighbour on each face
    int firstFace = gel->FirstSide(gel->Dimension() - 1);
    int lastFace = gel->FirstSide(gel->Dimension()) - 1;

    for (int side = firstFace; side <= lastFace; side++) {
      TPZGeoElSide gelside(gel, side);

      TPZGeoElSide neighbour = gelside.Neighbour();
      if (neighbour == gelside) {
        std::cout << "Element " << gel->Index() << " side " << side
                  << " has itself as neighbour!" << std::endl;
        gel->SetMaterialId(1002);
        hasErrors = true;
        continue;
      } else {
        int count = 0;
        while (neighbour != gelside) {
          neighbour = neighbour.Neighbour();
          count++;
        }
        if (count > 1) {
          std::cout << "Element " << gel->Index() << " side " << side
                    << " has more than one neighbour!" << std::endl;
          gel->SetMaterialId(1003);
          hasErrors = true;
          continue;
        }
      }
    }
  }

  // Check neighbors of 2D elements
  // All 2D elements should have at least one 3D neighbour on their face
  for (int iel = 0; iel < gmesh->NElements(); iel++) {
    TPZGeoEl *gel = gmesh->ElementVec()[iel];
    if (!gel) DebugStop();
    if (gel->HasSubElement()) continue;
    if (gel->Dimension() != 2) continue;

    int faceSide = gel->FirstSide(2);
    TPZGeoElSide gelside(gel, faceSide);
    TPZGeoElSide neighbour = gelside.Neighbour();
    if (neighbour == gelside) {
      std::cout << "2D Element " << gel->Index() << " has no neighbour on its face!" << std::endl;
      gel->SetMaterialId(1004);
      hasErrors = true;
    } else {
      if (neighbour.Element()->Dimension() != 3) {
        std::cout << "2D Element " << gel->Index() << " has a neighbour with wrong dimension: " << neighbour.Element()->Dimension() << std::endl;
        gel->SetMaterialId(1005);
        hasErrors = true;
      }
    }
  }

  return hasErrors;
}
