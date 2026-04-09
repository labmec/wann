#include "TPZWannGeometryTools.h"
#include "TPZWannAdaptivityTools.h"
#include "TPZRefPatternTools.h"

TPZGeoMesh* TPZWannGeometryTools::CreateGeoMesh(ProblemData* simData) {

  TPZGeoMesh* gmesh = ReadMeshFromGmsh(simData);
  
  // Divide pyramids if present in the OG mesh
  DividePyramids(gmesh);

  if (1) {
    std::ofstream out("gmeshorig.vtk");
    TPZVTKGeoMesh::PrintGMeshVTK(gmesh, out);
  }

  if (simData->m_Mesh.ToCylindrical) {
    ModifyGeometricMeshToCylWell(gmesh, simData);
  }

  if (simData->m_Mesh.customRefinement != 0) {
    std::string file = simData->m_Mesh.file;
    std::string baseName = file.substr(0, file.find_last_of('.'));
    std::string refProcessFile = baseName + "_refProcess.txt";
    std::string path(std::string(INPUTDIR) + "/" + refProcessFile);
    TPZWannGeometryTools::RefineFromFile(gmesh, path);
    if (simData->m_PostProc.verbosityLevel) {
      std::ofstream out("gmesh_customref.vtk");
      TPZVTKGeoMesh::PrintGMeshVTK(gmesh, out);
    }
  
  // Only perform uniform and directional refinements if no custom refinement is specified
  } else {
    if (simData->m_Mesh.NumUniformRef) {
      TPZCheckGeom checkgeom(gmesh);
      checkgeom.UniformRefine(simData->m_Mesh.NumUniformRef);

      if (simData->m_PostProc.verbosityLevel) {
        std::ofstream out("gmeshnonlin.vtk");
        TPZVTKGeoMesh::PrintGMeshVTK(gmesh, out);
      }
    }

    if (simData->m_Mesh.NumDirRef) {
      gRefDBase.InitializeRefPatterns(gmesh->Dimension());
      for (int i = 0; i < simData->m_Mesh.NumDirRef; i++) {
        TPZRefPatternTools::RefineDirectional(gmesh, {simData->ECurveHeel, simData->ECurveToe});
      }
      if (simData->m_PostProc.verbosityLevel) {
        std::ofstream out("gmeshnonlin_ref.vtk");
        TPZVTKGeoMesh::PrintGMeshVTK(gmesh, out);
      }
    }
  }

  // Create GeoelBCs in same location as geoels with ESurfWell
  CreatePressure2DEls(gmesh, simData);

  // Order nodes Id in the well according to the x-coordinate
  OrderIds(gmesh, simData);

  return gmesh;
}

TPZGeoMesh* TPZWannGeometryTools::ReadMeshFromGmsh(ProblemData* simData){

  std::string file = simData->m_Mesh.file;
  std::string path(std::string(INPUTDIR) + "/" + file);
  TPZGeoMesh* gmesh = new TPZGeoMesh();
  {
    TPZGmshReader reader;
    TPZManVector<std::map<std::string, int>, 4> stringtoint(4);
    stringtoint[3]["volume_reservoir"] = simData->EDomain;
    simData->m_Reservoir.matid = simData->EDomain;

    stringtoint[2]["surface_wellbore_cylinder"] = simData->ESurfWellCyl;
    stringtoint[2]["surface_wellbore_heel"] = simData->ESurfHeel;
    SetBC(simData, "surface_wellbore_heel", simData->ESurfHeel);
    stringtoint[2]["surface_wellbore_toe"] = simData->ESurfToe;
    SetBC(simData, "surface_wellbore_toe", simData->ESurfToe);
    stringtoint[2]["surface_farfield"] = simData->EFarField;
    SetBC(simData, "surface_farfield", simData->EFarField);
    stringtoint[2]["surface_cap_rock"] = simData->ECapRock;
    SetBC(simData, "surface_cap_rock", simData->ECapRock);
    
    stringtoint[1]["curve_wellbore"] = simData->ECurveWell;
    simData->m_Wellbore.matid = simData->ECurveWell;
    stringtoint[1]["curve_heel"] = simData->ECurveHeel;
    stringtoint[1]["curve_toe"] = simData->ECurveToe;

    stringtoint[0]["point_heel"] = simData->EPointHeel;
    SetBC(simData, "point_heel", simData->EPointHeel);
    stringtoint[0]["point_toe"] = simData->EPointToe;
    SetBC(simData, "point_toe", simData->EPointToe);
    
    reader.SetDimNamePhysical(stringtoint);
    reader.GeometricGmshMesh(path, gmesh);

    // Remove gmsh boundary elements and create GeoElBC so normals are consistent
    // int64_t nel = gmesh->NElements();
    // for(int64_t el = 0; el < nel; el++){
    //     TPZGeoEl *gel = gmesh->Element(el);
    //     if(!gel || gel->Dimension() != gmesh->Dimension()-1) continue;
    //     TPZGeoElSide gelside(gel);
    //     TPZGeoElSide neigh = gelside.Neighbour();
    //     gel->RemoveConnectivities();
    //     int matid = gel->MaterialId();
    //     delete gel;
    //     TPZGeoElBC gbc(neigh, matid);
    // }

    //remember to modify the matids in ProblemData according to the map
  }
  return gmesh;
}

bool TPZWannGeometryTools::SetBC(ProblemData* simData, const std::string& bcName, int matid) {
  auto it = simData->m_Reservoir.BCs.find(bcName);
  if (it != simData->m_Reservoir.BCs.end()) {
    it->second.matid = matid;
    return true;
  }

  it = simData->m_Wellbore.BCs.find(bcName);
  if (it != simData->m_Wellbore.BCs.end()) {
    it->second.matid = matid;
    return true;
  }

  return false;
}

void TPZWannGeometryTools::ModifyGeometricMeshToCylWell(TPZGeoMesh* gmesh, ProblemData* SimData) {
  const REAL cylradius = SimData->m_Wellbore.radius;
  const TPZManVector<REAL,3> cylcenter = {0.,0.,0.}, cylaxis = {1.,0.,0.};
  int64_t nel = gmesh->NElements();
  for(int64_t iel = 0; iel < nel ; iel++) {
    TPZGeoEl* gel = gmesh->Element(iel);
    if(!gel) continue;
    if(gel->HasSubElement()) DebugStop();
    if(gel->MaterialId() != SimData->ESurfWellCyl) continue;
    
    TPZManVector<int64_t, 4> nodeindices;
    gel->GetNodeIndices(nodeindices);
    const int nnodes = gel->NCornerNodes();
    //Moving the nodes to the cylinder surface
    TPZManVector<REAL,3> xnode(3,0);
    for (int in = 0; in < nnodes; in++)
    {
      gel->NodePtr(in)->GetCoordinates(xnode);
      // component of xnode that is orthogonal to cyl axis
      TPZManVector<REAL, 3> x_orth = xnode - cylcenter;
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
        for (int ix = 0; ix < 3; ix++)
        {
          xnode[ix] = cylcenter[ix] + (x_orth[ix] / computed_radius) * cylradius + dax * cylaxis[ix];
        }
        x_orth = xnode - cylcenter;
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

    TPZChangeEl::ChangeToCylinder(gmesh, iel, cylcenter, cylaxis, cylradius);
  }

  gmesh->BuildConnectivity();

  nel = gmesh->NElements();
  for(int64_t iel = 0; iel < nel ; iel++) {
    TPZGeoEl* gel = gmesh->Element(iel);
    if(!gel) continue;
    if(gel->HasSubElement()) DebugStop();
    // if(gel->MaterialId() != ESurfWellCyl && gel->MaterialId() != ESurfHeel && gel->MaterialId() != ESurfToe) continue;
    if(gel->MaterialId() != SimData->ESurfWellCyl) continue;
    int nsides = gel->NSides();
    // for (int iside = gel->NCornerNodes(); iside < nsides; iside++) {
    for (int iside = gel->FirstSide(1); iside < nsides; iside++) {
      TPZGeoElSide gelSide(gel,iside);
      TPZStack<TPZGeoElSide> allneigh;
      for(auto neigh = gelSide.Neighbour(); neigh != gelSide ; neigh++) {      
        if(neigh.Element()->IsGeoBlendEl()) {
          continue;
        }
        if (neigh.Element()->MaterialId() == SimData->ESurfWellCyl) {
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

void TPZWannGeometryTools::CreatePressure2DEls(TPZGeoMesh* gmesh, ProblemData* SimData) {
  const int nel = gmesh->NElements();
  for (int64_t iel = 0; iel < nel; iel++) {
    TPZGeoEl* gel = gmesh->Element(iel);
    if (gel->MaterialId() != SimData->ESurfWellCyl) continue;  
    if (gel->HasSubElement()) continue;
    // pressure2Dels.insert(iel);
    TPZGeoElBC bc(gel,gel->NSides()-1,SimData->EPressure2DSkin); 
    TPZGeoElBC bc2(gel,gel->NSides()-1,SimData->EPressureInterface);
    TPZGeoElBC bc3(gel,gel->NSides()-1,SimData->EHDivBoundInterface);
  }
}

void TPZWannGeometryTools::OrderIds(TPZGeoMesh* gmesh, ProblemData* SimData) {
  const int nel = gmesh->NElements();
  const REAL tol = 1.e-4;
  std::set<int64_t> pressure2Dels;
  std::set<REAL> nodeCoordsX;
  for (int64_t iel = 0; iel < nel; iel++) {
    TPZGeoEl* gel = gmesh->Element(iel);
    if (gel->MaterialId() != SimData->ESurfWellCyl) continue;  
    if (gel->HasSubElement()) continue;
    pressure2Dels.insert(iel);
    for (int i = 0; i < gel->NCornerNodes(); i++) {
      TPZManVector<REAL,3> coor(3);
      gel->NodePtr(i)->GetCoordinates(coor);
      InsertXCoorInSet(coor[0], nodeCoordsX, tol);
    }
  }

  std::map<REAL,std::set<int64_t>> xToNodes;
  for (auto iel : pressure2Dels) {
    TPZGeoEl* gel = gmesh->Element(iel);
    if (gel->MaterialId() != SimData->ESurfWellCyl) DebugStop();
    for (int i = 0; i < gel->NCornerNodes(); i++) {
      TPZManVector<REAL,3> coor(3);
      gel->NodePtr(i)->GetCoordinates(coor);
      REAL closestX = FindClosestX(coor[0], nodeCoordsX, tol);
      xToNodes[closestX].insert(gel->NodeIndex(i));
    }
  }

  for (auto& it : xToNodes) {
    const REAL x = it.first;
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
    REAL ref2 = *(--it);
    if (fabs(x - ref2) <= tol) {
      return;
    }
    nodeCoordsX.insert(x);
  }
}

REAL TPZWannGeometryTools::FindClosestX(const REAL x, const std::set<REAL>& nodeCoordsX, const REAL tol) {
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
  ref = *(--it);
  if (fabs(x - ref) <= tol) {
    closestX = ref;
    return closestX;
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

void TPZWannGeometryTools::hRefinement(TPZGeoMesh* gmesh, TPZVec<int>& refinementIndicator) {
  if (gmesh->NElements() != refinementIndicator.size()) {
    std::cout << "Refinement vector size " << refinementIndicator.size() 
              << " is different from number of elements in the mesh " << gmesh->NElements() << std::endl;
    DebugStop();
  }

  for (int64_t i = 0; i < refinementIndicator.size(); ++i) {
    if (refinementIndicator[i] == 0) continue;
    TPZVec<TPZGeoEl *> pv;
    TPZGeoEl* gel = gmesh->Element(i);
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
    TPZVec<int> indicatorVec(vecSize, 0);
    for (int i = 0; i < vecSize; ++i) {
      if (!(iss >> indicatorVec[i])) {
        std::cerr << "Error: Not enough entries for vector of size " << vecSize << " in line: '" << line << "'\n";
        DebugStop();
      }
    }

    hRefinement(og_gmesh, indicatorVec);
  }
}

void TPZWannGeometryTools::DividePyramids(TPZGeoMesh *gmesh) {
  gRefDBase.InitializeRefPatterns(EPiramide);
  auto refpatter = gRefDBase.FindRefPattern("PyrTwoTets");

  if (!refpatter) {
    std::cout << "Refinement pattern for pyramids not found!" << std::endl;
    DebugStop();
  }
  int64_t nelements = gmesh->NElements();
  for (int64_t el = 0; el < nelements; el++) {
    TPZGeoEl *geoel = gmesh->Element(el);
    if (!geoel)
      continue;
    if (geoel->Type() == EPiramide) {
      geoel->SetRefPattern(refpatter);
      TPZManVector<TPZGeoEl *> el(0);
      geoel->Divide(el);
    }
  }
}