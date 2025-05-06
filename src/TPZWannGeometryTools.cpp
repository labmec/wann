#include "TPZWannGeometryTools.h"

TPZGeoMesh* TPZWannGeometryTools::CreateGeoMesh(ProblemData* simData) {

  TPZGeoMesh* gmesh = ReadMeshFromGmsh(simData);

  if (simData->m_VerbosityLevel) {
    std::ofstream out("gmeshorig.vtk");
    TPZVTKGeoMesh::PrintGMeshVTK(gmesh, out);
  }

  if (simData->m_Mesh.ToCylindrical) {
    ModifyGeometricMeshToCylWell(gmesh, simData);
  }

  if (simData->m_Mesh.NumUniformRef) {
    TPZCheckGeom checkgeom(gmesh);
    checkgeom.UniformRefine(simData->m_Mesh.NumUniformRef);

    if (simData->m_VerbosityLevel) {
      std::ofstream out("gmeshnonlin.vtk");
      TPZVTKGeoMesh::PrintGMeshVTK(gmesh, out);
    }
  }

  if (simData->m_Mesh.NumDirRef) {
    gRefDBase.InitializeRefPatterns(gmesh->Dimension());
    for (int i = 0; i < simData->m_Mesh.NumDirRef; i++) {
      TPZRefPatternTools::RefineDirectional(gmesh, {simData->ECurveHeel, simData->ECurveToe});
    }
    if (simData->m_VerbosityLevel) {
      std::ofstream out("gmeshnonlin_ref.vtk");
      TPZVTKGeoMesh::PrintGMeshVTK(gmesh, out);
    }
  }

  // Create GeoelBCs in same location as geoels with ESurfWell
  CreatePressure2DElsAndOrderIds(gmesh, simData);

  return gmesh;
}

TPZGeoMesh* TPZWannGeometryTools::ReadMeshFromGmsh(ProblemData* simData){

  std::string file_name = simData->m_Mesh.file;
  TPZGeoMesh* gmesh = new TPZGeoMesh();
  {
    TPZGmshReader reader;
    TPZManVector<std::map<std::string, int>, 4> stringtoint(4);
    stringtoint[3]["volume_reservoir"] = simData->EDomain;
    simData->m_Reservoir.matid = simData->EDomain;

    stringtoint[2]["surface_wellbore_cylinder"] = simData->ESurfWellCyl;
    stringtoint[2]["surface_wellbore_heel"] = simData->ESurfHeel;
    simData->m_Reservoir.BCs["surface_wellbore_heel"].matid = simData->ESurfHeel;
    stringtoint[2]["surface_wellbore_toe"] = simData->ESurfToe;
    simData->m_Reservoir.BCs["surface_wellbore_toe"].matid = simData->ESurfToe;
    stringtoint[2]["surface_farfield"] = simData->EFarField;
    simData->m_Reservoir.BCs["surface_farfield"].matid = simData->EFarField;
    stringtoint[2]["surface_cap_rock"] = simData->EFarField;
    simData->m_Reservoir.BCs["surface_cap_rock"].matid = simData->EFarField;
    
    stringtoint[1]["curve_wellbore"] = simData->ECurveWell;
    simData->m_Wellbore.matid = simData->ECurveWell;
    stringtoint[1]["curve_heel"] = simData->ECurveHeel;
    stringtoint[1]["curve_toe"] = simData->ECurveToe;

    stringtoint[0]["point_heel"] = simData->EPointHeel;
    simData->m_Wellbore.BCs["point_heel"].matid = simData->EPointHeel;
    stringtoint[0]["point_toe"] = simData->EPointToe;
    simData->m_Wellbore.BCs["point_toe"].matid = simData->EPointToe;
    
    reader.SetDimNamePhysical(stringtoint);
    reader.GeometricGmshMesh(file_name, gmesh);

    //remember to modify the matids in ProblemData according to the map
  }
  return gmesh;
}

void TPZWannGeometryTools::ModifyGeometricMeshToCylWell(TPZGeoMesh* gmesh, ProblemData* SimData) {
  const REAL cylradius = 0.1;
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

void TPZWannGeometryTools::CreatePressure2DElsAndOrderIds(TPZGeoMesh* gmesh, ProblemData* SimData) {
  const int nel = gmesh->NElements();
  const REAL tol = 1.e-4;
  std::set<int64_t> pressure2Dels;
  std::set<REAL> nodeCoordsX;
  for (int64_t iel = 0; iel < nel; iel++) {
    TPZGeoEl* gel = gmesh->Element(iel);
    if (gel->MaterialId() != SimData->ESurfWellCyl) continue;  
    if (gel->HasSubElement()) continue;
    pressure2Dels.insert(iel);
    TPZGeoElBC bc(gel,gel->NSides()-1,SimData->EPressure2DSkin); 
    TPZGeoElBC bc2(gel,gel->NSides()-1,SimData->EPressureInterface);
    TPZGeoElBC bc3(gel,gel->NSides()-1,SimData->EHDivBoundInterface);
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

  int64_t maxid = gmesh->CreateUniqueNodeId();
  for (auto& it : xToNodes) {
    const REAL x = it.first;
    const auto& nodes = it.second;
    const int64_t nnodes = nodes.size();
    if (nnodes < 2) DebugStop();          
    for (auto& node : nodes) {
      gmesh->NodeVec()[node].SetNodeId(maxid++);
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