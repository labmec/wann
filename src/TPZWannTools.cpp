#include "TPZWannTools.h"

TPZGeoMesh* TPZWannTools::CreateGeoMesh(ProblemData* simData) {

  std::string file_name = simData->m_MeshName;
  TPZGeoMesh* gmesh = new TPZGeoMesh();
  {
    TPZGmshReader reader;
    TPZManVector<std::map<std::string, int>, 4> stringtoint(4);
    stringtoint[3]["volume_reservoir"] = simData->EDomain;

    stringtoint[2]["surface_wellbore_cylinder"] = simData->ESurfWellCyl;
    stringtoint[2]["surface_wellbore_heel"] = simData->ESurfHeel;
    stringtoint[2]["surface_wellbore_toe"] = simData->ESurfToe;
    stringtoint[2]["surface_farfield"] = simData->EFarField;
    stringtoint[2]["surface_cap_rock"] = simData->EFarField;
    
    stringtoint[1]["curve_wellbore"] = simData->ECurveWell;
    stringtoint[1]["curve_heel"] = simData->ECurveHeel;
    stringtoint[1]["curve_toe"] = simData->ECurveToe;

    stringtoint[0]["point_heel"] = simData->EPointHeel;
    stringtoint[0]["point_toe"] = simData->EPointToe;
    
    reader.SetDimNamePhysical(stringtoint);
    reader.GeometricGmshMesh(file_name, gmesh);

    //remember to modify the matids in ProblemData according to the map
  }

  if (simData->m_Verbose) {
    std::ofstream out("gmeshorig.vtk");
    TPZVTKGeoMesh::PrintGMeshVTK(gmesh, out);
  }

  if (simData->m_ToCylindrical) {
    ModifyGeometricMeshToCylWell(gmesh, simData);
  }

  if (simData->m_NumUniformRef) {
    TPZCheckGeom checkgeom(gmesh);
    checkgeom.UniformRefine(simData->m_NumUniformRef);

    if (simData->m_Verbose) {
      std::ofstream out("gmeshnonlin.vtk");
      TPZVTKGeoMesh::PrintGMeshVTK(gmesh, out);
    }
  }

  if (simData->m_NumDirRef) {
    gRefDBase.InitializeRefPatterns(gmesh->Dimension());
    for (int i = 0; i < simData->m_NumDirRef; i++) {
      TPZRefPatternTools::RefineDirectional(gmesh, {simData->ECurveHeel, simData->ECurveToe});
    }
    if (simData->m_Verbose) {
      std::ofstream out("gmeshnonlin_ref.vtk");
      TPZVTKGeoMesh::PrintGMeshVTK(gmesh, out);
    }
  }

  // Create GeoelBCs in same location as geoels with ESurfWell
  CreatePressure2DElsAndOrderIds(gmesh, simData);

  return gmesh;
}

void TPZWannTools::ModifyGeometricMeshToCylWell(TPZGeoMesh* gmesh, ProblemData* SimData) {
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

void TPZWannTools::CreatePressure2DElsAndOrderIds(TPZGeoMesh* gmesh, ProblemData* SimData) {
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

TPZMultiphysicsCompMesh* TPZWannTools::CreateMultiphysicsCompMesh(TPZGeoMesh* gmesh, ProblemData* SimData) {

  const int dim = gmesh->Dimension();
  
  TPZHDivApproxCreator hdivCreator(gmesh);
  hdivCreator.ProbType() = ProblemType::EDarcy;
  hdivCreator.SetDefaultOrder(SimData->m_Reservoir.pOrder);
  hdivCreator.SetShouldCondense(false);
  
  //Reservoir material
  auto ReservoirData = SimData->m_Reservoir;
  {

    TPZMixedDarcyFlow *reservoirMat = new TPZMixedDarcyFlow(SimData->EDomain, dim);
    reservoirMat->SetConstantPermeability(ReservoirData.perm);
    hdivCreator.InsertMaterialObject(reservoirMat);

    for (auto &bc : ReservoirData.BCs)
    {
      TPZFMatrix<STATE> val1(1, 1, 0.);
      TPZManVector<STATE> val2(1, 0);
      val2[0] = bc.value;
      TPZBndCondT<STATE> *BCond = reservoirMat->CreateBC(reservoirMat, bc.matid, bc.type, val1, val2);
      hdivCreator.InsertMaterialObject(BCond);
    }

    // TPZFMatrix<STATE> val1(1, 1, 0.);
    // TPZManVector<STATE> val2(1, pff);
    // TPZBndCondT<STATE> *BCond1 = matdarcy->CreateBC(matdarcy, EFarField, diri, val1, val2);
    // hdivCreator.InsertMaterialObject(BCond1);

    // val2[0] = 0.;
    // // TPZBndCondT<STATE>* BCond2 = matdarcy->CreateBC(matdarcy, ESurfWellCyl, diri, val1, val2);
    // // hdivCreator.InsertMaterialObject(BCond2);

    // TPZBndCondT<STATE> *BCond3 = matdarcy->CreateBC(matdarcy, ESurfHeel, neu, val1, val2);
    // hdivCreator.InsertMaterialObject(BCond3);

    // TPZBndCondT<STATE> *BCond4 = matdarcy->CreateBC(matdarcy, ESurfToe, neu, val1, val2);
    // hdivCreator.InsertMaterialObject(BCond4);
  }

  // Pressure skin material
  {
    TPZNullMaterialCS<STATE>* matl2proj = new TPZNullMaterialCS<STATE>(SimData->EPressure2DSkin, dim-1, 1/*nstate*/);
    hdivCreator.InsertMaterialObject(matl2proj);
  }

  // Wellbore material
  auto WellboreData = SimData->m_Wellbore;
  {
    const int dimwell = 1;
    TPZMixedDarcyFlow *wellboreMat = new TPZMixedDarcyFlow(SimData->ECurveWell, dimwell);
    wellboreMat->SetConstantPermeability(WellboreData.perm);
    hdivCreator.InsertMaterialObject(wellboreMat); // This will only be used in the creation of the multiphysics mesh since the dimension is smaller than the dimension of the geometric mesh
    
    for (auto &bc : WellboreData.BCs)
    {
      TPZFMatrix<STATE> val1(1, 1, 0.);
      TPZManVector<STATE> val2(1, 0);
      val2[0] = bc.value;
      TPZBndCondT<STATE> *BCond = wellboreMat->CreateBC(wellboreMat, bc.matid, bc.type, val1, val2);
      hdivCreator.InsertMaterialObject(BCond);
    }
    // val2[0] = 2.;
    // TPZBndCondT<STATE> *BCond5 = matdarcyWell->CreateBC(matdarcyWell, EPointHeel, diri, val1, val2);
    // hdivCreator.InsertMaterialObject(BCond5);
    // val2[0] = 0.;
    // TPZBndCondT<STATE> *BCond6 = matdarcyWell->CreateBC(matdarcyWell, EPointToe, neu, val1, val2);
    // hdivCreator.InsertMaterialObject(BCond6);
  }

  // Material for hdivbound elements in multiphysics mesh
  {
    TPZNullMaterialCS<STATE>* matnull = new TPZNullMaterialCS<STATE>(SimData->EHDivBoundInterface, dim-1, 1);
    hdivCreator.InsertMaterialObject(matnull);
  }

  int lagmultilevel = 1;
  TPZManVector<TPZCompMesh*,7> meshvec(hdivCreator.NumMeshes());
  hdivCreator.CreateAtomicMeshes(meshvec, lagmultilevel); // This method increments the lagmultilevel
  AddPressureSkinElements(meshvec[1], WellboreData.pOrder, lagmultilevel); // lagmultilevel is 2 here
  AddWellboreElements(meshvec,WellboreData.pOrder, lagmultilevel);
  EqualizePressureConnects(meshvec[1]);
  AddHDivBoundInterfaceElements(meshvec[0], ReservoirData.pOrder);
  TPZMultiphysicsCompMesh* cmesh = nullptr;
  hdivCreator.CreateMultiPhysicsMesh(meshvec, lagmultilevel, cmesh);
  
  // Add material for interface elements (has to be done after autobuild so it does not create the interface elements automatically)
  TPZLagrangeMultiplierCS<STATE> *matinterface = new TPZLagrangeMultiplierCS<STATE>(SimData->EPressureInterface, dim-1, 1);
  cmesh->InsertMaterialObject(matinterface);
  AddInterfaceElements(cmesh, SimData->EPressure2DSkin, SimData->EPressureInterface, lagmultilevel);

  return cmesh;
}





//Private methods
void TPZWannTools::InsertXCoorInSet(const REAL x, std::set<REAL>& nodeCoordsX, const REAL tol) {
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

REAL TPZWannTools::FindClosestX(const REAL x, const std::set<REAL>& nodeCoordsX, const REAL tol) {
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