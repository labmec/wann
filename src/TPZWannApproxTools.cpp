#include "TPZWannApproxTools.h"

TPZMultiphysicsCompMesh *TPZWannApproxTools::CreateMultiphysicsCompMesh(TPZGeoMesh *gmesh, ProblemData *SimData)
{

  const int dim = gmesh->Dimension();

  TPZHDivApproxCreator hdivCreator(gmesh);
  hdivCreator.ProbType() = ProblemType::EDarcy;
  hdivCreator.SetDefaultOrder(SimData->m_Reservoir.pOrder);
  hdivCreator.SetShouldCondense(false);

  // Reservoir material
  auto &ReservoirData = SimData->m_Reservoir;
  {

    TPZMixedDarcyFlow *reservoirMat = new TPZMixedDarcyFlow(SimData->EDomain, dim);
    reservoirMat->SetConstantPermeability(ReservoirData.perm);
    hdivCreator.InsertMaterialObject(reservoirMat);

    for (auto &bcpair : ReservoirData.BCs)
    {
      auto &bc = bcpair.second;
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
    TPZNullMaterialCS<STATE> *matl2proj = new TPZNullMaterialCS<STATE>(SimData->EPressure2DSkin, dim - 1, 1 /*nstate*/);
    hdivCreator.InsertMaterialObject(matl2proj);
  }

  // Wellbore material
  auto &WellboreData = SimData->m_Wellbore;
  {
    const int dimwell = 1;
    TPZMixedDarcyFlow *wellboreMat = new TPZMixedDarcyFlow(SimData->ECurveWell, dimwell);
    wellboreMat->SetConstantPermeability(WellboreData.perm);
    hdivCreator.InsertMaterialObject(wellboreMat); // This will only be used in the creation of the multiphysics mesh since the dimension is smaller than the dimension of the geometric mesh

    for (auto &bcpair : WellboreData.BCs)
    {
      auto &bc = bcpair.second;
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
    TPZNullMaterialCS<STATE> *matnull = new TPZNullMaterialCS<STATE>(SimData->EHDivBoundInterface, dim - 1, 1);
    hdivCreator.InsertMaterialObject(matnull);
  }

  int lagmultilevel = 1;
  TPZManVector<TPZCompMesh *, 7> meshvec(hdivCreator.NumMeshes());
  hdivCreator.CreateAtomicMeshes(meshvec, lagmultilevel);      // This method increments the lagmultilevel
  AddPressureSkinElements(meshvec[1], SimData, lagmultilevel); // lagmultilevel is 2 here
  AddWellboreElements(meshvec, SimData, lagmultilevel);
  EqualizePressureConnects(meshvec[1], SimData);
  AddHDivBoundInterfaceElements(meshvec[0], SimData);
  TPZMultiphysicsCompMesh *cmesh = nullptr;
  hdivCreator.CreateMultiPhysicsMesh(meshvec, lagmultilevel, cmesh);

  // Add material for interface elements (has to be done after autobuild so it does not create the interface elements automatically)
  TPZLagrangeMultiplierCS<STATE> *matinterface = new TPZLagrangeMultiplierCS<STATE>(SimData->EPressureInterface, dim - 1, 1);
  cmesh->InsertMaterialObject(matinterface);
  AddInterfaceElements(cmesh, SimData, lagmultilevel);

  return cmesh;
}

void TPZWannApproxTools::AddPressureSkinElements(TPZCompMesh *cmesh, ProblemData *SimData, const int laglevel)
{

  const int dim = cmesh->Dimension() - 1; // 2D for these pressure elements living in the boundary of the 3d well
  const int matid = SimData->EPressure2DSkin;
  TPZNullMaterial<STATE> *mat = new TPZNullMaterial<>(matid, dim);
  cmesh->SetAllCreateFunctionsContinuous();
  cmesh->ApproxSpace().CreateDisconnectedElements(true);
  cmesh->InsertMaterialObject(mat);
  auto &WellboreData = SimData->m_Wellbore;
  cmesh->SetDefaultOrder(WellboreData.pOrder);

  std::set<int> matidset = {matid};
  cmesh->AutoBuild(matidset);

  const bool isContinuousConnects = false; // the pressure skin elements should not be continuous if HDiv elements in the wellbore.
  if (isContinuousConnects)
  { // This leads to continuous elements
    std::set<REAL> nodeCoordsX;
    const int64_t nel = cmesh->NElements();
    std::set<int64_t> pressure2Dels;
    for (int64_t iel = 0; iel < nel; iel++)
    {
      TPZCompEl *cel = cmesh->Element(iel);
      if (!cel)
        continue;
      if (cel->Material()->Id() != matid)
        continue;
      TPZGeoEl *gel = cel->Reference();
      if (gel->HasSubElement())
        continue;
      if (gel->NNodes() != 4)
        DebugStop();
      TPZInterpolatedElement *intEl = dynamic_cast<TPZInterpolatedElement *>(cel);
      if (!intEl)
        DebugStop();
      // I am at a pressure 2D element at the 3d well boundary
      pressure2Dels.insert(iel);
      for (int i = 4; i < 8; i++)
      {
        TPZConnect &c = cel->Connect(i);
        REAL x0 = gel->NodePtr(i % 4)->Coord(0), x1 = gel->NodePtr((i + 1) % 4)->Coord(0);
        TPZWannGeometryTools::InsertXCoorInSet(x0, nodeCoordsX, 1.e-6);
        TPZWannGeometryTools::InsertXCoorInSet(x1, nodeCoordsX, 1.e-6);
        TPZWannGeometryTools::InsertXCoorInSet((x0 + x1) / 2., nodeCoordsX, 1.e-6);
        if (fabs(x0 - x1) < 1.e-6)
        {
          // Both nodes are on the same x coordinate
          c.SetOrder(1);
          c.SetNShape(0);
        }
      }

      // Setting the area connect to order 1
      TPZConnect &c = cel->Connect(8);
      c.SetOrder(1);
      c.SetNShape(0);

      for (int i = 0; i < 9; i++)
      {
        TPZConnect &c = cel->Connect(i);
        c.SetLagrangeMultiplier(laglevel);
      }
    }

    std::map<REAL, std::set<int64_t>> xToConnects;
    const REAL tol = 1.e-6;
    for (auto iel : pressure2Dels)
    {
      TPZCompEl *cel = cmesh->Element(iel);
      TPZGeoEl *gel = cel->Reference();
      if (gel->HasSubElement())
        DebugStop();
      for (int i = 0; i < gel->NCornerNodes(); i++)
      {
        TPZManVector<REAL, 3> coor(3);
        gel->NodePtr(i)->GetCoordinates(coor);
        REAL closestX = TPZWannGeometryTools::FindClosestX(coor[0], nodeCoordsX, tol);
        xToConnects[closestX].insert(cel->ConnectIndex(i));
      }
      for (int i = 4; i < 8; i++)
      {
        REAL x0 = gel->NodePtr(i % 4)->Coord(0), x1 = gel->NodePtr((i + 1) % 4)->Coord(0);
        if (fabs(x0 - x1) < 1.e-6)
        {
          continue;
        }
        REAL closestX = TPZWannGeometryTools::FindClosestX((x0 + x1) / 2., nodeCoordsX, tol);
        xToConnects[closestX].insert(cel->ConnectIndex(i));
      }
    }

    for (auto iel : pressure2Dels)
    {
      TPZCompEl *cel = cmesh->Element(iel);
      TPZGeoEl *gel = cel->Reference();

      for (int i = 0; i < gel->NCornerNodes(); i++)
      {
        TPZManVector<REAL, 3> coor(3);
        gel->NodePtr(i)->GetCoordinates(coor);
        REAL closestX = TPZWannGeometryTools::FindClosestX(coor[0], nodeCoordsX, tol);
        int64_t cindex = *xToConnects[closestX].begin();
        cel->SetConnectIndex(i, cindex);
      }
      for (int i = 4; i < 8; i++)
      {
        REAL x0 = gel->NodePtr(i % 4)->Coord(0), x1 = gel->NodePtr((i + 1) % 4)->Coord(0);
        if (fabs(x0 - x1) < 1.e-6)
        {
          continue;
        }
        REAL closestX = TPZWannGeometryTools::FindClosestX((x0 + x1) / 2., nodeCoordsX, tol);
        int64_t cindex = *xToConnects[closestX].begin();
        cel->SetConnectIndex(i, cindex);
      }
    }
  }
  else
  {
    const int64_t nel = cmesh->NElements();
    std::set<int64_t> pressure2Dels;
    std::map<REAL, TPZManVector<REAL, 3>> elCentToConnects;
    std::set<REAL> elCentX;
    const REAL tol = 1.e-6;
    for (int64_t iel = 0; iel < nel; iel++)
    {
      TPZCompEl *cel = cmesh->Element(iel);
      if (!cel)
        continue;
      if (cel->Material()->Id() != matid)
        continue;
      TPZGeoEl *gel = cel->Reference();
      if (gel->HasSubElement())
        DebugStop();
      if (gel->NNodes() != 4)
        DebugStop();

      TPZGeoElSide gelside(gel);
      TPZManVector<REAL, 3> cent(3);
      gelside.CenterX(cent);
      TPZWannGeometryTools::InsertXCoorInSet(cent[0], elCentX, 1.e-6);
      REAL closestX = TPZWannGeometryTools::FindClosestX(cent[0], elCentX, tol);
      if (elCentToConnects.find(closestX) == elCentToConnects.end())
      {
        TPZManVector<REAL> newconnects(3);
        newconnects[0] = cmesh->AllocateNewConnect(1, 1, WellboreData.pOrder);
        newconnects[1] = cmesh->AllocateNewConnect(1, 1, WellboreData.pOrder);
        newconnects[2] = cmesh->AllocateNewConnect(WellboreData.pOrder - 1, 1, WellboreData.pOrder);
        elCentToConnects[closestX] = newconnects;
      }

      TPZInterpolatedElement *intEl = dynamic_cast<TPZInterpolatedElement *>(cel);
      if (!intEl)
        DebugStop();
      // I am at a pressure 2D element at the 3d well boundary
      pressure2Dels.insert(iel);
      for (int i = 4; i < 8; i++)
      {
        TPZConnect &c = cel->Connect(i);
        REAL x0 = gel->NodePtr(i % 4)->Coord(0), x1 = gel->NodePtr((i + 1) % 4)->Coord(0);
        if (fabs(x0 - x1) < 1.e-6)
        {
          // Both nodes are on the same x coordinate
          c.SetOrder(1);
          c.SetNShape(0);
        }
      }

      // Setting the area connect to order 1
      TPZConnect &c = cel->Connect(8);
      c.SetOrder(1);
      c.SetNShape(0);

      for (int i = 0; i < 9; i++)
      {
        TPZConnect &c = cel->Connect(i);
        c.SetLagrangeMultiplier(laglevel);
      }
    }

    for (auto iel : pressure2Dels)
    {
      TPZCompEl *cel = cmesh->Element(iel);
      TPZGeoEl *gel = cel->Reference();
      TPZGeoElSide gelside(gel);
      TPZManVector<REAL, 3> cent(3);
      gelside.CenterX(cent);
      REAL closestX = TPZWannGeometryTools::FindClosestX(cent[0], elCentX, tol);
      auto it = elCentToConnects.find(closestX);
      TPZManVector<REAL, 3> newconnects;
      if (it != elCentToConnects.end())
      {
        newconnects = it->second;
      }
      else
      {
        DebugStop(); // closestX not found in elCentToConnects
      }

      for (int i = 0; i < gel->NCornerNodes(); i++)
      {
        TPZManVector<REAL, 3> coor(3);
        gel->NodePtr(i)->GetCoordinates(coor);
        if (fabs(coor[0] - closestX) < tol)
          DebugStop();
        if (coor[0] < closestX)
        {
          cel->SetConnectIndex(i, newconnects[0]);
        }
        else
        {
          cel->SetConnectIndex(i, newconnects[1]);
        }
      }
      for (int i = 4; i < 8; i++)
      {
        REAL x0 = gel->NodePtr(i % 4)->Coord(0), x1 = gel->NodePtr((i + 1) % 4)->Coord(0);
        if (fabs(x0 - x1) < 1.e-6)
        {
          continue;
        }
        const REAL midpoint = (x0 + x1) / 2.;
        if (fabs(midpoint - closestX) > tol)
          DebugStop(); // Has to be in the middle!
        cel->SetConnectIndex(i, newconnects[2]);
      }
    }
  }

  cmesh->ComputeNodElCon();
  cmesh->CleanUpUnconnectedNodes();
  for (int64_t i = 0; i < cmesh->NConnects(); i++)
  {
    TPZConnect &c = cmesh->ConnectVec()[i];
    if (c.NElConnected() == 0)
      continue;
    cmesh->Block().Set(c.SequenceNumber(), c.NShape() * c.NState());
  }
  cmesh->InitializeBlock();

  if (SimData->m_VerbosityLevel)
  {
    std::ofstream out("cmesh.txt");
    cmesh->Print(out);
  }
}

void TPZWannApproxTools::AddWellboreElements(TPZVec<TPZCompMesh *> &meshvec, ProblemData *SimData, const int laglevel)
{
  const int dimwell = 1;
  const int matid = SimData->ECurveWell;

  // First create the pressure elements
  TPZNullMaterial<STATE> *mat = new TPZNullMaterial<>(matid, dimwell);
  meshvec[1]->SetAllCreateFunctionsContinuous();
  meshvec[1]->ApproxSpace().CreateDisconnectedElements(true);
  meshvec[1]->InsertMaterialObject(mat);
  auto &WellboreData = SimData->m_Wellbore;
  meshvec[1]->SetDefaultOrder(WellboreData.pOrder);

  std::set<int> matidset = {matid};
  meshvec[1]->AutoBuild(matidset);

  const int64_t nel = meshvec[1]->NElements();
  for (int64_t iel = 0; iel < nel; iel++)
  {
    TPZCompEl *cel = meshvec[1]->Element(iel);
    if (!cel)
      continue;
    if (cel->Material()->Id() != matid)
      continue;
    if (cel->Reference()->HasSubElement())
      continue;
    TPZGeoEl *gel = cel->Reference();
    if (gel->NNodes() != 2)
      DebugStop();
    TPZInterpolatedElement *intEl = dynamic_cast<TPZInterpolatedElement *>(cel);
    if (!intEl)
      DebugStop();
    // I am at a wellbore element
    for (int i = 0; i < intEl->NConnects(); i++)
    {
      TPZConnect &c = cel->Connect(i);
      c.SetLagrangeMultiplier(laglevel);
    }
  }

  // Now create the flux elements in meshvec[0]
  TPZNullMaterial<STATE> *matFlux = new TPZNullMaterial<>(matid, dimwell);
  TPZNullMaterial<STATE> *matFluxBcHeel = new TPZNullMaterial<>(SimData->EPointHeel, dimwell - 1);
  TPZNullMaterial<STATE> *matFluxBcToe = new TPZNullMaterial<>(SimData->EPointToe, dimwell - 1);
  meshvec[0]->InsertMaterialObject(matFlux);
  meshvec[0]->InsertMaterialObject(matFluxBcHeel);
  meshvec[0]->InsertMaterialObject(matFluxBcToe);

  meshvec[0]->SetDimModel(dimwell);
  meshvec[0]->ApproxSpace().SetAllCreateFunctionsHDiv(dimwell);
  meshvec[0]->SetDefaultOrder(WellboreData.pOrder);

  std::set<int> matidsetFlux = {matid, SimData->EPointHeel, SimData->EPointToe};
  meshvec[0]->AutoBuild(matidsetFlux);
}

void TPZWannApproxTools::EqualizePressureConnects(TPZCompMesh *cmesh, ProblemData *SimData)
{
  cmesh->Reference()->ResetReference();
  cmesh->LoadReferences();
  TPZGeoMesh *gmesh = cmesh->Reference();
  const int dim = gmesh->Dimension();
  const int64_t nel = gmesh->NElements();
  for (auto &gel : gmesh->ElementVec())
  {
    if (!gel)
      continue;
    if (gel->MaterialId() != SimData->ECurveWell)
      continue;
    if (gel->HasSubElement())
      continue;
    TPZGeoElSide gelside(gel);
    TPZGeoElSide surfwellside = gelside.Neighbour();
    while (surfwellside.Element()->MaterialId() != SimData->EPressure2DSkin)
    {
      if (surfwellside.Element()->HasSubElement())
        DebugStop();
      if (surfwellside == gelside)
        DebugStop();
      surfwellside = surfwellside.Neighbour();
    }
    TPZCompEl *cel = gel->Reference();
    TPZCompEl *celneigh = surfwellside.Element()->Reference();

    // Set the internal connect of the wellbore as the edge of a 2d pressure skin element
    const int64_t cindex = celneigh->ConnectIndex(surfwellside.Side());
    cel->SetConnectIndex(2, cindex);

    // Set the nodal connect of the 1d wellbore elements equal to the 2d pressure skin element nodal connects
    const int64_t gindex0 = gel->NodeIndex(0), gindex1 = gel->NodeIndex(1);
    TPZGeoEl *gelsurf = surfwellside.Element();
    const int64_t nindex0 = gelsurf->SideNodeIndex(surfwellside.Side(), 0), nindex1 = gelsurf->SideNodeIndex(surfwellside.Side(), 1);
    const int sidenodelocindex0 = gelsurf->SideNodeLocIndex(surfwellside.Side(), 0), sidenodelocindex1 = gelsurf->SideNodeLocIndex(surfwellside.Side(), 1);
    if (nindex0 == gindex0 && nindex1 == gindex1)
    {
      cel->SetConnectIndex(0, celneigh->ConnectIndex(sidenodelocindex0));
      cel->SetConnectIndex(1, celneigh->ConnectIndex(sidenodelocindex1));
    }
    else if (nindex0 == gindex1 && nindex1 == gindex0)
    {
      cel->SetConnectIndex(0, celneigh->ConnectIndex(sidenodelocindex1));
      cel->SetConnectIndex(1, celneigh->ConnectIndex(sidenodelocindex0));
    }
    else
    {
      DebugStop();
    }

    std::cout << "Connect indexes of cel: ";
    for (int i = 0; i < cel->NConnects(); i++)
    {
      std::cout << cel->ConnectIndex(i) << " ";
    }
    std::cout << std::endl;
  }
  cmesh->CleanUpUnconnectedNodes();
}

void TPZWannApproxTools::AddHDivBoundInterfaceElements(TPZCompMesh *cmesh, ProblemData *SimData)
{
  cmesh->Reference()->ResetReference();
  cmesh->LoadReferences();
  const int dim = cmesh->Reference()->Dimension();
  const int matid = SimData->EHDivBoundInterface;
  TPZNullMaterial<STATE> *mat = new TPZNullMaterial<>(matid, dim - 1);
  cmesh->SetDimModel(dim);
  cmesh->SetAllCreateFunctionsHDiv();
  cmesh->InsertMaterialObject(mat);
  auto &ReservoirData = SimData->m_Reservoir;
  cmesh->SetDefaultOrder(ReservoirData.pOrder);
  std::set<int> matidset = {matid};
  cmesh->AutoBuild(matidset);
}

void TPZWannApproxTools::AddInterfaceElements(TPZMultiphysicsCompMesh *cmesh, ProblemData *SimData, const int laglevel)
{

  const int matidpressure = SimData->EPressure2DSkin;
  const int matidinterface = SimData->EPressureInterface;
  cmesh->Reference()->ResetReference();
  cmesh->LoadReferences();
  TPZGeoMesh *gmesh = cmesh->Reference();
  const int dim = gmesh->Dimension();
  const int64_t nel = gmesh->NElements();
  for (int64_t iel = 0; iel < nel; iel++)
  {
    TPZGeoEl *gel = gmesh->Element(iel);
    if (!gel)
      continue;
    if (gel->MaterialId() != matidinterface)
      continue;
    if (gel->HasSubElement())
      DebugStop();
    if (gel->Dimension() != dim - 1)
      DebugStop();
    TPZCompElSide comp_pressureSide, comp_hdivSide;
    TPZGeoElSide gelSide(gel);
    for (auto neigh = gelSide.Neighbour(); neigh != gelSide; neigh++)
    {
      if (neigh.Element()->MaterialId() == matidpressure)
      {
        comp_pressureSide = neigh.Reference();
      }
      if (neigh.Element()->MaterialId() == SimData->EHDivBoundInterface)
      {
        comp_hdivSide = neigh.Reference();
      }
    }
    if (!comp_pressureSide || !comp_hdivSide)
      DebugStop();
    TPZMultiphysicsInterfaceElement *interfaceel = new TPZMultiphysicsInterfaceElement(*cmesh, gel, comp_pressureSide, comp_hdivSide);
  }
}