#include "TPZWannApproxTools.h"
#include "TPZNonlinearWell.h"
#include "TPZNonLinearWellH1.h"
#include "TPZWannMixedDarcyNL.h"
#include "TPZWannDarcyNL.h"

TPZMultiphysicsCompMesh *TPZWannApproxTools::CreateMultiphysicsCompMesh(TPZGeoMesh *gmesh, ProblemData *SimData, TPZAnalyticSolution *exact, bool isDualProblem)
{

  const int dim = gmesh->Dimension();

  TPZHDivApproxCreator hdivCreator(gmesh);
  hdivCreator.ProbType() = ProblemType::EDarcy;
  hdivCreator.SetDefaultOrder(SimData->m_Numerics.reservoirPorder);
  hdivCreator.SetShouldCondense(false);

  TLaplaceExample1* exactsol = dynamic_cast<TLaplaceExample1 *>(exact);
  bool hasAnalyticSol = (exactsol != nullptr && exactsol->fExact != TLaplaceExample1::ENone);

  // Reservoir material
  auto &ReservoirData = SimData->m_Reservoir;
  auto &FluidData = SimData->m_Fluid;
  {
    // At this point, we only pass the absolute permeability to the material.
    // The mobility lambda is incorporated in the FastCondensedCompel.
    TPZWannMixedDarcyNL *reservoirMat = new TPZWannMixedDarcyNL(ReservoirData.matid, dim);
    TPZFMatrix<STATE> perm(3, 3, 0.);
    for (int i = 0; i < 3; i++) {
      perm(i, i) = ReservoirData.perm[i] / FluidData[0].viscosity; // TODO: assuming single phase flow for now
    }
    reservoirMat->SetConstantPermeability(perm);

    hdivCreator.InsertMaterialObject(reservoirMat);

    for (auto &bcpair : ReservoirData.BCs)
    {
      auto &bc = bcpair.second;
      TPZFMatrix<STATE> val1(1, 1, 0.);
      TPZManVector<STATE> val2(1, 0);
      val2[0] = bc.value;
      TPZBndCondT<STATE> *BCond = reservoirMat->CreateBC(reservoirMat, bc.matid, bc.type, val1, val2);
      if (hasAnalyticSol) BCond->SetForcingFunctionBC(exact->ExactSolution(), 3);
      hdivCreator.InsertMaterialObject(BCond);
    }
  }

  // Pressure skin material
  {
    TPZNullMaterialCS<STATE> *matl2proj = new TPZNullMaterialCS<STATE>(SimData->EPressure2DSkin, dim - 1, 1 /*nstate*/);
    hdivCreator.InsertMaterialObject(matl2proj);
  }

  // Wellbore material
  for (auto &WellboreData : SimData->m_Wellbore) {
    const int dimwell = 1;
    
    // TODO: assuming single phase flow for now
    TPZNonlinearWell *wellboreMat =
        new TPZNonlinearWell(WellboreData.matid, 2 * WellboreData.radius,
                             FluidData[0].viscosity, FluidData[0].density, 0.0, 0.0);
    if (hasAnalyticSol) {
      wellboreMat->SetExactSol(exact->ExactSolution(), 3);
      wellboreMat->SetForcingFunction(exact->ForceFunc(), 3);
    }

    // Insert material for the multiphysics mesh
    hdivCreator.InsertMaterialObject(wellboreMat);

    for (auto &bcpair : WellboreData.BCs) {
      auto &bc = bcpair.second;
      TPZFMatrix<STATE> val1(1, 1, 0.);
      TPZManVector<STATE> val2(1, 0);
      val2[0] = bc.value;
      TPZBndCondT<STATE> *BCond = wellboreMat->CreateBC(wellboreMat, bc.matid, bc.type, val1, val2);
      if (hasAnalyticSol) BCond->SetForcingFunctionBC(exact->ExactSolution(), 3);
      hdivCreator.InsertMaterialObject(BCond);
    }
  }

  // Material for HDivBound elements in multiphysics mesh
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

  if (SimData->m_PostProc.verbosityLevel) {
    std::ofstream out("cmesh.txt");
    cmesh->Print(out);
  }

  return cmesh;
}

TPZCompMesh *TPZWannApproxTools::CreateH1CompMesh(TPZGeoMesh *gmesh, ProblemData *SimData, TPZAnalyticSolution *exact) 
{
  const int dim = gmesh->Dimension();
  auto &ReservoirData = SimData->m_Reservoir;
  auto &WellboreData = SimData->m_Wellbore[0]; // TODO: generalize for multiple wellbores
  auto &FluidData = SimData->m_Fluid;
  std::set<int> reservoirMatIdSet;
  std::set<int> wellboreMatIdSet;
  int wellPorder = SimData->m_Numerics.wellPorder;

  TLaplaceExample1* exactsol = dynamic_cast<TLaplaceExample1 *>(exact);
  bool hasAnalyticSol = (exactsol != nullptr && exactsol->fExact != TLaplaceExample1::ENone);

  TPZCompMesh *cmesh = new TPZCompMesh(gmesh);
  cmesh->SetDimModel(dim);
  cmesh->SetDefaultOrder(SimData->m_Numerics.reservoirPorder);
  cmesh->SetAllCreateFunctionsContinuous();

  // Pressure skin material (as null material)
  {
    TPZNullMaterial<STATE> *mat = new TPZNullMaterial<>(SimData->EPressure2DSkin, dim - 1);
    cmesh->InsertMaterialObject(mat);
  }

  // Reservoir material and 3D boundary conditions
  {
    TPZWannDarcyNL *reservoirMat = new TPZWannDarcyNL(SimData->m_Reservoir.matid, dim);
    TPZFMatrix<STATE> perm(3, 3, 0.);
    for (int i = 0; i < 3; i++) {
      perm(i, i) = ReservoirData.perm[i] / FluidData[0].viscosity; // TODO: Assuming single phase flow for now
    }
    reservoirMat->SetConstantPermeability(perm);
    cmesh->InsertMaterialObject(reservoirMat);

    for (auto &bcpair : ReservoirData.BCs)
    {
      auto &bc = bcpair.second;
      reservoirMatIdSet.insert(bc.matid);
      TPZFMatrix<STATE> val1(1, 1, 0.);
      TPZManVector<STATE> val2(1, 0);
      val2[0] = bc.value;
      TPZBndCondT<STATE> *BCond = reservoirMat->CreateBC(reservoirMat, bc.matid, bc.type, val1, val2);
      cmesh->InsertMaterialObject(BCond);
    }
  }

  // Wellbore material and 1D boundary conditions
  {
    // TODO: generalize for multiple wellbores; assuming single phase flow for now
    TPZNonLinearWellH1 *wellboreMat = new TPZNonLinearWellH1(SimData->m_Wellbore[0].matid, 2 * WellboreData.radius,
                             FluidData[0].viscosity, FluidData[0].density, 0.0, 0.0);
    if (hasAnalyticSol) {
      wellboreMat->SetExactSol(exact->ExactSolution(), 3);
      wellboreMat->SetForcingFunction(exact->ForceFunc(), 3);
    }

    cmesh->InsertMaterialObject(wellboreMat);

    for (auto &bcpair : WellboreData.BCs)
    {
      auto &bc = bcpair.second;
      wellboreMatIdSet.insert(bc.matid);
      TPZFMatrix<STATE> val1(1, 1, 0.);
      TPZManVector<STATE> val2(1, 0);
      val2[0] = bc.value;
      TPZBndCondT<STATE> *BCond = wellboreMat->CreateBC(wellboreMat, bc.matid, bc.type, val1, val2);
      cmesh->InsertMaterialObject(BCond);
    }
  }

  cmesh->AutoBuild();

  // Ajust well polynomial order
  // Warning: block sizes will only be ajusted after colapsing the connects in EqualizeH1Connects
  for (int64_t iel = 0; iel < cmesh->NElements(); iel++) {
    TPZCompEl *cel = cmesh->Element(iel);
    if (!cel) continue;
    if (cel->Material()->Id() != SimData->m_Wellbore[0].matid) continue; // TODO: generalize for multiple wellbores
    if (cel->NConnects() > 3) DebugStop(); // Wellbore H1 elements should have only 3 connects
    TPZConnect &c = cel->Connect(2);
    c.SetOrder(wellPorder); // TODO: generalize for multiple wellbores
    c.SetNShape(wellPorder-1);
  }

  if (SimData->m_PostProc.verbosityLevel)
  {
    std::ofstream out("cmeshH1Init.txt");
    cmesh->Print(out);
  }

  // Set dependencies for surface well connects and ajust block sizes
  EqualizeH1Connects(cmesh, SimData);

  if (SimData->m_PostProc.verbosityLevel)
  {
    std::ofstream out("cmeshH1.txt");
    cmesh->Print(out);
  }

  return cmesh;
}

void TPZWannApproxTools::AddPressureSkinElements(TPZCompMesh *cmesh, ProblemData *SimData, const int laglevel)
{
  const int dim = cmesh->Dimension() - 1; // 2D for these pressure elements living in the boundary of the 3d well
  const int matid = SimData->EPressure2DSkin;
  int wellPorder = SimData->m_Numerics.wellPorder;
  int reservoirPorder = SimData->m_Numerics.reservoirPorder;
  TPZNullMaterial<STATE> *mat = new TPZNullMaterial<>(matid, dim);
  cmesh->SetAllCreateFunctionsContinuous();
  cmesh->ApproxSpace().CreateDisconnectedElements(true);
  cmesh->InsertMaterialObject(mat);
  cmesh->SetDefaultOrder(1);
  cmesh->AutoBuild(std::set<int>{matid});

  TPZGeoMesh *gmesh = cmesh->Reference();
  const size_t nWells = SimData->m_Wellbore.size();
  std::map<int, size_t> surfaceToWell;
  std::vector<std::vector<int64_t>> skinElements(nWells);
  std::vector<TPZManVector<REAL, 3>> heels(nWells), toes(nWells);
  std::vector<int> heelCount(nWells, 0), toeCount(nWells, 0);
  for (size_t iw = 0; iw < nWells; iw++) {
    const auto& well = SimData->m_Wellbore[iw];
    heels[iw].Resize(3);
    toes[iw].Resize(3);
    if (wellPorder < 1 || !surfaceToWell.emplace(well.matidSurf, iw).second) DebugStop();
    for (auto gel : gmesh->ElementVec()) {
      if (!gel || gel->HasSubElement() || gel->Dimension() != 0) continue;
      if (gel->MaterialId() == well.BCs.at("point_heel").matid) {
        gel->NodePtr(0)->GetCoordinates(heels[iw]);
        heelCount[iw]++;
      }
      if (gel->MaterialId() == well.BCs.at("point_toe").matid) {
        gel->NodePtr(0)->GetCoordinates(toes[iw]);
        toeCount[iw]++;
      }
    }
    if (heelCount[iw] != 1 || toeCount[iw] != 1) DebugStop();
  }

  // Skin and cylinder elements occupy the same geometric side. Use that
  // neighbor ring to identify ownership without assigning per-well skin IDs.
  for (int64_t iel = 0; iel < cmesh->NElements(); iel++) {
    TPZCompEl *cel = cmesh->Element(iel);
    if (!cel || cel->Material()->Id() != matid) continue;
    TPZGeoEl *gel = cel->Reference();
    if (gel->HasSubElement() || gel->NCornerNodes() != 4) DebugStop();
    TPZGeoElSide side(gel);
    TPZGeoElSide neighbor = side.Neighbour();
    std::set<size_t> owners;
    while (neighbor != side) {
      auto it = surfaceToWell.find(neighbor.Element()->MaterialId());
      if (it != surfaceToWell.end()) owners.insert(it->second);
      neighbor = neighbor.Neighbour();
    }
    if (owners.size() != 1) DebugStop();
    skinElements[*owners.begin()].push_back(iel);
  }

  const REAL tol = 1.e-6;
  for (size_t iw = 0; iw < nWells; iw++) {
    const auto& well = SimData->m_Wellbore[iw];
    if (skinElements[iw].empty()) DebugStop();
    const TPZManVector<REAL, 3> axis = toes[iw] - heels[iw];
    auto axialCoordinate = [&](const TPZManVector<REAL, 3>& point) {
      return TPZWannGeometryTools::ComputeAxialCoordinate(point, heels[iw], axis);
    };
    auto axialNode = [&](TPZGeoEl *gel, int node) {
      TPZManVector<REAL, 3> point(3);
      gel->NodePtr(node)->GetCoordinates(point);
      return axialCoordinate(point);
    };

    // Sharing is circumferential within one axial segment of one well only.
    std::set<REAL> axialCenters;
    std::map<REAL, TPZManVector<int64_t, 3>> centerToConnects;
    for (auto iel : skinElements[iw]) {
      TPZCompEl *cel = cmesh->Element(iel);
      auto *intEl = dynamic_cast<TPZInterpolatedElement *>(cel);
      if (!intEl) DebugStop();
      intEl->PRefine(wellPorder);
      TPZGeoEl *gel = cel->Reference();
      TPZManVector<REAL, 3> center(3);
      TPZGeoElSide(gel).CenterX(center);
      const REAL axial = axialCoordinate(center);
      TPZWannGeometryTools::InsertXCoorInSet(axial, axialCenters, tol);
      const REAL closest = TPZWannGeometryTools::FindClosestX(axial, axialCenters, tol);
      if (centerToConnects.find(closest) == centerToConnects.end()) {
        TPZManVector<int64_t, 3> connects(3);
        connects[0] = cmesh->AllocateNewConnect(1, 1, wellPorder);
        connects[1] = cmesh->AllocateNewConnect(1, 1, wellPorder);
        connects[2] = cmesh->AllocateNewConnect(wellPorder - 1, 1, wellPorder);
        centerToConnects.emplace(closest, connects);
      }
      for (int i = 4; i < 8; i++) {
        if (fabs(axialNode(gel, i % 4) - axialNode(gel, (i + 1) % 4)) < tol) {
          TPZConnect &connect = cel->Connect(i);
          connect.SetOrder(1);
          connect.SetNShape(0);
        }
      }
      cel->Connect(8).SetOrder(1);
      cel->Connect(8).SetNShape(0);
    }

    for (auto iel : skinElements[iw]) {
      TPZCompEl *cel = cmesh->Element(iel);
      TPZGeoEl *gel = cel->Reference();
      TPZManVector<REAL, 3> center(3);
      TPZGeoElSide(gel).CenterX(center);
      const REAL closest = TPZWannGeometryTools::FindClosestX(axialCoordinate(center), axialCenters, tol);
      const auto& connects = centerToConnects.at(closest);
      for (int i = 0; i < gel->NCornerNodes(); i++) {
        const REAL axial = axialNode(gel, i);
        if (fabs(axial - closest) < tol) DebugStop();
        cel->SetConnectIndex(i, connects[axial < closest ? 0 : 1]);
      }
      for (int i = 4; i < 8; i++) {
        const REAL a0 = axialNode(gel, i % 4), a1 = axialNode(gel, (i + 1) % 4);
        if (fabs(a0 - a1) < tol) continue;
        if (fabs((a0 + a1)/2. - closest) > tol) DebugStop();
        cel->SetConnectIndex(i, connects[2]);
      }
      for (int i = 0; i < cel->NConnects(); i++) {
        cel->Connect(i).SetLagrangeMultiplier(laglevel);
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

  if (SimData->m_PostProc.verbosityLevel)
  {
    std::ofstream out("cmesh.txt");
    cmesh->Print(out);
  }
}

void TPZWannApproxTools::EqualizeH1Connects(TPZCompMesh *cmesh, ProblemData *SimData) {
  std::set<REAL> nodeCoordsX;
  std::map<REAL,std::set<int64_t>> xToNodes;
  std::set<int64_t> pressure2Dels;
  REAL tol = 1e-6; // Tolerance for rounding
  int wellPorder = SimData->m_Numerics.wellPorder;

  // Ensure that references are updated
  cmesh->Reference()->ResetReference();
  cmesh->LoadReferences();

  const int dim = cmesh->Dimension();

  // Get x-coordinates of pressure skin nodes (fill pressure2Dels and nodeCoordsX)
  for (auto cel : cmesh->ElementVec()) {
    if (!cel) DebugStop();
    if (cel->Material()->Id() != SimData->EPressure2DSkin) continue;

    TPZGeoEl *gel = cel->Reference();
    if (!gel || gel->Dimension() != dim-1) DebugStop();
    if (gel->HasSubElement()) continue;

    pressure2Dels.insert(cel->Index());

    for (int inodes = 0; inodes < gel->NCornerNodes(); ++inodes) {
      // Get x-coord of corner connects
      REAL x0 = gel->NodePtr(inodes)->Coord(0);
      TPZWannGeometryTools::InsertXCoorInSet(x0, nodeCoordsX, tol);

      // Get x-coord of edge connects
      REAL x1 = gel->NodePtr((inodes + 1) % 4)->Coord(0);
      TPZWannGeometryTools::InsertXCoorInSet((x0 + x1) / 2., nodeCoordsX, tol);
    }
  }

  // Group connects by x-coordinate (fill xToNodes)
  for (int64_t iel : pressure2Dels) {
    TPZCompEl *cel = cmesh->Element(iel);
    if (!cel || cel->Material()->Id() != SimData->EPressure2DSkin) DebugStop();

    TPZGeoEl *gel = cel->Reference();
    if (!gel || gel->Dimension() != dim-1 || gel->HasSubElement()) DebugStop();

    // Corner connects
    for (int inodes = 0; inodes < gel->NCornerNodes(); inodes++) {
      REAL x0 = gel->NodePtr(inodes)->Coord(0);
      REAL closestX = TPZWannGeometryTools::FindClosestX(x0, nodeCoordsX, tol);
      xToNodes[closestX].insert(cel->ConnectIndex(inodes));
    }

    // Edge connects
    for (int inodes = 4; inodes < 8; inodes++) {
      REAL x0 = gel->NodePtr(inodes % 4)->Coord(0);
      REAL x1 = gel->NodePtr((inodes + 1) % 4)->Coord(0);
      TPZConnect &c = cel->Connect(inodes);
      if (fabs(x0 - x1) < tol) {
        // Remove edge connect if x-coordinates are equal
        c.SetOrder(1);
        c.SetNShape(0);
      } else {
        REAL closestX = TPZWannGeometryTools::FindClosestX((x0 + x1) / 2., nodeCoordsX, tol);
        xToNodes[closestX].insert(cel->ConnectIndex(inodes));
        c.SetOrder(wellPorder);
        c.SetNShape(wellPorder - 1); 
      }
    }

    // Remove face connect
    TPZConnect &c = cel->Connect(8);
    c.SetOrder(1);
    c.SetNShape(0);
  }

  // Print xToNodes and nodeCoordsX in a file
  if (SimData->m_PostProc.verbosityLevel) {
    std::ofstream out("xToNodes.txt");
    for (const auto& pair : xToNodes) {
      out << "x: " << pair.first << " -> nodes: ";
      for (const auto& node : pair.second) {
        out << node << " ";
      }
      out << "\n";
    }
    out.close();

    std::ofstream out2("nodeCoordsX.txt");
    for (const auto& coord : nodeCoordsX) {
      out2 << coord << "\n";
    }
    out2.close();
  }

  // Set dependency for connects in the same x-coordinate
  for (auto pair : xToNodes) {
    const REAL x = pair.first;
    const std::set<int64_t>& nodes = pair.second;

    auto it = *nodes.begin();

    // No need to set dependency for connects with no shape functions
    if (cmesh->ConnectVec()[it].NShape() == 0) continue;

    for (auto it2 : nodes) {
      if (it2 == it) continue; // Do not add dependency to itself
      TPZConnect &c = cmesh->ConnectVec()[it2];
      int nshape = c.NShape();

      TPZFNMatrix<1,STATE> val(nshape,nshape,0.);
      for (int i = 0; i < nshape; i++) {
        val(i,i) = 1.;
      }

      if (c.HasDependency()) {
        c.RemoveDepend();
      }
      c.AddDependency(it2 ,it,val,0,0,nshape,nshape);
    }
  }

  cmesh->ExpandSolution();
  cmesh->ComputeNodElCon();
  cmesh->CleanUpUnconnectedNodes();

  for (int64_t i = 0; i < cmesh->NConnects(); i++) {
    TPZConnect &c = cmesh->ConnectVec()[i];
    if (c.NElConnected() == 0)
      continue;
    cmesh->Block().Set(c.SequenceNumber(), c.NShape() * c.NState());
  }
  cmesh->InitializeBlock();
}

void TPZWannApproxTools::AddWellboreElements(TPZVec<TPZCompMesh *> &meshvec, ProblemData *SimData, const int laglevel) {
  const int dimwell = 1;
  const int wellPorder = SimData->m_Numerics.wellPorder;

  for (auto &WellboreData : SimData->m_Wellbore) {
    const int matid = WellboreData.matid;

    // TODO: I'm not sure if picking by names is the best approach
    const int matidHeel = WellboreData.BCs.at("point_heel").matid;
    const int matidToe = WellboreData.BCs.at("point_toe").matid;

    if (matidHeel == -1 || matidToe == -1) {
      std::cout << "Error finding point heel/toe material ids in AddWellboreElements."<< std::endl;
      DebugStop();
    }

    // First create the pressure elements
    TPZNullMaterial<STATE> *mat = new TPZNullMaterial<>(matid, dimwell);
    meshvec[1]->SetAllCreateFunctionsContinuous();
    meshvec[1]->ApproxSpace().CreateDisconnectedElements(true);
    meshvec[1]->InsertMaterialObject(mat);
    meshvec[1]->SetDefaultOrder(wellPorder);

    std::set<int> matidset = {matid};
    meshvec[1]->AutoBuild(matidset);

    // Set the lagrange level for the pressure elements of the wellbore
    const int64_t nel = meshvec[1]->NElements();
    for (int64_t iel = 0; iel < nel; iel++) {
      TPZCompEl *cel = meshvec[1]->Element(iel);
      if (!cel) continue;
      if (cel->Material()->Id() != matid) continue;
      if (cel->Reference()->HasSubElement()) DebugStop();
      TPZGeoEl *gel = cel->Reference();
      if (gel->NNodes() != 2) DebugStop();
      TPZInterpolatedElement *intEl = dynamic_cast<TPZInterpolatedElement *>(cel);
      if (!intEl) DebugStop();
      // I am at a wellbore element
      for (int i = 0; i < intEl->NConnects(); i++) {
        TPZConnect &c = cel->Connect(i);
        c.SetLagrangeMultiplier(laglevel);
      }
    }

    // Now create the flux elements in meshvec[0]
    TPZNullMaterial<STATE> *matFlux = new TPZNullMaterial<>(matid, dimwell);
    TPZNullMaterial<STATE> *matFluxBcHeel = new TPZNullMaterial<>(matidHeel, dimwell - 1);
    TPZNullMaterial<STATE> *matFluxBcToe = new TPZNullMaterial<>(matidToe, dimwell - 1);
    meshvec[0]->InsertMaterialObject(matFlux);
    meshvec[0]->InsertMaterialObject(matFluxBcHeel);
    meshvec[0]->InsertMaterialObject(matFluxBcToe);

    meshvec[0]->SetDimModel(dimwell);
    meshvec[0]->ApproxSpace().SetAllCreateFunctionsHDiv(dimwell);
    meshvec[0]->SetDefaultOrder(wellPorder);

    std::set<int> matidsetFlux = {matid, matidHeel, matidToe};
    meshvec[0]->AutoBuild(matidsetFlux);
  }
}

void TPZWannApproxTools::EqualizePressureConnects(TPZCompMesh *cmesh, ProblemData *SimData)
{
  cmesh->Reference()->ResetReference();
  cmesh->LoadReferences();
  TPZGeoMesh *gmesh = cmesh->Reference();
  const int dim = gmesh->Dimension();

  // Getter all the curve wellbore matids
  std::set<int> wellboreMatIds;
  for (auto &WellboreData : SimData->m_Wellbore) {
    wellboreMatIds.insert(WellboreData.matid);
  }

  const int64_t nel = gmesh->NElements();
  for (auto &gel : gmesh->ElementVec()) {
    if (!gel) continue;
    if (wellboreMatIds.find(gel->MaterialId()) == wellboreMatIds.end()) continue;
    if (gel->HasSubElement()) continue;
    TPZGeoElSide gelside(gel);
    TPZGeoElSide surfwellside = gelside.Neighbour();
    while (surfwellside.Element()->MaterialId() != SimData->EPressure2DSkin) {
      if (surfwellside == gelside) DebugStop();
      if (surfwellside.Element()->HasSubElement()) DebugStop();
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
  }
  cmesh->CleanUpUnconnectedNodes();

  if (SimData->m_PostProc.verbosityLevel)
  {
    std::ofstream out("cmesh.txt");
    cmesh->Print(out);
  }
}

void TPZWannApproxTools::AddHDivBoundInterfaceElements(TPZCompMesh *cmesh, ProblemData *SimData) {
  cmesh->Reference()->ResetReference();
  cmesh->LoadReferences();
  const int dim = cmesh->Reference()->Dimension();
  const int matid = SimData->EHDivBoundInterface;
  TPZNullMaterial<STATE> *mat = new TPZNullMaterial<>(matid, dim - 1);
  cmesh->SetDimModel(dim);
  cmesh->SetAllCreateFunctionsHDiv();
  cmesh->InsertMaterialObject(mat);
  auto &ReservoirData = SimData->m_Reservoir;
  cmesh->SetDefaultOrder(SimData->m_Numerics.reservoirPorder);
  std::set<int> matidset = {matid};
  cmesh->AutoBuild(matidset);
}

void TPZWannApproxTools::AddInterfaceElements(TPZMultiphysicsCompMesh *cmesh, ProblemData *SimData, const int laglevel) {
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
      continue;
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