#ifdef HAVE_CONFIG_H
#include <pz_config.h>
#endif

#include <DarcyFlow/TPZDarcyFlow.h>
#include <DarcyFlow/TPZMixedDarcyFlow.h>
#include <TPZGmshReader.h>
#include <TPZLinearAnalysis.h>
#include <TPZNullMaterial.h>
#include <TPZRefPatternTools.h>
#include <TPZSSpStructMatrix.h>  //symmetric sparse matrix storage
#include <TPZSimpleTimer.h>
#include <pzbuildmultiphysicsmesh.h>
#include <pzskylstrmatrix.h>
#include <pzstepsolver.h>

#include <iostream>

#include "TPZAnalyticSolution.h"
#include "TPZCompElH1.h"
#include "TPZElementMatrixT.h"
#include "TPZGenGrid2D.h"
#include "TPZGeoMeshTools.h"
#include "TPZRefPatternDataBase.h"
#include "TPZRefPatternTools.h"
#include "TPZSYSMPMatrix.h"
#include "TPZVTKGenerator.h"
#include "TPZVTKGeoMesh.h"
#include "pzcmesh.h"
#include "pzcompelwithmem.h"
#include "pzgmesh.h"
#include "pzlog.h"
#include "pzshapequad.h"
#include "pzvisualmatrix.h"
#include "tpzchangeel.h"
#include "pzvec_extras.h"
#include "TPZCylinderMap.h"
#include "tpzgeoelrefpattern.h"
#include "tpzchangeel.h"
#include "TPZHdivApproxCreator.h"
#include "TPZLinearAnalysis.h"
#include "TPZVTKGenerator.h"
#include "TPZMultiphysicsCompMesh.h"
#include "Projection/TPZL2ProjectionCS.h"

using namespace std;

enum EMatid { ENone,
              EDomain,              
              EFarField,
              ESurfWellCyl,
              ESurfHeel,
              ESurfToe,
              ECurveWell,
              ECurveHeel,
              ECurveToe,
              ESurfWellCylNonLin,
              ECapRock,
              EPressure2DSkin,
              EPressureInterface,
              EPointHeel,
              EPointToe };

const int global_nthread = 8;

TPZGeoMesh* ReadMeshFromGmsh(std::string file_name);
void ModifyGeometricMeshToCylWell(TPZGeoMesh* gmesh);
void CreatePressure2DElsAndOrderIds(TPZGeoMesh* gmesh);
void AddPressureSkinElements(TPZCompMesh* cmesh, const int pordWell, const int laglevel);
void AddWellboreElements(TPZVec<TPZCompMesh*>& meshvec, const int pordWell, const int laglevel);

int main() {
  std::cout << "--------- Starting simulation ---------" << std::endl;
#ifdef PZ_LOG
  TPZLogger::InitializePZLOG();
#endif

  // Mesh parameters
  const int nrefdirectional = 0;
  const int nref = 0;
  const int pOrder = 1, pOrdWell = 2;

  // Physical parameters
  const REAL reservoirPerm = 1.;
  const REAL wellPerm = reservoirPerm*10.;


  // TPZGeoMesh* gmesh = ReadMeshFromGmsh("../../geo/mesh3D_rev04.msh");
//   TPZGeoMesh* gmesh = ReadMeshFromGmsh("../../geo/mesh3D_rev05.msh");
  TPZGeoMesh* gmesh = ReadMeshFromGmsh("../../geo/mesh3D_rev06.msh");
  const int dim = gmesh->Dimension();
  {
    std::ofstream out("gmeshorig.vtk");
    TPZVTKGeoMesh::PrintGMeshVTK(gmesh, out);
  }
  
  ModifyGeometricMeshToCylWell(gmesh);

  {
    TPZCheckGeom checkgeom(gmesh);
    checkgeom.UniformRefine(nref);
    // const int nref = 3;
    // for (int i = 0; i < nref; i++) {
    //   const int64_t nel = gmesh->NElements();
    //   for (int64_t iel = 0; iel < nel; iel++) {
    //     TPZGeoEl* gel = gmesh->Element(iel);
    //     if(gel->MaterialId() != ESurfWellCyl) continue;
    //     TPZManVector<TPZGeoEl*,4> subels;
    //     gel->Divide(subels);
    //   }      
    // }
    std::ofstream out("gmeshnonlin.vtk");
    TPZVTKGeoMesh::PrintGMeshVTK(gmesh, out);
  }


  // refine directional towards heel and dedao
  if (nrefdirectional) {
    gRefDBase.InitializeRefPatterns(gmesh->Dimension());
    for (int i = 0; i < nrefdirectional; i++) {
      TPZRefPatternTools::RefineDirectional(gmesh, {ECurveHeel, ECurveToe});
    }
  }

  {
    std::ofstream out("gmeshnonlin_ref.vtk");
    TPZVTKGeoMesh::PrintGMeshVTK(gmesh, out);    
  }

  // Create GeoelBCs in same location as geoels with ESurfWell
  CreatePressure2DElsAndOrderIds(gmesh);
  
  // Create computational mesh for Darcy problem
  TPZHDivApproxCreator hdivCreator(gmesh);
  hdivCreator.ProbType() = ProblemType::EDarcy;
  // hdivCreator.HdivFamily() = hdivfam;
  // hdivCreator.IsRigidBodySpaces() = false;
  hdivCreator.SetDefaultOrder(pOrder);
  // hdivCreator.SetExtraInternalOrder(0);
  hdivCreator.SetShouldCondense(true);
  // hdivCreator.HybridType() = HybridizationType::ENone;

  // Insert Materials
  const int diri = 0, neu = 1, mixed = 2;
  TPZMixedDarcyFlow* matdarcy = new TPZMixedDarcyFlow(EDomain, dim);
  matdarcy->SetConstantPermeability(reservoirPerm);
  hdivCreator.InsertMaterialObject(matdarcy);

  TPZFMatrix<STATE> val1(1, 1, 0.);
  TPZManVector<STATE> val2(1, 1.);
  TPZBndCondT<STATE>* BCond1 = matdarcy->CreateBC(matdarcy, EFarField, diri, val1, val2);
  hdivCreator.InsertMaterialObject(BCond1);

  val2[0] = 0.;
  TPZBndCondT<STATE>* BCond2 = matdarcy->CreateBC(matdarcy, ESurfWellCyl, diri, val1, val2);
  hdivCreator.InsertMaterialObject(BCond2);

  TPZBndCondT<STATE>* BCond3 = matdarcy->CreateBC(matdarcy, ESurfHeel, neu, val1, val2);
  hdivCreator.InsertMaterialObject(BCond3);
  
  TPZBndCondT<STATE>* BCond4 = matdarcy->CreateBC(matdarcy, ESurfToe, neu, val1, val2);
  hdivCreator.InsertMaterialObject(BCond4);

  // Creating material for pressure skin
  TPZL2ProjectionCS<STATE>* matl2proj = new TPZL2ProjectionCS<STATE>(EPressure2DSkin, dim-1, 1/*nstate*/);
  hdivCreator.InsertMaterialObject(matl2proj);

  // Creating mateterial for 1D wellbore and bcs
  int dimwell = 1;
  TPZMixedDarcyFlow* matdarcyWell = new TPZMixedDarcyFlow(ECurveWell, dimwell);
  matdarcyWell->SetConstantPermeability(wellPerm);
  hdivCreator.InsertMaterialObject(matdarcyWell); // This will only be used in the creation of the multiphysics mesh since the dimension is smaller than the dimension of the geometric mesh
  val2[0] = 10.;
  TPZBndCondT<STATE>* BCond5 = matdarcyWell->CreateBC(matdarcyWell, EPointHeel, diri, val1, val2);
  hdivCreator.InsertMaterialObject(BCond5);
  val2[0] = 0.;
  TPZBndCondT<STATE>* BCond6 = matdarcyWell->CreateBC(matdarcyWell, EPointToe, neu, val1, val2);
  hdivCreator.InsertMaterialObject(BCond6);

  int lagmultilevel = 1;
  TPZManVector<TPZCompMesh*,7> meshvec(hdivCreator.NumMeshes());
  hdivCreator.CreateAtomicMeshes(meshvec, lagmultilevel); // This method increments the lagmultilevel
  AddPressureSkinElements(meshvec[1],pOrdWell,lagmultilevel); // lagmultilevel is 2 here
  AddWellboreElements(meshvec,pOrdWell,lagmultilevel);
  TPZMultiphysicsCompMesh* cmesh = nullptr;
  hdivCreator.CreateMultiPhysicsMesh(meshvec, lagmultilevel, cmesh);

  if (0) {
    std::cout << "Number of connects in multiphysics mesh: " << cmesh->NConnects() << std::endl;
    for (int i = 0; i < meshvec.size(); i++) {
      std::cout << "Number of connects in atomic mesh " << i << ": " << meshvec[i]->NConnects() << std::endl;
    }
    std::ofstream out("multiphysics_mesh.txt");
    cmesh->Print(out);
  }

  TPZLinearAnalysis analysis(cmesh);
  TPZSSpStructMatrix<STATE> skylstr(cmesh);
  analysis.SetStructuralMatrix(skylstr);

  TPZStepSolver<STATE> step;
  step.SetDirect(ELDLt);
  analysis.SetSolver(step);

  analysis.Run();

  std::string filename = "solution", filename1d = "solwell";
  TPZStack<std::string> fieldnames; 
  fieldnames.Push("Pressure");
  fieldnames.Push("Flux");
  TPZVTKGenerator vtkGen(cmesh, fieldnames, filename, 0);
  vtkGen.SetNThreads(0);
  vtkGen.Do();
  
  TPZVTKGenerator vtkGen1D(cmesh, fieldnames, filename1d, 0, 1);
  vtkGen1D.SetNThreads(0);
  vtkGen1D.Do();


  // delete cmesh;
  // delete gmesh;
  std::cout << "--------- Simulation finished ---------" << std::endl;
}

TPZGeoMesh*
ReadMeshFromGmsh(std::string file_name) {
  // read mesh from gmsh
  TPZGeoMesh* gmesh;
  gmesh = new TPZGeoMesh();
  {
    TPZGmshReader reader;
    TPZManVector<std::map<std::string, int>, 4> stringtoint(4);
    stringtoint[3]["volume_reservoir"] = EDomain;

    stringtoint[2]["surface_wellbore_cylinder"] = ESurfWellCyl;
    stringtoint[2]["surface_wellbore_heel"] = ESurfHeel;
    stringtoint[2]["surface_wellbore_toe"] = ESurfToe;
    stringtoint[2]["surface_farfield"] = EFarField;
    stringtoint[2]["surface_cap_rock"] = EFarField;
    
    stringtoint[1]["curve_wellbore"] = ECurveWell;
    stringtoint[1]["curve_heel"] = ECurveHeel;
    stringtoint[1]["curve_toe"] = ECurveToe;

    stringtoint[0]["point_heel"] = EPointHeel;
    stringtoint[0]["point_toe"] = EPointToe;
    
    reader.SetDimNamePhysical(stringtoint);
    reader.GeometricGmshMesh(file_name, gmesh);
  }

  return gmesh;
}

void ModifyGeometricMeshToCylWell(TPZGeoMesh* gmesh) {
  const REAL cylradius = 0.1;
  const TPZManVector<REAL,3> cylcenter = {0.,0.,0.}, cylaxis = {1.,0.,0.};
  int64_t nel = gmesh->NElements();
  for(int64_t iel = 0; iel < nel ; iel++) {
    TPZGeoEl* gel = gmesh->Element(iel);
    if(!gel) continue;
    if(gel->HasSubElement()) DebugStop();
    if(gel->MaterialId() != ESurfWellCyl) continue;
    
    TPZManVector<int64_t, 4> nodeindices;
    gel->GetNodeIndices(nodeindices);
    // AQUINATHAN: THIAGO WILL GENERATE THIS CYLINDER ALREADY ON GMSH
    // for(auto& it : nodeindices) {
    //   TPZManVector<REAL,3> coor(3);
    //   gmesh->NodeVec()[it].GetCoordinates(coor);
    //   REAL xcoortemp = coor[0];
    //   coor[0] = 0.;
    //   const REAL norm = Norm(coor)/cylradius;
    //   coor /= norm;
    //   coor[0] = xcoortemp;
    //   gmesh->NodeVec()[it].SetCoord(coor);
    // }
    // auto newgel = new TPZGeoElRefPattern<pzgeom::TPZCylinderMap<pzgeom::TPZGeoTriangle> >(nodeindices, ESurfWellCylNonLin, *gmesh);
    // pzgeom::TPZCylinderMap<pzgeom::TPZGeoTriangle>& cylmap = newgel->Geom();
    // cylmap.SetOrigin({0.,0.,0.});
    // cylmap.SetCylinderAxis({1.,0.,0.});
    // cylmap.Initialize(newgel);
    
    TPZChangeEl::ChangeToCylinder(gmesh, iel, cylcenter, cylaxis, cylradius);
  }

  gmesh->BuildConnectivity();

  nel = gmesh->NElements();
  for(int64_t iel = 0; iel < nel ; iel++) {
    TPZGeoEl* gel = gmesh->Element(iel);
    if(!gel) continue;
    if(gel->HasSubElement()) DebugStop();
    // if(gel->MaterialId() != ESurfWellCyl && gel->MaterialId() != ESurfHeel && gel->MaterialId() != ESurfToe) continue;
    if(gel->MaterialId() != ESurfWellCyl) continue;
    int nsides = gel->NSides();
    // for (int iside = gel->NCornerNodes(); iside < nsides; iside++) {
    for (int iside = gel->FirstSide(1); iside < nsides; iside++) {
      TPZGeoElSide gelSide(gel,iside);
      TPZStack<TPZGeoElSide> allneigh;
      for(auto neigh = gelSide.Neighbour(); neigh != gelSide ; neigh++) {      
        if(neigh.Element()->IsGeoBlendEl()) {
          continue;
        }
        if (neigh.Element()->MaterialId() == ESurfWellCyl) {
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

void InsertXCoorInSet(const REAL x, std::set<REAL>& nodeCoordsX, const REAL tol) {
  if (nodeCoordsX.size() == 0) {
    nodeCoordsX.insert(x);
    return;
  }
  auto it = std::lower_bound(nodeCoordsX.begin(), nodeCoordsX.end(), x);
  if (it == nodeCoordsX.end()) {
    REAL ref = *nodeCoordsX.begin();
    if (fabs(x - ref) > tol) {
      nodeCoordsX.insert(x);
    }
  } else {
    REAL ref1 = *it;
    if (fabs(x - ref1) <= tol) {
      return;
    }
    it++;
    if (it != nodeCoordsX.end()) {
      REAL ref2 = *it;
      if (fabs(x - ref2) <= tol) {
        return;
      }
    }
    nodeCoordsX.insert(x);
  }
}

REAL FindClosestX(const REAL x, const std::set<REAL>& nodeCoordsX, const REAL tol) {
  REAL closestX = -1000;
  auto it = std::lower_bound(nodeCoordsX.begin(), nodeCoordsX.end(), x);
  if (it == nodeCoordsX.end()) {
    closestX = *nodeCoordsX.begin();
    if(fabs(x - closestX) >= tol) DebugStop();
    return closestX;
  }
  REAL ref = *it;
  if (fabs(x - ref) <= tol) {
    closestX = ref;
    return closestX;
  }
  it++;
  ref = *it;
  if (fabs(x - ref) <= tol) {
    closestX = ref;
    return closestX;
  }
  DebugStop();
  return -1;
}

void CreatePressure2DElsAndOrderIds(TPZGeoMesh* gmesh) {
  const int nel = gmesh->NElements();
  const REAL tol = 1.e-6;
  std::set<int64_t> pressure2Dels;
  std::set<REAL> nodeCoordsX;
  for (int64_t iel = 0; iel < nel; iel++) {
    TPZGeoEl* gel = gmesh->Element(iel);
    if (gel->MaterialId() != ESurfWellCyl) continue;    
    pressure2Dels.insert(iel);
    TPZGeoElBC bc(gel,gel->NSides()-1,EPressure2DSkin); 
    TPZGeoElBC bc2(gel,gel->NSides()-1,EPressureInterface);
    for (int i = 0; i < gel->NCornerNodes(); i++) {
      TPZManVector<REAL,3> coor(3);
      gel->NodePtr(i)->GetCoordinates(coor);
      InsertXCoorInSet(coor[0], nodeCoordsX, tol);
    }
  }

  std::map<REAL,std::set<int64_t>> xToNodes;
  for (auto iel : pressure2Dels) {
    TPZGeoEl* gel = gmesh->Element(iel);
    if (gel->MaterialId() != ESurfWellCyl) DebugStop();
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

  // Print nodes for debugging
  // std::multimap<REAL,int64_t> xToNodeIds;
  // for (auto& it : xToNodes) {    
  //   for (auto& node : it.second) {
  //     const int64_t id = gmesh->NodeVec()[node].Id();
  //     const REAL x = gmesh->NodeVec()[node].Coord(0);
  //     xToNodeIds.insert({x,id});
  //   }        
  // }  
  // for (const auto& it : xToNodeIds) {
  //   std::cout << "X: " << it.first << " Node ID: " << it.second << std::endl;
  // }
}

void AddPressureSkinElements(TPZCompMesh* cmesh, const int pordWell, const int laglevel) {
  std::set<REAL> nodeCoordsX;

  const int dim = cmesh->Dimension()-1; // 2D for these pressure elements living in the boundary of the 3d well
  const int matid = EPressure2DSkin;
  TPZNullMaterial<STATE>* mat = new TPZNullMaterial<>(matid,dim);
  cmesh->SetAllCreateFunctionsContinuous();
  cmesh->ApproxSpace().CreateDisconnectedElements(true);
  cmesh->InsertMaterialObject(mat);
  cmesh->SetDefaultOrder(pordWell);

  std::set<int> matidset = {matid};
  cmesh->AutoBuild(matidset);

  const int64_t nel = cmesh->NElements();
  std::set<int64_t> pressure2Dels;
  for (int64_t iel = 0; iel < nel; iel++) {
    TPZCompEl* cel = cmesh->Element(iel);    
    if (!cel) continue;
    if (cel->Material()->Id() != matid) continue;
    TPZGeoEl* gel = cel->Reference();
    if(gel->NNodes() != 4) DebugStop();
    TPZInterpolatedElement* intEl = dynamic_cast<TPZInterpolatedElement*>(cel);
    if (!intEl) DebugStop();    
    // I am at a pressure 2D element at the 3d well boundary
    pressure2Dels.insert(iel);
    for(int i = 4 ; i < 8 ; i++){
      TPZConnect& c = cel->Connect(i);
      REAL x0 = gel->NodePtr(i%4)->Coord(0), x1 = gel->NodePtr((i+1)%4)->Coord(0);
      InsertXCoorInSet(x0, nodeCoordsX, 1.e-6);
      InsertXCoorInSet(x1, nodeCoordsX, 1.e-6);
      InsertXCoorInSet((x0+x1)/2., nodeCoordsX, 1.e-6);
      if (fabs(x0 - x1) < 1.e-6) {
        // Both nodes are on the same x coordinate
        c.SetOrder(1);
        c.SetNShape(0);
      }
    }

    // Setting the area connect to order 1
    TPZConnect& c = cel->Connect(8);
    c.SetOrder(1);
    c.SetNShape(0);

    for (int i = 0; i < 9; i++) {
      TPZConnect& c = cel->Connect(i);
      c.SetLagrangeMultiplier(laglevel);
    }
    
  }

  std::map<REAL,std::set<int64_t>> xToConnects;
  const REAL tol = 1.e-6;
  for (auto iel : pressure2Dels) {
    TPZCompEl* cel = cmesh->Element(iel);
    TPZGeoEl* gel = cel->Reference();
    for (int i = 0; i < gel->NCornerNodes(); i++) {
      TPZManVector<REAL,3> coor(3);
      gel->NodePtr(i)->GetCoordinates(coor);
      REAL closestX = FindClosestX(coor[0], nodeCoordsX, tol);
      xToConnects[closestX].insert(cel->ConnectIndex(i));
    }
    for (int i = 4; i < 8; i++) {
      REAL x0 = gel->NodePtr(i%4)->Coord(0), x1 = gel->NodePtr((i+1)%4)->Coord(0);
      if (fabs(x0 - x1) < 1.e-6) {
        continue;
      }
      REAL closestX = FindClosestX((x0+x1)/2., nodeCoordsX, tol);
      xToConnects[closestX].insert(cel->ConnectIndex(i));
    }    
  }

  for (auto iel : pressure2Dels) {
    TPZCompEl* cel = cmesh->Element(iel);
    TPZGeoEl* gel = cel->Reference();
    for (int i = 0; i < gel->NCornerNodes(); i++) {
      TPZManVector<REAL,3> coor(3);
      gel->NodePtr(i)->GetCoordinates(coor);
      REAL closestX = FindClosestX(coor[0], nodeCoordsX, tol);
      int64_t cindex = *xToConnects[closestX].begin();
      cel->SetConnectIndex(i,cindex);
    }
    for (int i = 4; i < 8; i++) {
      REAL x0 = gel->NodePtr(i%4)->Coord(0), x1 = gel->NodePtr((i+1)%4)->Coord(0);
      if (fabs(x0 - x1) < 1.e-6) {
        continue;
      }
      REAL closestX = FindClosestX((x0+x1)/2., nodeCoordsX, tol);
      int64_t cindex = *xToConnects[closestX].begin();
      cel->SetConnectIndex(i,cindex);
    }    
  }
  cmesh->ComputeNodElCon();
  cmesh->CleanUpUnconnectedNodes();  
  for (int64_t i = 0; i < cmesh->NConnects(); i++) {
    TPZConnect &c = cmesh->ConnectVec()[i];
    if (c.NElConnected() == 0) continue;
    cmesh->Block().Set(c.SequenceNumber(),c.NShape()*c.NState());
  }
  cmesh->InitializeBlock();

  {
    std::ofstream out("cmesh.txt");
    cmesh->Print(out);
  }
}

void AddWellboreElements(TPZVec<TPZCompMesh*>& meshvec, const int pordWell, const int laglevel) {
  const int dimwell = 1;
  const int matid = ECurveWell;

  // First create the pressure elements
  TPZNullMaterial<STATE>* mat = new TPZNullMaterial<>(matid,dimwell);
  meshvec[1]->SetAllCreateFunctionsContinuous();
  meshvec[1]->ApproxSpace().CreateDisconnectedElements(true);
  meshvec[1]->InsertMaterialObject(mat);
  meshvec[1]->SetDefaultOrder(pordWell);

  std::set<int> matidset = {matid};
  meshvec[1]->AutoBuild(matidset);

  const int64_t nel = meshvec[1]->NElements();
  for (int64_t iel = 0; iel < nel; iel++) {
    TPZCompEl* cel = meshvec[1]->Element(iel);    
    if (!cel) continue;
    if (cel->Material()->Id() != matid) continue;
    TPZGeoEl* gel = cel->Reference();
    if(gel->NNodes() != 2) DebugStop();
    TPZInterpolatedElement* intEl = dynamic_cast<TPZInterpolatedElement*>(cel);
    if (!intEl) DebugStop();    
    // I am at a wellbore element
    for (int i = 0; i < intEl->NConnects(); i++) {
      TPZConnect& c = cel->Connect(i);
      c.SetLagrangeMultiplier(laglevel);
    }
  }

  // Now create the flux elements in meshvec[0]
  TPZNullMaterial<STATE>* matFlux = new TPZNullMaterial<>(matid,dimwell);
  TPZNullMaterial<STATE>* matFluxBcHeel = new TPZNullMaterial<>(EPointHeel,dimwell-1);
  TPZNullMaterial<STATE>* matFluxBcToe = new TPZNullMaterial<>(EPointToe,dimwell-1);
  meshvec[0]->InsertMaterialObject(matFlux);
  meshvec[0]->InsertMaterialObject(matFluxBcHeel);
  meshvec[0]->InsertMaterialObject(matFluxBcToe);

  meshvec[0]->SetDimModel(dimwell);
  meshvec[0]->ApproxSpace().SetAllCreateFunctionsHDiv(dimwell);
  meshvec[0]->SetDefaultOrder(pordWell);
  
  std::set<int> matidsetFlux = {matid, EPointHeel, EPointToe};
  meshvec[0]->AutoBuild(matidsetFlux);

}