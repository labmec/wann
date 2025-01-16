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
#include "TPZHDivApproxCreator.h"
#include "TPZLinearAnalysis.h"
#include "TPZVTKGenerator.h"
#include "TPZMultiphysicsCompMesh.h"
#include "Projection/TPZL2ProjectionCS.h"
#include "TPZNullMaterialCS.h"
#include "TPZLagrangeMultiplierCS.h"

using namespace std;

enum EMatid { ENone,
              EDomain,              
              EFarField,
              ESurfWellCyl,
              ESurfHeel,
              ESurfToe, //5
              ECurveWell,
              ECurveHeel,
              ECurveToe,
              ESurfWellCylNonLin,
              ECapRock, // 10
              EPressure2DSkin,
              EPressureInterface,
              EPointHeel,
              EPointToe,
              EHDivBoundInterface };

const int global_nthread = 8;

TPZGeoMesh* ReadMeshFromGmsh(std::string file_name);
void ModifyGeometricMeshToCylWell(TPZGeoMesh* gmesh);
void CreatePressure2DElsAndOrderIds(TPZGeoMesh* gmesh);
void AddPressureSkinElements(TPZCompMesh* cmesh, const int pordWell, const int laglevel);
void AddWellboreElements(TPZVec<TPZCompMesh*>& meshvec, const int pordWell, const int laglevel);
void AddInterfaceElements(TPZMultiphysicsCompMesh* cmesh, const int matidpressure, const int matidinterface, const int laglevel);
void AddHDivBoundInterfaceElements(TPZCompMesh* cmesh, const int porder);
void EqualizePressureConnects(TPZCompMesh* cmesh);

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
  // TPZGeoMesh* gmesh = ReadMeshFromGmsh("../../geo/mesh3D_rev06.msh");
  // TPZGeoMesh* gmesh = ReadMeshFromGmsh("../../geo/mesh3D_rev07.msh");
  TPZGeoMesh* gmesh = ReadMeshFromGmsh("../../geo/mesh3D_rev08.msh");

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
  hdivCreator.SetShouldCondense(false);
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
  // TPZBndCondT<STATE>* BCond2 = matdarcy->CreateBC(matdarcy, ESurfWellCyl, diri, val1, val2);
  // hdivCreator.InsertMaterialObject(BCond2);

  TPZBndCondT<STATE>* BCond3 = matdarcy->CreateBC(matdarcy, ESurfHeel, neu, val1, val2);
  hdivCreator.InsertMaterialObject(BCond3);
  
  TPZBndCondT<STATE>* BCond4 = matdarcy->CreateBC(matdarcy, ESurfToe, neu, val1, val2);
  hdivCreator.InsertMaterialObject(BCond4);

  // Creating material for pressure skin
  // TPZL2ProjectionCS<STATE>* matl2proj = new TPZL2ProjectionCS<STATE>(EPressure2DSkin, dim-1, 1/*nstate*/);
  // matl2proj->SetScaleFactor(0.);
  // matl2proj->SetSol({0.});
  TPZNullMaterialCS<STATE>* matl2proj = new TPZNullMaterialCS<STATE>(EPressure2DSkin, dim-1, 1/*nstate*/);
  hdivCreator.InsertMaterialObject(matl2proj);

  // Creating mateterial for 1D wellbore and bcs
  int dimwell = 1;
  TPZMixedDarcyFlow* matdarcyWell = new TPZMixedDarcyFlow(ECurveWell, dimwell);
  matdarcyWell->SetConstantPermeability(wellPerm);
  hdivCreator.InsertMaterialObject(matdarcyWell); // This will only be used in the creation of the multiphysics mesh since the dimension is smaller than the dimension of the geometric mesh
  val2[0] = 2.;
  TPZBndCondT<STATE>* BCond5 = matdarcyWell->CreateBC(matdarcyWell, EPointHeel, diri, val1, val2);
  hdivCreator.InsertMaterialObject(BCond5);
  val2[0] = 0.;
  TPZBndCondT<STATE>* BCond6 = matdarcyWell->CreateBC(matdarcyWell, EPointToe, neu, val1, val2);
  hdivCreator.InsertMaterialObject(BCond6);

  // Create null material for hdivbound elements in multiphysics mesh
  TPZNullMaterialCS<STATE>* matnull = new TPZNullMaterialCS<STATE>(EHDivBoundInterface,2,1);
  hdivCreator.InsertMaterialObject(matnull);

  int lagmultilevel = 1;
  TPZManVector<TPZCompMesh*,7> meshvec(hdivCreator.NumMeshes());
  hdivCreator.CreateAtomicMeshes(meshvec, lagmultilevel); // This method increments the lagmultilevel
  AddPressureSkinElements(meshvec[1],pOrdWell,lagmultilevel); // lagmultilevel is 2 here
  AddWellboreElements(meshvec,pOrdWell,lagmultilevel);
  EqualizePressureConnects(meshvec[1]);
  AddHDivBoundInterfaceElements(meshvec[0], pOrder);
  TPZMultiphysicsCompMesh* cmesh = nullptr;
  hdivCreator.CreateMultiPhysicsMesh(meshvec, lagmultilevel, cmesh);
  
  // Add material for interface elements (has to be done after autobuild so it does not create the interface elements automatically)
  TPZLagrangeMultiplierCS<STATE> *matinterface = new TPZLagrangeMultiplierCS<STATE>(EPressureInterface, dim-1, 1);
  cmesh->InsertMaterialObject(matinterface);
  AddInterfaceElements(cmesh, EPressure2DSkin, EPressureInterface, lagmultilevel);



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
  skylstr.SetNumThreads(global_nthread);
  analysis.SetStructuralMatrix(skylstr);

  TPZStepSolver<STATE> step;
  step.SetDirect(ELDLt);
  analysis.SetSolver(step);

  analysis.Run();

  std::string filename = "solution", filename1d = "solwell", filenamePressureSkin = "solpressureskin";
  TPZStack<std::string> fieldnames; 
  fieldnames.Push("Pressure");
  fieldnames.Push("Flux");
  TPZVTKGenerator vtkGen(cmesh, fieldnames, filename, 0, 3, true);
  vtkGen.SetNThreads(0);
  vtkGen.Do();

  // matl2proj->SetScaleFactor(1.);
  // std::set<int> mats = {EPressure2DSkin};
  // TPZVTKGenerator vtkGen2D(cmesh, mats, {"Solution"}, filenamePressureSkin, 0);
  // vtkGen2D.SetNThreads(0);
  // vtkGen2D.Do();


  fieldnames.Push("Divergence");
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

REAL FindClosestX(const REAL x, const std::set<REAL>& nodeCoordsX, const REAL tol) {
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

void CreatePressure2DElsAndOrderIds(TPZGeoMesh* gmesh) {
  const int nel = gmesh->NElements();
  const REAL tol = 1.e-4;
  std::set<int64_t> pressure2Dels;
  std::set<REAL> nodeCoordsX;
  for (int64_t iel = 0; iel < nel; iel++) {
    TPZGeoEl* gel = gmesh->Element(iel);
    if (gel->MaterialId() != ESurfWellCyl) continue;  
    if (gel->HasSubElement()) continue;
    pressure2Dels.insert(iel);
    TPZGeoElBC bc(gel,gel->NSides()-1,EPressure2DSkin); 
    TPZGeoElBC bc2(gel,gel->NSides()-1,EPressureInterface);
    TPZGeoElBC bc3(gel,gel->NSides()-1,EHDivBoundInterface);
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
  const int dim = cmesh->Dimension()-1; // 2D for these pressure elements living in the boundary of the 3d well
  const int matid = EPressure2DSkin;
  TPZNullMaterial<STATE>* mat = new TPZNullMaterial<>(matid,dim);
  cmesh->SetAllCreateFunctionsContinuous();
  cmesh->ApproxSpace().CreateDisconnectedElements(true);
  cmesh->InsertMaterialObject(mat);
  cmesh->SetDefaultOrder(pordWell);

  std::set<int> matidset = {matid};
  cmesh->AutoBuild(matidset);

  const bool isContinuousConnects = false; // the pressure skin elements should not be continuous if HDiv elements in the wellbore.
  if(isContinuousConnects){ // This leads to continuous elements
    std::set<REAL> nodeCoordsX;
    const int64_t nel = cmesh->NElements();
    std::set<int64_t> pressure2Dels;
    for (int64_t iel = 0; iel < nel; iel++) {
      TPZCompEl* cel = cmesh->Element(iel);    
      if (!cel) continue;
      if (cel->Material()->Id() != matid) continue;
      TPZGeoEl* gel = cel->Reference();
      if(gel->HasSubElement()) continue;
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
      if(gel->HasSubElement()) DebugStop();
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
  }
  else {    
    const int64_t nel = cmesh->NElements();
    std::set<int64_t> pressure2Dels;
    std::map<REAL,TPZManVector<REAL,3>> elCentToConnects;
    std::set<REAL> elCentX;
    const REAL tol = 1.e-6;
    for (int64_t iel = 0; iel < nel; iel++) {
      TPZCompEl* cel = cmesh->Element(iel);    
      if (!cel) continue;
      if (cel->Material()->Id() != matid) continue;
      TPZGeoEl* gel = cel->Reference();
      if(gel->HasSubElement()) DebugStop();
      if(gel->NNodes() != 4) DebugStop();

      TPZGeoElSide gelside(gel);
      TPZManVector<REAL,3> cent(3);
      gelside.CenterX(cent);
      InsertXCoorInSet(cent[0], elCentX, 1.e-6);
      REAL closestX = FindClosestX(cent[0], elCentX, tol);      
      if (elCentToConnects.find(closestX) == elCentToConnects.end()) {
        TPZManVector<REAL> newconnects(3);
        newconnects[0] = cmesh->AllocateNewConnect(1, 1, pordWell);
        newconnects[1] = cmesh->AllocateNewConnect(1, 1, pordWell);
        newconnects[2] = cmesh->AllocateNewConnect(pordWell-1, 1, pordWell);
        elCentToConnects[closestX] = newconnects;        
      }


      TPZInterpolatedElement* intEl = dynamic_cast<TPZInterpolatedElement*>(cel);
      if (!intEl) DebugStop();    
      // I am at a pressure 2D element at the 3d well boundary
      pressure2Dels.insert(iel);
      for(int i = 4 ; i < 8 ; i++){
        TPZConnect& c = cel->Connect(i);
        REAL x0 = gel->NodePtr(i%4)->Coord(0), x1 = gel->NodePtr((i+1)%4)->Coord(0);
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

    for (auto iel : pressure2Dels) {
      TPZCompEl* cel = cmesh->Element(iel);
      TPZGeoEl* gel = cel->Reference();
      TPZGeoElSide gelside(gel);
      TPZManVector<REAL,3> cent(3);
      gelside.CenterX(cent);
      REAL closestX = FindClosestX(cent[0], elCentX, tol);
      auto it = elCentToConnects.find(closestX);
      TPZManVector<REAL,3> newconnects;
      if (it != elCentToConnects.end()) {
          newconnects = it->second;
      } else {
          DebugStop(); // closestX not found in elCentToConnects
      }

      for (int i = 0; i < gel->NCornerNodes(); i++) {
        TPZManVector<REAL,3> coor(3);
        gel->NodePtr(i)->GetCoordinates(coor);
        if (fabs(coor[0] - closestX) < tol) DebugStop();
        if (coor[0] < closestX) {
          cel->SetConnectIndex(i,newconnects[0]);
        } else {
          cel->SetConnectIndex(i,newconnects[1]);
        }
      }
      for (int i = 4; i < 8; i++) {
        REAL x0 = gel->NodePtr(i%4)->Coord(0), x1 = gel->NodePtr((i+1)%4)->Coord(0);
        if (fabs(x0 - x1) < 1.e-6) {
          continue;
        }
        const REAL midpoint = (x0+x1)/2.;
        if (fabs(midpoint - closestX) > tol) DebugStop(); // Has to be in the middle!
        cel->SetConnectIndex(i,newconnects[2]);        
      }    
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

void AddHDivBoundInterfaceElements(TPZCompMesh* cmesh, const int porder) {
  cmesh->Reference()->ResetReference();
  cmesh->LoadReferences();
  const int dim = cmesh->Reference()->Dimension();
  const int matid = EHDivBoundInterface;
  TPZNullMaterial<STATE>* mat = new TPZNullMaterial<>(matid, dim - 1);
  cmesh->SetDimModel(dim);
  cmesh->SetAllCreateFunctionsHDiv();
  cmesh->InsertMaterialObject(mat);
  cmesh->SetDefaultOrder(porder);

  std::set<int> matidset = {matid};
  cmesh->AutoBuild(matidset);
}

void AddInterfaceElements(TPZMultiphysicsCompMesh* cmesh, const int matidpressure, const int matidinterface, const int laglevel) {
  cmesh->Reference()->ResetReference();
  cmesh->LoadReferences();
  TPZGeoMesh* gmesh = cmesh->Reference();  
  const int dim = gmesh->Dimension();
  const int64_t nel = gmesh->NElements();
  for (int64_t iel = 0; iel < nel; iel++) {
    TPZGeoEl* gel = gmesh->Element(iel);
    if (!gel) continue;
    if (gel->MaterialId() != EPressureInterface) continue;
    if (gel->Dimension() != dim-1) DebugStop();
    TPZCompElSide comp_pressureSide, comp_hdivSide;
    TPZGeoElSide gelSide(gel);
    for(auto neigh = gelSide.Neighbour(); neigh != gelSide ; neigh++) {
      if(neigh.Element()->MaterialId() == EPressure2DSkin) {
        comp_pressureSide = neigh.Reference();
      }
      if(neigh.Element()->MaterialId() == EHDivBoundInterface) {
        comp_hdivSide = neigh.Reference();
      }
    }
    if (!comp_pressureSide || !comp_hdivSide) DebugStop();
    TPZMultiphysicsInterfaceElement* interfaceel = new TPZMultiphysicsInterfaceElement(*cmesh,gel,comp_pressureSide,comp_hdivSide);
  }
}

void EqualizePressureConnects(TPZCompMesh* cmesh) {
  cmesh->Reference()->ResetReference();
  cmesh->LoadReferences();
  TPZGeoMesh* gmesh = cmesh->Reference();
  const int dim = gmesh->Dimension();
  const int64_t nel = gmesh->NElements();
  for(auto& gel: gmesh->ElementVec()) {
    if (!gel) continue;
    if (gel->MaterialId() != ECurveWell) continue;
    if(gel->HasSubElement()) continue;
    TPZGeoElSide gelside(gel);
    TPZGeoElSide surfwellside = gelside.Neighbour();
    while (surfwellside.Element()->MaterialId() != EPressure2DSkin) {
      if (surfwellside.Element()->HasSubElement()) DebugStop();
      if(surfwellside == gelside) DebugStop();      
      surfwellside = surfwellside.Neighbour();
    }
    TPZCompEl* cel = gel->Reference();
    TPZCompEl* celneigh = surfwellside.Element()->Reference();

    // Set the internal connect of the wellbore as the edge of a 2d pressure skin element
    const int64_t cindex = celneigh->ConnectIndex(surfwellside.Side());
    cel->SetConnectIndex(2,cindex);

    // Set the nodal connect of the 1d wellbore elements equal to the 2d pressure skin element nodal connects
    const int64_t gindex0 = gel->NodeIndex(0), gindex1 = gel->NodeIndex(1);
    TPZGeoEl* gelsurf = surfwellside.Element();
    const int64_t nindex0 = gelsurf->SideNodeIndex(surfwellside.Side(),0), nindex1 = gelsurf->SideNodeIndex(surfwellside.Side(),1);
    const int sidenodelocindex0 = gelsurf->SideNodeLocIndex(surfwellside.Side(),0), sidenodelocindex1 = gelsurf->SideNodeLocIndex(surfwellside.Side(),1);
    if(nindex0 == gindex0 && nindex1 == gindex1){
      cel->SetConnectIndex(0,celneigh->ConnectIndex(sidenodelocindex0));
      cel->SetConnectIndex(1,celneigh->ConnectIndex(sidenodelocindex1));
    } else if(nindex0 == gindex1 && nindex1 == gindex0){
      cel->SetConnectIndex(0,celneigh->ConnectIndex(sidenodelocindex1));
      cel->SetConnectIndex(1,celneigh->ConnectIndex(sidenodelocindex0));
    } else {
      DebugStop();
    }

    std::cout << "Connect indexes of cel: ";
    for (int i = 0; i < cel->NConnects(); i++) {
      std::cout << cel->ConnectIndex(i) << " ";
    }
    std::cout << std::endl;
  }
  cmesh->CleanUpUnconnectedNodes();
}