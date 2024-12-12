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
              EPressureInterface };

const int global_nthread = 8;

TPZGeoMesh* ReadMeshFromGmsh(std::string file_name);
void ModifyGeometricMeshToCylWell(TPZGeoMesh* gmesh);
void CreatePressure2DElsAndOrderIds(TPZGeoMesh* gmesh);

int main() {
  std::cout << "--------- Starting simulation ---------" << std::endl;
#ifdef PZ_LOG
  TPZLogger::InitializePZLOG();
#endif

  const int nrefdirectional = 0;
  const int nref = 0;
  const int pOrder = 1;


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

  // Create GeoelBCs in same location as geoels with ESurfWell
  CreatePressure2DElsAndOrderIds(gmesh);

  //
  
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
  matdarcy->SetConstantPermeability(1.);
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

  TPZMultiphysicsCompMesh* cmesh = hdivCreator.CreateApproximationSpace();

  TPZLinearAnalysis analysis(cmesh);
  TPZSSpStructMatrix<STATE> skylstr(cmesh);
  analysis.SetStructuralMatrix(skylstr);

  TPZStepSolver<STATE> step;
  step.SetDirect(ELDLt);
  analysis.SetSolver(step);

  analysis.Run();

  std::string filename = "solution";
  TPZStack<std::string> fieldnames; 
  fieldnames.Push("Pressure");
  fieldnames.Push("Flux");
  TPZVTKGenerator vtkGen(cmesh, fieldnames, filename, 0);
  vtkGen.SetNThreads(0);
  vtkGen.Do();
  

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
    TPZManVector<std::map<std::string, int>, 4> stringtoint(8);
    stringtoint[3]["volume_reservoir"] = EDomain;

    stringtoint[2]["surface_wellbore_cylinder"] = ESurfWellCyl;
    stringtoint[2]["surface_wellbore_heel"] = ESurfHeel;
    stringtoint[2]["surface_wellbore_toe"] = ESurfToe;
    stringtoint[2]["surface_farfield"] = EFarField;
    stringtoint[2]["surface_cap_rock"] = EFarField;
    
    stringtoint[1]["curve_wellbore"] = ECurveWell;
    stringtoint[1]["curve_heel"] = ECurveHeel;
    stringtoint[1]["curve_toe"] = ECurveToe;
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
}

void CreatePressure2DElsAndOrderIds(TPZGeoMesh* gmesh) {
  const int nel = gmesh->NElements();
  const REAL tol = 1.e-6;
  std::set<int64_t> pressure2Dels;
  std::set<REAL> nodeCoordsX;
  for (int64_t iel = 0; iel < nel; iel++) {
    TPZGeoEl* gel = gmesh->Element(iel);
    if (gel->MaterialId() != ESurfWellCyl) continue;
    TPZGeoElBC bc(gel,gel->NSides()-1,EPressureInterface);
    TPZGeoElBC bc(gel,gel->NSides()-1,EPressure2DSkin);    
    pressure2Dels.insert(bc.CreatedElement()->Index());
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

}