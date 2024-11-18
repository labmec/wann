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
              ESurfWellCylNonLin };

const int global_nthread = 8;

TPZGeoMesh* ReadMeshFromGmsh(std::string file_name);
void ModifyGeometricMeshToCylWell(TPZGeoMesh* gmesh);

int main() {
  std::cout << "--------- Starting simulation ---------" << std::endl;
#ifdef PZ_LOG
  TPZLogger::InitializePZLOG();
#endif

  const int nrefdirectional = 0;

  TPZGeoMesh* gmesh = ReadMeshFromGmsh("../../geo/mesh3D_rev04.msh");

  {
    std::ofstream out("gmeshorig.vtk");
    TPZVTKGeoMesh::PrintGMeshVTK(gmesh, out);
  }
  
  ModifyGeometricMeshToCylWell(gmesh);

  {
    TPZCheckGeom checkgeom(gmesh);
    checkgeom.UniformRefine(2);
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

  

  delete gmesh;
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
    
    stringtoint[1]["curve_wellbore"] = ECurveWell;
    stringtoint[1]["curve_heel"] = ECurveHeel;
    stringtoint[1]["curve_toe"] = ECurveToe;
    reader.SetDimNamePhysical(stringtoint);
    reader.GeometricGmshMesh(file_name, gmesh);
  }

  return gmesh;
}

void ModifyGeometricMeshToCylWell(TPZGeoMesh* gmesh) {
  const REAL cylradius = 0.05;
  const TPZManVector<REAL,3> cylcenter = {0.,0.,0.}, cylaxis = {1.,0.,0.};
  int64_t nel = gmesh->NElements();
  for(int64_t iel = 0; iel < nel ; iel++) {
    TPZGeoEl* gel = gmesh->Element(iel);
    if(!gel) continue;
    if(gel->HasSubElement()) DebugStop();
    if(gel->MaterialId() != ESurfWellCyl) continue;
    
    TPZManVector<int64_t, 4> nodeindices;
    gel->GetNodeIndices(nodeindices);
    for(auto& it : nodeindices) {
      TPZManVector<REAL,3> coor(3);
      gmesh->NodeVec()[it].GetCoordinates(coor);
      REAL xcoortemp = coor[0];
      coor[0] = 0.;
      const REAL norm = Norm(coor)/cylradius;
      coor /= norm;
      coor[0] = xcoortemp;
      gmesh->NodeVec()[it].SetCoord(coor);
    }
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
    if(gel->MaterialId() != ESurfWellCyl) continue;
    int nsides = gel->NSides();
    for (int iside = gel->NCornerNodes(); iside < nsides; iside++) {
      TPZGeoElSide gelSide(gel,iside);
      TPZStack<TPZGeoElSide> allneigh;
      for(auto neigh = gelSide.Neighbour(); neigh != gelSide ; neigh++) {
        if(neigh.Element()->IsGeoBlendEl()) {
          continue;
        }
        allneigh.Push(neigh);
      }
      for(auto it : allneigh){
        TPZChangeEl::ChangeToGeoBlend(gmesh, it.Element()->Index());
      }
    }
  }  
}