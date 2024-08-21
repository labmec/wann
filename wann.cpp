#ifdef HAVE_CONFIG_H
#include <pz_config.h>
#endif

#include <DarcyFlow/TPZDarcyFlow.h>
#include <DarcyFlow/TPZMixedDarcyFlow.h>
#include <TPZGmshReader.h>
#include <TPZLinearAnalysis.h>
#include <TPZNullMaterial.h>
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
#include <TPZRefPatternTools.h>

using namespace std;

enum EMatid { ENone,
              EDomain,
              EWell,
              EHeel,
              EDedao,
              EFarField };

const int global_nthread = 8;

TPZGeoMesh* ReadMeshFromGmsh(std::string file_name);
TPZCompMesh* CreateCompMeshH1(TPZGeoMesh* gmesh, const int porder);
void PrintResults(TPZLinearAnalysis& an, TPZCompMesh* cmesh);

int main() {
  std::cout << "--------- Starting simulation ---------" << std::endl;
#ifdef PZ_LOG
  TPZLogger::InitializePZLOG();
#endif

  const int pord = 3;
  const int nrefdirectional = 3;
  // 1) Create gmesh
  TPZGeoMesh* gmesh = ReadMeshFromGmsh("../geo/mesh_rev01.msh");

  std::ofstream out("gmesh.vtk");

  // refine directional towards heel and dedao
  if(nrefdirectional){
    gRefDBase.InitializeRefPatterns(gmesh->Dimension());
    for (int i = 0; i < nrefdirectional; i++) {
      TPZRefPatternTools::RefineDirectional(gmesh, {EHeel,EDedao});
    }    
  } 

  TPZVTKGeoMesh::PrintGMeshVTK(gmesh, out);

  // 2) Create cmesh
  TPZCompMesh* cmeshH1 = CreateCompMeshH1(gmesh, pord);
  std::ofstream out2("cmeshH1.vtk");
  TPZVTKGeoMesh::PrintCMeshVTK(cmeshH1, out2);

  // 3) Solve problem
  TPZLinearAnalysis an(cmeshH1);
  TPZStepSolver<STATE> step;
  TPZSSpStructMatrix<STATE> strmat(cmeshH1);
  strmat.SetNumThreads(global_nthread);
  an.SetStructuralMatrix(strmat);
  step.SetDirect(ECholesky);
  an.SetSolver(step);  
  an.Run();

  // 4) Post process solution using vtk
  PrintResults(an, cmeshH1);

  delete cmeshH1;
  delete gmesh;

  std::cout << "--------- Simulation finished ---------" << std::endl;
}

TPZCompMesh* CreateCompMeshH1(TPZGeoMesh* gmesh, const int porder) {

  // create computational mesh
  TPZCompMesh* cmesh = new TPZCompMesh(gmesh);
  cmesh->SetDefaultOrder(porder);
  const int dim = gmesh->Dimension();
  cmesh->SetDimModel(dim);
  cmesh->SetAllCreateFunctionsContinuous();

  // create domain material
  TPZDarcyFlow* mat = new TPZDarcyFlow(EDomain, dim);
  mat->SetConstantPermeability(1.0);
  cmesh->InsertMaterialObject(mat);

  // create boundary conditions
  const int diri = 0, neu = 1, mixed = 2;
  // Far field
  TPZFMatrix<STATE> val1(1, 1, 0.0);
  TPZManVector<STATE, 1> val2(1, 1.0);
  TPZBndCond* bcFF = mat->CreateBC(mat, EFarField, diri, val1, val2);
  cmesh->InsertMaterialObject(bcFF);

  // create well material
  TPZDarcyFlow* matWell = new TPZDarcyFlow(EWell, dim-1);
  matWell->SetConstantPermeability(1.0e4);
  cmesh->InsertMaterialObject(matWell);

  // create heel boundary condition for well
  val2[0] = 0.0;
  TPZBndCond* bcHeel = matWell->CreateBC(matWell, EHeel, diri, val1, val2);
  cmesh->InsertMaterialObject(bcHeel);

  cmesh->AutoBuild();

  return cmesh;
}

TPZGeoMesh*
ReadMeshFromGmsh(std::string file_name) {
  // read mesh from gmsh
  TPZGeoMesh* gmesh;
  gmesh = new TPZGeoMesh();
  {
    TPZGmshReader reader;
    TPZManVector<std::map<std::string, int>, 4> stringtoint(5);
    stringtoint[2]["dom"] = EDomain;

    stringtoint[1]["wellbore"] = EWell;
    stringtoint[0]["heel"] = EHeel;
    stringtoint[0]["dedao"] = EDedao;
    stringtoint[1]["far_field"] = EFarField;
    reader.SetDimNamePhysical(stringtoint);
    reader.GeometricGmshMesh(file_name, gmesh);
  }

  return gmesh;
}

void PrintResults(TPZLinearAnalysis& an, TPZCompMesh* cmesh) {
  std::cout << "--------- Post Process ---------" << std::endl;
  TPZSimpleTimer postProc("Post processing time");

  const std::string plotfile = "postprocess";
  constexpr int vtkRes{1};
  TPZVec<std::string> fields = {
      "Pressure",
      "Flux"};
  static auto vtk = TPZVTKGenerator(cmesh, fields, plotfile, vtkRes);
  vtk.SetNThreads(0);
  vtk.Do();
  std::cout << "Total time = " << postProc.ReturnTimeDouble() / 1000. << " s" << std::endl;
}