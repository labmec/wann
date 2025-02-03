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
#include "pzfstrmatrix.h"
#include "pzvec_extras.h"

using namespace std;

enum EMatid { ENone,
              EDomain,
              EWell,
              EHeel,
              EDedao,
              EFarField };

struct PostProcElData {
  TPZGeoEl* gel;
  TPZManVector<REAL,1> qsi = {-10}; // a value that is not valid
  TPZManVector<REAL,3> x = {0.,0.,0.};
};


const int global_nthread = 8;

TPZGeoMesh* ReadMeshFromGmsh(std::string file_name);
TPZCompMesh* CreateCompMeshH1(TPZGeoMesh* gmesh, const int porder, const REAL pff, const REAL pheel, const REAL Kres, const REAL Kwell);
void PrintResults(TPZLinearAnalysis& an, TPZCompMesh* cmesh);
void FindElementsAndPtsToPostProc(TPZGeoMesh* gmesh, TPZStack<PostProcElData>& postProcData, const int matid, const int npts, const REAL x0, const REAL xf);
void PostProcDataForANN(TPZGeoMesh* gmesh, TPZStack<PostProcElData>& postProcData, const REAL pff, const REAL wellR);

int main() {
  std::cout << "--------- Starting simulation ---------" << std::endl;
#ifdef PZ_LOG
  TPZLogger::InitializePZLOG();
#endif

  const int pord = 5;
  const int nrefdirectional = 3;
  REAL pff = 22064967.11, pheel = 11767982.46;
  const REAL Kres = 1.e-13;
  // const REAL wellR = 0.0005;
  // const REAL wellR = 0.001;
  const REAL wellR = 0.005;
  const REAL mu = 0.005;
  REAL Kwell = M_PI * wellR * wellR * wellR * wellR / (8.0 * mu);
//   Kwell = 1000*Kres;
  // Kwell = 1e-6;
  // 1) Create gmesh
  TPZGeoMesh* gmesh = ReadMeshFromGmsh("../geo/mesh_rev01.msh");

  std::ofstream out("gmesh.vtk");

  // refine directional towards heel and dedao
  if (nrefdirectional) {
    gRefDBase.InitializeRefPatterns(gmesh->Dimension());
    for (int i = 0; i < nrefdirectional; i++) {
      TPZRefPatternTools::RefineDirectional(gmesh, {EHeel, EDedao});
    }
  }

  TPZVTKGeoMesh::PrintGMeshVTK(gmesh, out);

  // 2) Create cmesh
  TPZCompMesh* cmeshH1 = CreateCompMeshH1(gmesh, pord, pff, pheel, Kres, Kwell);
  std::ofstream out2("cmeshH1.vtk");
  TPZVTKGeoMesh::PrintCMeshVTK(cmeshH1, out2);

  // 3) Solve problem
  TPZLinearAnalysis an(cmeshH1);
  TPZStepSolver<STATE> step;
  // TPZSSpStructMatrix<STATE> strmat(cmeshH1);
  // TPZFStructMatrix<STATE> strmat(cmeshH1);
  TPZSkylineStructMatrix<STATE> strmat(cmeshH1);
  strmat.SetNumThreads(global_nthread);
  an.SetStructuralMatrix(strmat);
  step.SetDirect(ECholesky);
  an.SetSolver(step);
  an.Run();

  // 4) Post process solution using vtk
  PrintResults(an, cmeshH1);

  // 5) Post process pressure at well on several points using FindElementByMatId
  TPZStack<PostProcElData> postProcData;
  const REAL x0 = 0.0, xf = 400.0;
  const int npts = 4001;
  FindElementsAndPtsToPostProc(gmesh, postProcData, EWell, npts, x0, xf);

  // 6) Post process pressure and divergence of q at well on several points using postprocdata
  PostProcDataForANN(gmesh, postProcData, pff, wellR);

  delete cmeshH1;
  delete gmesh;

  std::cout << "--------- Simulation finished ---------" << std::endl;
}

void PostProcDataForANN(TPZGeoMesh* gmesh, TPZStack<PostProcElData>& postProcData, const REAL pff, const REAL wellR) {
  std::ofstream out("postprocdata.txt");
  out << std::setprecision(15);
  for (auto& data : postProcData) {
    TPZVec<REAL> qsi(3, 0.0);
    TPZCompEl* cel = data.gel->Reference();
    if(!cel) DebugStop();
    TPZMaterial* mat = cel->Material();
    TPZDarcyFlow* darcy = dynamic_cast<TPZDarcyFlow*>(mat);
    if(!darcy) DebugStop();
    const int pind = darcy->VariableIndex("Pressure");
    const int qind = darcy->VariableIndex("Flux");
    TPZManVector<STATE, 3> output(1);
    cel->Solution(data.qsi,pind,output);
    const REAL pressure = output[0];
    cel->Solution(data.qsi, qind, output);
    const REAL q = output[0];

    std::cout << "x = " << data.x << " p = " << pressure << " q = " << q << std::endl;
    out << data.x[0] << " " << pressure << " " << q << " " << wellR << std::endl;
  }
}

void FindElementsAndPtsToPostProc(TPZGeoMesh* gmesh, TPZStack<PostProcElData>& postProcData, const int matid, const int npts, const REAL x0, const REAL xf) {
  const REAL dx = (xf - x0) / (npts - 1);

  int64_t InitialElIndex = -1;
  for (auto gel : gmesh->ElementVec()) {
    if (gel && gel->MaterialId() == matid) {
      InitialElIndex = gel->Index();
      break;
    }
  }


  for (int i = 0; i < npts; i++) {
    const REAL x = x0 + i * dx;
    TPZManVector<REAL,3> qsi(3, 0.0);    
    TPZManVector<REAL,3> xvec(3, 0.);
    xvec[0] = x;
    TPZGeoEl* gel = TPZGeoMeshTools::FindElementByMatId(gmesh, xvec, qsi, InitialElIndex, {matid});
    if (!gel) {
      std::cout << "Element not found for x = " << x << std::endl;
      DebugStop();
    }
#ifdef PZDEBUG
    TPZManVector<REAL, 3> xcheck(3);
    gel->X(qsi, xcheck);
    REAL distance = dist(xvec, xcheck);
    if(distance > 1.e-8){
      DebugStop(); // check if the element found is the correct one
    }
#endif
    postProcData.Push({gel, qsi, xvec}); // Creates a struct entry with gel and qsi.
  }

}

TPZCompMesh* CreateCompMeshH1(TPZGeoMesh* gmesh, const int porder, const REAL pff, const REAL pheel, const REAL Kres, const REAL Kwell) {
  // create computational mesh
  TPZCompMesh* cmesh = new TPZCompMesh(gmesh);
  cmesh->SetDefaultOrder(porder);
  const int dim = gmesh->Dimension();
  cmesh->SetDimModel(dim);
  cmesh->SetAllCreateFunctionsContinuous();

  // create domain material
  TPZDarcyFlow* mat = new TPZDarcyFlow(EDomain, dim);
  mat->SetConstantPermeability(Kres);
  cmesh->InsertMaterialObject(mat);

  // create boundary conditions
  const int diri = 0, neu = 1, mixed = 2;
  // Far field
  TPZFMatrix<STATE> val1(1, 1, 0.0);
  TPZManVector<STATE, 1> val2(1, pff);
  TPZBndCond* bcFF = mat->CreateBC(mat, EFarField, diri, val1, val2);
  cmesh->InsertMaterialObject(bcFF);

  // create well material
  TPZDarcyFlow* matWell = new TPZDarcyFlow(EWell, dim - 1);
  matWell->SetConstantPermeability(Kwell);
  cmesh->InsertMaterialObject(matWell);

  // create heel boundary condition for well
  val2[0] = pheel;
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