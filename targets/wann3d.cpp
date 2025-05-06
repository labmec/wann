#ifdef HAVE_CONFIG_H
#include <pz_config.h>
#endif

#include <iostream>
#include "TPZWannGeometryTools.h"
#include "TPZWannApproxTools.h"
#include <TPZLinearAnalysis.h>
#include <TPZSSpStructMatrix.h>  //symmetric sparse matrix storage
#include <pzstepsolver.h>
#include "TPZVTKGenerator.h"

using namespace std;

const int global_nthread = 32;

void PostProceDataForANN(TPZGeoMesh* gmesh, ProblemData* SimData, std::string& filename, const int npts);

int main() {
  std::cout << "--------- Starting simulation ---------" << std::endl;
#ifdef PZ_LOG
  TPZLogger::InitializePZLOG();
#endif

  // Problem data
  ProblemData SimData;
  SimData.ReadJson("wann3d.json");
  
  TPZGeoMesh* gmesh = TPZWannGeometryTools::CreateGeoMesh(&SimData);
  
  TPZMultiphysicsCompMesh* cmesh = TPZWannApproxTools::CreateMultiphysicsCompMesh(gmesh, &SimData);
  
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

  fieldnames.Push("Divergence");
  TPZVTKGenerator vtkGen1D(cmesh, fieldnames, filename1d, 0, 1);
  vtkGen1D.SetNThreads(0);
  vtkGen1D.Do();

  std::string anndatafile = "anndata.txt";
  const int npts = 401;
  PostProceDataForANN(gmesh, &SimData, anndatafile, npts);

  delete cmesh;
  delete gmesh;
  std::cout << "--------- Simulation finished ---------" << std::endl;
}

void PostProceDataForANN(TPZGeoMesh* gmesh, ProblemData* SimData, std::string& filename, const int npts) {
  
  std::ofstream out(filename);

  auto& WellboreData = SimData->m_Wellbore;
  const REAL xf = 1., x0 = 0.;
  const REAL dx = (xf - x0) / (npts - 1);
  const REAL ywell = -WellboreData.radius /sqrt(2.), zwell = WellboreData.radius/sqrt(2.);

  int64_t InitialElIndex = -1;
  for (auto gel : gmesh->ElementVec()) {
    if (gel && gel->MaterialId() == SimData->ECurveWell) {
      InitialElIndex = gel->Index();
      break;
    }
  }

  TPZManVector<REAL,3> xvec(3, 0.);
  xvec[1] = ywell;
  xvec[2] = zwell;
  for (int i = 0; i < npts; i++) {
    const REAL x = x0 + i * dx;
    TPZManVector<REAL,3> qsi(1, 0.0);    
    xvec[0] = x;
    TPZGeoEl* gel = TPZGeoMeshTools::FindElementByMatId(gmesh, xvec, qsi, InitialElIndex, {SimData->ECurveWell});
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
    // Now postprocess the pressure and div q of the cel reference by gel
    TPZCompEl* cel = gel->Reference();
    TPZMaterial* mat = cel->Material();
    TPZMixedDarcyFlow* darcy = dynamic_cast<TPZMixedDarcyFlow*>(mat);
    if(!darcy) DebugStop();
    
    const int pind = darcy->VariableIndex("Pressure");
    const int qind = darcy->VariableIndex("Divergence");
    TPZManVector<STATE, 3> output(1);
    cel->Solution(qsi,pind,output);
    const REAL pressure = output[0];
    cel->Solution(qsi, qind, output);
    const REAL divq = output[0];

    auto& WellboreData = SimData->m_Wellbore;
    const REAL wellRad = WellboreData.radius;
    const REAL pff = WellboreData.BCs["surface_farfield"].value;
    const REAL K = divq / (pff - pressure);

    std::cout << "x = " << x << " p = " << pressure << " divq = " << divq << " K = " << K << std::endl;
    out << x << " " << K << " " << wellRad << std::endl;

  }
}