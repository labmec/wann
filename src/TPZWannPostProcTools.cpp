#include "TPZWannPostProcTools.h"

void TPZWannPostProcTools::GenerateTrainingData(TPZGeoMesh* gmesh, ProblemData* SimData) {

  auto& PostProcData = SimData->m_PostProc;
  std::string filename = PostProcData.training_data;
  std::string path = std::string(TRAININGDIR) + "/" + filename;
  std::ofstream out(path, std::ios::app);

  const int npts = PostProcData.training_resolution + 1;

  auto& WellboreData = SimData->m_Wellbore;
  const REAL dx = WellboreData.length / (npts - 1);
  const REAL ywell = -WellboreData.radius /sqrt(2.), zwell = WellboreData.radius/sqrt(2.);

  int64_t InitialElIndex = -1;
  for (auto gel : gmesh->ElementVec()) {
    if (gel && gel->MaterialId() == SimData->ECurveWell) {
      InitialElIndex = gel->Index();
      break;
    }
  }

  TPZManVector<REAL,3> xvec(3, 0.);
  xvec[1] = ywell + WellboreData.eccentricity[1];
  xvec[2] = zwell + WellboreData.eccentricity[2];
  for (int i = 0; i < npts; i++) {
    const REAL x = WellboreData.eccentricity[0] + i * dx;
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

void TPZWannPostProcTools::WriteWellboreVTK(TPZCompMesh* cmesh, ProblemData* SimData) {

  auto& PostProcData = SimData->m_PostProc;
  std::string filename = PostProcData.wellbore_vtk;

  TPZStack<std::string> fieldnames; 
  fieldnames.Push("Pressure");
  fieldnames.Push("Flux");
  fieldnames.Push("Divergence");

  TPZVTKGenerator vtk(cmesh, fieldnames, filename, PostProcData.vtk_resolution, 1);
  vtk.SetNThreads(PostProcData.nthreads);
  vtk.Do();
}

void TPZWannPostProcTools::WriteReservoirVTK(TPZCompMesh* cmesh, ProblemData* SimData) {

  auto& PostProcData = SimData->m_PostProc;
  std::string filename = PostProcData.reservoir_vtk;

  TPZStack<std::string> fieldnames; 
  fieldnames.Push("Pressure");
  fieldnames.Push("Flux");

  bool safe_check = SimData->m_Mesh.ToCylindrical? true : false;

  TPZVTKGenerator vtk(cmesh, fieldnames, filename, PostProcData.vtk_resolution, 3, safe_check);
  vtk.SetNThreads(PostProcData.nthreads);
  vtk.Do();
}

void TPZWannPostProcTools::WriteVTKs(TPZCompMesh* cmesh, ProblemData* SimData) {

  WriteReservoirVTK(cmesh, SimData);
  WriteWellboreVTK(cmesh, SimData);
}

void TPZWannPostProcTools::PostProcessAllData(TPZCompMesh* cmesh, TPZGeoMesh* gmesh, ProblemData* SimData) {

  WriteReservoirVTK(cmesh, SimData);
  WriteWellboreVTK(cmesh, SimData);
  GenerateTrainingData(gmesh, SimData);
}