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

  cmesh->Reference()->ResetReference();
  cmesh->LoadReferences();

  WriteReservoirVTK(cmesh, SimData);
  WriteWellboreVTK(cmesh, SimData);
  GenerateTrainingData(gmesh, SimData);
}

// Works for H(div) meshes.
// TPZVec<REAL> TPZWannPostProcTools::ComputeWellFluxes(TPZCompMesh* cmesh, ProblemData* SimData) {
//   TPZVec<REAL> fluxes(9,0.);

//   int64_t ncel = cmesh->NElements();
//   for (int64_t iel = 0; iel < ncel; iel++) {
//     TPZCompEl* cel = cmesh->Element(iel);
//     if (!cel) continue;
//     if (cel->Material()->Id() != SimData->EHDivBoundInterface) continue;

//     TPZGeoEl* gel = cel->Reference();
//     if (!gel) DebugStop();
//     TPZManVector<REAL,3> qsiCenter(gel->Dimension()); // parametric coords (size = gel dimension)
//     TPZManVector<REAL,3> xCenter(3,0.);               // physical coords (3D vector)

//     gel->CenterPoint(gel->NSides()-1, qsiCenter); // get parametric center of the volume side
//     gel->X(qsiCenter, xCenter);
//     REAL wellLength = SimData->m_Wellbore.length;

//     int segment;
//     for (int j = 0; j < 9; j++) {
//       if (xCenter[0] > j*wellLength/9.) {
//         segment = j;
//       }
//     }

//     // Loop on integration points
//     const TPZIntPoints &intrule = cel->GetIntegrationRule();
//     for (int ip = 0; ip < intrule.NPoints(); ip++) {
//       TPZManVector<REAL, 3> ptInElement(gel->Dimension());
//       REAL weight, detjac;
//       intrule.Point(ip, ptInElement, weight);
//       TPZFNMatrix<9, REAL> jacobian, axes, jacinv;
//       gel->Jacobian(ptInElement, jacobian, axes, detjac, jacinv);
//       weight *= fabs(detjac);

//       // get normal flux
//       REAL normalFLux = 0.;

//       fluxes[segment] += segment * normalFLux * weight;
//     }
//   }

//   return fluxes;
// }

// Works for H(div) meshes.
TPZVec<REAL> TPZWannPostProcTools::ComputeWellFluxes(TPZCompMesh* cmesh, ProblemData* SimData) {
  TPZVec<REAL> fluxes(9,0.);

  int64_t ncel = cmesh->NElements();
  for (int64_t iel = 0; iel < ncel; iel++) {
    TPZCompEl* cel = cmesh->Element(iel);
    if (!cel) continue;
    if (cel->Material()->Id() != SimData->EHDivBoundInterface) continue;

    TPZGeoEl* gel = cel->Reference();
    if (!gel) DebugStop();
    TPZManVector<REAL,3> qsiCenter(gel->Dimension());
    TPZManVector<REAL,3> xCenter(3,0.);

    int segment = 0; // Make it depend on of xCenter[0]

    int nconnects = cel->NConnects();
    if (nconnects != 1) DebugStop(); // expecting only one connect for boundary element

    int64_t connectIndex = cel->ConnectIndex(0);
    TPZConnect &connect = cel->Connect(0);
    int blockPos = connect.SequenceNumber();
    int blockSize = cmesh->Block().Size(blockPos);
    int solPos = cmesh->Block().Position(blockPos);

    TPZFMatrix<STATE> &sol = cmesh->Solution(); 

    for (int i = 0; i < blockSize; i++) {
      STATE value = sol(solPos + i,0);
      fluxes[segment] += value;
    }
  }

  return fluxes;
}

TPZVec<REAL> TPZWannPostProcTools::ComputeWellFluxesH1(TPZCompMesh* cmesh, ProblemData* SimData) {
  TPZVec<REAL> fluxes(9,0.);

  // Ensure references point to H1 mesh
  cmesh->Reference()->ResetReference();
  cmesh->LoadReferences();

  int64_t ncel = cmesh->NElements();
  for (int64_t iel = 0; iel < ncel; iel++) {
    TPZCompEl* cel = cmesh->Element(iel);
    if (!cel) continue;
    if (cel->Material()->Id() != SimData->EPressure2DSkin) continue;

    TPZGeoEl* gel = cel->Reference();
    if (!gel) DebugStop();

    int segment = 0; // Make it depend on of xCenter[0]

    // Get 3D neighbor element
    TPZGeoElSide gelside(gel, gel->NSides()-1);
    TPZGeoElSide neighside = gelside.HasNeighbour(SimData->EDomain);
    TPZGeoEl* gel3D = neighside.Element();
    if(!gel3D) DebugStop();
    TPZCompEl* cel3D = gel3D->Reference();
    if (!cel3D) DebugStop();
    if (cel3D->Material()->Id() != SimData->EDomain) DebugStop();

    int sideWell = neighside.Side();
    const TPZIntPoints* intrule = nullptr;
    intrule = gel3D->CreateSideIntegrationRule(sideWell, 4);
    for (int ip = 0; ip < intrule->NPoints(); ip++) {
      TPZManVector<REAL, 3> ptOnSide(gel->Dimension());
      TPZFNMatrix<9, REAL> jacobian, axes, jacinv;
      REAL weight, detjac;
      gel->Jacobian(ptOnSide, jacobian, axes, detjac, jacinv);

      // Compute normal
      // Must done in the integration point to account for curved sides
      TPZManVector<REAL, 3> v1(3), v2(3), normal(3);
      v1[0] = axes(0, 0);
      v1[1] = axes(0, 1);
      v1[2] = axes(0, 2);
      v2[0] = axes(1, 0);
      v2[1] = axes(1, 1);
      v2[2] = axes(1, 2);

      normal[0] = v1[1] * v2[2] - v1[2] * v2[1];
      normal[1] = v1[2] * v2[0] - v1[0] * v2[2];
      normal[2] = v1[0] * v2[1] - v1[1] * v2[0];

      REAL norm = sqrt(normal[0] * normal[0] + normal[1] * normal[1] +
                       normal[2] * normal[2]);
      if (norm > 1e-12) {
        normal[0] /= norm;
        normal[1] /= norm;
        normal[2] /= norm;
      }

      intrule->Point(ip, ptOnSide, weight);
      TPZManVector<REAL, 3> ptInElement(gel3D->Dimension());
      TPZTransform<> trans = gel3D->SideToSideTransform(sideWell, gel3D->NSides()-1);
      trans.Apply(ptOnSide, ptInElement);
      weight *= fabs(detjac);
      TPZManVector<REAL, 3> sigh(gel3D->Dimension(), 0.0);
      cel3D->Solution(ptInElement, 7, sigh);

      // get normal flux
      REAL normalFLux = sigh[0] * normal[0] + sigh[1] * normal[1] + sigh[2] * normal[2];
      fluxes[segment] += normalFLux * weight;
    }
  }

  return fluxes;
}