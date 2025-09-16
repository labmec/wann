
#include <iostream>
#include "TPZGenGrid2D.h"
#include "TPZGenGrid3D.h"
#include "TPZRefPatternDataBase.h"
#include "TPZVTKGeoMesh.h"
#include "TPZVTKGenerator.h"
#include "TPZNullMaterial.h"
#include "DarcyFlow/TPZDarcyFlow.h"
#include "DarcyFlow/TPZMixedDarcyFlow.h"
#include "TPZMultiphysicsCompMesh.h"
#include "pzmultiphysicselement.h"
#include "TPZLinearAnalysis.h"
#include "TPZSSpStructMatrix.h"
#include "pzstepsolver.h"
#include "TPZAnalyticSolution.h"
#include "pzlog.h"
#include "pzintel.h"

// ================
// Global variables
// ================

// Exact solution
TLaplaceExample1 gexact;

// Material IDs for domain and boundaries
enum EnumMatIds {
  EMatId = 1, 
  EBottom = 2,
  ERight = 3, 
  ETop = 4,   
  ELeft = 5,
  EFront = 6,
  EBack = 7 
};

// ===================
// Function prototypes
// ===================

// Creates a geometric mesh using TPZGenGrid3D
TPZGeoMesh* createGeoMesh3D(
  const TPZManVector<int, 3> &nelDiv, 
  const TPZManVector<REAL, 3> &minX, 
  const TPZManVector<REAL, 3> &maxX);

// Creates a 2D geometric mesh using TPZGenGrid2D
TPZGeoMesh* createGeoMesh2D(
  const TPZManVector<int, 2> &nelDiv, 
  const TPZManVector<REAL, 2> &minX, 
  const TPZManVector<REAL, 2> &maxX);

// Creates a computational mesh for H1 approximation
TPZCompMesh* createCompMeshH1(TPZGeoMesh *gmesh, int order = 1);

// Creates a computational mesh for mixed approximation
TPZMultiphysicsCompMesh* createCompMeshMixed(TPZGeoMesh *gmesh, int order = 1);

// Error estimation function - computes flux difference between H1 and Mixed methods
int64_t ErrorEstimation(TPZCompMesh* cmesh, TPZMultiphysicsCompMesh* cmeshMixed);

void RefineInteriorElement(TPZCompMesh* cmesh);

void RefineBoundaryElement(TPZCompMesh* cmesh);

// =============
// Main function
// =============

int main (int argc, char * const argv[]) {

  // --- Set up ---

  // Initialize logger
  #ifdef PZ_LOG
  TPZLogger::InitializePZLOG("logpz.txt");
  #endif

  // Initialize uniform refinements for 1D and 2D elements
  gRefDBase.InitializeUniformRefPattern(EOned);
  gRefDBase.InitializeUniformRefPattern(EQuadrilateral);
  gRefDBase.InitializeUniformRefPattern(ETriangle);

  const int nthreads = 0;
  int iteration = 0;

  // --- Solve darcy problem ---

  int order = 1; // Polynomial order

  // Initial geometric mesh
  TPZGeoMesh* gmesh = createGeoMesh3D({3, 3, 3}, {0., 0., 0.}, {1., 1., 1.});
  // TPZGeoMesh* gmesh = createGeoMesh2D({3, 3}, {0., 0.}, {1., 1.});
  std::ofstream outGmesh("gmesh.txt");
  gmesh->Print(outGmesh);

  while (iteration < 2) {
    iteration++;

    TPZMultiphysicsCompMesh* cmeshMixed = createCompMeshMixed(gmesh, order);
    TPZCompMesh* cmeshH1 = createCompMeshH1(gmesh, order);

    // Mixed solver
    TPZLinearAnalysis anMixed(cmeshMixed, RenumType::EMetis);
    TPZSSpStructMatrix<STATE> matMixed(cmeshMixed); 
    matMixed.SetNumThreads(nthreads); 
    anMixed.SetStructuralMatrix(matMixed);
    TPZStepSolver<STATE> stepMixed;
    stepMixed.SetDirect(ELDLt);
    anMixed.SetSolver(stepMixed);
    anMixed.Run();

    // H1 solver
    TPZLinearAnalysis anH1(cmeshH1);
    TPZSSpStructMatrix<STATE> matH1(cmeshH1);
    matH1.SetNumThreads(nthreads);
    anH1.SetStructuralMatrix(matH1);
    TPZStepSolver<STATE> stepH1;
    stepH1.SetDirect(ECholesky);
    anH1.SetSolver(stepH1);
    anH1.Run();

    // ---- Plotting ---

    {
      const std::string plotfile = "darcy_mixed_" + std::to_string(iteration);
      constexpr int vtkRes{0};
      TPZManVector<std::string, 2> fields = {"Flux", "Pressure"};
      auto vtk = TPZVTKGenerator(cmeshMixed, fields, plotfile, vtkRes);
      vtk.Do();
    }

    {
      const std::string plotfile = "darcy_h1_" + std::to_string(iteration);
      constexpr int vtkRes{0};
      TPZManVector<std::string, 2> fields = {"Flux", "Pressure"};
      auto vtk = TPZVTKGenerator(cmeshH1, fields, plotfile, vtkRes);
      vtk.Do();
    }

    // --- Fake refinement ---

    RefineBoundaryElement(cmeshMixed);

    // --- Clean up ---

    delete cmeshMixed;
    delete cmeshH1;
  }
}

// =========
// Functions
// =========

TPZGeoMesh* createGeoMesh3D(
  const TPZManVector<int, 3> &nelDiv, 
  const TPZManVector<REAL, 3> &minX, 
  const TPZManVector<REAL, 3> &maxX) {

  TPZGenGrid3D generator(minX, maxX, nelDiv, MMeshType::EHexahedral);

  generator.BuildVolumetricElements(EMatId);
  TPZGeoMesh *gmesh = generator.BuildBoundaryElements(EBottom, ELeft, EFront, ERight, EBack, ETop);

  return gmesh;
}

TPZGeoMesh* createGeoMesh2D(
  const TPZManVector<int, 2> &nelDiv, 
  const TPZManVector<REAL, 2> &minX, 
  const TPZManVector<REAL, 2> &maxX) {

  TPZGeoMesh *gmesh = new TPZGeoMesh;
  TPZGenGrid2D generator(nelDiv, minX, maxX);
  generator.SetElementType(MMeshType::EQuadrilateral);
  generator.Read(gmesh, EMatId); 
  generator.SetBC(gmesh, 4, EFront);  
  generator.SetBC(gmesh, 5, ERight);
  generator.SetBC(gmesh, 6, EBack);
  generator.SetBC(gmesh, 7, ELeft);

  return gmesh;
}

TPZCompMesh* createCompMeshH1(TPZGeoMesh *gmesh, int order) {
  TPZCompMesh *cmesh = new TPZCompMesh(gmesh);
  cmesh->SetDimModel(gmesh->Dimension());
  cmesh->SetDefaultOrder(order); // Polynomial order
  cmesh->SetAllCreateFunctionsContinuous(); // H1 Elements

  // Add materials (weak formulation)
  TPZDarcyFlow *mat = new TPZDarcyFlow(EMatId, gmesh->Dimension());
  mat->SetConstantPermeability(1.0); // Set constant permeability
  cmesh->InsertMaterialObject(mat);

  // Add boundary conditions
  TPZManVector<REAL,1> val2(1,0.); // Part that goes to the RHS vector
  TPZFMatrix<REAL> val1(1,1,0.); // Part that goes to the Stiffnes matrix
  TPZBndCondT<REAL> *bcond;


  val2[0] = 0.;
  bcond = mat->CreateBC(mat, ERight, 1, val1, val2);
  cmesh->InsertMaterialObject(bcond);

  bcond = mat->CreateBC(mat, ELeft, 1, val1, val2);
  cmesh->InsertMaterialObject(bcond);

  val2[0] = 0.;
  bcond = mat->CreateBC(mat, EFront, 0, val1, val2);
  cmesh->InsertMaterialObject(bcond);

  val2[0] = 1.;
  bcond = mat->CreateBC(mat, EBack, 0, val1, val2);
  cmesh->InsertMaterialObject(bcond);

  if (gmesh->Dimension() == 3) {
    val2[0] = 0.;
    bcond = mat->CreateBC(mat, EBottom, 1, val1, val2);
    cmesh->InsertMaterialObject(bcond);

    bcond = mat->CreateBC(mat, ETop, 1, val1, val2);
    cmesh->InsertMaterialObject(bcond);
  }

  cmesh->AutoBuild();

  {
    std::ofstream out("cmeshH1.txt");
    cmesh->Print(out);
  }

  return cmesh;
}

TPZMultiphysicsCompMesh* createCompMeshMixed(TPZGeoMesh *gmesh, int order) {

  // --- Flux atomic cmesh ----

  TPZCompMesh *cmeshFlux = new TPZCompMesh(gmesh);
  cmeshFlux->SetDimModel(gmesh->Dimension());
  cmeshFlux->SetDefaultOrder(order);
  if (order < 1) {
    cmeshFlux->ApproxSpace().SetHDivFamily(HDivFamily::EHDivConstant);
  }
  cmeshFlux->SetAllCreateFunctionsHDiv();

  // Add materials (weak formulation)
  TPZNullMaterial<STATE> *mat = new TPZNullMaterial(EMatId, gmesh->Dimension());  
  cmeshFlux->InsertMaterialObject(mat);

  // Create boundary conditions
  TPZManVector<REAL,1> val2(1,0.); // Part that goes to the RHS vector
  TPZFMatrix<REAL> val1(1,1,0.); // Part that goes to the Stiffnes matrix
  TPZBndCondT<REAL> *bcond;

  bcond = mat->CreateBC(mat, ERight, 0, val1, val2);
  cmeshFlux->InsertMaterialObject(bcond);
  
  bcond = mat->CreateBC(mat, ELeft, 0, val1, val2);
  cmeshFlux->InsertMaterialObject(bcond);

  bcond = mat->CreateBC(mat, EFront, 0, val1, val2);
  cmeshFlux->InsertMaterialObject(bcond);

  bcond = mat->CreateBC(mat, EBack, 0, val1, val2);
  cmeshFlux->InsertMaterialObject(bcond);

  if (gmesh->Dimension() == 3) {
    bcond = mat->CreateBC(mat, EBottom, 0, val1, val2);
    cmeshFlux->InsertMaterialObject(bcond);

    bcond = mat->CreateBC(mat, ETop, 0, val1, val2);
    cmeshFlux->InsertMaterialObject(bcond);
  }

  cmeshFlux->AutoBuild();

  {
    std::ofstream out("cmeshFlux.txt");
    cmeshFlux->Print(out);
  }

  // --- Pressure atomic cmesh ---

  TPZCompMesh *cmeshPressure = new TPZCompMesh(gmesh);
  cmeshPressure->SetDimModel(gmesh->Dimension());
  cmeshPressure->SetDefaultOrder(order);
  if (order < 1) {
    cmeshPressure->SetAllCreateFunctionsDiscontinuous();
  } else {
    cmeshPressure->SetAllCreateFunctionsContinuous();
    cmeshPressure->ApproxSpace().CreateDisconnectedElements(true);
  }

  // Add materials (weak formulation)
  cmeshPressure->InsertMaterialObject(mat);

  // Set up the computational mesh
  cmeshPressure->AutoBuild();

  int ncon = cmeshPressure->NConnects();
  const int lagLevel = 1; // Lagrange multiplier level
  for(int i=0; i<ncon; i++)
  {
      TPZConnect &newnod = cmeshPressure->ConnectVec()[i]; 
      newnod.SetLagrangeMultiplier(lagLevel);
  }

  {
    std::ofstream out("cmeshPressure.txt");
    cmeshPressure->Print(out);
  }

  // --- Multiphysics mesh ---

  TPZMultiphysicsCompMesh *cmesh = new TPZMultiphysicsCompMesh(gmesh);
  cmesh->SetDimModel(gmesh->Dimension());
  cmesh->SetDefaultOrder(1);
  cmesh->ApproxSpace().Style() = TPZCreateApproximationSpace::EMultiphysics;

  // Add materials (weak formulation)
  TPZMixedDarcyFlow *matDarcy = new TPZMixedDarcyFlow(EMatId, gmesh->Dimension());  
  matDarcy->SetConstantPermeability(1.0);
  cmesh->InsertMaterialObject(matDarcy);

  // Create, set and add boundary conditions
  val2[0] = 0.;
  bcond = matDarcy->CreateBC(matDarcy, ERight, 1, val1, val2);
  cmesh->InsertMaterialObject(bcond);
  
  bcond = matDarcy->CreateBC(matDarcy, ELeft, 1, val1, val2);
  cmesh->InsertMaterialObject(bcond);

  val2[0] = 0.;
  bcond = matDarcy->CreateBC(matDarcy, EFront, 0, val1, val2);
  cmesh->InsertMaterialObject(bcond);

  val2[0] = 1.;
  bcond = matDarcy->CreateBC(matDarcy, EBack, 0, val1, val2);
  cmesh->InsertMaterialObject(bcond);

  if (gmesh->Dimension() == 3) {
    val2[0] = 0.;
    bcond = matDarcy->CreateBC(matDarcy, EBottom, 1, val1, val2);
    cmesh->InsertMaterialObject(bcond);

    bcond = matDarcy->CreateBC(matDarcy, ETop, 1, val1, val2);
    cmesh->InsertMaterialObject(bcond);
  }

  // Incorporate the atomic meshes into the multiphysics mesh
  TPZManVector<TPZCompMesh *,2> cmeshes(2);
  cmeshes[0] = cmeshFlux;
  cmeshes[1] = cmeshPressure;

  TPZManVector<int> active(cmeshes.size(),1);    
  cmesh->BuildMultiphysicsSpace(active, cmeshes);

  {
    std::ofstream out("cmeshMP.txt");
    cmesh->Print(out);
  }

  return cmesh;
}

void RefineInteriorElement(TPZCompMesh* cmesh) {
  TPZGeoMesh *gmesh = cmesh->Reference();
  if (!gmesh) DebugStop();
  int dim  = gmesh->Dimension();

  int64_t nels = gmesh->NElements();
  for (int64_t el = 0; el < nels; ++el) {
    bool flag = true;
    TPZGeoEl *gel = gmesh->Element(el);
    if (!gel || gel->Dimension() != dim || gel->HasSubElement()) continue;
    int firstside = gel->FirstSide(dim-1);
    int lastside = gel->FirstSide(dim);
    for (int side = firstside; side < lastside; ++side){
      TPZGeoElSide gelSide(gel,side);
      TPZGeoElSide neigh = gelSide.Neighbour();
      if (neigh.Element()->Dimension() == (dim-1)) {
        flag = false;
      }
    }
    TPZVec<TPZGeoEl *> pv;
    if (flag == true) gel->Divide(pv);
  }
}

void RefineBoundaryElement(TPZCompMesh* cmesh) {
  TPZGeoMesh *gmesh = cmesh->Reference();
  if (!gmesh) DebugStop();
  int dim  = gmesh->Dimension();

  int64_t nels = gmesh->NElements();
  for (int64_t el = 0; el < nels; ++el) {
    bool flag = false;
    TPZGeoEl *gel = gmesh->Element(el);
    if (!gel || gel->Dimension() != dim || gel->HasSubElement()) continue;
    int firstside = gel->FirstSide(dim-1);
    int lastside = gel->FirstSide(dim);
    for (int side = firstside; side < lastside; ++side){
      TPZGeoElSide gelSide(gel,side);
      TPZGeoElSide neigh = gelSide.Neighbour();
      if (gelSide.HasNeighbour(EBack)) {
        flag = true;
      }
    }
    TPZVec<TPZGeoEl *> pv;
    if (flag == true) {
      gel->Divide(pv);
      break;
    }
  }
}

int64_t ErrorEstimation(TPZCompMesh* cmesh, TPZMultiphysicsCompMesh* cmeshMixed) {
  cmesh->LoadReferences();

  int64_t ncel = cmeshMixed->NElements();
  TPZAdmChunkVector<TPZCompEl *> &elementvec_m = cmeshMixed->ElementVec();

  // Initialize error storage
  REAL totalError = 0.0;
  int64_t imax = 0;
  REAL maxError = 0.0;
  
  // Number of integration points for flux evaluation
  const int numIntPts = 5; // Can be adjusted based on accuracy needs
  
  for (int64_t cel = 0; cel < ncel; ++cel) {

    TPZCompEl *celMixed = elementvec_m[cel];
    TPZGeoEl *gel = celMixed->Reference();
    TPZCompEl *celH1 = gel->Reference();

    if (!gel || gel->Dimension() != 2 || gel->HasSubElement()) continue;
    if (!celMixed || !celH1) continue;

    REAL hk = 1.0; // ElementDiameter(gel); // Element size
    
    TPZIntPoints *intrule = gel->CreateSideIntegrationRule(gel->NSides()-1, numIntPts);
    if (!intrule) continue;
      
    REAL fluxError = 0.0;
    REAL balanceError = 0.0;
    REAL elementError = 0.0;
    int npts = intrule->NPoints();
    
    for (int ip = 0; ip < npts; ++ip) {
      TPZManVector<REAL,3> ptInElement(gel->Dimension());
      REAL weight;
      intrule->Point(ip, ptInElement, weight);
      
      // Get Jacobian for integration
      TPZFMatrix<REAL> jacobian, axes, jacinv;
      REAL detjac;
      gel->Jacobian(ptInElement, jacobian, axes, detjac, jacinv);
      weight *= fabs(detjac);
      
      TPZManVector<REAL,3> x(3,0.0);
      gel->X(ptInElement, x); // Real coordinates for x

      TPZManVector<REAL,3> force(1,0.0);
      gexact.ForceFunc()(x, force); // Exact divergence (forcing function)

      // Compute flux from H1 method
      TPZVec<REAL> fluxH1(3,0.0);
      celH1->Solution(ptInElement, 7, fluxH1);
      
      // Get flux from Mixed method
      TPZVec<REAL> fluxMixed(3,0.0);
      TPZVec<REAL> divFluxMixed(1,0.0);
      TPZMultiphysicsElement *celMulti = dynamic_cast<TPZMultiphysicsElement*>(celMixed);
      if (celMulti) {
        celMulti->Solution(ptInElement, 1, fluxMixed);
        celMulti->Solution(ptInElement, 5, divFluxMixed);
      }
      
      // Difference between fluxes
      REAL fluxDiffNorm = 0.0;
      for (int d = 0; d < gel->Dimension(); ++d) {
        REAL diffFlux = fluxH1[d] - fluxMixed[d];
        fluxDiffNorm += diffFlux * diffFlux;
      }

      // Balance of flux
      REAL balanceDiffNorm = pow(divFluxMixed[0] - force[0], 2);

      // Add weighted contribution to element error
      fluxError += fluxDiffNorm * weight;
      balanceError += balanceDiffNorm * weight;
    }
    
    // For debugging purposes
    // std::cout << "Element " << el << ": Flux error = " << fluxError 
    //           << ", Balance error = " << balanceError << std::endl;

    elementError = pow((sqrt(fluxError) + (hk/M_PI)*sqrt(balanceError)), 2);
    totalError += elementError;

    if (elementError > maxError) {
      maxError = elementError;
      imax = gel->Index();
    }
    delete intrule;
  }
  
  totalError = sqrt(totalError);
  
  std::cout << "\nError estimation in energy norm: " << totalError << std::endl;

  return imax;
}