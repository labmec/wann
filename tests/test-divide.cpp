
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
#include "pzstepsolver.h"
#include "TPZAnalyticSolution.h"
#include "pzlog.h"
#include "pzintel.h"

// ================
// Global variables
// ================

std::map<MMeshType, std::string> gMeshTypeNames = {
    {MMeshType::EQuadrilateral, "Quadrilateral"},
    {MMeshType::ETriangular, "Triangular"},
    {MMeshType::EHexahedral, "Hexahedral"},
    {MMeshType::ETetrahedral, "Tetrahedral"}
};

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
template <MMeshType elType>
TPZGeoMesh* createGeoMesh3D(
  const TPZManVector<int, 3> &nelDiv, 
  const TPZManVector<REAL, 3> &minX, 
  const TPZManVector<REAL, 3> &maxX);

// Creates a 2D geometric mesh using TPZGenGrid2D
template <MMeshType elType>
TPZGeoMesh* createGeoMesh2D(
  const TPZManVector<int, 2> &nelDiv, 
  const TPZManVector<REAL, 2> &minX, 
  const TPZManVector<REAL, 2> &maxX);

void Refine(TPZGeoMesh* gmesh);

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
  gRefDBase.InitializeUniformRefPattern(ETetraedro);
  gRefDBase.InitializeUniformRefPattern(ECube);
  
  const int nthreads = 0;
  int nref = 6;
  int dim = 2;
  int ndiv = 50;
  const auto eltype = MMeshType::EQuadrilateral;

  // Initial geometric mesh
  TPZGeoMesh* gmesh = nullptr;
  if (dim == 3) {
    gmesh = createGeoMesh3D<eltype>({ndiv, ndiv, ndiv}, {0., 0., 0.}, {1., 1., 1.});
  } else {
    gmesh = createGeoMesh2D<eltype>({ndiv, ndiv}, {0., 0.}, {1., 1.});
  }
  {
    std::ofstream out("gmesh0.vtk");
    TPZVTKGeoMesh::PrintGMeshVTK(gmesh, out);
  }

  std::cout << "Mesh dimension: " << gmesh->Dimension() << "\n";
  std::cout << "Number of elements: " << gmesh->NElements() << "\n";
  std::cout << "Element type: " << gMeshTypeNames[eltype] << "\n";
  for(int iref = 0; iref < nref; iref++) {
    
    auto start = std::chrono::high_resolution_clock::now();
    Refine(gmesh);
    auto end = std::chrono::high_resolution_clock::now();
    std::chrono::duration<REAL> diff = end - start;
    std::cout << "Refinement " << iref+1 << " completed in " << diff.count() << " s\n";
    {
      std::string filename = "gmesh" + std::to_string(iref+1) + ".vtk";
      std::ofstream out(filename);
      TPZVTKGeoMesh::PrintGMeshVTK(gmesh, out);
    }
  }
}

// =========
// Functions
// =========

template <MMeshType elType>
TPZGeoMesh* createGeoMesh3D(
  const TPZManVector<int, 3> &nelDiv, 
  const TPZManVector<REAL, 3> &minX, 
  const TPZManVector<REAL, 3> &maxX) {

  TPZGenGrid3D generator(minX, maxX, nelDiv, MMeshType::EHexahedral);

  generator.BuildVolumetricElements(EMatId);
  TPZGeoMesh *gmesh = generator.BuildBoundaryElements(EBottom, ELeft, EFront, ERight, EBack, ETop);

  return gmesh;
}

template <MMeshType elType>
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

void Refine(TPZGeoMesh* gmesh) {
  if (!gmesh) DebugStop();
  int dim  = gmesh->Dimension();

  int64_t nels = gmesh->NElements();
  for (int64_t el = 0; el < nels; ++el) {
    TPZGeoEl *gel = gmesh->Element(el);
    if (!gel || gel->Dimension() != dim || gel->HasSubElement()) continue;
    TPZVec<TPZGeoEl *> subels;
    gel->Divide(subels);
  }
}