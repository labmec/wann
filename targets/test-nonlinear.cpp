
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
#include "pzfstrmatrix.h"
#include "pzstepsolver.h"
#include "TPZAnalyticSolution.h"
#include "pzlog.h"
#include "pzintel.h"
#include "pzvec_extras.h"
#include "TPZWannMixedDarcyNL.h"

// ================
// Global variables
// ================

// Exact solution
TLaplaceExample1 gexact;

// Material IDs for domain and boundaries
enum EnumMatIds
{
    EMatId = 1,
    EBottom = 2,
    ERight = 3,
    ETop = 4,
    ELeft = 5,
    EFront = 6,
    EBack = 7
};

int nthreads = 0;

auto SetBoundaryCondition = [](TPZVec<REAL> x, TPZVec<REAL> u, TPZFMatrix<REAL> &du) {
    u[0] = x[1];
};

// ===================
// Function prototypes
// ===================

// Creates a geometric mesh using TPZGenGrid3D
TPZGeoMesh *createGeoMesh3D(
    const TPZManVector<int, 3> &nelDiv,
    const TPZManVector<REAL, 3> &minX,
    const TPZManVector<REAL, 3> &maxX);

// Creates a 2D geometric mesh using TPZGenGrid2D
TPZGeoMesh *createGeoMesh2D(
    const TPZManVector<int, 2> &nelDiv,
    const TPZManVector<REAL, 2> &minX,
    const TPZManVector<REAL, 2> &maxX);

// Creates a computational mesh for mixed approximation
TPZMultiphysicsCompMesh *createCompMeshMixed(TPZGeoMesh *gmesh, int order = 1, bool isNL = false);

void SolveLinear(int order, TPZGeoMesh *gmesh);

void SolveNonLinear(int order, TPZGeoMesh *gmesh);

// =============
// Main function
// =============

int main(int argc, char *const argv[])
{

// --- Set up ---

// Initialize logger
#ifdef PZ_LOG
    TPZLogger::InitializePZLOG();
#endif

    // Initialize uniform refinements for 1D and 2D elements
    gRefDBase.InitializeUniformRefPattern(EOned);
    gRefDBase.InitializeUniformRefPattern(EQuadrilateral);
    gRefDBase.InitializeUniformRefPattern(ETriangle);

    // --- Solve darcy problem ---

    int order = 1; // Polynomial order

    // Initial geometric mesh
    bool isNL = true; // Non-linear flag
    TPZGeoMesh *gmesh = createGeoMesh3D({3, 3, 3}, {0., 0., 0.}, {1., 1., 1.});
    // TPZGeoMesh* gmesh = createGeoMesh2D({3, 3}, {0., 0.}, {1., 1.});

    if (isNL)
    {
        std::cout << "Solving non-linear Darcy problem..." << std::endl;
        SolveNonLinear(order, gmesh);
    }
    else
    {
        std::cout << "Solving linear Darcy problem..." << std::endl;
        SolveLinear(order, gmesh);
    }
}

// =========
// Functions
// =========

TPZGeoMesh *createGeoMesh3D(
    const TPZManVector<int, 3> &nelDiv,
    const TPZManVector<REAL, 3> &minX,
    const TPZManVector<REAL, 3> &maxX)
{

    TPZGenGrid3D generator(minX, maxX, nelDiv, MMeshType::EHexahedral);

    generator.BuildVolumetricElements(EMatId);
    TPZGeoMesh *gmesh = generator.BuildBoundaryElements(EBottom, ELeft, EFront, ERight, EBack, ETop);

    return gmesh;
}

TPZGeoMesh *createGeoMesh2D(
    const TPZManVector<int, 2> &nelDiv,
    const TPZManVector<REAL, 2> &minX,
    const TPZManVector<REAL, 2> &maxX)
{

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

TPZMultiphysicsCompMesh *createCompMeshMixed(TPZGeoMesh *gmesh, int order, bool isNL)
{

    // --- Flux atomic cmesh ----

    TPZCompMesh *cmeshFlux = new TPZCompMesh(gmesh);
    cmeshFlux->SetDimModel(gmesh->Dimension());
    cmeshFlux->SetDefaultOrder(order);
    if (order < 1)
    {
        cmeshFlux->ApproxSpace().SetHDivFamily(HDivFamily::EHDivConstant);
    }
    cmeshFlux->SetAllCreateFunctionsHDiv();

    // Add materials (weak formulation)
    TPZNullMaterial<STATE> *mat = new TPZNullMaterial(EMatId, gmesh->Dimension());
    cmeshFlux->InsertMaterialObject(mat);

    // Create boundary conditions
    TPZManVector<REAL, 1> val2(1, 0.); // Part that goes to the RHS vector
    TPZFMatrix<REAL> val1(1, 1, 0.);   // Part that goes to the Stiffnes matrix
    TPZBndCondT<REAL> *bcond;

    bcond = mat->CreateBC(mat, ERight, 0, val1, val2);
    cmeshFlux->InsertMaterialObject(bcond);

    bcond = mat->CreateBC(mat, ELeft, 0, val1, val2);
    cmeshFlux->InsertMaterialObject(bcond);

    bcond = mat->CreateBC(mat, EFront, 0, val1, val2);
    cmeshFlux->InsertMaterialObject(bcond);

    bcond = mat->CreateBC(mat, EBack, 0, val1, val2);
    cmeshFlux->InsertMaterialObject(bcond);

    if (gmesh->Dimension() == 3)
    {
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
    if (order < 1)
    {
        cmeshPressure->SetAllCreateFunctionsDiscontinuous();
    }
    else
    {
        cmeshPressure->SetAllCreateFunctionsContinuous();
        cmeshPressure->ApproxSpace().CreateDisconnectedElements(true);
    }

    // Add materials (weak formulation)
    cmeshPressure->InsertMaterialObject(mat);

    // Set up the computational mesh
    cmeshPressure->AutoBuild();

    int ncon = cmeshPressure->NConnects();
    const int lagLevel = 1; // Lagrange multiplier level
    for (int i = 0; i < ncon; i++)
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
    TPZMixedDarcyFlow *matDarcy = nullptr;
    if (isNL)
    {
        TPZWannMixedDarcyNL *temp_mat = new TPZWannMixedDarcyNL(EMatId, gmesh->Dimension());
        temp_mat->SetConstantPermeability(1.0);
        matDarcy = temp_mat;
    }
    else
    {
        matDarcy->SetConstantPermeability(1.0);
    }
    cmesh->InsertMaterialObject(matDarcy);

    // Create, set and add boundary conditions
    val2[0] = 0.;
    bcond = matDarcy->CreateBC(matDarcy, ERight, 1, val1, val2);
    cmesh->InsertMaterialObject(bcond);
    // bcond->SetForcingFunctionBC(SetBoundaryCondition, 1);

    bcond = matDarcy->CreateBC(matDarcy, ELeft, 1, val1, val2);
    cmesh->InsertMaterialObject(bcond);
    // bcond->SetForcingFunctionBC(SetBoundaryCondition, 1);

    val2[0] = 0.;
    bcond = matDarcy->CreateBC(matDarcy, EFront, 0, val1, val2);
    cmesh->InsertMaterialObject(bcond);

    val2[0] = -1.;
    bcond = matDarcy->CreateBC(matDarcy, EBack, 1, val1, val2);
    cmesh->InsertMaterialObject(bcond);

    if (gmesh->Dimension() == 3)
    {
        val2[0] = 0.;
        bcond = matDarcy->CreateBC(matDarcy, EBottom, 1, val1, val2);
        cmesh->InsertMaterialObject(bcond);
        // bcond->SetForcingFunctionBC(SetBoundaryCondition, 1);

        bcond = matDarcy->CreateBC(matDarcy, ETop, 1, val1, val2);
        cmesh->InsertMaterialObject(bcond);
        // bcond->SetForcingFunctionBC(SetBoundaryCondition, 1);
    }

    // Incorporate the atomic meshes into the multiphysics mesh
    TPZManVector<TPZCompMesh *, 2> cmeshes(2);
    cmeshes[0] = cmeshFlux;
    cmeshes[1] = cmeshPressure;

    TPZManVector<int> active(cmeshes.size(), 1);
    cmesh->BuildMultiphysicsSpace(active, cmeshes);

    {
        std::ofstream out("cmeshMP.txt");
        cmesh->Print(out);
    }

    return cmesh;
}

void SolveNonLinear(int order, TPZGeoMesh *gmesh)
{

    TPZMultiphysicsCompMesh *cmeshMixed = createCompMeshMixed(gmesh, order, true);

    // Mixed solver
    TPZLinearAnalysis anMixed(cmeshMixed, RenumType::EMetis);
#ifdef PZ_USING_MKL
    TPZSSpStructMatrix<STATE> matMixed(cmeshMixed);
#else
    TPZFStructMatrix<STATE> matMixed(cmeshMixed);
#endif
    matMixed.SetNumThreads(nthreads);
    anMixed.SetStructuralMatrix(matMixed);
    TPZStepSolver<STATE> stepMixed;
    stepMixed.SetDirect(ELDLt);
    anMixed.SetSolver(stepMixed);

    const int maxIter = 10;
    const REAL tol = 1.e-6;
    TPZFMatrix<STATE> Sol(cmeshMixed->NEquations(), 1, 0.);
    for (int iteration = 0; iteration < maxIter; iteration++)
    {
        std::cout << "Non-linear iteration " << iteration << std::endl;
        anMixed.Assemble();
        // check convergence
        if (iteration > 0)
        {
            TPZMatrix<STATE> &rhs = anMixed.Rhs();
            REAL norm = 0.;
            for (int i = 0; i < rhs.Rows(); i++)
            {
                norm += rhs.GetVal(i, 0) * rhs.GetVal(i, 0);
            }
            norm = sqrt(norm);
            std::cout << "norm = " << norm << std::endl;
            if (norm < tol)
            {
                std::cout << "Converged!" << std::endl;
                break;
            }
        }
        anMixed.Solve();
        TPZMatrix<STATE> &dsol = anMixed.Solution();
        Sol += dsol;
        // anMixed.LoadSolution(Sol);
        cmeshMixed->LoadSolution(Sol);
        cmeshMixed->TransferMultiphysicsSolution();
    }

    // ---- Plotting ---

    {
        const std::string plotfile = "darcy_mixed";
        constexpr int vtkRes{0};
        TPZManVector<std::string, 2> fields = {"Flux", "Pressure"};
        auto vtk = TPZVTKGenerator(cmeshMixed, fields, plotfile, vtkRes);
        vtk.Do();
    }

    // --- Clean up ---

    delete cmeshMixed;
}

void SolveLinear(int order, TPZGeoMesh *gmesh)
{
    TPZMultiphysicsCompMesh *cmeshMixed = createCompMeshMixed(gmesh, order, false);

    // Mixed solver
    TPZLinearAnalysis anMixed(cmeshMixed, RenumType::EMetis);
#ifdef PZ_USING_MKL
    TPZSSpStructMatrix<STATE> matMixed(cmeshMixed);
#else
    TPZFStructMatrix<STATE> matMixed(cmeshMixed);
#endif
    matMixed.SetNumThreads(nthreads);
    anMixed.SetStructuralMatrix(matMixed);
    TPZStepSolver<STATE> stepMixed;
    stepMixed.SetDirect(ELDLt);
    anMixed.SetSolver(stepMixed);
    anMixed.Run();

    // ---- Plotting ---

    {
        const std::string plotfile = "darcy_mixed";
        constexpr int vtkRes{0};
        TPZManVector<std::string, 2> fields = {"Flux", "Pressure"};
        auto vtk = TPZVTKGenerator(cmeshMixed, fields, plotfile, vtkRes);
        vtk.Do();
    }

    // --- Clean up ---

    delete cmeshMixed;
}