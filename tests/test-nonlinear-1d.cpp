
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
#include "TPZNonlinearWell.h"
#include "TPZNonLinearWellH1.h"
#include "TPZWannAnalysis.h"
#include "tpzgeoelrefpattern.h"

// ================
// Global variables
// ================

enum EnumMatIds
{
    EMatId = 1,
    ERight = 2,
    ELeft = 3,
    ENone = -1
};

int nthreads = 0;

// ===================
// Function prototypes
// ===================

// Creates a 1D geometric mesh
TPZGeoMesh *Create1DMesh(int nel, REAL x0, REAL x1);

// Creates a computational mesh for mixed approximation
TPZMultiphysicsCompMesh *createCompMeshMixed(TPZGeoMesh *gmesh, ProblemData *simData);

// Creates a computational mesh for H1 approximation
TPZCompMesh *createCompMeshH1(TPZGeoMesh *gmesh, ProblemData *simData);

// Nonlinear solver for the well problem
void SolveNonLinear(TPZGeoMesh *gmesh, ProblemData *simData);

// =============
// Main function
// =============

int main(int argc, char *const argv[])
{

// Initialize logger
#ifdef PZ_LOG
    TPZLogger::InitializePZLOG();
#endif

    // Set up the problem data
    ProblemData simData;
    simData.m_Numerics.maxIterations = 10;
    simData.m_Numerics.res_tol = 1.e-6;
    simData.m_Numerics.corr_tol = 1.e-6;
    simData.m_Numerics.nthreads = nthreads;

    simData.m_Wellbore.radius = 0.05;
    simData.m_Wellbore.length = 400.;
    simData.m_Wellbore.BCs["left"] = {ELeft, 1, 0.05};
    simData.m_Wellbore.BCs["right"] = {ERight, 1, 0.0};
    simData.m_Wellbore.pOrder = 2;

    simData.m_Fluid.viscosity = 0.001;
    simData.m_Fluid.density = 800;

    simData.m_Reservoir.BCs["pres"] = {ENone, 0, 2.2e7};
    simData.m_Reservoir.perm = 1.e-11; // Pseudo resistivity

    // Non-linear solver
    TPZGeoMesh *gmesh = Create1DMesh(20, 0., simData.m_Wellbore.length);
    SolveNonLinear(gmesh, &simData);
}

// =========
// Functions
// =========

TPZMultiphysicsCompMesh *createCompMeshMixed(TPZGeoMesh *gmesh, ProblemData *simData)
{

    int order = simData->m_Wellbore.pOrder;

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

    bcond = mat->CreateBC(mat, ERight, simData->m_Wellbore.BCs["right"].type, val1, val2);
    cmeshFlux->InsertMaterialObject(bcond);

    bcond = mat->CreateBC(mat, ELeft, simData->m_Wellbore.BCs["left"].type, val1, val2);
    cmeshFlux->InsertMaterialObject(bcond);

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
    cmesh->SetDefaultOrder(order);
    cmesh->ApproxSpace().Style() = TPZCreateApproximationSpace::EMultiphysics;

    // Problem parameters
    REAL Dw = simData->m_Wellbore.radius * 2;
    REAL mu = simData->m_Fluid.viscosity;
    REAL rho = simData->m_Fluid.density;
    REAL pres = simData->m_Reservoir.BCs["pres"].value;
    REAL Kvw = simData->m_Reservoir.perm;

    TPZNonlinearWell *matWell = new TPZNonlinearWell(EMatId, Dw, mu, rho, pres, Kvw);
    cmesh->InsertMaterialObject(matWell);

    // Create, set and add boundary conditions
    val2[0] = -simData->m_Wellbore.BCs["left"].value;
    bcond = matWell->CreateBC(matWell, ELeft, simData->m_Wellbore.BCs["left"].type, val1, val2);
    cmesh->InsertMaterialObject(bcond);

    val2[0] = simData->m_Wellbore.BCs["right"].value;
    bcond = matWell->CreateBC(matWell, ERight, simData->m_Wellbore.BCs["right"].type, val1, val2);
    cmesh->InsertMaterialObject(bcond);

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

TPZCompMesh *createCompMeshH1(TPZGeoMesh *gmesh, ProblemData *simData)
{
    int order = simData->m_Wellbore.pOrder;

    TPZCompMesh *cmeshH1 = new TPZCompMesh(gmesh);
    cmeshH1->SetDimModel(gmesh->Dimension());
    cmeshH1->SetDefaultOrder(order);
    cmeshH1->SetAllCreateFunctionsContinuous();

    // Problem parameters
    REAL Dw = simData->m_Wellbore.radius * 2;
    REAL mu = simData->m_Fluid.viscosity;
    REAL rho = simData->m_Fluid.density;
    REAL pres = simData->m_Reservoir.BCs["pres"].value;
    REAL Kvw = simData->m_Reservoir.perm;

    TPZNonLinearWellH1 *matWell = new TPZNonLinearWellH1(EMatId, Dw, mu, rho, pres, Kvw);
    cmeshH1->InsertMaterialObject(matWell);

    // Create boundary conditions
    TPZManVector<REAL, 1> val2(1, 0.); // Part that goes to the RHS vector
    TPZFMatrix<REAL> val1(1, 1, 0.);   // Part that goes to the Stiffnes matrix
    TPZBndCondT<REAL> *bcond;

    // For some reason its necessary to change the sign. TODO: look this up.
    val2[0] = -simData->m_Wellbore.BCs["left"].value;
    bcond = matWell->CreateBC(matWell, ELeft, simData->m_Wellbore.BCs["left"].type, val1, val2);
    cmeshH1->InsertMaterialObject(bcond);

    val2[0] = simData->m_Wellbore.BCs["right"].value;
    bcond = matWell->CreateBC(matWell, ERight, simData->m_Wellbore.BCs["right"].type, val1, val2);
    cmeshH1->InsertMaterialObject(bcond);

    cmeshH1->AutoBuild();
    return cmeshH1;
}

void SolveNonLinear(TPZGeoMesh *gmesh, ProblemData *simData)
{
    TPZMultiphysicsCompMesh *cmeshMixed = createCompMeshMixed(gmesh, simData);
    TPZCompMesh *cmeshH1 = createCompMeshH1(gmesh, simData);
    TPZWannAnalysis anMixed(cmeshMixed, RenumType::EMetis);
    TPZWannAnalysis anH1(cmeshH1, RenumType::EMetis);

    std::cout << "\nSolving non-linear well (Mixed H(div))" << std::endl;

    anMixed.SetProblemData(simData);
    anMixed.Initialize();
    anMixed.NewtonIteration();

    std::cout << "\nSolving non-linear well (H1)" << std::endl;

    anH1.SetProblemData(simData);
    anH1.Initialize();
    anH1.NewtonIteration();

    // --- Plotting ---

    {
        const std::string plotfile = "testWell_mixed";
        constexpr int vtkRes{0};
        TPZManVector<std::string, 2> fields = {"Flux", "Pressure", "Divergence"};
        auto vtk = TPZVTKGenerator(cmeshMixed, fields, plotfile, vtkRes);
        vtk.Do();
    }

    {
        const std::string plotfile = "testWell_H1";
        constexpr int vtkRes{0};
        TPZManVector<std::string, 2> fields = {"Flux", "Pressure"};
        auto vtk = TPZVTKGenerator(cmeshH1, fields, plotfile, vtkRes);
        vtk.Do();
    }

    // --- Clean up ---

    delete cmeshMixed;
    delete cmeshH1;
}

TPZGeoMesh *Create1DMesh(int nel, REAL x0, REAL x1)
{
    TPZGeoMesh *gmesh = new TPZGeoMesh;

    // Create nodes
    int nnodes = nel + 1;
    gmesh->NodeVec().Resize(nnodes);
    REAL dx = (x1 - x0) / nel;

    for (int i = 0; i < nnodes; i++)
    {
        TPZManVector<REAL, 3> coord(3, 0.0);
        coord[0] = x0 + i * dx; // x coordinate
        coord[1] = 0.0;         // y coordinate
        coord[2] = 0.0;         // z coordinate
        gmesh->NodeVec()[i].Initialize(coord, *gmesh);
    }

    // Create 1D elements
    for (int i = 0; i < nel; i++)
    {
        TPZManVector<int64_t, 2> nodes(2);
        nodes[0] = i;
        nodes[1] = i + 1;
        int64_t index;
        new TPZGeoElRefPattern<pzgeom::TPZGeoLinear>(nodes, EMatId, *gmesh, index);
    }

    // Create boundary conditions (points at the ends)
    TPZManVector<int64_t, 1> nodeBC(1);

    // Left boundary (x = x0)
    nodeBC[0] = 0;
    new TPZGeoElRefPattern<pzgeom::TPZGeoPoint>(nodeBC, ELeft, *gmesh); // BC material ID = -1

    // Right boundary (x = x1)
    nodeBC[0] = nnodes - 1;
    new TPZGeoElRefPattern<pzgeom::TPZGeoPoint>(nodeBC, ERight, *gmesh); // BC material ID = -2

    gmesh->SetDimension(1);
    gmesh->BuildConnectivity();
    return gmesh;
}

// =========================================
// Old versions of solver, keep just in case
// =========================================

// void SolveNonLinearOld(int order, TPZGeoMesh *gmesh)
// {

//     TPZMultiphysicsCompMesh *cmeshMixed = createCompMeshMixed(gmesh, simData, order);

//     // Mixed solver
//     TPZLinearAnalysis anMixed(cmeshMixed, RenumType::EMetis);
// #ifdef PZ_USING_MKL
//     TPZSSpStructMatrix<STATE> matMixed(cmeshMixed);
// #else
//     TPZFStructMatrix<STATE> matMixed(cmeshMixed);
// #endif
//     matMixed.SetNumThreads(nthreads);
//     anMixed.SetStructuralMatrix(matMixed);
//     TPZStepSolver<STATE> stepMixed;
//     stepMixed.SetDirect(ELDLt);
//     anMixed.SetSolver(stepMixed);

//     const int maxIter = 10;
//     const REAL tol = 1.e-6;
//     TPZFMatrix<STATE> Sol(cmeshMixed->NEquations(), 1, 0.);
//     for (int iteration = 0; iteration < maxIter; iteration++)
//     {
//         std::cout << "Non-linear iteration " << iteration << std::endl;
//         anMixed.Assemble();
//         // check convergence
//         if (iteration > 0)
//         {
//             TPZMatrix<STATE> &rhs = anMixed.Rhs();
//             REAL norm = 0.;
//             for (int i = 0; i < rhs.Rows(); i++)
//             {
//                 norm += rhs.GetVal(i, 0) * rhs.GetVal(i, 0);
//             }
//             norm = sqrt(norm);
//             std::cout << "norm = " << norm << std::endl;
//             if (norm < tol)
//             {
//                 std::cout << "Converged!" << std::endl;
//                 break;
//             }
//         }
//         anMixed.Solve();
//         TPZMatrix<STATE> &dsol = anMixed.Solution();
//         Sol += dsol;
//         // anMixed.LoadSolution(Sol);
//         cmeshMixed->LoadSolution(Sol);
//         cmeshMixed->TransferMultiphysicsSolution();
//     }

//     // ---- Plotting ---

//     {
//         const std::string plotfile = "darcy_mixed";
//         constexpr int vtkRes{0};
//         TPZManVector<std::string, 2> fields = {"Flux", "Pressure"};
//         auto vtk = TPZVTKGenerator(cmeshMixed, fields, plotfile, vtkRes);
//         vtk.Do();
//     }

//     // --- Clean up ---

//     delete cmeshMixed;
// }