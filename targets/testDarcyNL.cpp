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

// Wann includes
#include "TPZWannMixedDarcyNL.h"
#include "TPZWannDarcyNL.h"
#include "TPZWannAnalysis.h"
#include "TPZWannGeometryTools.h"
#include "ProblemData.h"

// ===================
// Function prototypes
// ===================

// Creates a computational mesh for mixed approximation
TPZMultiphysicsCompMesh *createCompMeshMixed(TPZGeoMesh *gmesh, ProblemData *SimData, bool isNL = false);

// Creates a computational mesh for H1 approximation
TPZCompMesh *createCompMeshH1(TPZGeoMesh *gmesh, ProblemData *SimData, bool isNL = false);

// Non-linear solvers
void SolveNonLinear(ProblemData *SimData, TPZGeoMesh *gmesh);

// =============
// Main function
// =============

int main(int argc, char *const argv[]) {

// Initialize logger
#ifdef PZ_LOG
  TPZLogger::InitializePZLOG();
#endif

  // Problem data
  ProblemData simData;
  simData.m_Numerics.maxIterations = 10;
  simData.m_Numerics.res_tol = 1.e-6; 
  simData.m_Numerics.corr_tol = 1.e-6;
  simData.m_Numerics.nthreads = 0;

  simData.m_Reservoir.BCs["right"] = {ERight, 1, 0.0};
  simData.m_Reservoir.BCs["left"] = {ELeft, 1, 0.0};
  simData.m_Reservoir.BCs["front"] = {EFront, 0, 3.0};
  simData.m_Reservoir.BCs["back"] = {EBack, 1, 1.0};

  TPZGeoMesh *gmesh = TPZWannGeometryTools::readGeoMesh(&simData);

  int order = simData.m_Reservoir.pOrder;
  SolveNonLinear(&simData, gmesh);
}

// =========
// Functions
// =========

TPZMultiphysicsCompMesh *createCompMeshMixed(TPZGeoMesh *gmesh, ProblemData *SimData, bool isNL)
{

    // --- Flux atomic cmesh ----

    TPZCompMesh *cmeshFlux = new TPZCompMesh(gmesh);
    cmeshFlux->SetDimModel(gmesh->Dimension());
    cmeshFlux->SetDefaultOrder(SimData->m_Reservoir.pOrder);
    if (SimData->m_Reservoir.pOrder < 1)
    {
        cmeshFlux->ApproxSpace().SetHDivFamily(HDivFamily::EHDivConstant);
    }
    cmeshFlux->SetAllCreateFunctionsHDiv();

    // Add materials (weak formulation)
    TPZNullMaterial<STATE> *mat = new TPZNullMaterial(SimData->EDomain, gmesh->Dimension());
    cmeshFlux->InsertMaterialObject(mat);

    // Create boundary conditions
    TPZManVector<REAL, 1> val2(1, 0.); // Part that goes to the RHS vector
    TPZFMatrix<REAL> val1(1, 1, 0.);   // Part that goes to the Stiffnes matrix
    TPZBndCondT<REAL> *bcond;

    bcond = mat->CreateBC(mat, SimData->EFarfield, 0, val1, val2);
    cmeshFlux->InsertMaterialObject(bcond);

    bcond = mat->CreateBC(mat, SimData->ESurfHeel, 0, val1, val2);
    cmeshFlux->InsertMaterialObject(bcond);

    bcond = mat->CreateBC(mat, SimData->ESurfToe, 0, val1, val2);
    cmeshFlux->InsertMaterialObject(bcond);

    bcond = mat->CreateBC(mat, SimData->ESurfWellCyl, 0, val1, val2);
    cmeshFlux->InsertMaterialObject(bcond);

    cmeshFlux->AutoBuild();

    {
        std::ofstream out("cmeshFlux.txt");
        cmeshFlux->Print(out);
    }

    // --- Pressure atomic cmesh ---

    TPZCompMesh *cmeshPressure = new TPZCompMesh(gmesh);
    cmeshPressure->SetDimModel(gmesh->Dimension());
    cmeshPressure->SetDefaultOrder(SimData->m_Reservoir.pOrder);
    if (SimData->m_Reservoir.pOrder < 1)
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

    // --- Multiphysics mesh ---

    TPZMultiphysicsCompMesh *cmesh = new TPZMultiphysicsCompMesh(gmesh);
    cmesh->SetDimModel(gmesh->Dimension());
    cmesh->SetDefaultOrder(1);
    cmesh->ApproxSpace().Style() = TPZCreateApproximationSpace::EMultiphysics;

    // Add materials (weak formulation)
    TPZMixedDarcyFlow *matDarcy = nullptr;
    if (isNL)
    {
        TPZWannMixedDarcyNL *temp_mat = new TPZWannMixedDarcyNL(SimData->EDomain, gmesh->Dimension());
        matDarcy = temp_mat;
    }
    else
    {
        TPZMixedDarcyFlow *temp_mat = new TPZMixedDarcyFlow(SimData->EDomain, gmesh->Dimension());
        matDarcy = temp_mat;
    }
    matDarcy->SetConstantPermeability(1.0);
    cmesh->InsertMaterialObject(matDarcy);

    // Create, set and add boundary conditions
    val2[0] = 3.;
    bcond = matDarcy->CreateBC(matDarcy, SimData->EFarField, 1, val1, val2);
    cmesh->InsertMaterialObject(bcond);

    val2[0] = 0.;
    bcond = matDarcy->CreateBC(matDarcy, SimData->ESurfHeel, 1, val1, val2);
    cmesh->InsertMaterialObject(bcond);

    val2[0] = 0.;
    bcond = matDarcy->CreateBC(matDarcy, SimData->ESurfToe, 1, val1, val2);
    cmesh->InsertMaterialObject(bcond);

    val2[0] = 1.0;
    bcond = matDarcy->CreateBC(matDarcy, SimData->ESurfWellCyl, 0, val1, val2);
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

TPZCompMesh *createCompMeshH1(TPZGeoMesh *gmesh, ProblemData *SimData, bool isNL) 
{
    
    TPZCompMesh *cmesh = new TPZCompMesh(gmesh);
    cmesh->SetDimModel(gmesh->Dimension());
    cmesh->SetDefaultOrder(SimData->m_Reservoir.pOrder);
    cmesh->SetAllCreateFunctionsContinuous();

    // Add materials (weak formulation)
    TPZDarcyFlow *matDarcy = nullptr;
    if (isNL)
    {
        TPZWannDarcyNL *temp_mat = new TPZWannDarcyNL(SimData->EDomain, gmesh->Dimension());
        matDarcy = temp_mat;
    }
    else
    {
        TPZDarcyFlow *temp_mat = new TPZDarcyFlow(SimData->EDomain, gmesh->Dimension());
        matDarcy = temp_mat;
    }
    matDarcy->SetConstantPermeability(1.0);
    cmesh->InsertMaterialObject(matDarcy);

    // Create, set and add boundary conditions
    TPZManVector<REAL, 1> val2(1, 0.); // Part that goes to the RHS vector
    TPZFMatrix<REAL> val1(1, 1, 0.);   // Part that goes to the Stiffnes matrix
    TPZBndCondT<REAL> *bcond;

    val2[0] = 3.;
    bcond = matDarcy->CreateBC(matDarcy, SimData->EFarField, 0, val1, val2);
    cmesh->InsertMaterialObject(bcond);

    val2[0] = 1.;
    bcond = matDarcy->CreateBC(matDarcy, SimData->ESurfWellCyl, 0, val1, val2);
    cmesh->InsertMaterialObject(bcond);

    val2[0] = 0.;
    bcond = matDarcy->CreateBC(matDarcy, SimData->ESurfHeel, 1, val1, val2);
    cmesh->InsertMaterialObject(bcond);

    val2[0] = 0.0;
    bcond = matDarcy->CreateBC(matDarcy, SimData->ESurfToe, 1, val1, val2);
    cmesh->InsertMaterialObject(bcond);

    cmesh->AutoBuild();
    return cmesh;
}

void SolveLinear(int order, TPZGeoMesh *gmesh)
{
    TPZMultiphysicsCompMesh *cmeshMixed = createCompMeshMixed(gmesh, order, false);
    TPZCompMesh *cmeshH1 = createCompMeshH1(gmesh, order, false);

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

    // H1 solver
    TPZLinearAnalysis anH1(cmeshH1, RenumType::EMetis);
#ifdef PZ_USING_MKL
    TPZSSpStructMatrix<STATE> matH1(cmeshH1);
#else
    TPZFStructMatrix<STATE> matH1(cmeshH1);
#endif
    matH1.SetNumThreads(nthreads);
    anH1.SetStructuralMatrix(matH1);
    TPZStepSolver<STATE> stepH1;
    stepH1.SetDirect(ECholesky);
    anH1.SetSolver(stepH1);
    anH1.Run();

    // ---- Plotting ---

    {
        const std::string plotfile = "Darcy_Mixed";
        constexpr int vtkRes{0};
        TPZManVector<std::string, 2> fields = {"Flux", "Pressure"};
        auto vtk = TPZVTKGenerator(cmeshMixed, fields, plotfile, vtkRes);
        vtk.Do();
    }

    {
        const std::string plotfile = "Darcy_H1";
        constexpr int vtkRes{0};
        TPZManVector<std::string, 2> fields = {"Flux", "Pressure"};
        auto vtk = TPZVTKGenerator(cmeshH1, fields, plotfile, vtkRes);
        vtk.Do();
    }

    // --- Clean up ---

    delete cmeshMixed;
    delete cmeshH1;
}

void SolveNonLinear(ProblemData *SimData, TPZGeoMesh *gmesh)
{
    TPZMultiphysicsCompMesh *cmeshMixed = createCompMeshMixed(gmesh, SimData, true);
    TPZCompMesh *cmeshH1 = createCompMeshH1(gmesh, SimData, true);

    {
        std::ofstream out("cmeshH1.txt");
        cmeshH1->Print(out);
    }

    // Mixed solver
    TPZWannAnalysis anMixed(cmeshMixed, RenumType::EMetis);
    anMixed.SetProblemData(SimData);
    anMixed.Initialize();
    anMixed.NewtonIteration();

    // H1 solver
    TPZWannAnalysis anH1(cmeshH1, RenumType::EMetis);
    anH1.SetProblemData(SimData);
    anH1.Initialize();

    {
        std::ofstream out("cmeshH1.txt");
        cmeshH1->Print(out);
    }

    anH1.NewtonIteration();

    // ---- Plotting ---

    {
        const std::string plotfile = "DarcyNL_Mixed";
        constexpr int vtkRes{0};
        TPZManVector<std::string, 2> fields = {"Flux", "Pressure"};
        auto vtk = TPZVTKGenerator(cmeshMixed, fields, plotfile, vtkRes);
        vtk.Do();
    }

    {
        const std::string plotfile = "DarcyNL_H1";
        constexpr int vtkRes{0};
        TPZManVector<std::string, 2> fields = {"Flux", "Pressure"};
        auto vtk = TPZVTKGenerator(cmeshH1, fields, plotfile, vtkRes);
        vtk.Do();
    }

    // --- Clean up ---

    delete cmeshMixed;
    delete cmeshH1;
}

// =========
// Old stuff
// =========

void SolveNonLinearOld(int order, TPZGeoMesh *gmesh)
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