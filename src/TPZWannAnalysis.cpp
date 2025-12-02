//
//  Created by Giovane Avancini and Giovanni Taraschi on 17/11/25.
//

#include "TPZWannAnalysis.h"
#include "TPZNonlinearWell.h"

using namespace std;

TPZWannAnalysis::TPZWannAnalysis() : TPZLinearAnalysis() {}

TPZWannAnalysis::TPZWannAnalysis(TPZMultiphysicsCompMesh *cmesh,
                                 const RenumType &renumtype)
    : TPZLinearAnalysis(cmesh, renumtype) {}

TPZWannAnalysis::~TPZWannAnalysis() {}

void TPZWannAnalysis::SetProblemData(ProblemData *simData)
{
    fSimData = simData;
}

ProblemData *TPZWannAnalysis::GetProblemData() { return fSimData; }

int TPZWannAnalysis::GetNumberOfIterations() { return fKiteration; }

void TPZWannAnalysis::Initialize()
{
#ifdef PZ_USING_MKL
    TPZSSpStructMatrix<STATE> matrix(fCompMesh);
#else
    TPZSkylineStructMatrix<STATE> matrix(fCompMesh);
#endif
    int n_threads = fSimData->m_Numerics.nthreads;
    matrix.SetNumThreads(n_threads);
    SetStructuralMatrix(matrix);
    std::set<int> neumannMatids;
    FillNeumannBCMatids(neumannMatids);
    SetInitialSolution(neumannMatids);
    ApplyEquationFilter(neumannMatids);
    int nreducedeq = fStructMatrix->NReducedEquations();
    TPZStepSolver<STATE> step;
    step.SetDirect(ELDLt);
    SetSolver(step);

    std::cout << "Number of equations: " << fCompMesh->NEquations() << std::endl;
    std::cout << "Number of elements: " << fCompMesh->NElements() << std::endl;
}

void TPZWannAnalysis::NewtonIteration()
{
    TPZMultiphysicsCompMesh *cmesh =
        dynamic_cast<TPZMultiphysicsCompMesh *>(Mesh());
    if (!cmesh)
        DebugStop();

    int matIter = fSimData->m_Numerics.maxIterations;
    REAL res_norm = 1.0;
    REAL corr_norm = 1.0;
    REAL res_tol = fSimData->m_Numerics.res_tol;
    REAL corr_tol = fSimData->m_Numerics.corr_tol;
    bool converged = false;

    TPZFMatrix<STATE> sol = Solution();
    for (fKiteration = 0; fKiteration < matIter; fKiteration++)
    {
        Assemble();

        // Check residual convergence
        if (fKiteration > 0)
        {
            TPZFMatrix<STATE> rhs = Rhs();
            res_norm = Norm(rhs);
            std::cout << "------Newton iteration: " << fKiteration << std::endl;
            std::cout << "---------Residual norm: " << res_norm << std::endl;
            std::cout << "---------Correction norm: " << corr_norm << std::endl;
            if (res_norm < res_tol || corr_norm < corr_tol)
            {
                std::cout << "------Iterative method converged with res_norm: "
                          << res_norm << std::endl;
                std::cout << "------Number of iterations = " << fKiteration
                          << std::endl;
                converged = true;
                fSolution = sol;
                break;
            }
        }
        Solve();
        TPZFMatrix<STATE> dsol = Solution();
        REAL eta = LineSearchStep(sol, dsol, ENONE);
        corr_norm = eta * Norm(dsol);
        sol += eta * dsol;
        cmesh->LoadSolution(sol);
        cmesh->TransferMultiphysicsSolution();

        PostProcessIteration(cmesh->Dimension(), fKiteration);
    }

    if (!converged)
    {
        std::cout << "------Iterative method did not converge. res_norm: "
                  << res_norm << " corr_norm: " << corr_norm << std::endl;
    }
}

REAL TPZWannAnalysis::LineSearchStep(TPZFMatrix<STATE> &sol, TPZFMatrix<STATE> &dsol, LineSearchMethod method)
{
    REAL eta = 1.0;
    switch (method)
    {
    case ENONE:
        eta = 1.0;
        break;
    case EBISSECT:
        // To be implemented
        eta = BissectionMethod(sol, dsol);
        break;
    case EGOLDEN:
        // To be implemented
        DebugStop();
        break;
    case ESECANT:
        // To be implemented
        DebugStop();
        break;
    }
    return eta;
}

REAL TPZWannAnalysis::BissectionMethod(TPZFMatrix<STATE> &sol_n, TPZFMatrix<STATE> &dsol)
{
    TPZNonlinearWell::fAssembleRHSOnly = true;

    TPZFMatrix<STATE> &res_L = Rhs(); // lower point
    TPZFMatrix<STATE> aux(1, 1, 0.0);
    dsol.MultAdd(res_L, res_L, aux, 1.0, 0.0, 1);
    REAL s_L = aux(0, 0);

    TPZFMatrix<STATE> sol_n1 = sol_n;
    sol_n1 += dsol; // upper point
    Mesh()->LoadSolution(sol_n1);
    Mesh()->TransferMultiphysicsSolution();
    Assemble(); // Here, assemble only computes the RHS (residual)
    TPZFMatrix<STATE> &res_U = Rhs();
    dsol.MultAdd(res_U, res_U, aux, 1.0, 0.0, 1);
    REAL s_U = aux(0, 0);

    REAL eta = 1.0;
    if (std::abs(s_U) < 1e-6 || s_L * s_U > 0)
    {
        return eta; // upper point is root
    }

    eta = 0.5;
    REAL deta = 0.25;
    REAL s_M = 1.0;
    while (std::abs(s_M) > 1e-6 && deta > 1e-3)
    {
        sol_n1 = sol_n + eta * dsol; // mid point
        Mesh()->LoadSolution(sol_n1);
        Mesh()->TransferMultiphysicsSolution();
        Assemble(); // Here, assemble only computes the RHS (residual)
        TPZFMatrix<STATE> &res_M = Rhs();
        dsol.MultAdd(res_M, res_M, aux, 1.0, 0.0, 1);
        s_M = aux(0, 0);

        if (s_L * s_M < 0)
        {
            s_U = s_M;
            eta -= deta;
        }
        else if (s_U * s_M < 0)
        {
            s_L = s_M;
            eta += deta;
        }
        else
        {
            DebugStop(); // Something is wrong
        }

        std::cout << "Bissection step eta: " << eta << " s_M: " << s_M << std::endl;

        deta *= 0.5;
    }
    TPZNonlinearWell::fAssembleRHSOnly = false;
    return eta;
}

void TPZWannAnalysis::Assemble()
{
    auto start_time_ass = std::chrono::steady_clock::now();

    TPZLinearAnalysis::Assemble();

    auto total_time_ass = std::chrono::duration_cast<std::chrono::milliseconds>(
                              std::chrono::steady_clock::now() - start_time_ass)
                              .count() /
                          1000.;
    std::cout << "---------Time to assemble: " << total_time_ass << " seconds"
              << std::endl;
}

void TPZWannAnalysis::Solve()
{
    auto start_time_solve = std::chrono::steady_clock::now();

    TPZLinearAnalysis::Solve();

    auto total_time_solve =
        std::chrono::duration_cast<std::chrono::milliseconds>(
            std::chrono::steady_clock::now() - start_time_solve)
            .count() /
        1000.;
    std::cout << "---------Time to solve: " << total_time_solve << " seconds"
              << std::endl;
}

void TPZWannAnalysis::FillNeumannBCMatids(std::set<int> &neumannMatids)
{
    for (auto it = fSimData->m_Reservoir.BCs.begin(); it != fSimData->m_Reservoir.BCs.end(); ++it)
    {
        int matid = it->second.matid;
        int bc_type = it->second.type;
        if (bc_type == 1)
            neumannMatids.insert(matid);
    }
    for (auto it = fSimData->m_Wellbore.BCs.begin(); it != fSimData->m_Wellbore.BCs.end(); ++it)
    {
        int matid = it->second.matid;
        int bc_type = it->second.type;
        if (bc_type == 1)
            neumannMatids.insert(matid);
    }
}

void TPZWannAnalysis::SetInitialSolution(std::set<int> &neumannMatids)
{
    fCompMesh->LoadReferences();
    TPZGeoMesh *gmesh = fCompMesh->Reference();

    TPZFMatrix<STATE> &cmesh_sol = fCompMesh->Solution();
    for (auto el : gmesh->ElementVec())
    {
        int elMatID = el->MaterialId();

        if (neumannMatids.find(elMatID) == neumannMatids.end())
            continue;

        TPZCompEl *compEl = el->Reference();

        int64_t nConnects = compEl->NConnects();

        if (nConnects != 1)
            DebugStop();

        int64_t seq = compEl->Connect(0).SequenceNumber();
        int ncorner = el->NCornerNodes();
        REAL volume = el->Dimension() == 0 ? 1.0 : el->Volume();

        auto firstEq = fCompMesh->Block().Position(seq);

        int64_t blockSize = fCompMesh->Block().Size(seq);

        std::string bc_name;
        REAL val;
        bool found = false;
        for (auto it = fSimData->m_Reservoir.BCs.begin(); it != fSimData->m_Reservoir.BCs.end(); ++it)
        {
            if (it->second.matid == elMatID)
            {
                bc_name = it->first;
                val = fSimData->m_Reservoir.BCs[bc_name].value;
                found = true;
                break;
            }
        }
        if (!found)
        {
            for (auto it = fSimData->m_Wellbore.BCs.begin(); it != fSimData->m_Wellbore.BCs.end(); ++it)
            {
                if (it->second.matid == elMatID)
                {
                    bc_name = it->first;
                    val = fSimData->m_Wellbore.BCs[bc_name].value;
                    found = true;
                    break;
                }
            }
        }
        if (!found)
            DebugStop();
        val *= volume / ncorner;
        for (int64_t eq = firstEq; eq < firstEq + blockSize; eq++)
        {
            if (eq - firstEq < ncorner)
                cmesh_sol.PutVal(eq, 0, val);
        }
    }

    fCompMesh->TransferMultiphysicsSolution();

    // When the internal dofs are condensed, the analysis solution size is
    // different from the cmesh solution size Analysis only holds the independent
    // equations, while cmesh holds all equations The independent equations are
    // stored first in the cmesh solution vector, so we just need to copy them to
    // the analysis solution
    int cmesh_neq = fCompMesh->NEquations();
    TPZFMatrix<STATE> &sol = Solution();
    for (int i = 0; i < cmesh_neq; i++)
    {
        sol.PutVal(i, 0, cmesh_sol.GetVal(i, 0));
    }
}

void TPZWannAnalysis::ApplyEquationFilter(std::set<int> &neumannMatids)
{
    fCompMesh->LoadReferences();
    std::set<int64_t> removeEquations;
    TPZGeoMesh *gmesh = fCompMesh->Reference();
    TPZFMatrix<STATE> sol = fCompMesh->Solution();
    for (auto el : gmesh->ElementVec())
    {
        int elMatID = el->MaterialId();

        if (neumannMatids.find(elMatID) == neumannMatids.end())
            continue;

        TPZCompEl *compEl = el->Reference();

        int64_t nConnects = compEl->NConnects();

        if (nConnects != 1)
            DebugStop();

        int64_t seq = compEl->Connect(0).SequenceNumber();
        int ncorner = el->NCornerNodes();

        auto firstEq = fCompMesh->Block().Position(seq);

        int64_t blockSize = fCompMesh->Block().Size(seq);

        for (int64_t eq = firstEq; eq < firstEq + blockSize; eq++)
        {
            removeEquations.insert(eq);
        }
    }

    TPZEquationFilter filter(fCompMesh->NEquations());
    filter.ExcludeEquations(removeEquations);
    fStructMatrix->EquationFilter() = filter;
}

void TPZWannAnalysis::PostProcessIteration(int dimToPost, int it)
{
    auto start_time_pp = std::chrono::steady_clock::now();

    TPZStack<std::string, 10> scalnames;

    scalnames.Push("Pressure");
    scalnames.Push("Flux");

    int div = 0;
    if (dimToPost < 0)
    {
        dimToPost = Mesh()->Reference()->Dimension();
    }

    std::string file = "solution_it" + std::to_string(it) + ".vtk";
    const int vtkRes = fSimData->m_PostProc.vtk_resolution;

    const std::string plotfile = file.substr(0, file.find(".")); // sem o .vtk no final
    auto vtk = TPZVTKGenerator(fCompMesh, scalnames, plotfile, vtkRes, dimToPost);
    vtk.SetStep(it);
    int nthreads = fSimData->m_PostProc.nthreads;
    vtk.SetNThreads(nthreads);
    vtk.Do();

    auto total_time_pp = std::chrono::duration_cast<std::chrono::milliseconds>(std::chrono::steady_clock::now() - start_time_pp).count() / 1000.;
    cout << "Total time post process = " << total_time_pp << " seconds" << endl;
}