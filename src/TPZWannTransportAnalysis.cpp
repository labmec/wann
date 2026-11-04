// Adapted from subFlow/src/TSFTransportAnalysis.cpp (LabMeC, MIT).
#include "TPZWannTransportAnalysis.h"
#include "TPZWannTransportDataTransfer.h"
#include "TPZSpStructMatrix.h"
#include "TPZSkylineNSymStructMatrix.h"
#include "TPZVTKGenerator.h"
#include "pzstepsolver.h"
#include "pzcmesh.h"
#include <cmath>
#include <iostream>
#include <stdexcept>

TPZWannTransportAnalysis::TPZWannTransportAnalysis(TPZCompMesh* mesh)
    : TPZLinearAnalysis(mesh, RenumType::ENone) {}

void TPZWannTransportAnalysis::Initialize(TPZMultiphysicsCompMesh* flow,
                                         const ProblemData& data) {
    TPZWannTransportDataTransfer::Initialize(flow, Mesh(), data, fAlgebraicTransport);
    fData = &data;
#ifdef PZ_USING_MKL
    TPZSpStructMatrix<STATE> structure(Mesh());
#else
    TPZSkylineNSymStructMatrix<STATE> structure(Mesh());
#endif
    structure.SetNumThreads(0);
    SetStructuralMatrix(structure);
    TPZStepSolver<STATE> solver;
    solver.SetDirect(ELU);
    solver.SetMatrix(fStructMatrix->Create());
    SetSolver(solver);
    const int64_t neq = Mesh()->NEquations();
    Rhs().Redim(neq, 1);
    fMassMatrix.Redim(neq, 1);
    auto& cells = fAlgebraicTransport.fCellsData;
    Solution().Redim(neq, 1);
    TPZFMatrix<STATE>& solution = Solution();
    for (size_t i = 0; i < cells.fVolume.size(); i++) {
        const int64_t equation = cells.fEqNumber[i];
        solution(equation, 0) = cells.fSaturation[i];
        fMassMatrix(equation, 0) = cells.fVolume[i] * cells.fPorosity[i];
    }
    LoadSolution();
}

void TPZWannTransportAnalysis::Assemble() {
    if (!fData) throw std::logic_error("Initialize transport before assembly");
    auto matrix = MatrixSolver<STATE>().Matrix();
    matrix->Zero();
    // Sparse Zero()/PutVal() do not invalidate NeoPZ's previous factorization.
    matrix->SetIsDecomposed(ENoDecompose);
    TPZFMatrix<STATE>& rhs = Rhs();
    rhs.Zero();
    // Unlike subFlow's raw sparse-array access, this also supports skyline LU.
    auto add = [&](int64_t row, int64_t col, STATE value) {
        if (matrix->PutVal(row, col, matrix->GetVal(row, col) + value) == 0)
            throw std::runtime_error("Missing transport matrix connection");
    };
    const auto& cells = fAlgebraicTransport.fCellsData;
    for (size_t cell = 0; cell < cells.fVolume.size(); cell++) {
        const int64_t eq = cells.fEqNumber[cell];
        TPZFNMatrix<1, STATE> ek(1, 1, 0.), ef(1, 1, 0.);
        fAlgebraicTransport.Contribute(cell, ek, ef);
        rhs(eq, 0) += ef(0, 0);
        add(eq, eq, fMassMatrix(eq, 0));
    }
    const REAL dt = fData->m_Transport->timeStep;
    const auto boundaries = fData->TransportBoundarySaturations();
    for (const auto& group : fAlgebraicTransport.fInterfaceData) {
        const int matid = group.first;
        const auto& faces = group.second;
        const REAL exteriorFractionalFlow = matid == fData->ETransportInterface
            ? 0. : fData->FractionalFlow(boundaries.at(matid)).first;
        for (size_t face = 0; face < faces.fIntegralFlux.size(); face++) {
            const auto lr = faces.fLeftRightVolIndex[face];
            const int64_t left = cells.fEqNumber[lr.first];
            if (lr.second >= 0) {
                const int64_t right = cells.fEqNumber[lr.second];
                TPZFNMatrix<4, STATE> ek(2, 2, 0.);
                TPZFNMatrix<2, STATE> ef(2, 1, 0.);
                fAlgebraicTransport.ContributeInterface(face, matid, dt, ek, ef);
                rhs(left, 0) += ef(0, 0);
                rhs(right, 0) += ef(1, 0);
                add(left, left, ek(0, 0)); add(left, right, ek(0, 1));
                add(right, left, ek(1, 0)); add(right, right, ek(1, 1));
            } else {
                TPZFNMatrix<1, STATE> ek(1, 1, 0.), ef(1, 1, 0.);
                fAlgebraicTransport.ContributeBC(face, matid, dt, exteriorFractionalFlow, ek, ef);
                rhs(left, 0) += ef(0, 0);
                add(left, left, ek(0, 0));
            }
        }
    }
}

void TPZWannTransportAnalysis::RunTimeStep() {
    if (!fData) throw std::logic_error("Initialize transport before stepping");
    fData->ValidateTransport();
    const auto& newton = fData->m_Numerics.NewtonTransport;
    auto& cells = fAlgebraicTransport.fCellsData;
    const TPZFMatrix<STATE> accepted = Solution();
    TPZFMatrix<STATE> iterate = accepted;
    auto load = [&](const TPZFMatrix<STATE>& state) {
        fSolution = state;
        LoadSolution();
        cells.UpdateSaturations(state);
        cells.UpdateFractionalFlows(*fData);
    };
    try {
        for (int iteration = 0; iteration <= newton.maxIterations; iteration++) {
            Assemble();
            const REAL residual = Norm(Rhs());
            if (newton.printNewtonProgress)
                std::cout << "Transport Newton iteration " << iteration
                          << ", residual norm = " << residual << '\n';
            if (std::isfinite(residual) && residual < newton.res_tol) {
                cells.fSaturationLastState = cells.fSaturation;
                if (newton.printNewtonProgress)
                    std::cout << "Transport Newton converged after " << iteration
                              << " corrections\n";
                return;
            }
            if (iteration == newton.maxIterations || !std::isfinite(residual)) break;
            Solve();
            const TPZFMatrix<STATE> correction = Solution();
            // Full corrections in subFlow are supplemented with bounded backtracking.
            // Previous-time saturation remains fixed during every trial/iteration.
            bool found = false;
            for (REAL alpha = 1.; alpha >= 1.e-8; alpha *= 0.5) {
                TPZFMatrix<STATE> trial = iterate + alpha * correction;
                bool bounded = true;
                for (int64_t i = 0; i < trial.Rows(); i++)
                    if (!std::isfinite(trial(i, 0)) || trial(i, 0) < -1.e-12 || trial(i, 0) > 1. + 1.e-12)
                        bounded = false;
                if (!bounded) continue;
                load(trial);
                Assemble();
                const REAL trialResidual = Norm(Rhs());
                if (std::isfinite(trialResidual) && trialResidual < residual) {
                    iterate = trial;
                    found = true;
                    if (newton.printNewtonProgress)
                        std::cout << "  Accepted correction norm = " << alpha * Norm(correction)
                                  << ", alpha = " << alpha << '\n';
                    break;
                }
            }
            if (!found) break;
        }
    } catch (...) {
        load(accepted);
        throw;
    }
    load(accepted);
    if (newton.printNewtonProgress)
        std::cout << "Transport Newton failed to converge. Previous state restored.\n";
    throw std::runtime_error("Transport Newton iteration failed; previous state restored");
}

void TPZWannTransportAnalysis::PostProcessTimeStep(const std::string& filename, int step, bool safeCheck) {
    TPZStack<std::string> fields;
    fields.Push("Sw"); fields.Push("So");
    TPZVTKGenerator vtk(Mesh(), fields, filename, 0, Mesh()->Dimension(), safeCheck);
    vtk.SetStep(step);
    vtk.SetNThreads(0);
    vtk.Do();
}
