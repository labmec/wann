/**
 * @file
 * @brief Class and method contracts for TPZWannAnalysis.
 */

//
//  Created by Giovane Avancini and Giovanni Taraschi on 17/11/25.
//

#pragma once

#include "ProblemData.h"
#include "TPZLinearAnalysis.h"
#include "TPZMultiphysicsCompMesh.h"
#include "TPZSSpStructMatrix.h"
#include "pzskylstrmatrix.h"
#include "pzstepsolver.h"
#include <stdio.h>
#include "TPZVTKGenerator.h"
#include "TPZSYSMPMatrix.h"

/**
 * @brief Newton driver and boundary-equation filtering for WANN computational meshes.
 * Requires a mesh and a borrowed ProblemData pointer before Initialize().
 * Mixed meshes filter prescribed normal-flux equations; H1 meshes filter
 * prescribed pressure equations. NewtonIteration() runs the entire iteration
 * loop, loading corrections into the mesh and transferring atomic solutions.
 * Line-search routines change mesh/RHS state and use shared well-material flags.
 */
class TPZWannAnalysis : public TPZLinearAnalysis
{

public:
    /**
     * @brief Correction scaling methods; ESECANT is declared but unimplemented.
     */
    enum LineSearchMethod
    {
        ENONE,
        EBISSECTION,
        EGOLDENRATIO,
        ESECANT
    };

protected:
    /**
     * @brief Borrowed configuration pointer set by SetProblemData().
     */
    ProblemData *fSimData;

    /**
     * @brief Current Newton loop index, or iteration limit on exhaustion.
     */
    int fKiteration = 0;

    /**
     * @brief Reserved flag; currently not consulted by the implementation.
     */
    bool fIsFirstAssemble = true;

    /**
     * @brief Whether the attached mesh is TPZMultiphysicsCompMesh.
     */
    bool fIsMultiphysics = true;

    /**
     * @brief Method used to scale Newton corrections after the first step.
     */
    LineSearchMethod fLineSearchMethod = ENONE;

public:

    /**
     * @brief Construct an object; configure simulation/material data before use.
     */
    TPZWannAnalysis();

    /**
     * @brief Attach the computational mesh and detect mixed versus H1 operation.
     * @param[in] cmesh Mesh used by the NeoPZ analysis.
     * @param[in] renumtype Equation renumbering strategy.
     */
    TPZWannAnalysis(TPZCompMesh *cmesh,
                    const RenumType &renumtype = RenumType::EDefault);

    /**
     * @brief Destroy this object; borrowed simulation pointers are not deleted.
     */
    ~TPZWannAnalysis();

    /**
     * @brief Configure threaded assembly, initialize prescribed BC values, filter their
     * equations, and choose LDLt (mixed) or Cholesky (H1) direct solution.
     * Requires SetProblemData() and a valid computational mesh beforehand.
     */
    void Initialize();

    /**
     * @brief Store a borrowed simulation-data pointer; the object must outlive this analysis.
     */
    void SetProblemData(ProblemData *simData);

    /**
     * @brief Return the borrowed simulation-data pointer without transferring ownership.
     */
    ProblemData *GetProblemData();

    /**
     * @brief Return the Newton loop counter, rather than an independently counted solve total.
     */
    int GetNumberOfIterations();

    /**
     * @brief Run the full Newton loop up to maxIterations, testing residual and
     * scaled-correction norms against configured tolerances after the first step.
     * Loads each update into the mesh; reports nonconvergence without returning status.
     */
    void NewtonIteration();

    /**
     * @brief Choose a correction scale for sol + eta*dsol. ENONE returns 1;
     * bisection and golden ratio evaluate trial solutions; ESECANT triggers DebugStop().
     */
    REAL LineSearchStep(TPZFMatrix<STATE> &sol, TPZFMatrix<STATE> &dsol, LineSearchMethod method = ENONE);

    /**
     * @brief Search for a zero of dsol^T RHS along eta in [0,1].
     * Assumes a multiphysics mesh and changes mesh/RHS state. The early-return
     * branch leaves the mixed-well residual-only flag enabled.
     */
    REAL BissectionMethod(TPZFMatrix<STATE> &sol, TPZFMatrix<STATE> &dsol);

    /**
     * @brief Search eta in [0,1] to reduce the residual norm using golden sections.
     * Changes mesh/RHS state and toggles the mixed-well residual-only flag;
     * the H1 well residual-only flag is not toggled here.
     */
    REAL GoldenRatioMethod(TPZFMatrix<STATE> &sol, TPZFMatrix<STATE> &dsol);

    /**
     * @brief Print a finite-difference Jacobian consistency estimate using step and 2*step.
     * Requires current assembled Jacobian/RHS and no line search. Equation filtering
     * is not supported by this diagnostic; trial mesh solutions are not restored.
     */
    void TaylorTest(TPZFMatrix<STATE> &sol, REAL step = 1e-5);

    /**
     * @brief Delegate to NeoPZ assembly and measure elapsed time (printing is disabled).
     */
    void Assemble() override;

    /**
     * @brief Delegate to NeoPZ linear solution and print elapsed solve time.
     */
    void Solve() override;

    /**
     * @brief Insert reservoir/well boundary material IDs of the requested type into
     * bcMatids; existing entries are retained.
     */
    void FillBCMatids(std::set<int> &bcMatids, int type);

    /**
     * @brief Load prescribed values into boundary connects and copy independent
     * equations into the analysis solution. H1 initializes corner values; mixed
     * boundary coefficients include the boundary-measure/ncorner factor. Well
     * endpoints convert heel-to-toe flux using the adjacent line side orientation.
     */
    void SetInitialSolution(std::set<int> &bcMatids);

    /**
     * @brief Exclude all connect equations on the specified boundary material IDs
     * from the structural matrix, preserving the initialized values during corrections.
     */
    void ApplyEquationFilter(std::set<int> &bcMatids);

    /**
     * @brief Write Pressure/Flux to solution_itN VTK output (N = it) using configured
     * resolution/threads. Negative dimToPost selects the geometric mesh dimension.
     */
    void PostProcessIteration(int dimToPost = -1, int it = -1);

    /**
     * @brief Select the line-search method for Newton steps after the first iteration.
     */
    void SetLineSearchMethod(LineSearchMethod method) { fLineSearchMethod = method; }
};