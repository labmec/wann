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

class TPZWannAnalysis : public TPZLinearAnalysis
{

public:
    enum LineSearchMethod
    {
        ENONE,
        EBISSECT,
        EGOLDEN,
        ESECANT
    };

    /// Data transfer object
    ProblemData *fSimData;

    /// Number of iterations
    int fKiteration = 0;

    bool fIsFirstAssemble = true;

    bool fIsMultiphysics = true;

    /// Default constructor
    TPZWannAnalysis();

    /// Constructor based on a cmesh and optimization band directive
    TPZWannAnalysis(TPZCompMesh *cmesh,
                    const RenumType &renumtype = RenumType::EDefault);

    /// Default destructor
    ~TPZWannAnalysis();

    /// Configurates iternal members
    void Initialize();

    /// Set data transfer object
    void SetProblemData(ProblemData *simData);

    /// Get data transfer object
    ProblemData *GetProblemData();

    /// Get the number of iterations
    int GetNumberOfIterations();

    /// Perform a Newton iteration
    void NewtonIteration();

    REAL LineSearchStep(TPZFMatrix<STATE> &sol, TPZFMatrix<STATE> &dsol, LineSearchMethod method = ENONE);

    REAL BissectionMethod(TPZFMatrix<STATE> &sol, TPZFMatrix<STATE> &dsol);

    REAL GoldenRatioMethod(TPZFMatrix<STATE> &sol, TPZFMatrix<STATE> &dsol);

    /// Perform a Taylor test in the current Newton iteration
    /// Used only for debugging purposes (check correctness of Jacobian implementation)
    /// To get correct results, no line search must be performed
    void TaylorTest(TPZFMatrix<STATE> &sol, REAL step = 1e-5);

    /// override assemble to have timers
    void Assemble() override;

    /// override solve to have timers
    void Solve() override;

    /// Fill the set of material ids related to BCs of given type
    /// This method is used to apply BCs without using BigNumbers
    void FillBCMatids(std::set<int> &bcMatids, int type);

    /// Set the initial solution based on the boundary conditions
    /// This method is used to apply BCs without using BigNumbers
    void SetInitialSolution(std::set<int> &bcMatids);

    /// Remove equations related to specified BCs from the system
    /// This method is used to apply BCs without using BigNumbers
    void ApplyEquationFilter(std::set<int> &bcMatids);

    /// Render a vtk file with requested variables for a time step
    void PostProcessIteration(int dimToPost = -1, int it = -1);
};