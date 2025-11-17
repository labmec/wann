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

class TPZWannAnalysis : public TPZLinearAnalysis
{

public:
    /// Data transfer object
    ProblemData *fSimData;

    /// Number of iterations
    int fKiteration = 0;

    bool fIsFirstAssemble = true;

    /// Default constructor
    TPZWannAnalysis();

    /// Constructor based on a cmesh and optimization band directive
    TPZWannAnalysis(TPZMultiphysicsCompMesh *cmesh,
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

    /// override assemble to have timers
    void Assemble() override;

    /// override solve to have timers
    void Solve() override;

    /// Fill the set of material ids related to Neumann BCs
    /// This method is used to apply Neumann BCs without using BigNumbers
    void FillNeumannBCMatids(std::set<int> &neumannMatids);

    /// Set the initial solution based on the boundary conditions
    /// This method is used to apply Neumann BCs without using BigNumbers
    void SetInitialSolution(std::set<int> &neumannMatids);

    /// Remove equations related to Neumann BCs from the system
    /// This method is used to apply Neumann BCs without using BigNumbers
    void ApplyEquationFilter(std::set<int> &neumannMatids);
};