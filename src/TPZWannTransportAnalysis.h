#pragma once

#include "TPZLinearAnalysis.h"
#include "TPZWannAlgebraicTransport.h"

class TPZMultiphysicsCompMesh;

/** Transport Newton/assembly driver adapted from subFlow::TSFTransportAnalysis.
 * Borrows transport mesh and configuration; they must outlive the analysis.
 * Initialize extracts a frozen WANN mixed flow. Each successful RunTimeStep
 * commits one backward-Euler step. Failure restores the last accepted state.
 */
class TPZWannTransportAnalysis : public TPZLinearAnalysis {
    const ProblemData* fData = nullptr;
    TPZFMatrix<STATE> fMassMatrix;
public:
    TPZWannAlgebraicTransport fAlgebraicTransport;
    explicit TPZWannTransportAnalysis(TPZCompMesh* mesh);
    void Initialize(TPZMultiphysicsCompMesh* flow, const ProblemData& data);
    void Assemble() override;
    void RunTimeStep();
    void PostProcessTimeStep(const std::string& filename, int step, bool safeCheck = false);
};
