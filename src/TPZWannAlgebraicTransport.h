#pragma once

#include "ProblemData.h"
#include "pzfmatrix.h"
#include <cstdint>
#include <vector>

/** Incompressible, gravity-free subset of subFlow::TSFAlgebraicTransport.
 * Cell/face arrays and local residual/Jacobian assembly retain subFlow's design.
 * ef is -R and ek is dR/dS. See subFlow-LICENSE.txt for the original license.
 */
class TPZWannAlgebraicTransport {
public:
    struct TCellData {
        std::vector<int64_t> fEqNumber, fGeoIndex;
        std::vector<int> fMatId;
        std::vector<REAL> fVolume, fPorosity;
        std::vector<REAL> fSaturation, fSaturationLastState;
        std::vector<REAL> fWaterfractionalflow, fDerivativeWfractionalflow;
        void UpdateSaturations(const TPZFMatrix<STATE>& solution);
        void UpdateFractionalFlows(const ProblemData& data);
    };
    struct TInterfaceData {
        /// Left cell; right cell is -1 for external boundaries.
        std::vector<std::pair<int64_t, int64_t>> fLeftRightVolIndex;
        std::vector<int64_t> fGeoIndex;
        /// Positive left -> right, or outward from left at a boundary [m^3/s].
        std::vector<REAL> fIntegralFlux;
    };
    TCellData fCellsData;
    std::map<int, TInterfaceData> fInterfaceData;

    void Contribute(int64_t cell, TPZFMatrix<STATE>& ek, TPZFMatrix<STATE>& ef) const;
    void ContributeInterface(int64_t face, int matid, REAL dt,
                             TPZFMatrix<STATE>& ek, TPZFMatrix<STATE>& ef) const;
    void ContributeBC(int64_t face, int matid, REAL dt, REAL exteriorFractionalFlow,
                      TPZFMatrix<STATE>& ek, TPZFMatrix<STATE>& ef) const;
    REAL WaterVolume() const;
    /// Current outward water flow summed over all external boundaries.
    REAL BoundaryWaterRate(const ProblemData& data) const;
};
