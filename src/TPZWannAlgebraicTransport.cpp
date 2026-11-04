// Adapted from subFlow/src/TSFAlgebraicTransport.cpp (LabMeC, MIT).
// Retains the accumulation and advective upwind contributions, with b = 1.
#include "TPZWannAlgebraicTransport.h"

void TPZWannAlgebraicTransport::TCellData::UpdateSaturations(const TPZFMatrix<STATE>& sol) {
    for (size_t i = 0; i < fVolume.size(); i++) fSaturation[i] = sol(fEqNumber[i], 0);
}

void TPZWannAlgebraicTransport::TCellData::UpdateFractionalFlows(const ProblemData& data) {
    fWaterfractionalflow.resize(fSaturation.size());
    fDerivativeWfractionalflow.resize(fSaturation.size());
    for (size_t i = 0; i < fSaturation.size(); i++) {
        const auto flow = data.FractionalFlow(fSaturation[i]);
        fWaterfractionalflow[i] = flow.first;
        fDerivativeWfractionalflow[i] = flow.second;
    }
}

void TPZWannAlgebraicTransport::Contribute(int64_t cell, TPZFMatrix<STATE>& ek,
                                         TPZFMatrix<STATE>& ef) const {
    const auto& c = fCellsData;
    const REAL capacity = c.fVolume[cell] * c.fPorosity[cell];
    ef(0, 0) = -capacity * (c.fSaturation[cell] - c.fSaturationLastState[cell]);
    ek(0, 0) = capacity;
}

void TPZWannAlgebraicTransport::ContributeInterface(int64_t face, int matid,
    REAL dt, TPZFMatrix<STATE>& ek, TPZFMatrix<STATE>& ef) const {
    const auto& faces = fInterfaceData.at(matid);
    const auto lr = faces.fLeftRightVolIndex[face];
    const REAL flux = faces.fIntegralFlux[face];
    const REAL beta = flux > 0. ? 1. : 0.;
    const auto& c = fCellsData;
    const REAL fwL = c.fWaterfractionalflow[lr.first];
    const REAL fwR = c.fWaterfractionalflow[lr.second];
    const REAL dfL = c.fDerivativeWfractionalflow[lr.first];
    const REAL dfR = c.fDerivativeWfractionalflow[lr.second];
    const REAL factor = flux * dt;
    ef(0, 0) = -(beta * fwL + (1. - beta) * fwR) * factor;
    ef(1, 0) = -ef(0, 0);
    ek(0, 0) = dfL * beta * factor;
    ek(0, 1) = dfR * (1. - beta) * factor;
    ek(1, 0) = -ek(0, 0);
    ek(1, 1) = -ek(0, 1);
}

void TPZWannAlgebraicTransport::ContributeBC(int64_t face, int matid,
    REAL dt, REAL exteriorFractionalFlow, TPZFMatrix<STATE>& ek, TPZFMatrix<STATE>& ef) const {
    const auto& faces = fInterfaceData.at(matid);
    const auto left = faces.fLeftRightVolIndex[face].first;
    const REAL factor = faces.fIntegralFlux[face] * dt;
    if (factor < 0.) {
        ef(0, 0) = -exteriorFractionalFlow * factor;
    } else {
        ef(0, 0) = -fCellsData.fWaterfractionalflow[left] * factor;
        ek(0, 0) = fCellsData.fDerivativeWfractionalflow[left] * factor;
    }
}

REAL TPZWannAlgebraicTransport::WaterVolume() const {
    REAL water = 0.;
    const auto& c = fCellsData;
    for (size_t i = 0; i < c.fVolume.size(); i++)
        water += c.fVolume[i] * c.fPorosity[i] * c.fSaturation[i];
    return water;
}

REAL TPZWannAlgebraicTransport::BoundaryWaterRate(const ProblemData& data) const {
    REAL rate = 0.;
    for (const auto& boundary : data.TransportBoundarySaturations()) {
        const auto it = fInterfaceData.find(boundary.first);
        if (it == fInterfaceData.end()) continue;
        const auto& faces = it->second;
        const REAL exteriorFractionalFlow = data.FractionalFlow(boundary.second).first;
        for (size_t f = 0; f < faces.fIntegralFlux.size(); f++) {
            const REAL flux = faces.fIntegralFlux[f];
            const auto left = faces.fLeftRightVolIndex[f].first;
            const REAL fw = flux < 0. ? exteriorFractionalFlow
                                     : fCellsData.fWaterfractionalflow[left];
            rate += fw * flux;
        }
    }
    return rate;
}
