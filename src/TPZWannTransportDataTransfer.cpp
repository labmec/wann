// Cell/face correspondence follows subFlow::TSFDataTransfer (LabMeC, MIT).
#include "TPZWannTransportDataTransfer.h"
#include "TPZMultiphysicsCompMesh.h"
#include "TPZInterfaceEl.h"
#include "TPZMaterial.h"
#include "pzmultiphysicselement.h"
#include "pzinterpolationspace.h"
#include "pzquad.h"
#include <algorithm>
#include <cmath>
#include <memory>
#include <stdexcept>

void TPZWannTransportDataTransfer::Initialize(TPZMultiphysicsCompMesh* flow,
    TPZCompMesh* transport, const ProblemData& data, TPZWannAlgebraicTransport& algebraic) {
    data.ValidateTransport();
    if (!flow || !transport || flow->Reference() != transport->Reference())
        throw std::invalid_argument("Flow and transport must share geometry");
    algebraic = TPZWannAlgebraicTransport();
    auto& cells = algebraic.fCellsData;
    const int dim = transport->Dimension();
    std::map<int64_t, int64_t> geoToCell;
    std::map<int64_t, TPZInterpolationSpace*> geoToFlux;
    for (auto* cel : flow->ElementVec()) {
        if (!cel || cel->Dimension() != dim || cel->Material()->Id() != data.m_Reservoir.matid) continue;
        auto* mp = dynamic_cast<TPZMultiphysicsElement*>(cel);
        auto* flux = mp ? dynamic_cast<TPZInterpolationSpace*>(mp->Element(0)) : nullptr;
        if (!flux) throw std::invalid_argument("Transport currently requires uncondensed mixed reservoir elements");
        geoToFlux.emplace(cel->Reference()->Index(), flux);
    }
    for (auto* cel : transport->ElementVec()) {
        if (!cel || cel->Dimension() != dim || cel->Material()->Id() != data.m_Reservoir.matid) continue;
        if (cel->NConnects() != 1 || cel->Connect(0).NDof() != 1)
            throw std::invalid_argument("Transport requires one saturation DOF per cell");
        auto* gel = cel->Reference();
        const int64_t index = cells.fVolume.size();
        geoToCell.emplace(gel->Index(), index);
        cells.fGeoIndex.push_back(gel->Index());
        cells.fMatId.push_back(gel->MaterialId());
        cells.fEqNumber.push_back(transport->Block().Position(cel->Connect(0).SequenceNumber()));
        cells.fVolume.push_back(gel->Volume());
        cells.fPorosity.push_back(data.m_Reservoir.porosity);
        cells.fSaturation.push_back(data.m_Transport->initialWaterSaturation);
    }
    if (cells.fVolume.empty()) throw std::invalid_argument("No reservoir transport cells");
    cells.fSaturationLastState = cells.fSaturation;
    cells.UpdateFractionalFlows(data);

    for (auto* cel : transport->ElementVec()) {
        auto* face = dynamic_cast<TPZInterfaceElement*>(cel);
        if (!face) continue;
        const auto left = face->LeftElementSide();
        const auto right = face->RightElementSide();
        auto* leftGel = left.Element()->Reference();
        const int64_t leftCell = geoToCell.at(leftGel->Index());
        const auto rightIt = geoToCell.find(right.Element()->Reference()->Index());
        const int64_t rightCell = rightIt == geoToCell.end() ? -1 : rightIt->second;
        const int matid = cel->Material()->Id();
        if ((rightCell < 0) != (matid != data.ETransportInterface))
            throw std::invalid_argument("Inconsistent transport interface classification");
        auto* flux = geoToFlux.at(leftGel->Index());
        auto* gel = face->Reference();
        TPZGeoElSide faceSide(gel);
        TPZTransform<> faceToSide(dim - 1);
        faceSide.SideTransform3(left.Reference(), faceToSide);
        TPZTransform<> sideToVolume = leftGel->SideToSideTransform(left.Side(), leftGel->NSides() - 1);
        TPZTransform<> transform = sideToVolume.Multiply(faceToSide);
        std::unique_ptr<TPZIntPoints> rule(gel->CreateSideIntegrationRule(gel->NSides() - 1, 2 * flux->MaxOrder() + 4));
        REAL integrated = 0.;
        for (int ip = 0; ip < rule->NPoints(); ip++) {
            TPZManVector<REAL, 3> qsi(dim - 1, 0.), volumePoint(dim, 0.), normal(3, 0.);
            REAL weight, determinant;
            rule->Point(ip, qsi, weight);
            TPZFNMatrix<9, REAL> jac, axes, inverse;
            gel->Jacobian(qsi, jac, axes, determinant, inverse);
            faceSide.Normal(qsi, leftGel, right.Element()->Reference(), normal);
            transform.Apply(qsi, volumePoint);
            TPZMaterialDataT<STATE> point;
            flux->InitMaterialData(point);
            point.fNeedsSol = true;
            flux->ComputeRequiredData(point, volumePoint);
            REAL normalFlux = 0.;
            for (int d = 0; d < 3; d++) normalFlux += point.sol[0][d] * normal[d];
            integrated += weight * std::abs(determinant) * normalFlux;
        }
        auto& faces = algebraic.fInterfaceData[matid];
        faces.fLeftRightVolIndex.emplace_back(leftCell, rightCell);
        faces.fGeoIndex.push_back(gel->Index());
        faces.fIntegralFlux.push_back(integrated);
    }
}

REAL TPZWannTransportDataTransfer::MaxCellFluxImbalance(const TPZWannAlgebraicTransport& algebraic) {
    std::vector<REAL> balance(algebraic.fCellsData.fVolume.size(), 0.);
    for (const auto& group : algebraic.fInterfaceData) {
        const auto& faces = group.second;
        for (size_t i = 0; i < faces.fIntegralFlux.size(); i++) {
            const auto lr = faces.fLeftRightVolIndex[i];
            balance[lr.first] += faces.fIntegralFlux[i];
            if (lr.second >= 0) balance[lr.second] -= faces.fIntegralFlux[i];
        }
    }
    REAL maximum = 0.;
    for (REAL value : balance) maximum = std::max(maximum, std::abs(value));
    return maximum;
}
