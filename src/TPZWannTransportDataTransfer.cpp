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

namespace {
REAL ReservoirFaceRate(TPZInterfaceElement& face, const TPZCompElSide& left,
                       TPZGeoEl* rightGel, TPZInterpolationSpace* flux) {
    auto* gel = face.Reference();
    auto* leftGel = left.Element()->Reference();
    const int dim = leftGel->Dimension();
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
        faceSide.Normal(qsi, leftGel, rightGel, normal);
        transform.Apply(qsi, volumePoint);
        TPZMaterialDataT<STATE> point;
        flux->InitMaterialData(point);
        point.fNeedsSol = true;
        flux->ComputeRequiredData(point, volumePoint);
        REAL normalFlux = 0.;
        for (int d = 0; d < 3; d++) normalFlux += point.sol[0][d] * normal[d];
        integrated += weight * std::abs(determinant) * normalFlux;
    }
    return integrated;
}

REAL WellEndpointRate(const TPZCompElSide& endpoint, TPZInterpolationSpace& flux) {
    auto* line = endpoint.Element()->Reference();
    TPZManVector<REAL, 1> qsi(1);
    line->CenterPoint(endpoint.Side(), qsi);
    TPZFNMatrix<1, REAL> jac, inverse;
    TPZFNMatrix<3, REAL> axes;
    REAL determinant;
    line->Jacobian(qsi, jac, axes, determinant, inverse);
    TPZMaterialDataT<STATE> point;
    flux.InitMaterialData(point);
    point.fNeedsSol = true;
    flux.ComputeRequiredData(point, qsi);
    REAL rate = 0.;
    for (int d = 0; d < 3; ++d)
        rate += point.sol[0][d] * axes(0, d);
    return qsi[0] < 0. ? -rate : rate;
}
}

void TPZWannTransportDataTransfer::Initialize(TPZMultiphysicsCompMesh* flow,
    TPZCompMesh* transport, const ProblemData& data, TPZWannAlgebraicTransport& algebraic) {
    data.ValidateTransport();
    if (!flow || !transport || flow->Reference() != transport->Reference())
        throw std::invalid_argument("Flow and transport must share geometry");
    algebraic = TPZWannAlgebraicTransport();
    auto& cells = algebraic.fCellsData;
    const int dim = transport->Dimension();
    std::map<int, std::pair<REAL, REAL>> wellCells; // line ID -> (area, initial Sw)
    for (const auto& well : data.m_Wellbore) {
        const auto saturation = data.m_Transport->initialWellWaterSaturation.find(well.name);
        if (saturation != data.m_Transport->initialWellWaterSaturation.end())
            wellCells.emplace(well.matid, std::make_pair(std::acos(-1.) * well.radius * well.radius, saturation->second));
    }
    auto isCell = [&](TPZCompEl* cel) {
        if (!cel || !cel->Material()) return false;
        const int id = cel->Material()->Id();
        return (cel->Dimension() == dim && id == data.m_Reservoir.matid) ||
               (cel->Dimension() == 1 && wellCells.count(id));
    };
    std::map<int64_t, int64_t> connectToCell;
    std::map<int64_t, TPZInterpolationSpace*> geoToFlux;
    for (auto* cel : flow->ElementVec()) {
        if (!isCell(cel)) continue;
        auto* mp = dynamic_cast<TPZMultiphysicsElement*>(cel);
        auto* flux = mp ? dynamic_cast<TPZInterpolationSpace*>(mp->Element(0)) : nullptr;
        if (!flux) throw std::invalid_argument("Transport currently requires uncondensed mixed flow elements");
        geoToFlux.emplace(cel->Reference()->Index(), flux);
    }
    for (auto* cel : transport->ElementVec()) {
        if (!isCell(cel)) continue;
        if (cel->NConnects() != 1 || cel->Connect(0).NDof() != 1)
            throw std::invalid_argument("Transport requires one saturation DOF per cell");
        auto* gel = cel->Reference();
        const int64_t index = cells.fVolume.size();
        if (!connectToCell.emplace(cel->ConnectIndex(0), index).second)
            throw std::invalid_argument("Physical transport cells must have distinct saturation connects");
        cells.fGeoIndex.push_back(gel->Index());
        cells.fMatId.push_back(gel->MaterialId());
        cells.fEqNumber.push_back(transport->Block().Position(cel->Connect(0).SequenceNumber()));
        const auto well = wellCells.find(gel->MaterialId());
        const bool isWell = well != wellCells.end();
        cells.fVolume.push_back(gel->Volume() * (isWell ? well->second.first : 1.));
        cells.fPorosity.push_back(isWell ? 1. : data.m_Reservoir.porosity);
        cells.fSaturation.push_back(isWell ? well->second.second : data.m_Transport->initialWaterSaturation);
    }
    if (cells.fVolume.empty()) throw std::invalid_argument("No reservoir transport cells");
    cells.fSaturationLastState = cells.fSaturation;
    cells.UpdateFractionalFlows(data);

    const auto boundaries = data.TransportBoundarySaturations();
    for (auto* cel : transport->ElementVec()) {
        auto* face = dynamic_cast<TPZInterfaceElement*>(cel);
        if (!face) continue;
        const auto left = face->LeftElementSide();
        const auto right = face->RightElementSide();
        auto* leftGel = left.Element()->Reference();
        const int64_t leftCell = connectToCell.at(left.Element()->ConnectIndex(0));
        const int matid = cel->Material()->Id();
        const bool external = boundaries.count(matid);
        int64_t rightCell = -1;
        if (!external) {
            rightCell = connectToCell.at(right.Element()->ConnectIndex(0));
            if (leftCell == rightCell)
                throw std::invalid_argument("Transport interface connects a cell to itself");
        }
        auto* flux = geoToFlux.at(leftGel->Index());
        const REAL integrated = leftGel->Dimension() == 1
            ? WellEndpointRate(left, *flux)
            : ReservoirFaceRate(*face, left, right.Element()->Reference(), flux);
        if (leftGel->Dimension() == 1 && rightCell >= 0) {
            const auto rightRate = WellEndpointRate(right, *geoToFlux.at(right.Element()->Reference()->Index()));
            if (std::abs(integrated + rightRate) > 1.e-10 * (1. + std::abs(integrated)))
                throw std::invalid_argument("Inconsistent flow traces at a well transport point");
        }
        auto& faces = algebraic.fInterfaceData[matid];
        faces.fLeftRightVolIndex.emplace_back(leftCell, rightCell);
        faces.fGeoIndex.push_back(face->Reference()->Index());
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
