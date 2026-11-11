// Adapted from subFlow/src/TSFApproxCreator.cpp (LabMeC, MIT).
#include "TPZWannTransportMesh.h"
#include "TPZWannGeometryTools.h"
#include "TPZWannTransportMaterial.h"
#include "TPZInterfaceEl.h"
#include "TPZBndCondT.h"
#include "pzcmesh.h"
#include "pzgeoelbc.h"
#include <algorithm>
#include <limits>
#include <memory>
#include <set>
#include <stdexcept>

namespace {
using SkinLinks = std::vector<std::pair<TPZGeoEl*, TPZGeoEl*>>;

void AddWellSkin(TPZGeoMesh& geometry, int64_t nel, int lineId, int wallId, int skinId,
                 SkinLinks& links) {
    std::vector<TPZGeoEl*> segments, walls;
    for (int64_t i = 0; i < nel; i++) {
        auto* gel = geometry.Element(i);
        if (!gel) continue;
        if (gel->MaterialId() == lineId && gel->Dimension() == 1) segments.push_back(gel);
        if (gel->MaterialId() == wallId && gel->Dimension() == 2) walls.push_back(gel);
    }
    if (segments.empty() || walls.empty())
        throw std::invalid_argument("Well transport requires line and wall elements");

    TPZManVector<REAL, 3> origin(3), end(3), point(3);
    segments.front()->NodePtr(0)->GetCoordinates(origin);
    segments.front()->NodePtr(1)->GetCoordinates(end);
    const TPZManVector<REAL, 3> axis = end - origin;
    auto span = [&](TPZGeoEl* gel) {
        REAL low = std::numeric_limits<REAL>::infinity(), high = -low;
        for (int i = 0; i < gel->NCornerNodes(); ++i) {
            gel->NodePtr(i)->GetCoordinates(point);
            const REAL x = TPZWannGeometryTools::ComputeAxialCoordinate(point, origin, axis);
            low = std::min(low, x);
            high = std::max(high, x);
        }
        return std::make_pair(low, high);
    };
    std::map<REAL, std::pair<REAL, TPZGeoEl*>> axialSegments;
    for (auto* segment : segments) {
        const auto [low, high] = span(segment);
        if (!(low < high) || !axialSegments.emplace(low, std::make_pair(high, segment)).second)
            throw std::invalid_argument("Invalid well segment axial span");
    }
    std::set<TPZGeoEl*> covered;
    for (auto* wall : walls) {
        const auto [low, high] = span(wall);
        auto segment = axialSegments.upper_bound(low + (high - low) / 2.);
        constexpr REAL tolerance = 1.e-6;
        if (segment == axialSegments.begin())
            throw std::invalid_argument("Well wall face has no axial segment");
        --segment;
        if (!(high > low) || low < segment->first - tolerance || high > segment->second.first + tolerance)
            throw std::invalid_argument("Well wall face crosses an axial segment endpoint");
        TPZGeoElBC skin(wall, wall->NSides() - 1, skinId);
        links.emplace_back(skin.CreatedElement(), segment->second.second);
        covered.insert(segment->second.second);
    }
    if (covered.size() != segments.size())
        throw std::invalid_argument("Well segment has no saturation skin");
}
}

TPZCompMesh* TPZWannTransportMesh::Build(TPZGeoMesh* geometry, const ProblemData& data) {
    data.ValidateTransport();
    const auto boundaries = data.TransportBoundarySaturations();
    const int reservoirId = data.m_Reservoir.matid;
    const int dim = geometry->Dimension();
    if (dim != 2 && dim != 3) throw std::invalid_argument("Transport requires a 2D or 3D reservoir");
    const int64_t nel = geometry->NElements();
    for (auto* gel : geometry->ElementVec()) {
        if (!gel) continue;
        const int id = gel->MaterialId();
        if (id == data.ETransportInterface || id == data.ETransportWellSkin ||
            id == data.ETransportWellWall || id == data.ETransportWellPoint)
            throw std::invalid_argument("Transport auxiliary ID already exists in geometry");
        if (gel->HasSubElement()) throw std::invalid_argument("Refined geometry is not supported yet");
    }
    auto mesh = std::make_unique<TPZCompMesh>(geometry);
    mesh->SetDimModel(dim);
    mesh->SetDefaultOrder(0);
    mesh->SetAllCreateFunctionsDiscontinuous();
    auto* material = new TPZWannTransportMaterial(reservoirId, dim);
    mesh->InsertMaterialObject(material);
    mesh->InsertMaterialObject(new TPZWannTransportMaterial(data.ETransportInterface, dim - 1));

    SkinLinks links;
    std::set<int> wellIds;
    std::map<int, TPZWannTransportMaterial*> endpointMaterials;
    for (const auto& well : data.m_Wellbore) {
        if (!data.m_Transport->initialWellWaterSaturation.count(well.name)) continue;
        if (dim != 3) throw std::invalid_argument("Well transport requires a 3D reservoir");
        auto* wellMaterial = new TPZWannTransportMaterial(well.matid, 1);
        mesh->InsertMaterialObject(wellMaterial);
        wellIds.insert(well.matid);
        endpointMaterials.emplace(well.BCs.at("point_heel").matid, wellMaterial);
        endpointMaterials.emplace(well.BCs.at("point_toe").matid, wellMaterial);
        AddWellSkin(*geometry, nel, well.matid, well.matidSurf, data.ETransportWellSkin, links);
    }
    if (!wellIds.empty()) {
        mesh->InsertMaterialObject(new TPZWannTransportMaterial(data.ETransportWellSkin, 2));
        mesh->InsertMaterialObject(new TPZWannTransportMaterial(data.ETransportWellWall, 2));
        mesh->InsertMaterialObject(new TPZWannTransportMaterial(data.ETransportWellPoint, 0));
    }
    for (const auto& [id, saturation] : boundaries) {
        auto endpoint = endpointMaterials.find(id);
        auto* parent = endpoint == endpointMaterials.end() ? material : endpoint->second;
        TPZFMatrix<STATE> val1(1, 1, 0.);
        TPZManVector<STATE, 1> val2(1, saturation);
        mesh->InsertMaterialObject(parent->CreateBC(parent, id, 0, val1, val2));
    }

    geometry->ResetReference();
    mesh->AutoBuild();
    geometry->ResetReference();
    mesh->LoadReferences();
    for (const auto& [skin, segment] : links) {
        auto* skinCell = skin->Reference();
        auto* wellCell = segment->Reference();
        if (!skinCell || !wellCell)
            throw std::invalid_argument("Missing saturation skin or well segment");
        // NeoPZ's discontinuous builder leaves elements below the mesh
        // dimension without connects. Give each 1D segment one P0 connect.
        if (wellCell->NConnects() == 0)
            wellCell->SetConnectIndex(0, mesh->AllocateNewConnect(1, 1, 0));
        skinCell->SetConnectIndex(0, wellCell->ConnectIndex(0));
    }

    for (int64_t el = 0; el < nel; el++) {
        auto* gel = geometry->Element(el);
        if (!gel) continue;
        const int cellDim = gel->Dimension();
        const int id = gel->MaterialId();
        const bool reservoir = cellDim == dim && id == reservoirId;
        if (!reservoir && !(cellDim == 1 && wellIds.count(id))) continue;
        for (int side = 0; side < gel->NSides() - 1; side++) {
            if (gel->SideDimension(side) != cellDim - 1) continue;
            TPZGeoElSide left(gel, side), right;
            int interfaceId = -1;
            int neighbors = 0;
            for (auto neigh = left.Neighbour(); neigh != left; neigh = neigh.Neighbour()) {
                auto* other = neigh.Element();
                const int otherId = other->MaterialId();
                if (other->Dimension() == cellDim && otherId == id) {
                    interfaceId = reservoir ? data.ETransportInterface : data.ETransportWellPoint;
                } else if (reservoir && otherId == data.ETransportWellSkin) {
                    interfaceId = data.ETransportWellWall;
                } else if (other->Dimension() == cellDim - 1 && boundaries.count(otherId)) {
                    interfaceId = otherId;
                } else {
                    continue;
                }
                right = neigh;
                ++neighbors;
            }
            if (neighbors != 1)
                throw std::invalid_argument("Every transport side needs exactly one neighbor or boundary");
            if (right.Element()->Dimension() == cellDim && gel->Index() > right.Element()->Index()) continue;
            TPZGeoElBC gbc(left, interfaceId);
            auto leftComp = left.Reference(), rightComp = right.Reference();
            if (!leftComp.Element() || !rightComp.Element())
                throw std::invalid_argument("Missing computational element at transport interface");
            new TPZInterfaceElement(*mesh, gbc.CreatedElement(), leftComp, rightComp);
        }
    }
    geometry->BuildConnectivity();
    mesh->ComputeNodElCon();
    mesh->CleanUpUnconnectedNodes();
    mesh->ExpandSolution();
    return mesh.release();
}
