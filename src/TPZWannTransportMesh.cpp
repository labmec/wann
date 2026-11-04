// Adapted from subFlow/src/TSFApproxCreator.cpp (LabMeC, MIT).
#include "TPZWannTransportMesh.h"
#include "TPZWannTransportMaterial.h"
#include "TPZInterfaceEl.h"
#include "TPZBndCondT.h"
#include "pzcmesh.h"
#include "pzgeoelbc.h"
#include <memory>
#include <stdexcept>

TPZCompMesh* TPZWannTransportMesh::Build(TPZGeoMesh* geometry, const ProblemData& data) {
    data.ValidateTransport();
    const auto boundaries = data.TransportBoundarySaturations();
    const int reservoirId = data.m_Reservoir.matid;
    const int dim = geometry->Dimension();
    if (dim != 2 && dim != 3) throw std::invalid_argument("Transport requires a 2D or 3D reservoir");
    const int64_t nel = geometry->NElements();
    for (auto* gel : geometry->ElementVec()) {
        if (!gel) continue;
        if (gel->MaterialId() == data.ETransportInterface)
            throw std::invalid_argument("Transport interface ID already exists in geometry");
        if (gel->HasSubElement()) throw std::invalid_argument("Refined geometry is not supported yet");
    }
    auto mesh = std::make_unique<TPZCompMesh>(geometry);
    mesh->SetDimModel(dim);
    mesh->SetDefaultOrder(0);
    mesh->SetAllCreateFunctionsDiscontinuous();
    auto* material = new TPZWannTransportMaterial(reservoirId, dim);
    mesh->InsertMaterialObject(material);
    for (const auto& boundary : boundaries) {
        TPZFMatrix<STATE> val1(1, 1, 0.);
        TPZManVector<STATE, 1> val2(1, boundary.second);
        mesh->InsertMaterialObject(material->CreateBC(material, boundary.first, 0, val1, val2));
    }
    mesh->InsertMaterialObject(new TPZWannTransportMaterial(data.ETransportInterface, dim - 1));
    mesh->AutoBuild();
    geometry->ResetReference();
    mesh->LoadReferences();

    // Same P0/interface topology as subFlow, but explicitly select reservoir IDs.
    // Coincident WANN auxiliary surfaces are not transport cells or neighbors.
    for (int64_t el = 0; el < nel; el++) {
        auto* gel = geometry->Element(el);
        if (!gel || gel->Dimension() != dim || gel->MaterialId() != reservoirId) continue;
        for (int side = 0; side < gel->NSides() - 1; side++) {
            if (gel->SideDimension(side) != dim - 1) continue;
            TPZGeoElSide left(gel, side), right, boundary;
            int nvolume = 0, nboundary = 0;
            for (auto neigh = left.Neighbour(); neigh != left; neigh = neigh.Neighbour()) {
                auto* other = neigh.Element();
                if (other->Dimension() == dim && other->MaterialId() == reservoirId) {
                    right = neigh;
                    nvolume++;
                } else if (other->Dimension() == dim - 1 && boundaries.count(other->MaterialId())) {
                    boundary = neigh;
                    nboundary++;
                }
            }
            if (nvolume > 1 || (nvolume == 0 && nboundary != 1))
                throw std::invalid_argument("Every reservoir face needs one neighbor or one configured boundary");
            if (nvolume == 1) {
                if (gel->Index() > right.Element()->Index()) continue;
                TPZGeoElBC gbc(left, data.ETransportInterface);
                auto leftComp = left.Reference(), rightComp = right.Reference();
                new TPZInterfaceElement(*mesh, gbc.CreatedElement(), leftComp, rightComp);
            } else {
                TPZGeoElBC gbc(left, boundary.Element()->MaterialId());
                auto leftComp = left.Reference(), rightComp = boundary.Reference();
                new TPZInterfaceElement(*mesh, gbc.CreatedElement(), leftComp, rightComp);
            }
        }
    }
    geometry->BuildConnectivity();
    mesh->ComputeNodElCon();
    mesh->CleanUpUnconnectedNodes();
    mesh->ExpandSolution();
    return mesh.release();
}
