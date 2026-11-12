#include <cmath>
#include <iostream>
#include <memory>
#include <vector>

#include "TPZHDivApproxCreator.h"
#include "TPZLinearAnalysis.h"
#include "TPZWannStructMatrix.h"
#include "pzstepsolver.h"
#include "DarcyFlow/TPZMixedDarcyFlow.h"
#include "TPZWannGeometryTools.h"
#include "TPZWannPostProcTools.h"

int main(int argc, char* argv[]) {
    if (argc > 2) {
        std::cerr << "Usage: " << argv[0] << " [input.json]\n";
        return 1;
    }
#ifdef PZ_LOG
    TPZLogger::InitializePZLOG();
#endif
    std::string input = argc == 2 ? argv[1] : "test-sideOrient-1d.json";
    if (input.find(".json") == std::string::npos) input += ".json";
    ProblemData data;
    data.ReadJson(input);
    if (data.m_Wellbore.size() != 1) {
        std::cerr << "The orientation test requires one well.\n";
        return 1;
    }
    const auto& well = data.m_Wellbore[0];
    const auto& heel = well.BCs.at("point_heel");
    const auto& toe = well.BCs.at("point_toe");
    if (heel.type != 0 || toe.type != 0) {
        std::cerr << "The orientation test requires pressures at both endpoints.\n";
        return 1;
    }

    std::unique_ptr<TPZGeoMesh> geometry(TPZWannGeometryTools::ReadMeshFromGmsh(&data));
    geometry->SetDimension(1);
    TPZManVector<REAL, 3> heelPoint(3, 0.), toePoint(3, 0.);
    bool foundHeel = false, foundToe = false;
    for (auto* element : geometry->ElementVec()) {
        if (!element || element->Dimension() != 0) continue;
        if (element->MaterialId() == heel.matid) {
            element->NodePtr(0)->GetCoordinates(heelPoint);
            foundHeel = true;
        }
        if (element->MaterialId() == toe.matid) {
            element->NodePtr(0)->GetCoordinates(toePoint);
            foundToe = true;
        }
    }
    const TPZManVector<REAL, 3> axis = toePoint - heelPoint;
    const REAL length = Norm(axis);
    if (!foundHeel || !foundToe || !(length > 0.)) {
        std::cerr << "Could not locate a valid heel and toe in the mesh.\n";
        return 1;
    }

    TPZHDivApproxCreator creator(geometry.get());
    creator.SetProbType(ProblemType::EDarcy);
    creator.SetDefaultOrder(data.m_Numerics.wellPorder);
    creator.SetShouldCondense(false);
    const REAL mobility = 1. / data.m_Numerics.pressureScale;
    auto* material = new TPZMixedDarcyFlow(well.matid, 1);
    material->SetConstantPermeability(mobility);
    creator.InsertMaterialObject(material);
    for (const auto& [name, bc] : well.BCs) {
        TPZFMatrix<STATE> val1(1, 1, 0.);
        TPZManVector<STATE, 1> val2(1, bc.value);
        creator.InsertMaterialObject(material->CreateBC(material, bc.matid, bc.type, val1, val2));
    }

    std::vector<std::unique_ptr<TPZCompMesh>> atomicMeshes;
    std::unique_ptr<TPZMultiphysicsCompMesh> mesh(creator.CreateApproximationSpace());
    for (auto* atomic : mesh->MeshVector()) atomicMeshes.emplace_back(atomic);

    TPZLinearAnalysis analysis(mesh.get());
    TPZWannStructMatrix matrix(mesh.get());
    matrix.SetNumThreads(data.m_Numerics.nthreads);
    analysis.SetStructuralMatrix(matrix);
    TPZStepSolver<STATE> solver;
    solver.SetDirect(ELDLt);
    analysis.SetSolver(solver);
    analysis.Run();

    const REAL expectedRate = mobility * (heel.value - toe.value) / length;
    int segments = 0;
    for (auto* element : mesh->ElementVec()) {
        if (!element || !element->Reference() || element->Reference()->MaterialId() != well.matid) continue;
        ++segments;
        for (REAL endpoint : {-1., 1.}) {
            TPZManVector<REAL, 1> point(1, endpoint);
            TPZManVector<STATE, 3> flux(3, 0.);
            element->Solution(point, material->VariableIndex("Flux"), flux);
            const REAL rate = Dot(flux, axis) / length;
            if (!(std::abs(rate - expectedRate) < 1.e-8 * (1. + std::abs(expectedRate)))) {
                std::cerr << "Well flux " << rate << " differs from expected " << expectedRate << '\n';
                return 1;
            }
        }
    }
    if (segments == 0) return 1;
    std::cout << "Verified heel-to-toe flux " << expectedRate << " on " << segments << " segments.\n";
    TPZWannPostProcTools::WriteWellboreVTK(mesh.get(), &data);
}
