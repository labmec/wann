#include <memory>
#include <vector>

#include "TPZHDivApproxCreator.h"
#include "TPZVTKGenerator.h"
#include "pzgeoelbc.h"
#include "pzlog.h"
#include "TPZWannAnalysis.h"
#include "TPZWannH1Well.h"
#include "TPZWannMixedWell.h"

namespace {
constexpr int wellId = 1;
constexpr int heelId = 3;
constexpr int toeId = 2;
constexpr int segmentCount = 20;
constexpr REAL reservoirPressure = 2.2e7;
constexpr REAL exchangeCoefficient = 1.e-11;

void BuildWellLine(TPZGeoMesh& geometry, REAL length) {
    geometry.SetDimension(1);
    geometry.NodeVec().Resize(segmentCount + 1);
    for (int i = 0; i <= segmentCount; ++i) {
        TPZManVector<REAL, 3> point = {length * i / segmentCount, 0., 0.};
        geometry.NodeVec()[i].Initialize(i, point, geometry);
    }

    TPZGeoEl *first = nullptr, *last = nullptr;
    for (int i = 0; i < segmentCount; ++i) {
        TPZManVector<int64_t, 2> nodes = {i, i + 1};
        int64_t index;
        last = geometry.CreateGeoElement(EOned, nodes, wellId, index);
        if (i == 0) first = last;
    }
    geometry.BuildConnectivity();
    TPZGeoElBC heel(TPZGeoElSide(first, 0), heelId);
    TPZGeoElBC toe(TPZGeoElSide(last, 1), toeId);
}

template<class MeshBuilder, class Material>
void InsertWellMaterial(MeshBuilder& builder, Material* material, const ProblemData& data,
                        bool mixed) {
    builder.InsertMaterialObject(material);
    for (const auto& [name, bc] : data.m_Wellbore[0].BCs) {
        TPZFMatrix<STATE> val1(1, 1, 0.);
        // Input rates are positive heel to toe. The H1 load uses inward flux.
        TPZManVector<STATE, 1> val2(1, !mixed && name == "point_toe" ? -bc.value : bc.value);
        builder.InsertMaterialObject(material->CreateBC(material, bc.matid, bc.type, val1, val2));
    }
}
}

int main() {
#ifdef PZ_LOG
    TPZLogger::InitializePZLOG();
#endif
    ProblemData data;
    data.m_Numerics.maxIterations = 10;
    data.m_Numerics.res_tol = 1.e-6;
    data.m_Numerics.corr_tol = 1.e-6;
    data.m_Numerics.nthreads = 0;
    data.m_Numerics.wellPorder = 2;
    data.m_Wellbore.Resize(1);
    auto& well = data.m_Wellbore[0];
    well.matid = wellId;
    well.radius = 0.05;
    well.length = 400.;
    well.BCs["point_heel"] = {heelId, 1, 0.05};
    well.BCs["point_toe"] = {toeId, 1, 0.};
    data.m_Fluid.Resize(1);
    data.m_Fluid[0] = {"oil", 0.001, 800.};
    const auto& fluid = data.m_Fluid[0];
    const REAL diameter = 2. * well.radius;

    TPZGeoMesh geometry;
    BuildWellLine(geometry, well.length);
    TPZHDivApproxCreator creator(&geometry);
    creator.SetProbType(ProblemType::EDarcy);
    creator.SetDefaultOrder(data.m_Numerics.wellPorder);
    creator.SetShouldCondense(false);
    InsertWellMaterial(creator, new TPZWannMixedWell(wellId, diameter, fluid.viscosity,
        fluid.density, reservoirPressure, exchangeCoefficient), data, true);
    std::vector<std::unique_ptr<TPZCompMesh>> atomicMeshes;
    std::unique_ptr<TPZMultiphysicsCompMesh> mixed(creator.CreateApproximationSpace());
    for (auto* mesh : mixed->MeshVector()) atomicMeshes.emplace_back(mesh);

    geometry.ResetReference();
    TPZCompMesh h1(&geometry);
    h1.SetDimModel(1);
    h1.SetDefaultOrder(data.m_Numerics.wellPorder);
    h1.SetAllCreateFunctionsContinuous();
    InsertWellMaterial(h1, new TPZWannH1Well(wellId, diameter, fluid.viscosity,
        fluid.density, reservoirPressure, exchangeCoefficient), data, false);
    h1.AutoBuild();

    for (auto* mesh : {static_cast<TPZCompMesh*>(mixed.get()), &h1}) {
        const bool isMixed = mesh == mixed.get();
        geometry.ResetReference();
        TPZWannAnalysis analysis(mesh, RenumType::EMetis);
        analysis.SetProblemData(&data);
        analysis.Initialize();
        analysis.NewtonIteration();
        if (analysis.GetNumberOfIterations() >= data.m_Numerics.maxIterations) return 1;

        TPZStack<std::string> fields;
        fields.Push("Flux");
        fields.Push("Pressure");
        if (isMixed) fields.Push("Divergence");
        const char* filename = isMixed ? "testWell_mixed" : "testWell_H1";
        TPZVTKGenerator vtk(mesh, fields, filename, 0, 1);
        vtk.Do();
    }
}
