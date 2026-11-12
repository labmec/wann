#include <memory>
#include <vector>

#include "TPZVTKGenerator.h"
#include "TPZWannAnalysis.h"
#include "TPZWannApproxTools.h"
#include "TPZWannGeometryTools.h"

int main() {
#ifdef PZ_LOG
    TPZLogger::InitializePZLOG();
#endif
    ProblemData data;
    data.m_Mesh.file = "../geo/old/mesh3D_rev08.msh";
    data.m_PostProc.verbosityLevel = 0;
    data.m_Numerics.nthreads = 0;
    data.m_Numerics.reservoirPorder = 1;
    data.m_Numerics.wellPorder = 2;
    data.m_Numerics.maxIterations = 10;
    data.m_Numerics.res_tol = 1.e-6;
    data.m_Numerics.corr_tol = 1.e-6;
    data.m_Fluid.Resize(1);
    data.m_Fluid[0] = {"fluid", 1., 1.};
    data.m_Reservoir.matid = 110;
    data.m_Reservoir.perm.Fill(1.);
    // Material IDs from the historical reservoir mesh. The wall is a
    // prescribed-pressure reservoir boundary in this standalone Darcy test.
    data.m_Reservoir.BCs["surface_wellbore_cylinder"] = {103, 0, 1.};
    data.m_Reservoir.BCs["surface_wellbore_toe"] = {104, 1, 0.};
    data.m_Reservoir.BCs["surface_wellbore_heel"] = {105, 1, 0.};
    data.m_Reservoir.BCs["surface_farfield"] = {106, 0, 3.};
    data.m_Reservoir.BCs["surface_cap_rock"] = {107, 1, 0.};
    data.UpdateAuxiliaryMaterialIds();

    std::unique_ptr<TPZGeoMesh> geometry(TPZWannGeometryTools::ReadMeshFromGmsh(&data));
    std::vector<std::unique_ptr<TPZCompMesh>> atomicMeshes;
    std::unique_ptr<TPZMultiphysicsCompMesh> mixed(
        TPZWannApproxTools::CreateMixedCompMesh(geometry.get(), &data));
    for (auto* mesh : mixed->MeshVector()) atomicMeshes.emplace_back(mesh);
    std::unique_ptr<TPZCompMesh> h1(TPZWannApproxTools::CreateH1CompMesh(geometry.get(), &data));

    TPZManVector<std::string, 2> fields = {"Flux", "Pressure"};
    for (auto* mesh : {static_cast<TPZCompMesh*>(mixed.get()), h1.get()}) {
        geometry->ResetReference();
        TPZWannAnalysis analysis(mesh, RenumType::EMetis);
        analysis.SetProblemData(&data);
        analysis.Initialize();
        analysis.NewtonIteration();
        if (analysis.GetNumberOfIterations() >= data.m_Numerics.maxIterations) return 1;

        const char* filename = mesh == mixed.get() ? "DarcyNL_Mixed" : "DarcyNL_H1";
        TPZVTKGenerator vtk(mesh, fields, filename, 0, geometry->Dimension());
        vtk.Do();
    }
}
