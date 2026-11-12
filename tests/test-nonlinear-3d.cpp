#include <cmath>
#include <iostream>
#include <memory>
#include <vector>

#include "TPZGenGrid3D.h"
#include "TPZMaterial.h"
#include "TPZVTKGenerator.h"
#include "TPZWannAnalysis.h"
#include "TPZWannApproxTools.h"

int main() {
#ifdef PZ_LOG
    TPZLogger::InitializePZLOG();
#endif
    enum MaterialId { Domain = 1, Bottom, Right, Top, Left, Front, Back };
    TPZGenGrid3D generator(TPZManVector<REAL, 3>{0., 0., 0.},
                           TPZManVector<REAL, 3>{1., 1., 1.},
                           TPZManVector<int, 3>{1, 1, 1}, MMeshType::EHexahedral);
    generator.BuildVolumetricElements(Domain);
    std::unique_ptr<TPZGeoMesh> geometry(generator.BuildBoundaryElements(Bottom, Left, Front, Right, Back, Top));

    ProblemData data;
    data.m_PostProc.verbosityLevel = 0;
    data.m_Numerics.nthreads = 0;
    data.m_Numerics.reservoirPorder = 1;
    data.m_Numerics.wellPorder = 2;
    data.m_Numerics.maxIterations = 10;
    data.m_Numerics.res_tol = 1.e-6;
    data.m_Numerics.corr_tol = 1.e-6;
    data.m_Fluid.Resize(1);
    data.m_Fluid[0] = {"fluid", 1., 1.};
    data.m_Reservoir.matid = Domain;
    data.m_Reservoir.perm.Fill(1.);
    data.m_Reservoir.BCs["bottom"] = {Bottom, 1, 0.};
    data.m_Reservoir.BCs["right"] = {Right, 1, 0.};
    data.m_Reservoir.BCs["top"] = {Top, 1, 0.};
    data.m_Reservoir.BCs["left"] = {Left, 1, 0.};
    data.m_Reservoir.BCs["front"] = {Front, 0, 3.};
    data.m_Reservoir.BCs["back"] = {Back, 1, 1.};
    data.UpdateAuxiliaryMaterialIds();

    TPZManVector<std::string, 2> fields = {"Flux", "Pressure"};
    for (bool isMixed : {true, false}) {
        geometry->ResetReference();
        std::vector<std::unique_ptr<TPZCompMesh>> atomicMeshes;
        std::unique_ptr<TPZCompMesh> mesh;
        if (isMixed) {
            auto* mixed = TPZWannApproxTools::CreateMixedCompMesh(geometry.get(), &data);
            mesh.reset(mixed);
            for (auto* atomic : mixed->MeshVector()) atomicMeshes.emplace_back(atomic);
        } else {
            mesh.reset(TPZWannApproxTools::CreateH1CompMesh(geometry.get(), &data));
        }
        TPZWannAnalysis analysis(mesh.get(), RenumType::EMetis);
        analysis.SetProblemData(&data);
        analysis.Initialize();
        analysis.NewtonIteration();
        if (analysis.GetNumberOfIterations() >= data.m_Numerics.maxIterations) return 1;

        // With unit mobility and the prescribed boundary data, p = 3 - y and q = (0, 1, 0).
        for (auto* element : mesh->ElementVec()) {
            if (!element || !element->Material() || element->Material()->Id() != Domain) continue;
            TPZManVector<REAL, 3> center(3, 0.);
            TPZManVector<STATE, 3> flux(3, 0.);
            TPZManVector<STATE, 1> pressure(1, 0.);
            element->Solution(center, element->Material()->VariableIndex("Flux"), flux);
            element->Solution(center, element->Material()->VariableIndex("Pressure"), pressure);
            if (!(std::abs(pressure[0] - 2.5) < 1.e-8 && std::abs(flux[0]) < 1.e-8 &&
                  std::abs(flux[1] - 1.) < 1.e-8 && std::abs(flux[2]) < 1.e-8)) {
                std::cerr << "Unexpected pressure or flux at the cube center.\n";
                return 1;
            }
        }

        const char* filename = isMixed ? "DarcyNL_Mixed" : "DarcyNL_H1";
        TPZVTKGenerator vtk(mesh.get(), fields, filename, 0, geometry->Dimension());
        vtk.Do();
    }
}
