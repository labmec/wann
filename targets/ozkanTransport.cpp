#include "TPZWannAnalysis.h"
#include "TPZWannApproxTools.h"
#include "TPZWannGeometryTools.h"
#include "TPZWannTransportAnalysis.h"
#include "TPZWannTransportDataTransfer.h"
#include "TPZWannTransportMesh.h"
#include "pzlog.h"

#include <cmath>
#include <iostream>
#include <memory>
#include <stdexcept>
#include <vector>

int main() {
    try {
        TPZLogger::InitializePZLOG();
        ProblemData problem;
        problem.ReadJson("ozkanTransport.json");
        const auto& well = problem.m_Wellbore[0];
        const REAL injectionRate = well.BCs.at("point_heel").value;

        std::unique_ptr<TPZGeoMesh> geometry(TPZWannGeometryTools::CreateGeoMesh(&problem));
        std::vector<std::unique_ptr<TPZCompMesh>> atomicMeshes;
        std::unique_ptr<TPZMultiphysicsCompMesh> flow(
            TPZWannApproxTools::CreateMixedCompMesh(geometry.get(), &problem));
        for (auto* mesh : flow->MeshVector()) atomicMeshes.emplace_back(mesh);
        TPZWannAnalysis flowAnalysis(flow.get(), RenumType::EMetis);
        flowAnalysis.SetProblemData(&problem);
        flowAnalysis.Initialize();
        flowAnalysis.NewtonIteration();
        if (flowAnalysis.GetNumberOfIterations() >= problem.m_Numerics.NewtonFlow.maxIterations)
            throw std::runtime_error("Ozkan coupled flow did not converge");

        std::unique_ptr<TPZCompMesh> transport(TPZWannTransportMesh::Build(geometry.get(), problem));
        TPZWannTransportAnalysis analysis(transport.get());
        analysis.Initialize(flow.get(), problem);
        const auto& algebraic = analysis.fAlgebraicTransport;
        REAL wellInflow = 0.;
        for (REAL flux : algebraic.fInterfaceData.at(well.matidSurf).fIntegralFlux)
            wellInflow -= flux;
        const REAL imbalance = TPZWannTransportDataTransfer::MaxCellFluxImbalance(algebraic);
        std::cout << "Well surface net inflow = " << wellInflow
                  << ", max cell flow imbalance = " << imbalance << '\n';
        if (!(std::abs(wellInflow - injectionRate) < 1.e-6 * injectionRate &&
              imbalance < 1.e-6 * injectionRate))
            throw std::runtime_error("Transferred Ozkan flow failed conservation checks");

        const bool safeVtk = problem.m_Mesh.ToCylindrical;
        analysis.PostProcessTimeStep("ozkan-transport", 0, safeVtk);
        const REAL initialWater = algebraic.WaterVolume();
        REAL boundaryWater = 0.;
        for (int step = 1; step <= problem.m_Transport->numSteps; ++step) {
            std::cout << "Ozkan transport step " << step << '\n';
            analysis.RunTimeStep();
            boundaryWater += problem.m_Transport->timeStep * algebraic.BoundaryWaterRate(problem);
            const REAL water = algebraic.WaterVolume();
            const REAL balance = water - initialWater + boundaryWater;
            if (!(std::abs(balance) < 1.e-4))
                throw std::runtime_error("Ozkan transport failed water conservation");
            std::cout << "Step " << step << ", water volume = " << water
                      << ", water balance error = " << balance << '\n';
            analysis.PostProcessTimeStep("ozkan-transport", step, safeVtk);
        }
    } catch (const std::exception& error) {
        std::cerr << "Ozkan transport failed: " << error.what() << '\n';
        return 1;
    }
    return 0;
}
