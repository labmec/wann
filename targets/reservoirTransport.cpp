/** Initial target to test transport part.
 * Solves WANN mixed Darcy on a unit 3D channel, freezes its flow, and advances
 * subFlow-derived P0 transport. No JSON/Gmsh input or well elements are needed.
 * Writes saturation fields at initialization and after each time step.
 */
#include "TPZWannMixedDarcyNL.h"
#include "TPZWannAnalysis.h"
#include "TPZWannTransportMesh.h"
#include "TPZWannTransportAnalysis.h"
#include "TPZWannTransportDataTransfer.h"
#include "TPZHDivApproxCreator.h"
#include "TPZGenGrid3D.h"
#include <algorithm>
#include <cmath>
#include <iostream>
#include <memory>
#include <stdexcept>
#include <string>

namespace {
enum { Reservoir = 1, Bottom = 2, Left = 3, Front = 4, Right = 5, Back = 6, Top = 7 };
constexpr int nx = 8;

void Require(bool condition, const std::string& message) {
    if (!condition) throw std::runtime_error(message);
}

void CheckDirection(bool reverse) {
    TPZGenGrid3D generator(TPZManVector<REAL, 3>{0., 0., 0.},
                          TPZManVector<REAL, 3>{1., 1., 1.},
                          TPZManVector<int, 3>{nx, 2, 2}, MMeshType::EHexahedral);
    generator.BuildVolumetricElements(Reservoir);
    std::unique_ptr<TPZGeoMesh> geometry(generator.BuildBoundaryElements(Bottom, Left, Front, Right, Back, Top));

    ProblemData problem;
    problem.m_Reservoir.matid = Reservoir;
    problem.m_Reservoir.porosity = 0.25;
    problem.m_Reservoir.BCs["left"] = {Left, 0, reverse ? 0. : 1.};
    problem.m_Reservoir.BCs["right"] = {Right, 0, reverse ? 1. : 0.};
    problem.m_Reservoir.BCs["bottom"] = {Bottom, 1, 0.};
    problem.m_Reservoir.BCs["front"] = {Front, 1, 0.};
    problem.m_Reservoir.BCs["back"] = {Back, 1, 0.};
    problem.m_Reservoir.BCs["top"] = {Top, 1, 0.};
    problem.m_Numerics.nthreads = 0;
    problem.m_Numerics.NewtonFlow.maxIterations = 5;
    problem.m_Numerics.NewtonFlow.res_tol = 1.e-10;
    problem.m_Numerics.NewtonFlow.corr_tol = 1.e-10;

    TPZHDivApproxCreator creator(geometry.get());
    creator.ProbType() = ProblemType::EDarcy;
    creator.SetDefaultOrder(1);
    creator.SetShouldCondense(false);
    auto* material = new TPZWannMixedDarcyNL(Reservoir, 3);
    material->SetConstantPermeability(1.);
    creator.InsertMaterialObject(material);
    for (const auto& entry : problem.m_Reservoir.BCs) {
        const auto& bc = entry.second;
        TPZFMatrix<STATE> val1(1, 1, 0.);
        TPZManVector<STATE, 1> val2(1, bc.value);
        creator.InsertMaterialObject(material->CreateBC(material, bc.matid, bc.type, val1, val2));
    }
    TPZManVector<TPZCompMesh*, 7> atomic(creator.NumMeshes());
    int level = 1;
    creator.CreateAtomicMeshes(atomic, level);
    // Atomic meshes must outlive the multiphysics mesh, which borrows them.
    std::vector<std::unique_ptr<TPZCompMesh>> atomicOwners;
    for (auto* mesh : atomic) atomicOwners.emplace_back(mesh);
    TPZMultiphysicsCompMesh* mixed = nullptr;
    creator.CreateMultiPhysicsMesh(atomic, level, mixed);
    std::unique_ptr<TPZMultiphysicsCompMesh> flow(mixed);
    {
        TPZWannAnalysis analysis(flow.get(), RenumType::ENone);
        analysis.SetProblemData(&problem);
        analysis.Initialize();
        analysis.NewtonIteration();
    }
    const TPZFMatrix<STATE> frozenSolution = flow->Solution();

    problem.UpdateAuxiliaryMaterialIds();
    problem.m_Fluid.push_back({"water", 1., 1.});
    problem.m_Fluid.push_back({"oil", 1., 1.});
    problem.m_Numerics.NewtonTransport.printNewtonProgress = true;
    auto& settings = problem.m_Transport.emplace();
    settings.timeStep = 0.01;
    settings.numSteps = 20;
    for (const auto& [name, bc] : problem.m_Reservoir.BCs)
        settings.boundaryWaterSaturation[name] = 1.;
    std::unique_ptr<TPZCompMesh> transport(TPZWannTransportMesh::Build(geometry.get(), problem));

    for (int scenario = 0; scenario < 3; scenario++) {
        const bool constant = scenario == 0;
        const bool quadratic = scenario == 2;
        const REAL initial = constant ? 0.37 : 0.;
        settings.initialWaterSaturation = initial;
        settings.relativePermeability = quadratic ? ProblemData::RelativePermeability::Quadratic
                                                  : ProblemData::RelativePermeability::Linear;
        problem.m_Fluid[1] = {"oil", quadratic ? 5. : 1., 1.};
        for (auto& boundary : settings.boundaryWaterSaturation)
            boundary.second = constant ? initial : 1.;
        TPZWannTransportAnalysis analysis(transport.get());
        analysis.Initialize(flow.get(), problem);
        auto& algebraic = analysis.fAlgebraicTransport;
        const auto& cells = algebraic.fCellsData;
        Require(TPZWannTransportDataTransfer::MaxCellFluxImbalance(algebraic) < 1.e-10,
                "Transferred flow is not cellwise conservative");
        const int inlet = reverse ? Right : Left, outlet = reverse ? Left : Right;
        auto total = [&](int id) {
            REAL sum = 0.;
            for (REAL value : algebraic.fInterfaceData.at(id).fIntegralFlux) sum += value;
            return sum;
        };
        Require(std::abs(total(inlet) + 1.) < 1.e-10 && std::abs(total(outlet) - 1.) < 1.e-10,
                "Transferred flow has wrong magnitude or orientation");
        const REAL initialWater = algebraic.WaterVolume();
        REAL integratedOutflow = 0.;
        const std::string output = "reservoir-transport-" + std::to_string(reverse) + "-" + std::to_string(scenario);
        analysis.PostProcessTimeStep(output, 0);
        for (int step = 1; step <= settings.numSteps; step++) {
            std::cout << "Transport direction " << (reverse ? "reverse" : "forward")
                      << ", scenario " << scenario << ", step " << step << '\n';
            analysis.RunTimeStep();
            integratedOutflow += settings.timeStep * algebraic.BoundaryWaterRate(problem);
            Require(std::abs(algebraic.WaterVolume() - initialWater + integratedOutflow) < 1.e-9,
                    "Water inventory does not match integrated boundary flow");
            for (size_t i = 0; i < cells.fVolume.size(); i++) {
                const REAL sw = cells.fSaturation[i];
                Require(sw >= -1.e-11 && sw <= 1. + 1.e-11, "Saturation outside [0,1]");
                if (constant) Require(std::abs(sw - initial) < 1.e-10, "Constant composition not preserved");
                if (scenario == 1 && step == 1) {
                    // Exact discrete backward-Euler upwind solution for this 1D flow.
                    auto* gel = geometry->Element(cells.fGeoIndex[i]);
                    TPZManVector<REAL, 3> qsi(3, 0.), x(3, 0.);
                    gel->X(qsi, x);
                    int column = std::min(nx - 1, int(x[0] * nx));
                    if (reverse) column = nx - 1 - column;
                    const REAL courant = settings.timeStep * nx / problem.m_Reservoir.porosity;
                    const REAL expected = std::pow(courant / (1. + courant), column + 1);
                    Require(std::abs(sw - expected) < 1.e-9, "First step differs from discrete analytical solution");
                }
            }
            analysis.PostProcessTimeStep(output, step);
        }
        const TPZFMatrix<STATE> currentFlow = flow->Solution();
        Require(Norm(currentFlow - frozenSolution) < 1.e-14, "Transport modified the frozen Darcy solution");
        std::cout << "PASS: " << (reverse ? "reverse" : "forward") << ", "
                  << (constant ? "constant composition" : quadratic ? "nonlinear oil-water" : "linear tracer")
                  << ", water-balance error = " << std::abs(algebraic.WaterVolume() - initialWater + integratedOutflow) << '\n';
    }
    // A rejected nonlinear step must neither advance history nor retain a trial state.
    settings.timeStep = 0.05;
    problem.m_Numerics.NewtonTransport.maxIterations = 1;
    TPZWannTransportAnalysis rejected(transport.get());
    rejected.Initialize(flow.get(), problem);
    const TPZFMatrix<STATE> previous = rejected.Solution();
    bool failed = false;
    std::cout << "Transport rollback check (expected Newton failure)\n";
    try { rejected.RunTimeStep(); } catch (const std::runtime_error&) { failed = true; }
    const TPZFMatrix<STATE> restored = rejected.Solution();
    Require(failed && Norm(restored - previous) < 1.e-14 &&
            rejected.fAlgebraicTransport.WaterVolume() == 0., "Rejected step did not restore state");
}
}

int main() {
    try {
        CheckDirection(false);
        CheckDirection(true);
        std::cout << "All reservoir transport checks passed.\n";
        return 0;
    } catch (const std::exception& error) {
        std::cerr << "Reservoir transport check failed: " << error.what() << '\n';
        return 1;
    }
}
