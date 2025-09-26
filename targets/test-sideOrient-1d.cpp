#ifdef HAVE_CONFIG_H
#include <pz_config.h>
#endif

#include <iostream>
#include "TPZWannGeometryTools.h"
#include "TPZWannApproxTools.h"
#include "TPZWannPostProcTools.h"
#include "TPZWannEstimationTools.h"
#include <TPZLinearAnalysis.h>
#include <TPZSSpStructMatrix.h> //symmetric sparse matrix storage>
#include <pzskylstrmatrix.h>
#include <pzstepsolver.h>
#include <TPZAnalyticSolution.h>

using namespace std;

const int global_nthread = 0;

TPZMultiphysicsCompMesh* CreateMultiphysicsCompMesh(TPZGeoMesh * gmesh, ProblemData *SimData, TPZAnalyticSolution *exact);

int main(int argc, char *argv[])
{

    TLaplaceExample1 exact; // Global variable to be used in the material objects
    exact.fDimension = 3;
    exact.fExact = TLaplaceExample1::EConst;
    TLaplaceExample1::gC = -1.0;

    std::string jsonfile = "test-sideOrient-1d.json";

    if (argc > 2)
    {
        std::cout << argv[0] << " being called with too many arguments." << std::endl;
        DebugStop();
    }
    else if (argc == 2)
    {
        jsonfile = argv[1];
    }

    if (jsonfile.find(".json") == std::string::npos)
    {
        jsonfile += ".json";
    }
    std::cout << "Using json file: " << jsonfile << std::endl;
    std::cout << "--------- Starting simulation ---------" << std::endl;
#ifdef PZ_LOG
    TPZLogger::InitializePZLOG();
#endif

    // Problem data
    ProblemData SimData;
    SimData.ReadJson(jsonfile);

    TPZGeoMesh *gmesh = TPZWannGeometryTools::CreateGeoMesh(&SimData);

    TPZMultiphysicsCompMesh *cmesh = CreateMultiphysicsCompMesh(gmesh, &SimData, &exact);

    // ---- Hdiv analysis ----
    TPZLinearAnalysis analysis(cmesh);
#ifdef PZ_USING_MKL
    TPZSSpStructMatrix<STATE> skylstr(cmesh);
#else
    TPZSkylStrMatrix skylstr(cmesh);
#endif
    skylstr.SetNumThreads(global_nthread);
    analysis.SetStructuralMatrix(skylstr);

    TPZStepSolver<STATE> step;
    step.SetDirect(ELDLt);
    analysis.SetSolver(step);
    analysis.Run();

    std::cout << "--------- Simulation finished ---------" << std::endl;

    TPZWannPostProcTools::PostProcessAllData(cmesh, gmesh, &SimData);

    delete cmesh;

    delete gmesh;
}

TPZMultiphysicsCompMesh* CreateMultiphysicsCompMesh(TPZGeoMesh * gmesh, ProblemData *SimData, TPZAnalyticSolution *exact)
{
    const int dim = 1;

    TPZHDivApproxCreator hdivCreator(gmesh);
    hdivCreator.ProbType() = ProblemType::EDarcy;
    hdivCreator.SetDefaultOrder(SimData->m_Wellbore.pOrder);
    hdivCreator.SetShouldCondense(false);

    TLaplaceExample1 *exactsol = dynamic_cast<TLaplaceExample1 *>(exact);
    bool hasAnalyticSol = (exactsol != nullptr && exactsol->fExact != TLaplaceExample1::ENone);

    // Wellbore material
    auto &WellboreData = SimData->m_Wellbore;
    {
        TPZMixedDarcyFlow *wellboreMat = new TPZMixedDarcyFlow(SimData->ECurveWell, dim);
        if (hasAnalyticSol)
        {
            wellboreMat->SetExactSol(exact->ExactSolution(), 3);
            wellboreMat->SetForcingFunction(exact->ForceFunc(), 3);
            wellboreMat->SetConstantPermeability(1.0);
        }
        else
        {
            wellboreMat->SetConstantPermeability(WellboreData.perm);
        }
        hdivCreator.InsertMaterialObject(wellboreMat);

        for (auto &bcpair : WellboreData.BCs)
        {
            auto &bc = bcpair.second;
            TPZFMatrix<STATE> val1(1, 1, 0.);
            TPZManVector<STATE> val2(1, 0);
            val2[0] = bc.value;
            TPZBndCondT<STATE> *BCond = wellboreMat->CreateBC(wellboreMat, bc.matid, bc.type, val1, val2);
            if (hasAnalyticSol)
                BCond->SetForcingFunctionBC(exact->ExactSolution(), 3);
            hdivCreator.InsertMaterialObject(BCond);
        }
    }

    int lagmultilevel = 1;
    TPZManVector<TPZCompMesh *, 7> meshvec(hdivCreator.NumMeshes());
    hdivCreator.CreateAtomicMeshes(meshvec, lagmultilevel);      // This method increments the lagmultilevel

    TPZMultiphysicsCompMesh *cmesh = nullptr;
    hdivCreator.CreateMultiPhysicsMesh(meshvec, lagmultilevel, cmesh);

    if (SimData->m_VerbosityLevel)
    {
        std::ofstream out("cmesh.txt");
        cmesh->Print(out);
    }

    return cmesh;
}