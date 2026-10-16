#include "TPZWannAnalysis.h"
#include "TPZHDivApproxCreator.h"
#include "DarcyFlow/TPZMixedDarcyFlow.h"
#include "pzgeoelbc.h"
#include "pzmultiphysicscompel.h"
#include "pzinterpolationspace.h"
#include <cmath>
#include <iostream>

int main() {
  for (REAL angle : {0., 45., 90., 135., 225., 315.}) {
    const REAL radians = angle * std::acos(-1.) / 180.;
    const REAL tx = std::cos(radians), ty = std::sin(radians);
    TPZGeoMesh geometry;
    geometry.SetDimension(1);
    geometry.NodeVec().Resize(4);
    ProblemData sim;
    sim.m_Wellbore.Resize(2);
    TPZGeoEl *lines[2];
    for (int w = 0; w < 2; ++w) {
      auto &well = sim.m_Wellbore[w];
      well.matid = w + 1;
      well.BCs["point_heel"] = {-10-2*w, 1, 0.};
      well.BCs["point_toe"] = {-11-2*w, 1, 0.};
      for (int end = 0; end < 2; ++end) {
        TPZManVector<REAL,3> x = {2.*end*tx, 2.*end*ty, REAL(w)};
        geometry.NodeVec()[2*w+end].Initialize(2*w+end, x, geometry);
      }
      // The second well has the opposite local node ordering.
      TPZManVector<int64_t,2> nodes = {2*w+w, 2*w+1-w};
      int64_t index;
      lines[w] = geometry.CreateGeoElement(EOned, nodes, well.matid, index);
    }
    geometry.BuildConnectivity();
    for (int w = 0; w < 2; ++w) {
      TPZGeoElBC heel(TPZGeoElSide(lines[w], w), -10-2*w);
      TPZGeoElBC toe(TPZGeoElSide(lines[w], 1-w), -11-2*w);
    }
    TPZHDivApproxCreator creator(&geometry);
    creator.ProbType() = ProblemType::EDarcy;
    creator.SetDefaultOrder(1);
    creator.SetShouldCondense(false);
    std::set<int> bcIds;
    for (auto &well : sim.m_Wellbore) {
      auto *mat = new TPZMixedDarcyFlow(well.matid, 1);
      mat->SetConstantPermeability(1.);
      creator.InsertMaterialObject(mat);
      for (auto &entry : well.BCs) {
        TPZFMatrix<STATE> v1(1,1,0.);
        TPZManVector<STATE,1> v2(1,0.);
        creator.InsertMaterialObject(mat->CreateBC(mat, entry.second.matid, 1, v1, v2));
        bcIds.insert(entry.second.matid);
      }
    }
    TPZManVector<TPZCompMesh*,7> meshes(creator.NumMeshes());
    int level = 1;
    creator.CreateAtomicMeshes(meshes, level);
    TPZMultiphysicsCompMesh *mesh = nullptr;
    creator.CreateMultiPhysicsMesh(meshes, level, mesh);
    {
      TPZWannAnalysis analysis(mesh, RenumType::ENone);
      analysis.SetProblemData(&sim);
      for (REAL heel : {-0.05, 0., 0.05}) for (REAL toe : {-0.07, 0., 0.07}) {
        for (int w = 0; w < 2; ++w) {
          sim.m_Wellbore[w].BCs.at("point_heel").value = (w+1)*heel;
          sim.m_Wellbore[w].BCs.at("point_toe").value = (w+1)*toe;
        }
        analysis.SetInitialSolution(bcIds);
        for (int w = 0; w < 2; ++w) {
          auto *mp = dynamic_cast<TPZMultiphysicsElement*>(lines[w]->Reference());
          auto *flux = dynamic_cast<TPZInterpolationSpace*>(mp->Element(0));
          for (int end = 0; end < 2; ++end) {
            TPZMaterialDataT<STATE> data;
            flux->InitMaterialData(data);
            data.fNeedsSol = true;
            TPZManVector<REAL,1> qsi(1, (end == w) ? -1. : 1.);
            flux->ComputeRequiredData(data, qsi);
            const REAL actual = data.sol[0][0]*tx + data.sol[0][1]*ty;
            const REAL expected = (w+1)*(end == 0 ? heel : toe);
            if (std::abs(actual-expected) > 1.e-10) {
              std::cerr << "Endpoint mismatch: angle=" << angle << " well=" << w
                        << " end=" << end << " actual=" << actual << " expected=" << expected << '\n';
              return 1;
            }
          }
        }
      }
    }
    delete mesh;
    for (auto *atomic : meshes) delete atomic;
  }
  std::cout << "Mixed endpoint flux checks passed.\n";
}
