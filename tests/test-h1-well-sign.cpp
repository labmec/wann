#include "TPZWannH1Well.h"
#include "TPZMaterialDataT.h"
#include <algorithm>
#include <cmath>
#include <iostream>

int main() {
    TPZWannH1Well material(1, 0.2, 0.001, 1000., 0., 0.);
    TPZMaterialDataT<STATE> data;
    data.fShapeType = TPZMaterialData::EScalarShape;
    data.phi.Redim(1, 1);
    data.phi(0, 0) = 1.;
    data.dphix.Redim(1, 1);
    data.dphix(0, 0) = 1.;
    data.sol.Resize(1);
    data.sol[0].Resize(1);
    data.sol[0][0] = 0.;
    data.dsol.Resize(1);
    data.dsol[0].Redim(1, 1);

    auto close = [](REAL a, REAL b) {
        return std::abs(a-b) <= 1.e-6 * std::max({REAL(1.e-10), std::abs(a), std::abs(b)});
    };
    auto residual = [&](REAL gradient) {
        data.dsol[0](0, 0) = gradient;
        TPZFMatrix<STATE> ef(1, 1, 0.);
        material.ContributeResidual(data, 1., ef);
        return ef(0, 0);
    };

    for (bool first : {true, false}) {
        TPZWannH1Well::fIsFirstIteration = first;
        for (REAL gradient : {-1., -1.e-8, 0., 1.e-8, 1.}) {
            data.dsol[0](0, 0) = gradient;
            TPZFMatrix<STATE> ek(1, 1, 0.), ef(1, 1, 0.);
            material.Contribute(data, 1., ek, ef);
            TPZManVector<STATE, 1> flux;
            material.Solution(data, material.VariableIndex("Flux"), flux);
            if (!close(ef(0, 0), flux[0]) ||
                !close(residual(gradient), ef(0, 0)) || ek(0, 0) <= 0. ||
                (gradient != 0. && flux[0]*gradient >= 0.)) {
                std::cerr << "Residual/flux mismatch at gradient " << gradient << '\n';
                return 1;
            }
            const REAL step = std::max(REAL(1.e-12), std::abs(gradient)*1.e-5);
            const REAL tangent = -(residual(gradient+step)-residual(gradient-step))/(2.*step);
            if (!close(tangent, ek(0, 0))) {
                std::cerr << "Tangent mismatch at gradient " << gradient << '\n';
                return 1;
            }
        }
    }
    TPZWannH1Well::fIsFirstIteration = true;
    std::cout << "H1 well: signed flux, both residual paths and finite-difference tangent passed\n";
}
