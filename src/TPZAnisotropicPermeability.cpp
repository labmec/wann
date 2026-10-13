//
// Created by Gi Taraschi 27/05/26
// Tensors are stored in global 3D coordinates; 1D/2D inputs use the x-axis/xy-plane.
//

#include "TPZAnisotropicPermeability.h"
#include <cmath>

namespace {
// Preserve three-component flux multiplication while accepting lower-dimensional tensors.
TPZFMatrix<STATE> NormalizePermeability(TPZFMatrix<STATE> tensor) {
    const auto size = tensor.Rows();
    if (tensor.Cols() != size || (size < 1 || size > 3)) {
        PZError << "Permeability tensor must be 1x1, 2x2 or 3x3.";
        DebugStop();
    }
    if (size == 1) {
        const STATE permeability = tensor(0, 0);
        if (permeability <= 0.) {
            PZError << "Invalid 1D permeability value.";
            DebugStop();
        }
        TPZFMatrix<STATE> embedded(3, 3, 0.);
        for (int i = 0; i < 3; i++) {
            embedded(i, i) = permeability;
        }
        return embedded;
    }
    if (size == 2) {
        const STATE determinant = tensor(0, 0) * tensor(1, 1)
                                - tensor(0, 1) * tensor(1, 0);
        if (determinant <= 0.) {
            PZError << "Invalid 2D permeability matrix determinant.";
            DebugStop();
        }
        TPZFMatrix<STATE> embedded(3, 3, 0.);
        for (int i = 0; i < 2; i++) {
            for (int j = 0; j < 2; j++) {
                embedded(i, j) = tensor(i, j);
            }
        }
        embedded(2, 2) = std::sqrt(determinant);
        return embedded;
    }
    return tensor;
}
}


void TPZAnisotropicPermeability::SetConstantPermeability(const STATE constant) {
    fConstantPermeabilityScalar = constant;
    fConstantPermeability.Redim(3, 3);
    fConstantInversePermeability.Redim(3, 3);
    fConstantPermeability = 0.;
    fConstantInversePermeability = 0.;
    for (int i = 0; i < 3; i++) {
        fConstantPermeability(i, i) = constant;
        fConstantInversePermeability(i, i) = 1. / constant;
    }
    fHomogeneous = true;
    fIsotropic = true;
}

void TPZAnisotropicPermeability::SetConstantPermeability(TPZFMatrix<STATE> constant) {
    constant = NormalizePermeability(constant);
    fConstantPermeability = constant;
    fHomogeneous = true;
    fIsotropic = false;
    STATE determinant = constant(0, 0) * (constant(1, 1) * constant(2, 2) - constant(1, 2) * constant(2, 1))
                       -constant(0, 1) * (constant(1, 0) * constant(2, 2) - constant(1, 2) * constant(2, 0))
                       +constant(0, 2) * (constant(1, 0) * constant(2, 1) - constant(1, 1) * constant(2, 0));
    if (determinant <= 0.) {
        PZError << "Invalid permeability matrix determinant.";
        DebugStop();
    }
    fConstantPermeabilityScalar = pow(determinant, 1./3.);
    TPZFMatrix<STATE> aux_matrix = constant;
    aux_matrix.Inverse(fConstantInversePermeability, ELU);
}

void TPZAnisotropicPermeability::SetPermeabilityFunction(IsotropicFunctionType &perm_function) {
    fIsotropicPermeabilityFunction = perm_function;
    fHomogeneous = false;
    fIsotropic = true;
}

void TPZAnisotropicPermeability::SetPermeabilityFunction(AnisotropicFunctionType &perm_function) {
    fAnisotropicPermeabilityFunction = perm_function;
    fHomogeneous = false;
    fIsotropic = false;
}

STATE TPZAnisotropicPermeability::GetPermeability(const TPZVec<REAL> &coord) {
    if (fIsotropic) {
        return fHomogeneous ? fConstantPermeabilityScalar : fIsotropicPermeabilityFunction(coord);
    } else {
        std::cout << "Anisotropic permeability can not be represented as a scalar." << std::endl;
        return 0.;
    }
}

void TPZAnisotropicPermeability::GetPermeability(const TPZVec<REAL> &coord, TPZFMatrix<STATE> &perm_matrix) {
    perm_matrix.Redim(3, 3);
    perm_matrix = 0;
    if (fIsotropic) {
        STATE perm_value = fHomogeneous ? fConstantPermeabilityScalar : fIsotropicPermeabilityFunction(coord);
        for (int i = 0; i < 3; i++) {
            perm_matrix(i, i) = perm_value;
        }
    } else {
        if (fHomogeneous) {
            perm_matrix = fConstantPermeability;
        } else {
            perm_matrix = NormalizePermeability(fAnisotropicPermeabilityFunction(coord));
        }
    }
}

void TPZAnisotropicPermeability::GetInversePermeability(const TPZVec<REAL> &coord, TPZFMatrix<STATE> &inv_perm_matrix) {
    if (fHomogeneous) {
        inv_perm_matrix = fConstantInversePermeability;
    } else {
        TPZFMatrix<STATE> perm_matrix;
        GetPermeability(coord, perm_matrix);
        perm_matrix.Inverse(inv_perm_matrix, ELU);
    }
}

int TPZAnisotropicPermeability::ClassId() const {
    return Hash("TPZAnisotropicPermeability");
}
