//
// Created by Gi Taraschi 27/05/26
//

#include "TPZAnisotropicPermeability.h"

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
    if (perm_matrix.Cols() != perm_matrix.Rows() && perm_matrix.Cols() != 3) {
        PZError << "Permeability matrix must be 3x3.";
        DebugStop();
    }

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
            perm_matrix = fAnisotropicPermeabilityFunction(coord);
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
