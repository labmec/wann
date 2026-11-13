/**
 * @file
 * @brief Class and method contracts for TPZAnisotropicPermeability.
 */

//
// Created by Gi Taraschi 27/05/26
//

#pragma once

#include <functional>
#include "pzreal.h"
#include "pzvec.h"
#include "pzfmatrix.h"

// Alias to improve readability of the permeability function type
/**
 * @brief Spatial callback returning scalar permeability at physical coordinates.
 */
using IsotropicFunctionType = std::function<STATE(const TPZVec<REAL> &coord)>;
/**
 * @brief Spatial callback returning a 1x1, 2x2, or 3x3 permeability tensor.
 */
using AnisotropicFunctionType = std::function<TPZFMatrix<STATE>(const TPZVec<REAL> &coord)>;

// Forward declaration of dummy BC interface class
class TPZAnisotropicPermeabilityBC;

/**
 * @brief Constant or spatially varying permeability provider for Darcy materials.
 * Stores and returns 3x3 tensors in global coordinates without a space-dimension
 * member. Scalar input gives k I; 1x1 input uses global x and is embedded as k I;
 * 2x2 input uses global xy and gets Kzz = sqrt(det(K2D)), with zero z coupling.
 * Call a setter before relying on the cached inverse. Invertible, physically
 * admissible inputs are the caller's responsibility; determinant checks do not
 * constitute a complete positive-definiteness test.
 * @see TPZDarcyAnisotropic
 * @see TPZMixedDarcyAnisotropic
 */
class TPZAnisotropicPermeability : public virtual TPZSavable {

public:
    using TInterfaceBC = TPZAnisotropicPermeabilityBC;

    /**
     * @brief Construct an object; configure simulation/material data before use.
     */
    TPZAnisotropicPermeability() : fConstantPermeability(1.), fConstantInversePermeability(1.), fIsotropicPermeabilityFunction(NULL), fAnisotropicPermeabilityFunction(NULL){
    }

    /**
     * @brief Construct a copy using the existing base/member copy semantics.
     */
    TPZAnisotropicPermeability(const TPZAnisotropicPermeability &copy) = default;
    /**
     * @brief Assign configuration using the existing base/member copy semantics.
     */
    TPZAnisotropicPermeability &operator=(const TPZAnisotropicPermeability &copy) = default;

    /**
     * @brief Configure homogeneous scalar permeability k I3 and cache its inverse.
     * @param[in] constant Positive scalar permeability; no positivity check is performed.
     */
    void SetConstantPermeability(STATE constant);

    /**
     * @brief Configure a constant tensor, normalize it to 3x3, and cache its LU inverse.
     * @param[in] constant Square tensor of size 1, 2, or 3 in the documented global frame.
     * Size and determinant checks can trigger DebugStop(). The field is marked
     * anisotropic even when the supplied tensor is a multiple of the identity.
     */
    void SetConstantPermeability(TPZFMatrix<STATE> constant);

    /**
     * @brief Copy a spatial tensor callback and select nonhomogeneous anisotropic mode.
     * @param[in] perm_function Callback returning a square 1x1, 2x2, or 3x3 tensor.
     * Normalization/inversion occur at evaluation; captured references must remain valid.
     */
    void SetPermeabilityFunction(AnisotropicFunctionType &perm_function);

    /**
     * @brief Copy a spatial scalar callback and select nonhomogeneous isotropic mode.
     * @param[in] perm_function Callback returning positive scalar permeability.
     * Captured references must remain valid for subsequent evaluation.
     */
    void SetPermeabilityFunction(IsotropicFunctionType &perm_function);

    /**
     * @brief Return scalar permeability for scalar-configured fields.
     * @param[in] coord Physical coordinates used for a spatial callback.
     * @return Constant/callback scalar, or zero with a diagnostic for tensor fields.
     * Does not return an average of an anisotropic tensor.
     */
    STATE GetPermeability(const TPZVec<REAL> &coord);

    /**
     * @brief Return the cached constant inverse or evaluate and invert a spatial tensor.
     * @param[in] coord Physical coordinates passed to the permeability callback.
     * @param[out] inv_perm_matrix Replaced with the 3x3 inverse tensor.
     */
    void GetInversePermeability(const TPZVec<REAL> &coord, TPZFMatrix<STATE> &inv_perm_matrix);

    /**
     * @brief Evaluate permeability as a global 3x3 tensor, with lower-dimensional embedding.
     * @param[in] coord Physical coordinates used for a spatial callback.
     * @param[out] perm_matrix Resized/replaced with the evaluated 3x3 tensor.
     */
    void GetPermeability(const TPZVec<REAL> &coord, TPZFMatrix<STATE> &perm_matrix);

    /**
     * @brief Return the NeoPZ runtime class identifier.
     */
    [[nodiscard]] int ClassId() const override;

    /**
     * @brief No-op serialization hook; this class currently restores no permeability state.
     */
    void Read(TPZStream &buf, void *context) override {};

    /**
     * @brief No-op serialization hook; this class currently writes no permeability state.
     */
    void Write(TPZStream &buf, int withclassid) const override {};

private:

    // Member variable to describe a constant permeability field
    STATE fConstantPermeabilityScalar = 1.;
    TPZFMatrix<STATE> fConstantPermeability;
    TPZFMatrix<STATE> fConstantInversePermeability;
    bool fHomogeneous = true;
    bool fIsotropic = true;

    // Member variables to describe a varying permeability field
    IsotropicFunctionType fIsotropicPermeabilityFunction = NULL;
    AnisotropicFunctionType fAnisotropicPermeabilityFunction = NULL;
};

// Dummy BC interface class
class TPZMaterial;
/**
 * @brief Empty NeoPZ boundary interface for permeability-bearing materials.
 * Does not check compatibility or add its own boundary assembly.
 */
class TPZAnisotropicPermeabilityBC : public TPZAnisotropicPermeability {
protected:
    // this method is your chance to verify if the material to which this
    // BC interface applies is compatible with this boundary interface
    // it is called in the method SetMaterial of class TPZBndCondBase
    /**
     * @brief No-op compatibility hook called by NeoPZ boundary-material binding.
     */
    static void SetMaterialImpl(TPZMaterial *) {}
};
