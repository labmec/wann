/**
 * @file
 * @brief Class and method contracts for TPZMixedDarcyAnisotropic.
 */

#pragma once

#include "TPZMatBase.h"
#include "TPZMatCombinedSpaces.h"
#include "TPZMatErrorCombinedSpaces.h"
#include "TPZAnisotropicPermeability.h"

/**
 * @brief Mixed H(div)-flux / scalar-pressure Darcy material.
 * Uses datavec[0] for flux and datavec[1] for pressure. Local unknowns are
 * ordered as flux followed by pressure. At volume quadrature points,
 * fDeformedDirections is 3 x nFlux, phi for pressure is nPressure x 1,
 * and divphi is nFlux x 1. Permeability inverses remain 3x3 in 1D/2D/3D.
 * @see TPZAnisotropicPermeability
 * @see TPZWannMixedDarcyNL
 */
class TPZMixedDarcyAnisotropic : public TPZMatBase<STATE, TPZMatCombinedSpacesT<STATE>,
        TPZMatErrorCombinedSpaces<STATE>, TPZAnisotropicPermeability> {

    using TBase = TPZMatBase<STATE, TPZMatCombinedSpacesT<STATE>,
            TPZMatErrorCombinedSpaces<STATE>, TPZAnisotropicPermeability>;

public:
    /**
     * @brief Construct an object; configure simulation/material data before use.
     */
    TPZMixedDarcyAnisotropic();
    /**
     * @brief Construct a material with the given identifier and dimension.
     */
    [[maybe_unused]] TPZMixedDarcyAnisotropic(int id, int dim);
    /**
     * @brief Construct a copy using the existing base/member copy semantics.
     */
    TPZMixedDarcyAnisotropic(const TPZMixedDarcyAnisotropic &copy);
    /**
     * @brief Assign configuration using the existing base/member copy semantics.
     */
    TPZMixedDarcyAnisotropic &operator=(const TPZMixedDarcyAnisotropic &copy);

    /**
     * @brief Return the concrete material name used in diagnostics.
     */
    [[nodiscard]] std::string Name() const override { return "TPZMixedDarcyAnisotropic"; }
    /**
     * @brief Return the configured topological problem dimension.
     */
    [[nodiscard]] int Dimension() const override { return this->fDim; }
    /**
     * @brief Return one scalar state per approximation space.
     */
    [[nodiscard]] int NStateVariables() const override { return 1; }
    /**
     * @brief Return the number of error slots declared by this material.
     */
    int NEvalErrors() const override { return 5; }

    /**
     * @brief Set the problem dimension; values outside 1 through 3 trigger DebugStop().
     */
    virtual void SetDimension(int dim);

    /**
     * @brief Add the mixed Darcy matrix and pressure-equation source load.
     * @param[in] datavec Flux at index 0, scalar pressure at index 1.
     * @param[in] weight Physical integration weight supplied by NeoPZ.
     * @param[in,out] ek Flux-first matrix; A_ij = v_i^T K^{-1} v_j and
     * off-diagonal entries are -phi_j div(v_i), multiplied by weight.
     * @param[in,out] ef Source load in pressure rows; flux rows precede them.
     */
    void Contribute(const TPZVec<TPZMaterialDataT<STATE>> &datavec, REAL weight, TPZFMatrix<STATE> &ek,
                    TPZFMatrix<STATE> &ef) override;

    /**
     * @brief Add pressure loads (type 0), normal-flux penalties (type 1), or
     * Robin/special terms (types 2/4/5). Boundary phi is the flux trace basis,
     * not the volume fDeformedDirections vector matrix.
     */
    void ContributeBC(const TPZVec<TPZMaterialDataT<STATE>> &datavec, REAL weight, TPZFMatrix<STATE> &ek,
                      TPZFMatrix<STATE> &ef, TPZBndCondT<STATE> &bc) override;

    /**
     * @brief Map a post-processing field name to its numeric variable index.
     * See docs/materials.md for supported fields and current implementation limits.
     */
    [[nodiscard]] int VariableIndex(const std::string &name) const override;
    /**
     * @brief Return the output component count for a numeric post-processing variable.
     */
    [[nodiscard]] int NSolutionVariables(int var) const override;
    /**
     * @brief Evaluate a requested post-processing field from the current point data.
     * Exact fields require an installed exact-solution callback. Registered names
     * do not all have a corresponding implemented output; see docs/materials.md.
     */
    void Solution(const TPZVec<TPZMaterialDataT<STATE>> &datavec, int var, TPZVec<STATE> &solOut) override;
    /**
     * @brief Resize/fill five pointwise squared errors: pressure, flux, divergence
     * residual, pressure gradient, and flux-plus-divergence (H(div)). Exact pressure
     * and gradient come from fExactSol; exact divergence comes from fForcingFunction.
     */
    void Errors(const TPZVec<TPZMaterialDataT<STATE>> &data, TPZVec<REAL> &errors) override;
    /**
     * @brief Clear optional volume requirements in every approximation-space data entry.
     */
    void FillDataRequirements(TPZVec<TPZMaterialDataT<STATE> > &datavec) const override;
    /**
     * @brief Clear optional boundary requirements and request the flux-space normal.
     */
    void FillBoundaryConditionDataRequirements(int type, TPZVec<TPZMaterialDataT<STATE> > &datavec) const override;
    /**
     * @brief Return the NeoPZ runtime class identifier.
     */
    [[nodiscard]] int ClassId() const override;
    /**
     * @brief Allocate a material copy using the concrete copy constructor.
     * @return A new material pointer for ownership by the caller/NeoPZ mesh.
     */
    [[nodiscard]] TPZMaterial *NewMaterial() const override;
    /**
     * @brief Write the stored configuration or material summary to out.
     */
    void Print(std::ostream & out) const override;

protected:
    /// Topological dimension, independent of 3x3 tensor storage.
    int fDim;
};
