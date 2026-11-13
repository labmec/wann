/**
 * @file
 * @brief Class and method contracts for TPZDarcyAnisotropic.
 */

#pragma once

#include "TPZMatBase.h"
#include "TPZMatSingleSpace.h"
#include "TPZMatErrorSingleSpace.h"
#include "TPZAnisotropicPermeability.h"

/**
 * @brief Pressure-only H1 Darcy material with tensor-valued permeability.
 * Assembles the linear weak form of -div(K grad(p)) = f. Flux is reconstructed
 * as -K grad(p). Assembly multiplies the 3x3 tensor directly by data.dphix;
 * the supplied derivative matrix must have compatible rows and coordinates.
 * @see TPZAnisotropicPermeability
 * @see TPZWannDarcyNL
 */
class TPZDarcyAnisotropic : public TPZMatBase<STATE, TPZMatSingleSpaceT<STATE>,
        TPZMatErrorSingleSpace<STATE>, TPZAnisotropicPermeability> {

    // type alias to improve constructor readability
    using TBase = TPZMatBase<STATE, TPZMatSingleSpaceT<STATE>, TPZMatErrorSingleSpace<STATE>, TPZAnisotropicPermeability>;

public:

    /**
     * @brief Construct an object; configure simulation/material data before use.
     */
    TPZDarcyAnisotropic();

    /**
     * @brief Construct a material with the given identifier and dimension.
     */
    TPZDarcyAnisotropic(int id, int dim);

            /**
             * @brief Construct a copy using the existing base/member copy semantics.
             */
            TPZDarcyAnisotropic(const TPZDarcyAnisotropic &copy);

            /**
             * @brief Assign configuration using the existing base/member copy semantics.
             */
            TPZDarcyAnisotropic& operator=(const TPZDarcyAnisotropic &copy);

    /**
     * @brief Return the concrete material name used in diagnostics.
     */
    [[nodiscard]] std::string Name() const override { return "TPZDarcyAnisotropic"; }

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
    int NEvalErrors() const override { return 3; }

    /**
     * @brief Set the problem dimension; values outside 1 through 3 trigger DebugStop().
     */
    virtual void SetDimension(int dim);

    /**
     * @brief Add the linear pressure stiffness and source load at one quadrature point.
     * @param[in] data Scalar pressure shape values and compatible derivatives.
     * @param[in] weight Physical integration weight supplied by NeoPZ.
     * @param[in,out] ek Accumulates grad(phi_i)^T K grad(phi_j) contributions.
     * @param[in,out] ef Accumulates the forcing-function load; zero if absent.
     */
    void Contribute(const TPZMaterialDataT<STATE> &data, STATE weight, TPZFMatrix<STATE> &ek,
                    TPZFMatrix<STATE> &ef) override;

    /**
     * @brief Add type-0 pressure penalty, type-1 load, or type-2 Robin terms.
     * BC callbacks supply pressure and derivatives used to form the boundary value.
     * ek and ef are accumulated, not cleared.
     */
    void ContributeBC(const TPZMaterialDataT<STATE> &data, STATE weight, TPZFMatrix<STATE> &ek, TPZFMatrix<STATE> &ef,
                      TPZBndCondT<STATE> &bc) override;

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
    void Solution(const TPZMaterialDataT<STATE> &data, int var, TPZVec<STATE> &solOut) override;

    /**
     * @brief Set u_len = 1, du_row = fDim, and du_col = 1 for pressure data.
     */
    void GetSolDimensions(uint64_t &u_len, uint64_t &du_row, uint64_t &du_col) const override;

    /**
     * @brief Fill three pointwise squared errors: pressure-plus-energy, pressure L2,
     * and permeability-weighted gradient error. Requires an exact callback and an
     * errors vector already sized to NEvalErrors(); integration is done by NeoPZ.
     */
    void Errors(const TPZMaterialDataT<STATE> &data, TPZVec<REAL> &errors) override;

    /**
     * @brief Resize names and label the three H1/L2/gradient error slots.
     */
    virtual void ErrorNames(TPZVec<std::string> &names) const override;

    /**
     * @brief Clear optional volume data requirements for linear stiffness/load assembly.
     */
    void FillDataRequirements(TPZMaterialData &data) const override;

    /**
     * @brief Set optional boundary requirements: solution for type 50, normals for
     * types 1/3 or when a volume forcing function is installed.
     */
    void FillBoundaryConditionDataRequirements(int type, TPZMaterialData &data) const override;

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
