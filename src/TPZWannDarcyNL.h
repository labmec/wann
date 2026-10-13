/**
 * @file
 * @brief Class and method contracts for TPZWannDarcyNL.
 */

#pragma once

#include "TPZDarcyAnisotropic.h"

/**
 * @brief H1 reservoir material assembled as a Newton tangent and negative residual.
 * The current permeability law is spatial, not pressure dependent; the NL
 * suffix refers to participation in the coupled nonlinear solve.
 * Pressure Dirichlet values must be imposed by the analysis before assembly.
 * @see TPZWannAnalysis
 * @see TPZDarcyAnisotropic
 */
class TPZWannDarcyNL : public TPZDarcyAnisotropic {

    using TBase = TPZDarcyAnisotropic;

public:

    /**
     * @brief Construct an object; configure simulation/material data before use.
     */
    TPZWannDarcyNL();

    /**
     * @brief Construct a material with the given identifier and dimension.
     */
    [[maybe_unused]] TPZWannDarcyNL(int id, int dim);

    /**
     * @brief Construct a copy using the existing base/member copy semantics.
     */
    TPZWannDarcyNL(const TPZWannDarcyNL &copy);

    /**
     * @brief Assign configuration using the existing base/member copy semantics.
     */
    TPZWannDarcyNL &operator=(const TPZWannDarcyNL &copy);

    /**
     * @brief Return the concrete material name used in diagnostics.
     */
    [[nodiscard]] std::string Name() const override { return "TPZWannDarcyNL"; }

    /**
     * @brief Return the NeoPZ runtime class identifier.
     */
    [[nodiscard]] int ClassId() const override { return Hash("TPZWannDarcyNL") ^ (TBase::ClassId() << 1); }

    /**
     * @brief Delegate stream restoration to the base material implementation.
     */
    void Read(TPZStream &buf, void *context) override { TBase::Read(buf, context); }

    /**
     * @brief Delegate stream serialization to the base material implementation.
     */
    void Write(TPZStream &buf, int withclassid) const override { TBase::Write(buf, withclassid); }

    /**
     * @brief Forward the constant scalar/tensor configuration to TPZAnisotropicPermeability.
     */
    void SetConstantPermeability(STATE constant)
    {
        TPZAnisotropicPermeability::SetConstantPermeability(constant);
    }

    /**
     * @brief Forward the constant scalar/tensor configuration to TPZAnisotropicPermeability.
     */
    void SetConstantPermeability(TPZFMatrix<STATE> constant)
    {
        TPZAnisotropicPermeability::SetConstantPermeability(constant);
    }

    /**
     * @brief Forward the spatial permeability callback to TPZAnisotropicPermeability.
     */
    void SetPermeabilityFunction(AnisotropicFunctionType &perm_function)
    {
        TPZAnisotropicPermeability::SetPermeabilityFunction(perm_function);
    }

    /**
     * @brief Forward the spatial permeability callback to TPZAnisotropicPermeability.
     */
    void SetPermeabilityFunction(IsotropicFunctionType &perm_function)
    {
        TPZAnisotropicPermeability::SetPermeabilityFunction(perm_function);
    }

    /**
     * @brief Forward the scalar or tensor query to TPZAnisotropicPermeability.
     */
    STATE GetPermeability(const TPZVec<REAL> &coord)
    {
        return TPZAnisotropicPermeability::GetPermeability(coord);
    }

    /**
     * @brief Forward the scalar or tensor query to TPZAnisotropicPermeability.
     */
    void GetPermeability(const TPZVec<REAL> &coord, TPZFMatrix<STATE> &perm_matrix)
    {
        TPZAnisotropicPermeability::GetPermeability(coord, perm_matrix);
    }

    /**
     * @brief Accumulate the pressure tangent and negative residual at the current solution.
     * data.sol/data.dsol are required. ek uses the linear Darcy stiffness; ef adds
     * -K grad(p) against test gradients and the optional volume source load.
     * @param[in] data Pressure shape/solution data at the integration point.
     * @param[in] weight Physical integration weight supplied by NeoPZ.
     * @param[in,out] ek Element tangent matrix receiving additive contributions.
     * @param[in,out] ef Element negative-residual vector receiving additive contributions.
     */
    void Contribute(const TPZMaterialDataT<STATE> &data, REAL weight, TPZFMatrix<STATE> &ek,
                    TPZFMatrix<STATE> &ef) override;

    /**
     * @brief Add the negative type-1 boundary load; type-0 pressure is imposed by
     * TPZWannAnalysis and contributes nothing here. Other types trigger DebugStop().
     * The callback branch uses the scalar permeability getter.
     * @param[in] data Pressure shape/solution data at the integration point.
     * @param[in] weight Physical integration weight supplied by NeoPZ.
     * @param[in,out] ek Element tangent matrix receiving additive contributions.
     * @param[in,out] ef Element negative-residual vector receiving additive contributions.
     * @param[in] bc Boundary type and prescribed values.
     */
    void ContributeBC(const TPZMaterialDataT<STATE> &data, REAL weight, TPZFMatrix<STATE> &ek,
                      TPZFMatrix<STATE> &ef, TPZBndCondT<STATE> &bc) override;

    /**
     * @brief Request the current solution for residual assembly.
     */
    void FillDataRequirements(TPZMaterialData &data) const override;

    /**
     * @brief Request boundary normals and clear the solution requirement.
     */
    void FillBoundaryConditionDataRequirements(int type, TPZMaterialData &data) const override;

    /**
     * @brief Set u_len = 1, du_row = fDim, and du_col = 1 for pressure data.
     */
    void GetSolDimensions(uint64_t &u_len, uint64_t &du_row, uint64_t &du_col) const override;

    /**
     * @brief Evaluate a requested post-processing field from the current point data.
     * Exact fields require an installed exact-solution callback. Registered names
     * do not all have a corresponding implemented output; see docs/materials.md.
     */
    void Solution(const TPZMaterialDataT<STATE> &data, int var, TPZVec<STATE> &solOut) override;

};
