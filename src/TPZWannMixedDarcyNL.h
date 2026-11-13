/**
 * @file
 * @brief Class and method contracts for TPZWannMixedDarcyNL.
 */

#pragma once

#include "TPZMixedDarcyAnisotropic.h"

/**
 * @brief Mixed reservoir material assembled as a Newton tangent and negative residual.
 * Uses the flux/pressure ordering of TPZMixedDarcyAnisotropic. The current
 * volume assembly has no forcing-function contribution and no pressure-dependent
 * permeability. Prescribed normal flux is imposed through analysis filtering.
 * @see TPZWannAnalysis
 */
class TPZWannMixedDarcyNL : public TPZMixedDarcyAnisotropic {

    using TBase = TPZMixedDarcyAnisotropic;

public:

    /**
     * @brief Construct an object; configure simulation/material data before use.
     */
    TPZWannMixedDarcyNL();

    /**
     * @brief Construct a material with the given identifier and dimension.
     */
    [[maybe_unused]] TPZWannMixedDarcyNL(int id, int dim);

    /**
     * @brief Construct a copy using the existing base/member copy semantics.
     */
    TPZWannMixedDarcyNL(const TPZWannMixedDarcyNL &copy);

    /**
     * @brief Assign configuration using the existing base/member copy semantics.
     */
    TPZWannMixedDarcyNL &operator=(const TPZWannMixedDarcyNL &copy);

    /**
     * @brief Return the concrete material name used in diagnostics.
     */
    [[nodiscard]] std::string Name() const override { return "TPZWannMixedDarcyNL"; }

    /**
     * @brief Return the NeoPZ runtime class identifier.
     */
    [[nodiscard]] int ClassId() const override { return Hash("TPZWannMixedDarcyNL") ^ (TBase::ClassId() << 1); }

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
     * @brief Accumulate the mixed Darcy tangent and negative current-solution residual.
     * Uses datavec[0].sol/divsol and datavec[1].sol with flux-first ordering.
     * The current implementation adds no volume forcing-function term.
     * @param[in] datavec Flux data at index 0 and pressure data at index 1.
     * @param[in] weight Physical integration weight supplied by NeoPZ.
     * @param[in,out] ek Element tangent matrix receiving additive contributions.
     * @param[in,out] ef Element negative-residual vector receiving additive contributions.
     */
    void Contribute(const TPZVec<TPZMaterialDataT<STATE>> &datavec, REAL weight, TPZFMatrix<STATE> &ek,
                    TPZFMatrix<STATE> &ef) override;

    /**
     * @brief Add the type-0 pressure load to flux rows. Type-1 prescribed normal
     * flux contributes nothing here because TPZWannAnalysis initializes and filters
     * its equations. The assembly switch has no Robin contribution.
     * @param[in] datavec Flux data at index 0 and pressure data at index 1.
     * @param[in] weight Physical integration weight supplied by NeoPZ.
     * @param[in,out] ek Element tangent matrix receiving additive contributions.
     * @param[in,out] ef Element negative-residual vector receiving additive contributions.
     * @param[in] bc Boundary type and prescribed values.
     */
    void ContributeBC(const TPZVec<TPZMaterialDataT<STATE>> &datavec, REAL weight, TPZFMatrix<STATE> &ek,
                      TPZFMatrix<STATE> &ef, TPZBndCondT<STATE> &bc) override;

    /**
     * @brief Request the current solution in each space for residual evaluation.
     */
    void FillDataRequirements(TPZVec<TPZMaterialDataT<STATE> > &datavec) const override;

    /**
     * @brief Request solutions and normals in each boundary data entry.
     */
    void FillBoundaryConditionDataRequirements(int type, TPZVec<TPZMaterialDataT<STATE> > &datavec) const override;

};
