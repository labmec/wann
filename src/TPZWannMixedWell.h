/**
 * @file
 * @brief Class and method contracts for TPZWannMixedWell.
 */

//
// Created by Giovane and Giovanni on 17/09/25
//

#pragma once

#include "TPZMatBase.h"
#include "TPZMatCombinedSpaces.h"
#include "TPZMatErrorCombinedSpaces.h"

/**
 * @brief Mixed one-dimensional well material for axial flow rate and pressure.
 * Uses datavec[0] for flow and datavec[1] for pressure. The constitutive law is
 * laminar resistance times Q or turbulent resistance times Q*abs(Q)^(3/4).
 * Current assembly extracts only the x component of flux vectors. The separate
 * residual-only path forces laminar flow even when full assembly is turbulent.
 * Errors() is not implemented.
 * @see TPZWannH1Well
 * @see TPZWannAnalysis
 */
class TPZWannMixedWell : public TPZMatBase<STATE, TPZMatCombinedSpacesT<STATE>,
                                           TPZMatErrorCombinedSpaces<STATE>>
{

    // type alias to improve constructor readability
    using TBase = TPZMatBase<STATE, TPZMatCombinedSpacesT<STATE>,
                             TPZMatErrorCombinedSpaces<STATE>>;

public:

    /**
     * @brief Construct an object; configure simulation/material data before use.
     */
    TPZWannMixedWell();

    /**
     * @brief Configure a 1D well and compute laminar/turbulent coefficients.
     * @param[in] id Material identifier.
     * @param[in] Dw Positive well diameter.
     * @param[in] mu Positive fluid viscosity.
     * @param[in] rho Positive fluid density.
     * @param[in] pres Reference reservoir pressure for the pseudo-resistivity term.
     * @param[in] Kvw Pseudo-resistivity coefficient; zero disables that term.
     */
    [[maybe_unused]] TPZWannMixedWell(int id, REAL Dw, REAL mu, REAL rho, REAL pres, REAL Kvw);

    /**
     * @brief Construct using the current copy implementation.
     * See docs/nonlinear-wells.md for well-coefficient copy limitations.
     */
    TPZWannMixedWell(const TPZWannMixedWell &copy);

    /**
     * @brief Assign using the current base/member copy implementation.
     * See docs/nonlinear-wells.md for well-coefficient copy limitations.
     */
    TPZWannMixedWell &operator=(const TPZWannMixedWell &copy);

    /**
     * @brief Return the concrete material name used in diagnostics.
     */
    [[nodiscard]] std::string Name() const override { return "TPZWannMixedWell"; }

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
     * @brief Add the axial flow/pressure tangent and negative residual, flux rows first.
     * Laminar first iteration; later regime selection uses velocity-based Reynolds
     * number. Includes pseudo-resistivity and an optional dual functional.
     * Delegates to ContributeResidual() when fAssembleRHSOnly is set.
     * @param[in] datavec Flux data at index 0 and pressure data at index 1.
     * @param[in] weight Physical integration weight supplied by NeoPZ.
     * @param[in,out] ek Element tangent matrix receiving additive contributions.
     * @param[in,out] ef Element negative-residual vector receiving additive contributions.
     */
    void Contribute(const TPZVec<TPZMaterialDataT<STATE>> &datavec, REAL weight, TPZFMatrix<STATE> &ek,
                    TPZFMatrix<STATE> &ef) override;

    /**
     * @brief Add the flow/pressure negative residual without a matrix.
     * @note This implementation explicitly forces the laminar branch; it is not
     * the turbulent residual used by the full Contribute() path.
     * @param[in] datavec Flux data at index 0 and pressure data at index 1.
     * @param[in] weight Physical integration weight supplied by NeoPZ.
     * @param[in,out] ef Element negative-residual vector receiving additive contributions.
     */
    void ContributeResidual(const TPZVec<TPZMaterialDataT<STATE>> &datavec, REAL weight,
                            TPZFMatrix<STATE> &ef);

    /**
     * @brief Add type-0 pressure endpoint loads to flux rows; type-1 flow is
     * initialized and filtered by TPZWannAnalysis, so no penalty is assembled.
     * @param[in] datavec Flux data at index 0 and pressure data at index 1.
     * @param[in] weight Physical integration weight supplied by NeoPZ.
     * @param[in,out] ek Element tangent matrix receiving additive contributions.
     * @param[in,out] ef Element negative-residual vector receiving additive contributions.
     * @param[in] bc Boundary type and prescribed values.
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
     * @brief Unimplemented error evaluation hook; calling it triggers DebugStop().
     */
    void Errors(const TPZVec<TPZMaterialDataT<STATE>> &data, TPZVec<REAL> &errors) override;

    /**
     * @brief Request current solution data in every approximation space.
     */
    void FillDataRequirements(TPZVec<TPZMaterialDataT<STATE>> &datavec) const override;

    /**
     * @brief Clear optional boundary requirements and request the flux-space normal.
     */
    void FillBoundaryConditionDataRequirements(int type, TPZVec<TPZMaterialDataT<STATE>> &datavec) const override;

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
    void Print(std::ostream &out) const override;

    /**
     * @brief Store the dual-functional callback, evaluated against div(flux test).
     * The order argument is currently unused.
     */
    void SetForcingFunctionDual(ForcingFunctionType<REAL> func, int order) {
        fForcingFunctionDual = func;
    }

    /**
     * @brief Shared switch selecting residual-only assembly for every instance.
     */
    static bool fAssembleRHSOnly;
    /**
     * @brief Shared switch forcing laminar assembly in the first Newton iteration.
     */
    static bool fIsFirstIteration;

protected:

    /// Topological dimension; well physical constructors set this to one.
    int fDim;

    /**
     * @brief Well diameter.
     */
    REAL fDw;

    /**
     * @brief Fluid viscosity (in the pressure-scaled units supplied by the builder).
     */
    REAL fMu;

    /**
     * @brief Fluid density (in the pressure-scaled units supplied by the builder).
     */
    REAL fRho;

    /**
     * @brief Reference pressure in the pseudo-resistivity source.
     */
    REAL fPres;

    /**
     * @brief Pseudo-resistivity coupling coefficient.
     */
    REAL fKvw;

    /**
     * @brief Turbulent coefficient; maps Q*abs(Q)^(3/4) to pressure gradient.
     */
    REAL fC; // this stores c^{-7/4}, where c is computed as (2.252610888 Dw^(19/7))/(mu^(1/7) rho^(3/7))

    /**
     * @brief Laminar resistance 128*mu/(pi*Dw^4).
     */
    REAL fCLin; // this stores 128 mu /(pi Dw^4)

    ForcingFunctionType<REAL> fForcingFunctionDual{nullptr};
};
