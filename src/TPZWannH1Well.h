/**
 * @file
 * @brief Class and method contracts for TPZWannH1Well.
 */

#pragma once

#include "TPZMatBase.h"
#include "TPZMatSingleSpace.h"
#include "TPZMatErrorSingleSpace.h"

/**
 * @brief Pressure-only one-dimensional well material with laminar/turbulent friction.
 * Uses pressure derivatives for the axial flow rate, a laminar first Newton
 * iteration, and a transition criterion based on Reynolds number 1187.38.
 * Optional pseudo-resistivity contributes Kvw*(pres-p). Physical assembly needs
 * the diameter/viscosity/density constructor; the dimension-only constructor
 * does not initialize those coefficients. Current flux output uses global x.
 * @see TPZWannMixedWell
 * @see TPZWannAnalysis
 */
class TPZWannH1Well : public TPZMatBase<STATE, TPZMatSingleSpaceT<STATE>,
        TPZMatErrorSingleSpace<STATE>> {

    // type alias to improve constructor readability
    using TBase = TPZMatBase<STATE, TPZMatSingleSpaceT<STATE>, TPZMatErrorSingleSpace<STATE>>;

public:

    /**
     * @brief Construct an object; configure simulation/material data before use.
     */
    TPZWannH1Well();

    /**
     * @brief Configure a 1D well and compute laminar/turbulent coefficients.
     * @param[in] id Material identifier.
     * @param[in] Dw Positive well diameter.
     * @param[in] mu Positive fluid viscosity.
     * @param[in] rho Positive fluid density.
     * @param[in] pres Reference reservoir pressure for the pseudo-resistivity term.
     * @param[in] Kvw Pseudo-resistivity coefficient; zero disables that term.
     */
    [[maybe_unused]] TPZWannH1Well(int id, REAL Dw, REAL mu, REAL rho, REAL pres, REAL Kvw);

    /**
     * @brief Construct a material with the given identifier and dimension.
     */
    TPZWannH1Well(int id, int dim);

            /**
             * @brief Construct using the current copy implementation.
             * See docs/nonlinear-wells.md for well-coefficient copy limitations.
             */
            TPZWannH1Well(const TPZWannH1Well &copy);

            /**
             * @brief Assign using the current base/member copy implementation.
             * See docs/nonlinear-wells.md for well-coefficient copy limitations.
             */
            TPZWannH1Well& operator=(const TPZWannH1Well &copy);

    /**
     * @brief Return the concrete material name used in diagnostics.
     */
    [[nodiscard]] std::string Name() const override { return "TPZWannH1Well"; }

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
     * @brief Add the axial well tangent and negative residual at the current pressure.
     * Uses fCLin in the first iteration and switches to the turbulent power law
     * above the critical pressure gradient. Includes pseudo-resistivity terms.
     * Delegates to ContributeResidual() when fAssembleRHSOnly is set.
     * @param[in] data Pressure shape/solution data at the integration point.
     * @param[in] weight Physical integration weight supplied by NeoPZ.
     * @param[in,out] ek Element tangent matrix receiving additive contributions.
     * @param[in,out] ef Element negative-residual vector receiving additive contributions.
     */
    void Contribute(const TPZMaterialDataT<STATE> &data, STATE weight, TPZFMatrix<STATE> &ek,
                    TPZFMatrix<STATE> &ef) override;

    /**
     * @brief Accumulate only the well residual; choose the regime from the laminar
     * flow estimate and recheck turbulent Reynolds consistency. Does not add a matrix.
     * @param[in] data Pressure shape/solution data at the integration point.
     * @param[in] weight Physical integration weight supplied by NeoPZ.
     * @param[in,out] ef Element negative-residual vector receiving additive contributions.
     */
    void ContributeResidual(const TPZMaterialDataT<STATE> &data, REAL weight,
                            TPZFMatrix<STATE> &ef);

    /**
     * @brief Leave type-0 pressure to analysis filtering; add type-1 flow loads or
     * type-2 Robin terms. Other types trigger DebugStop().
     * @param[in] data Pressure shape/solution data at the integration point.
     * @param[in] weight Physical integration weight supplied by NeoPZ.
     * @param[in,out] ek Element tangent matrix receiving additive contributions.
     * @param[in,out] ef Element negative-residual vector receiving additive contributions.
     * @param[in] bc Boundary type and prescribed values.
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
     * @brief Fill three pointwise squared pressure/gradient errors against an exact
     * callback, using unit gradient weighting rather than the well friction law.
     * Requires errors to be pre-sized to three.
     */
    void Errors(const TPZMaterialDataT<STATE> &data, TPZVec<REAL> &errors) override;

    /**
     * @brief Resize names and label the three H1/L2/gradient error slots.
     */
    virtual void ErrorNames(TPZVec<std::string> &names) const override;

    /**
     * @brief Request pressure and derivative solution data for nonlinear assembly.
     */
    void FillDataRequirements(TPZMaterialData &data) const override;

    /**
     * @brief Clear optional solution data; request normals for boundary types 1/3
     * or when a volume forcing function is installed.
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
     * @brief Turbulent coefficient; maps pressure-gradient magnitude to flow.
     */
    REAL fC; // this stores c = (2.252610888 Dw^(19/7))/(mu^(1/7) rho^(3/7))

    /**
     * @brief Laminar conductance pi*Dw^4/(128*mu).
     */
    REAL fCLin; // this stores (pi Dw^4) / (128 mu)
};

