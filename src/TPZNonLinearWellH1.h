//
// Created by Gustavo Batistela on 5/13/21.
//

#pragma once

#include "TPZMatBase.h"
#include "TPZMatSingleSpace.h"
#include "TPZMatErrorSingleSpace.h"

/**
 * @ingroup material
 * @brief This class implements an H1-conforming approximation for the Poiseulli flow.
 */

class TPZNonLinearWellH1 : public TPZMatBase<STATE, TPZMatSingleSpaceT<STATE>,
        TPZMatErrorSingleSpace<STATE>> {

    // type alias to improve constructor readability
    using TBase = TPZMatBase<STATE, TPZMatSingleSpaceT<STATE>, TPZMatErrorSingleSpace<STATE>>;

public:
    /**
     * @brief Default constructor
     */
    TPZNonLinearWellH1();
    
    /**
     * @brief Class constructor
     * @param [in] id material id
     * @param [in] Dw Well diameter
     * @param [in] mu Fluid viscosity
     * @param [in] rho Fluid density
     * @param [in] pres Reservoir pressure
     * @param [in] Kvw Pseudo resistivity
     */
    [[maybe_unused]] TPZNonLinearWellH1(int id, REAL Dw, REAL mu, REAL rho, REAL pres, REAL Kvw);


    /**
	 * @brief Class constructor
	 * @param [in] id material id
	 * @param [in] dim problem dimension
	 */
    TPZNonLinearWellH1(int id, int dim);

            TPZNonLinearWellH1(const TPZNonLinearWellH1 &copy);

            TPZNonLinearWellH1& operator=(const TPZNonLinearWellH1 &copy);

    /**
	 * @brief Returns a 'std::string' with the name of the material
	 */
    [[nodiscard]] std::string Name() const override { return "TPZNonLinearWellH1"; }

    /**
	 * @brief Returns the problem dimension
	 */
    [[nodiscard]] int Dimension() const override { return this->fDim; }

    /**
	 * @brief Returns the number of state variables
	 */
    [[nodiscard]] int NStateVariables() const override { return 1; }

    /**
	 * @brief Returns the number of errors to be evaluated
     *
     * Returns the number of errors to be evaluated, that is, the number of error norms associated
     * with the problem.
     */
    int NEvalErrors() const override { return 3; }

    /**
     * @brief Sets problem dimension
     */
    virtual void SetDimension(int dim);

    /**
     * @brief It computes a contribution to the stiffness matrix and load vector at one integration point
     * @param[in] data stores all input data
     * @param[in] weight is the weight of the integration rule
     * @param[out] ek is the element matrix
     * @param[out] ef is the rhs vector
     */
    void Contribute(const TPZMaterialDataT<STATE> &data, STATE weight, TPZFMatrix<STATE> &ek,
                    TPZFMatrix<STATE> &ef) override;

    void ContributeResidual(const TPZMaterialDataT<STATE> &data, REAL weight,
                            TPZFMatrix<STATE> &ef);

    /**
     * @brief It computes a contribution to the stiffness matrix and load vector at one BC integration point
     * @param[in] data stores all input data
     * @param[in] weight is the weight of the integration rule
     * @param[out] ek is the element matrix
     * @param[out] ef is the rhs vector
     * @param[in] bc is the boundary condition material
     */
    void ContributeBC(const TPZMaterialDataT<STATE> &data, STATE weight, TPZFMatrix<STATE> &ek, TPZFMatrix<STATE> &ef,
                      TPZBndCondT<STATE> &bc) override;

    /**
     * @brief Returns an integer associated with a post-processing variable name
     * @param [in] name string containing the name of the post-processing variable. Ex: "Pressure".
     */
    [[nodiscard]] int VariableIndex(const std::string &name) const override;

    /**
     * @brief Returns an integer with the dimension of a post-processing variable
     * @param [in] var index of the post-processing variable, according to TPZDarcyFlow::VariableIndex method.
     */
    [[nodiscard]] int NSolutionVariables(int var) const override;

    /**
     * @brief Returns the solution associated with the var index based on the
     * finite element approximation at a point
     * @param [in] data material data associated with a given integration point
     * @param [in] var index of the variable to be calculated
     * @param [out] solOut vector to store the solution
     */
    void Solution(const TPZMaterialDataT<STATE> &data, int var, TPZVec<STATE> &solOut) override;

    /**
     * @brief Get the dimensions of the solution for each state variable
     *
     * This will be used for initializing the corresponding TPZMaterialData
     * @param [out] u_len solution dimension
     * @param [out] du_row number of rows of the derivative
     * @param [out] du_col number of columns of the derivative
     */
    void GetSolDimensions(uint64_t &u_len, uint64_t &du_row, uint64_t &du_col) const override;

    /**
     * @brief Calculates the approximation error at a point
     * @param [in] data material data of the integration point
     * @param [out] errors calculated errors
     */
    void Errors(const TPZMaterialDataT<STATE> &data, TPZVec<REAL> &errors) override;

    /** @brief Fills in the name of the errors that are computed */
    virtual void ErrorNames(TPZVec<std::string> &names) const override;

    /*
     * @brief fill requirements for volumetric contribute
     */
    void FillDataRequirements(TPZMaterialData &data) const override;

    /*
     * @brief fill requirements for boundary contribute
     */
    void FillBoundaryConditionDataRequirements(int type, TPZMaterialData &data) const override;

    /**
     * @brief Returns an unique class identifier
     */
    [[nodiscard]] int ClassId() const override;

    /**
     * @brief Creates another material of the same type
     */
    [[nodiscard]] TPZMaterial *NewMaterial() const override;

    /**
     * @brief Prints data associated with the material.
     */
    void Print(std::ostream & out) const override;

    static bool fAssembleRHSOnly;

    static bool fIsFirstIteration;

protected:
    /**
     * @brief Problem dimension
     */
    int fDim;

    REAL fDw;

    REAL fMu;

    REAL fRho;

    REAL fPres;

    REAL fKvw;

    REAL fC; // this stores c = (2.252610888 Dw^(19/7))/(mu^(1/7) rho^(3/7))

    REAL fCLin; // this stores (pi Dw^4) / (128 mu)
};

