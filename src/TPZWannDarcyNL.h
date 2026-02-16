#ifndef TPZWannDarcyNL_H
#define TPZWannDarcyNL_H

#include "Material/DarcyFlow/TPZDarcyFlow.h"

class TPZWannDarcyNL : public TPZDarcyFlow {

    using TBase = TPZDarcyFlow;

public:
    /**
     * @brief Default constructor
     */
    TPZWannDarcyNL();

    /**
	 * @brief Class constructor
	 * @param [in] id material id
	 * @param [in] dim problem dimension
	 */
    [[maybe_unused]] TPZWannDarcyNL(int id, int dim);

    /**
             copy constructor
     */
    TPZWannDarcyNL(const TPZWannDarcyNL &copy);
    /**
             copy constructor
     */
    TPZWannDarcyNL &operator=(const TPZWannDarcyNL &copy);
    /**
	 * @brief Returns a 'std::string' with the name of the material
	 */
    [[nodiscard]] std::string Name() const override { return "TPZWannDarcyNL"; }

    /**
     * @brief It computes a contribution to the stiffness matrix and load vector at one integration point
     * @param[in] datavec stores all input data
     * @param[in] weight is the weight of the integration rule
     * @param[out] ek is the element matrix
     * @param[out] ef is the rhs vector
     */
    void Contribute(const TPZMaterialDataT<STATE> &data, REAL weight, TPZFMatrix<STATE> &ek,
                    TPZFMatrix<STATE> &ef) override;

    /**
     * @brief It computes a contribution to the stiffness matrix and load vector at one BC integration point
     * @param[in] datavec stores all input data
     * @param[in] weight is the weight of the integration rule
     * @param[out] ek is the element matrix
     * @param[out] ef is the rhs vector
     * @param[in] bc is the boundary condition material
     */
    void ContributeBC(const TPZMaterialDataT<STATE> &data, REAL weight, TPZFMatrix<STATE> &ek,
                      TPZFMatrix<STATE> &ef, TPZBndCondT<STATE> &bc) override;

    /*
     * @brief Fill requirements for volumetric contribute
     */
    void FillDataRequirements(TPZMaterialData &data) const override;

    /**
     * @brief Fill requirements for boundary contribute
     */
    void FillBoundaryConditionDataRequirements(int type, TPZMaterialData &data) const override;

    /**
     * @brief Get the dimensions of the solution for each state variable
     *
     * This will be used for initializing the corresponding TPZMaterialData
     * @param [out] u_len solution dimension
     * @param [out] du_row number of rows of the derivative
     * @param [out] du_col number of columns of the derivative
     */
    void GetSolDimensions(uint64_t &u_len, uint64_t &du_row, uint64_t &du_col) const override;

};

#endif //TPZWannDarcyNL_H
