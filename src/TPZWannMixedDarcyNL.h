#ifndef TPZWannMixedDarcyNL_H
#define TPZWannMixedDarcyNL_H

#include "Material/DarcyFlow/TPZMixedDarcyFlow.h"
#include "TPZAnisotropicPermeability.h"

class TPZWannMixedDarcyNL : public TPZMixedDarcyFlow, public TPZAnisotropicPermeability {

    using TBase = TPZMixedDarcyFlow;

public:
    /**
     * @brief Default constructor
     */
    TPZWannMixedDarcyNL();

    /**
	 * @brief Class constructor
	 * @param [in] id material id
	 * @param [in] dim problem dimension
	 */
    [[maybe_unused]] TPZWannMixedDarcyNL(int id, int dim);

    /**
             copy constructor
     */
    TPZWannMixedDarcyNL(const TPZWannMixedDarcyNL &copy);
    /**
             copy constructor
     */
    TPZWannMixedDarcyNL &operator=(const TPZWannMixedDarcyNL &copy);
    /**
	 * @brief Returns a 'std::string' with the name of the material
	 */
    [[nodiscard]] std::string Name() const override { return "TPZWannMixedDarcyNL"; }

    /**
     * @brief Returns a unique class identifier
     */
    [[nodiscard]] int ClassId() const override { return Hash("TPZWannMixedDarcyNL") ^ (TBase::ClassId() << 1); }

    void Read(TPZStream &buf, void *context) override { TBase::Read(buf, context); }

    void Write(TPZStream &buf, int withclassid) const override { TBase::Write(buf, withclassid); }

    void SetConstantPermeability(STATE constant)
    {
        TPZAnisotropicPermeability::SetConstantPermeability(constant);
    }

    void SetConstantPermeability(TPZFMatrix<STATE> constant)
    {
        TPZAnisotropicPermeability::SetConstantPermeability(constant);
    }

    void SetPermeabilityFunction(AnisotropicFunctionType &perm_function)
    {
        TPZAnisotropicPermeability::SetPermeabilityFunction(perm_function);
    }

    void SetPermeabilityFunction(IsotropicFunctionType &perm_function)
    {
        TPZAnisotropicPermeability::SetPermeabilityFunction(perm_function);
    }

    STATE GetPermeability(const TPZVec<REAL> &coord)
    {
        return TPZAnisotropicPermeability::GetPermeability(coord);
    }

    void GetPermeability(const TPZVec<REAL> &coord, TPZFMatrix<STATE> &perm_matrix)
    {
        TPZAnisotropicPermeability::GetPermeability(coord, perm_matrix);
    }

    /**
     * @brief It computes a contribution to the stiffness matrix and load vector at one integration point
     * @param[in] datavec stores all input data
     * @param[in] weight is the weight of the integration rule
     * @param[out] ek is the element matrix
     * @param[out] ef is the rhs vector
     */
    void Contribute(const TPZVec<TPZMaterialDataT<STATE>> &datavec, REAL weight, TPZFMatrix<STATE> &ek,
                    TPZFMatrix<STATE> &ef) override;

    /**
     * @brief It computes a contribution to the stiffness matrix and load vector at one BC integration point
     * @param[in] datavec stores all input data
     * @param[in] weight is the weight of the integration rule
     * @param[out] ek is the element matrix
     * @param[out] ef is the rhs vector
     * @param[in] bc is the boundary condition material
     */
    void ContributeBC(const TPZVec<TPZMaterialDataT<STATE>> &datavec, REAL weight, TPZFMatrix<STATE> &ek,
                      TPZFMatrix<STATE> &ef, TPZBndCondT<STATE> &bc) override;

    /*
     * @brief Fill requirements for volumetric contribute
     */
    void FillDataRequirements(TPZVec<TPZMaterialDataT<STATE> > &datavec) const override;

    /*
     * @brief Fill requirements for boundary contribute
     */
    void FillBoundaryConditionDataRequirements(int type, TPZVec<TPZMaterialDataT<STATE> > &datavec) const override;

};

#endif //TPZWannMixedDarcyNL_H
