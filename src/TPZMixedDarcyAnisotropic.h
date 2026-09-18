//
// Created by Gustavo Batistela on 5/13/21.
//

#ifndef TPZMIXEDDARCYANISOTROPIC_H
#define TPZMIXEDDARCYANISOTROPIC_H

#include "TPZMatBase.h"
#include "TPZMatCombinedSpaces.h"
#include "TPZMatErrorCombinedSpaces.h"
#include "TPZAnisotropicPermeability.h"

/**
 * @ingroup material
 * @brief Mixed Darcy flow material with tensor-valued permeability.
 */

class TPZMixedDarcyAnisotropic : public TPZMatBase<STATE, TPZMatCombinedSpacesT<STATE>,
        TPZMatErrorCombinedSpaces<STATE>, TPZAnisotropicPermeability> {

    using TBase = TPZMatBase<STATE, TPZMatCombinedSpacesT<STATE>,
            TPZMatErrorCombinedSpaces<STATE>, TPZAnisotropicPermeability>;

public:
    TPZMixedDarcyAnisotropic();
    [[maybe_unused]] TPZMixedDarcyAnisotropic(int id, int dim);
    TPZMixedDarcyAnisotropic(const TPZMixedDarcyAnisotropic &copy);
    TPZMixedDarcyAnisotropic &operator=(const TPZMixedDarcyAnisotropic &copy);

    [[nodiscard]] std::string Name() const override { return "TPZMixedDarcyAnisotropic"; }
    [[nodiscard]] int Dimension() const override { return this->fDim; }
    [[nodiscard]] int NStateVariables() const override { return 1; }
    int NEvalErrors() const override { return 5; }

    virtual void SetDimension(int dim);

    void Contribute(const TPZVec<TPZMaterialDataT<STATE>> &datavec, REAL weight, TPZFMatrix<STATE> &ek,
                    TPZFMatrix<STATE> &ef) override;

    void ContributeBC(const TPZVec<TPZMaterialDataT<STATE>> &datavec, REAL weight, TPZFMatrix<STATE> &ek,
                      TPZFMatrix<STATE> &ef, TPZBndCondT<STATE> &bc) override;

    [[nodiscard]] int VariableIndex(const std::string &name) const override;
    [[nodiscard]] int NSolutionVariables(int var) const override;
    void Solution(const TPZVec<TPZMaterialDataT<STATE>> &datavec, int var, TPZVec<STATE> &solOut) override;
    void Errors(const TPZVec<TPZMaterialDataT<STATE>> &data, TPZVec<REAL> &errors) override;
    void FillDataRequirements(TPZVec<TPZMaterialDataT<STATE> > &datavec) const override;
    void FillBoundaryConditionDataRequirements(int type, TPZVec<TPZMaterialDataT<STATE> > &datavec) const override;
    [[nodiscard]] int ClassId() const override;
    [[nodiscard]] TPZMaterial *NewMaterial() const override;
    void Print(std::ostream & out) const override;

protected:
    int fDim;
};

#endif //TPZMIXEDDARCYANISOTROPIC_H