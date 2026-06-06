//
// Created by Gi Taraschi 27/05/26
//

#ifndef TPZANISOTROPICPERMEABILITY_H
#define TPZANISOTROPICPERMEABILITY_H

#include <functional>
#include "pzreal.h"
#include "pzvec.h"
#include "pzfmatrix.h"

// Alias to improve readability of the permeability function type
using IsotropicFunctionType = std::function<STATE(const TPZVec<REAL> &coord)>;
using AnisotropicFunctionType = std::function<TPZFMatrix<STATE>(const TPZVec<REAL> &coord)>;

// Forward declaration of dummy BC interface class
class TPZAnisotropicPermeabilityBC;

/**
 * @brief  This class implements the interface with the methods required to
 * handle the permeability field of an anisotropic material.
 */
class TPZAnisotropicPermeability : public virtual TPZSavable {

public:
    using TInterfaceBC = TPZAnisotropicPermeabilityBC;

    TPZAnisotropicPermeability() : fConstantPermeability(1.), fConstantInversePermeability(1.), fIsotropicPermeabilityFunction(NULL), fAnisotropicPermeabilityFunction(NULL){
    }
    
    TPZAnisotropicPermeability(const TPZAnisotropicPermeability &copy) = default;
    TPZAnisotropicPermeability &operator=(const TPZAnisotropicPermeability &copy) = default;
    /**
     * @brief Set a constant isotropic permeability to the material
     * @param [in] constant permeability value
     */
    void SetConstantPermeability(STATE constant);

        /**
     * @brief Set a constant isotropic permeability to the material
     * @param [in] constant permeability value
     */
    void SetConstantPermeability(TPZFMatrix<STATE> constant);

    /**
     * @brief Set a varying anisotropic permeability field to the material
     * @param [in] perm_function function that describes the permeability field
     */
    void SetPermeabilityFunction(AnisotropicFunctionType &perm_function);

    /**
     * @brief Set a varying isotropic permeability field to the material
     * @param [in] perm_function function that describes the permeability field
     */
    void SetPermeabilityFunction(IsotropicFunctionType &perm_function);

    /**
     * @brief Return the permeability value at a coordinate. In the anisotropic case, returns an average value
     * @param [in] coord coordinate of interest
     */
    STATE GetPermeability(const TPZVec<REAL> &coord);

    /**
     * @brief Return the permeability value (or its average value) at a coordinate
     * @param [in] coord coordinate of interest
     * @param [out] perm_matrix permeability matrix at the given coordinate
     */
    void GetPermeability(const TPZVec<REAL> &coord, TPZFMatrix<STATE> &perm_matrix);

    void GetInversePermeability(const TPZVec<REAL> &coord, TPZFMatrix<STATE> &inv_perm_matrix);

    [[nodiscard]] int ClassId() const override;

    void Read(TPZStream &buf, void *context) override {};

    void Write(TPZStream &buf, int withclassid) const override {};

private:

    // Member variable to describe a constant permeability field
    STATE fConstantPermeabilityScalar = 1.;
    TPZFMatrix<STATE> fConstantPermeability;
    TPZFMatrix<STATE> fConstantInversePermeability;

    bool fHomogeneous = true;
    bool fIsotropic = true;

    // Member variables to describe a varying permeability field
    IsotropicFunctionType fIsotropicPermeabilityFunction = NULL;
    AnisotropicFunctionType fAnisotropicPermeabilityFunction = NULL;
};

// Dummy BC interface class
class TPZMaterial;
class TPZAnisotropicPermeabilityBC : public TPZAnisotropicPermeability {
protected:
    // this method is your chance to verify if the material to which this
    // BC interface applies is compatible with this boundary interface
    // it is called in the method SetMaterial of class TPZBndCondBase
    static void SetMaterialImpl(TPZMaterial *) {}
};

#endif //TPZANISOTROPICPERMEABILITY_H
