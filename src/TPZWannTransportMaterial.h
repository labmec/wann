#pragma once

#include "TPZMatBase.h"
#include "TPZMatSingleSpace.h"
#include "TPZMatInterfaceSingleSpace.h"

/** Mesh/interface and output subset of subFlow::TSFTransportMaterial.
 * Physics is assembled by TPZWannAlgebraicTransport, not these callbacks.
 */
class TPZWannTransportMaterial : public TPZMatBase<STATE,
    TPZMatSingleSpaceT<STATE>, TPZMatInterfaceSingleSpace<STATE>> {
    using TBase = TPZMatBase<STATE, TPZMatSingleSpaceT<STATE>, TPZMatInterfaceSingleSpace<STATE>>;
    int fDimension;
public:
    TPZWannTransportMaterial(int id, int dim) : TBase(id), fDimension(dim) {}
    int Dimension() const override { return fDimension; }
    int NStateVariables() const override { return 1; }
    void GetSolDimensions(uint64_t& length, uint64_t& rows, uint64_t& columns) const override {
        length = 1; rows = fDimension; columns = 1;
    }
    std::string Name() const override { return "TPZWannTransportMaterial"; }
    int ClassId() const override { return Hash("TPZWannTransportMaterial") ^ (TBase::ClassId() << 1); }
    TPZMaterial* NewMaterial() const override { return new TPZWannTransportMaterial(*this); }
    void FillDataRequirements(TPZMaterialData& data) const override { data.SetAllRequirements(false); data.fNeedsSol = true; }
    void FillDataRequirementsInterface(TPZMaterialData& data) const override { data.SetAllRequirements(false); }
    void Contribute(const TPZMaterialDataT<STATE>&, REAL, TPZFMatrix<STATE>&, TPZFMatrix<STATE>&) override {}
    void ContributeBC(const TPZMaterialDataT<STATE>&, REAL, TPZFMatrix<STATE>&, TPZFMatrix<STATE>&, TPZBndCondT<STATE>&) override {}
    void ContributeInterface(const TPZMaterialDataT<STATE>&, const TPZMaterialDataT<STATE>&,
        const TPZMaterialDataT<STATE>&, REAL, TPZFMatrix<STATE>&, TPZFMatrix<STATE>&) override {}
    void ContributeBCInterface(const TPZMaterialDataT<STATE>&, const TPZMaterialDataT<STATE>&,
        REAL, TPZFMatrix<STATE>&, TPZFMatrix<STATE>&, TPZBndCondT<STATE>&) override {}
    void SolutionInterface(const TPZMaterialDataT<STATE>&, const TPZMaterialDataT<STATE>&,
        const TPZMaterialDataT<STATE>&, int, TPZVec<STATE>&) override;
    int VariableIndex(const std::string& name) const override;
    int NSolutionVariables(int var) const override;
    void Solution(const TPZMaterialDataT<STATE>& data, int var, TPZVec<STATE>& sol) override;
};
