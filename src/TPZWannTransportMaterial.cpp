// Adapted from subFlow/src/TSFTransportMaterial.cpp (LabMeC, MIT).
#include "TPZWannTransportMaterial.h"
#include <stdexcept>

int TPZWannTransportMaterial::VariableIndex(const std::string& name) const {
    if (name == "Sw") return 0;
    if (name == "So") return 1;
    return TBase::VariableIndex(name);
}
int TPZWannTransportMaterial::NSolutionVariables(int var) const {
    if (var == 0 || var == 1) return 1;
    return TBase::NSolutionVariables(var);
}
void TPZWannTransportMaterial::Solution(const TPZMaterialDataT<STATE>& data, int var,
                                       TPZVec<STATE>& sol) {
    if (var != 0 && var != 1) throw std::invalid_argument("Unknown transport output");
    sol.Resize(1);
    sol[0] = var == 0 ? data.sol[0][0] : 1. - data.sol[0][0];
}
void TPZWannTransportMaterial::SolutionInterface(const TPZMaterialDataT<STATE>&,
    const TPZMaterialDataT<STATE>&, const TPZMaterialDataT<STATE>&, int, TPZVec<STATE>&) {
    throw std::logic_error("Transport output is defined on reservoir cells only");
}
