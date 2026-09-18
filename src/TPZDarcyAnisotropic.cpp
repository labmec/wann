//
// Created by Gustavo Batistela on 5/13/21.
//

#include "TPZDarcyAnisotropic.h"
#include "pzaxestools.h"

TPZDarcyAnisotropic::TPZDarcyAnisotropic() : TPZRegisterClassId(&TPZDarcyAnisotropic::ClassId),
                               TBase(), fDim(-1) {}

TPZDarcyAnisotropic::TPZDarcyAnisotropic(int id, int dim) : TPZRegisterClassId(&TPZDarcyAnisotropic::ClassId),
                                              TBase(id), fDim(dim) {
                                              }

TPZDarcyAnisotropic::TPZDarcyAnisotropic(const TPZDarcyAnisotropic &copy) : TPZMatBase(copy), fDim(copy.fDim)
{
    *this = copy;
}

TPZDarcyAnisotropic& TPZDarcyAnisotropic::operator=(const TPZDarcyAnisotropic &copy) {
    TPZMatBase::operator=(copy);
    fDim = copy.fDim;
    return *this;
}

void TPZDarcyAnisotropic::SetDimension(int dim) {
    if (dim > 3 || dim < 1) DebugStop();
    fDim = dim;
}

void TPZDarcyAnisotropic::Contribute(const TPZMaterialDataT<STATE> &data, STATE weight, TPZFMatrix<STATE> &ek,
                              TPZFMatrix<STATE> &ef) {
    const TPZFMatrix<REAL> &phi = data.phi;
    const TPZFMatrix<REAL> &dphi = data.dphix;
    const TPZVec<REAL> &x = data.x;
    TPZFNMatrix<1, REAL> Aux(1, 1, 1.);

    STATE source_term = 0;
    if (this->HasForcingFunction()) {
        TPZManVector<STATE, 1> res(1);
        fForcingFunction(x, res);
        source_term = -res[0];
    }

    // Stiffness matrix
    TPZFMatrix<STATE> perm(3, 3, 0.);
    GetPermeability(x, perm);
    TPZFMatrix<STATE> perm_dot_dphi;
    perm.Multiply(dphi, perm_dot_dphi);
    ek.AddContribution(0, 0, perm_dot_dphi, 1, dphi, 0, weight);

    // Source term
    ef.AddContribution(0, 0, phi, 0, Aux, 0, -source_term * weight);
}

void TPZDarcyAnisotropic::ContributeBC(const TPZMaterialDataT<STATE> &data, STATE weight, TPZFMatrix<STATE> &ek,
                                TPZFMatrix<STATE> &ef, TPZBndCondT<STATE> &bc) {

    const TPZFMatrix<REAL> &phi = data.phi;
    const TPZFMatrix<REAL> &axes = data.axes;
    int phr = phi.Rows();
    int in, jn;

    STATE v2 = bc.Val2()[0];

    if (bc.HasForcingFunctionBC()) {
        TPZManVector<STATE, 1> rhs_val(1);
        TPZFNMatrix<1, STATE> mat_val(fDim, 1);
        bc.ForcingFunctionBC()(data.x, rhs_val, mat_val);
        // MinusKGradU/Flux;
        TPZFMatrix<STATE> perm(3, 3, 0.);
        GetPermeability(data.x, perm);
        TPZFMatrix<STATE> Flux(fDim, 1, 0.);
        perm.MultAdd(mat_val, Flux, Flux, -1., 0.);
        TPZManVector<REAL,3> normal(3,0.);
        for (int i = 0; i < fDim; i++) {
            normal[i] = data.normal[i];
        }
        if(bc.Type() == 0) {
            v2 = rhs_val[0];
        } else if(bc.Type() == 1) {
            v2 = 0.;
            for (int i = 0; i < fDim; i++) {
                v2 += -Flux[i] * normal[i];
            }
        } else if(bc.Type() == 2) {
            v2 = 0.;
            for (int i = 0; i < fDim; i++) {
                v2 += -Flux[i] * normal[i];
            }
            v2 += bc.Val1()(0,0) * rhs_val[0];
        }
    }

    switch (bc.Type()) {
        case 0 : // Dirichlet condition
            for (in = 0; in < phr; in++) {
                ef(in, 0) += (STATE) (TPZMaterial::fBigNumber * phi(in, 0) * weight) * v2;
                for (jn = 0; jn < phr; jn++) {
                    ek(in, jn) += TPZMaterial::fBigNumber * phi(in, 0) * phi(jn, 0) * weight;
                }
            }
            break;
        case 1 : // Neumann condition
            for (in = 0; in < phi.Rows(); in++) {
                ef(in, 0) += v2 * (STATE) (phi(in, 0) * weight);
            }
            break;
        case 2 : // Robin condition
            for (in = 0; in < phi.Rows(); in++) {
                ef(in, 0) += v2 * (STATE) (phi(in, 0) * weight);
                for (jn = 0; jn < phi.Rows(); jn++) {
                    ek(in, jn) += bc.Val1()(0, 0) * (STATE) (phi(in, 0) * phi(jn, 0) * weight);
                }
            }
            break;
        default:
            PZError << __PRETTY_FUNCTION__
                    << "\nBoundary condition type not implemented. Please use one of the following:\n"
                    << "\t 0: Dirichlet\n"
                    << "\t 1: Neumann\n"
                    << "\t 2: Robin\n";
            DebugStop();
    }
}

int TPZDarcyAnisotropic::VariableIndex(const std::string &name) const {

    if (!strcmp("Solution", name.c_str())) return 1;
    if (!strcmp("Pressure", name.c_str())) return 1;
    if (!strcmp("Derivative", name.c_str())) return 2;
    if (!strcmp("GradU", name.c_str())) return 2;
    if (!strcmp("KDuDx", name.c_str())) return 3;
    if (!strcmp("KDuDy", name.c_str())) return 4;
    if (!strcmp("KDuDz", name.c_str())) return 5;
    if (!strcmp("NormKDu", name.c_str())) return 6;
    if (!strcmp("MinusKGradU", name.c_str())) return 7;
    if (!strcmp("Flux", name.c_str())) return 7;
    if (!strcmp("POrder", name.c_str())) return 8;
    if (!strcmp("ExactPressure", name.c_str())) return 9;
    if (!strcmp("ExactSolution", name.c_str())) return 9;
    if (!strcmp("ExactFlux", name.c_str())) return 10;
    if (!strcmp("Div", name.c_str())) return 11;
    if (!strcmp("Divergence", name.c_str())) return 11;
    if (!strcmp("ExactDiv", name.c_str())) return 12;
    if (!strcmp("ExactDivergence", name.c_str())) return 12;
    if (!strcmp("FluxL2", name.c_str())) return 13;
    if (!strcmp("EstimatedError", name.c_str())) return 100;
    if (!strcmp("TrueError", name.c_str())) return 101;
    if (!strcmp("EffectivityIndex", name.c_str())) return 102;
    if (!strcmp("ResidualError", name.c_str())) return 103;

    return TPZMatBase::VariableIndex(name);
}

int TPZDarcyAnisotropic::NSolutionVariables(int var) const {

    if (var == 1) return 1;      // Solution/Pressure
    if (var == 2) return fDim;   // Derivative/GradU
    if (var == 3) return 1;      // KDuDx;
    if (var == 4) return 1;      // KDuDy;
    if (var == 5) return 1;      // KDuDz;
    if (var == 6) return 1;      // NormKDu;
    if (var == 7) return 3;   // MinusKGradU/Flux;
    if (var == 8) return 1;      // POrder
    if (var == 9) return 1;      // ExactPressure/ExactSolution
    if (var == 10) return fDim;  // ExactFlux
    if (var == 11) return 1;     // Div/Divergence
    if (var == 12) return 1;     // ExactDiv/ExactDivergence
    if (var == 13) return fDim;  // FluxL2
    if (var == 100) return 1;  // EstimatedError
    if (var == 101) return 1;  // TrueError
    if (var == 102) return 1;  // EffectivityIndex
    if (var == 103) return 1;  // ResidualError


    return TPZMatBase::NSolutionVariables(var);
}

void TPZDarcyAnisotropic::Solution(const TPZMaterialDataT<STATE> &data, int var, TPZVec<STATE> &solOut) {

    if(data.fShapeType == TPZMaterialData::EEmpty) {
        solOut.Resize(NSolutionVariables(var));
        solOut.Fill(0.);
        return;
    }
    switch (var) {
        case 1: {
            // Solution/Pressure
            solOut[0] = data.sol[0][0];
            return;
        }
        case 2: {
            // Derivative/GradU
            TPZFNMatrix<9, STATE> dsoldx;
            TPZAxesTools<STATE>::Axes2XYZ(data.dsol[0], dsoldx, data.axes);
            for (int id = 0; id < fDim; id++) {
                solOut[id] = dsoldx(id, 0);
            }
            return;
        }
        case 3: {
            // KDuDx;
            TPZFNMatrix<9, STATE> dsoldx;
            TPZAxesTools<STATE>::Axes2XYZ(data.dsol[0], dsoldx, data.axes);
            TPZFMatrix<STATE> perm(3, 3, 0.);
            GetPermeability(data.x, perm);
            TPZFMatrix<STATE> flux(3, 1, 0.);
            perm.MultAdd(dsoldx, flux, flux);
            solOut[0] = flux(0, 0);
            return;
        }
        case 4: {
            // KDuDy;
            TPZFNMatrix<9, STATE> dsoldx;
            TPZAxesTools<STATE>::Axes2XYZ(data.dsol[0], dsoldx, data.axes);
            TPZFMatrix<STATE> perm(3, 3, 0.);
            GetPermeability(data.x, perm);
            TPZFMatrix<STATE> flux(3, 1, 0.);
            perm.MultAdd(dsoldx, flux, flux);
            solOut[0] = flux(1, 0);
            return;
        }
        case 5: {
            // KDuDz;
            TPZFNMatrix<9, STATE> dsoldx;
            TPZAxesTools<STATE>::Axes2XYZ(data.dsol[0], dsoldx, data.axes);
            TPZFMatrix<STATE> perm(3, 3, 0.);
            GetPermeability(data.x, perm);
            TPZFMatrix<STATE> flux(3, 1, 0.);
            perm.MultAdd(dsoldx, flux, flux);
            solOut[0] = flux(2, 0);
            return;
        }
        case 6: {
            // NormKDu;
            TPZFNMatrix<9, STATE> dsoldx;
            TPZAxesTools<STATE>::Axes2XYZ(data.dsol[0], dsoldx, data.axes);

            TPZFMatrix<STATE> perm(3, 3, 0.);
            GetPermeability(data.x, perm);
            TPZFMatrix<STATE> flux(3, 1, 0.);
            perm.MultAdd(dsoldx, flux, flux);
            STATE res = 0;
            for (int id = 0; id < fDim; id++) {
                res += flux(id, 0) * flux(id, 0);
            }
            solOut[0] = sqrt(res);
            return;
        }
        case 7: {
            // MinusKGradU/Flux;
            TPZFNMatrix<9, STATE> dsoldx;
            TPZAxesTools<STATE>::Axes2XYZ(data.dsol[0], dsoldx, data.axes);
            TPZFMatrix<STATE> perm(3, 3, 0.);
            GetPermeability(data.x, perm);
            TPZFMatrix<STATE> flux(3, 1, 0.);
            perm.MultAdd(dsoldx, flux, flux, -1., 0.);
            for (int id = 0; id < 3; id++) {
                solOut[id] = flux(id, 0);
            }
            return;
        }
        case 8: {
            // POrder
            solOut[0] = data.p;
            return;
        }
        case 9: {
            // ExactPressure/ExactSolution
            TPZVec<STATE> exact_pressure(1);
            TPZFMatrix<STATE> exact_flux(fDim, 1);
            fExactSol(data.x, exact_pressure, exact_flux);
            solOut[0] = exact_pressure[0];
            return;
        }
        case 10: {
            // ExactFlux
            TPZVec<STATE> exact_pressure(1);
            TPZFMatrix<STATE> exact_flux(fDim, 1);
            fExactSol(data.x, exact_pressure, exact_flux);
            for (int id = 0; id < fDim; id++) {
                solOut[id] = exact_flux[id];
            }
            return;
        }
        case 11: {
            // Div/Divergence
            TPZFNMatrix<9, STATE> dsoldx;
            TPZAxesTools<STATE>::Axes2XYZ(data.dsol[0], dsoldx, data.axes);
            TPZFMatrix<STATE> perm(3, 3, 0.);
            GetPermeability(data.x, perm);
            TPZFMatrix<STATE> flux(3, 1, 0.);
            perm.MultAdd(dsoldx, flux, flux, -1., 0.);
            STATE res = 0;
            for (int id = 0; id < fDim; id++) {
                res += flux(id, 0);
            }
            solOut[0] = res;
            return;
        }
        case 12: {
            // ExactDiv/ExactDivergence
            TPZVec<STATE> exact_pressure(1);
            TPZFMatrix<STATE> exact_flux(fDim, 1);
            fExactSol(data.x, exact_pressure, exact_flux);
            STATE res = 0;
            for (int id = 0; id < fDim; id++) {
                res += exact_flux(id, 0);
            }
            solOut[0] = res;

            return;
        }

    default: {
            PZError << __PRETTY_FUNCTION__ << "\n Post-processing variable index not implemented!\n";
            DebugStop();
        }
    }
}

void TPZDarcyAnisotropic::GetSolDimensions(uint64_t &u_len, uint64_t &du_row, uint64_t &du_col) const {
    u_len=1;
    du_row=fDim;
    du_col=1;
}

void TPZDarcyAnisotropic::ErrorNames(TPZVec<std::string> &names) const {
    int nerr = NEvalErrors();
    names.Resize(nerr);
    names[0] = "H1_Norm";
    names[1] = "L2_Norm";
    names[2] = "H1_Semi_Norm";
}

void TPZDarcyAnisotropic::Errors(const TPZMaterialDataT<STATE> &data,
                          TPZVec<REAL> &errors) {
    const TPZVec<REAL> &x = data.x;
    const TPZVec<STATE> &sol = data.sol[0];
    const TPZFMatrix<STATE> &dsol = data.dsol[0];
    const TPZFMatrix<REAL> &axes = data.axes;

#ifdef PZDEBUG
    if(!this->HasExactSol()){
        PZError<<__PRETTY_FUNCTION__;
        PZError<<"\nThe material has no associated exact solution. Aborting...\n";
        DebugStop(); 
    }
#endif
    if(errors.size() != NEvalErrors()) DebugStop();
//    errors.Resize(NEvalErrors(), 0.);

    TPZManVector<STATE,1> exact_pressure(1, 0);
    TPZFNMatrix<3,STATE> exact_flux(fDim, 1, 0);
    fExactSol(x, exact_pressure, exact_flux);

    TPZFNMatrix<3,STATE> gradu(3,1);
    TPZAxesTools<STATE>::Axes2XYZ(dsol,gradu,axes);

    // errors[1] - L2 norm error
    REAL diff = fabs(sol[0] - exact_pressure[0]);
    errors[1] = diff * diff;

    // errors[2] - H1 semi-norm using the permeability tensor
    TPZFMatrix<STATE> perm(3, 3, 0.);
    GetPermeability(data.x, perm);

    TPZManVector<STATE, 3> graduDiff(fDim, 0);
    for (int id = 0; id < fDim; id++) {
        graduDiff[id] = gradu(id) - exact_flux(id, 0);
    }
    diff = 0;
    for (int id = 0; id < fDim; id++) {
        for (int jd = 0; jd < fDim; jd++) {
            diff += graduDiff[id] * perm(id, jd) * graduDiff[jd];
        }
    }
    errors[2] = diff;

    // errors[0] - H1 norm
    errors[0] = errors[1] + errors[2];

    // TODO confirm with Phil is the following norms are correct
    // errors[3] - L2 norm of the x-component of the flux
    // errors[4] - L2 norm of the y-component of the flux, if applicable
    // errors[5] - L2 norm of the z-component of the flux, if applicable
    TPZFNMatrix<9, STATE> dsoldx;
    TPZAxesTools<STATE>::Axes2XYZ(dsol, dsoldx, axes);
    TPZFMatrix<STATE> flux_sol(fDim, 1, 0.);
    perm.MultAdd(dsoldx, flux_sol, flux_sol, -1., 0.);

    for (int id = 0; id < fDim && 3 + id < errors.size(); id++) {
        diff = fabs(exact_flux(id, 0) - flux_sol(id, 0));
        errors[3 + id] = diff * diff;
    }
}

int TPZDarcyAnisotropic::ClassId() const {
    return Hash("TPZDarcyAnisotropic") ^ TBase::ClassId() << 1;
}

TPZMaterial *TPZDarcyAnisotropic::NewMaterial() const {
    return new TPZDarcyAnisotropic(*this);
}

void TPZDarcyAnisotropic::Print(std::ostream &out) const {
    out << "Material Name: " << this->Name() << "\n";
    out << "Material Id: " << this->Id() << "\n";
    out << "Dimension: " << this->Dimension() << "\n\n";
}

void TPZDarcyAnisotropic::FillDataRequirements(TPZMaterialData &data) const {
    data.SetAllRequirements(false);
}


void TPZDarcyAnisotropic::FillBoundaryConditionDataRequirements(int type, TPZMaterialData &data) const {

    data.SetAllRequirements(false);
    if (type == 50) {
        data.fNeedsSol = true;
    }
    if (type == 3 || type == 1) {
        data.fNeedsNormal = true;
    }
    if (HasForcingFunction()) {
        data.fNeedsNormal = true;
    }
}
