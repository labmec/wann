//
// Created by Gustavo Batistela on 5/13/21.
//

#include "TPZNonLinearWellH1.h"
#include "pzaxestools.h"

bool TPZNonLinearWellH1::fAssembleRHSOnly = false;
bool TPZNonLinearWellH1::fIsFirstIteration = true;

TPZNonLinearWellH1::TPZNonLinearWellH1() : TPZRegisterClassId(&TPZNonLinearWellH1::ClassId),
                               TBase(), fDim(-1) {}

TPZNonLinearWellH1::TPZNonLinearWellH1(int id, int dim) : TPZRegisterClassId(&TPZNonLinearWellH1::ClassId),
                                              TBase(id), fDim(dim) {
                                              }

[[maybe_unused]] TPZNonLinearWellH1::TPZNonLinearWellH1(int id, REAL Dw, REAL mu, REAL rho, REAL pres, REAL Kvw) : TPZRegisterClassId(&TPZNonLinearWellH1::ClassId),
                                                                                                               TBase(id), fDim(1), fDw(Dw), fMu(mu), fRho(rho), fPres(pres), fKvw(Kvw)
{
    fC = (2.252610888 * pow(fDw, 19. / 7.)) / (pow(fMu, 1. / 7.) * pow(fRho, 3. / 7.));
    fCLin = (M_PI * pow(fDw, 4)) / (128. * fMu);
}

TPZNonLinearWellH1::TPZNonLinearWellH1(const TPZNonLinearWellH1 &copy) : TPZMatBase(copy), fDim(copy.fDim)
{
    *this = copy;
}

TPZNonLinearWellH1& TPZNonLinearWellH1::operator=(const TPZNonLinearWellH1 &copy) {
    TPZMatBase::operator=(copy);
    fDim = copy.fDim;
    return *this;
}

void TPZNonLinearWellH1::SetDimension(int dim) {
    if (dim > 3 || dim < 1) DebugStop();
    fDim = dim;
}

void TPZNonLinearWellH1::Contribute(const TPZMaterialDataT<STATE> &data, STATE weight, TPZFMatrix<STATE> &ek,
                              TPZFMatrix<STATE> &ef) {
    
    if (fAssembleRHSOnly)
    {
        ContributeResidual(data, weight, ef);
        return;
    }

    const TPZFMatrix<STATE> &dphix = data.dphix;
    const TPZFMatrix<STATE> &phi = data.phi;
    TPZFNMatrix<1, REAL> Aux(1, 1, 1.);

    REAL psol = data.sol[0][0];
    REAL dpsol = data.dsol[0][0];

    // Compute the Q and dp/dx in which the flow changes from laminar to turbulent
    REAL Qcrit = (1187.38 * fMu * M_PI * fDw) / (4 * fRho);
    REAL DeltaPcrit = Qcrit/fCLin;

    bool turbulent = false;
    if (!fIsFirstIteration && std::abs(dpsol) > DeltaPcrit) { 
      turbulent = true;
    }

    REAL signal = (dpsol >= 0.) ? 1. : -1.;
    REAL factor = fCLin * weight;
    if (turbulent) {
        factor = fC * (4./7.) * signal * (1.0/ std::pow(std::abs(dpsol), 3./7.)) * weight;
    }

    ek.AddContribution(0, 0, dphix, 1, dphix, 0, factor);
    ek.AddContribution(0, 0, phi, 0, phi, 1, fKvw * weight); // Pseudo-resistivity
    
    factor = -fCLin * weight * dpsol;
    if (turbulent) {
        factor = -fC * weight * std::pow(std::abs(dpsol), 4./7.);
    }

    ef.AddContribution(0, 0, dphix, 1, Aux, 0, factor);
    ef.AddContribution(0, 0, phi, 0, Aux, 0, fKvw * weight * (fPres - psol)); // Pseudo-resistivity
}

void TPZNonLinearWellH1::ContributeBC(const TPZMaterialDataT<STATE> &data, STATE weight, TPZFMatrix<STATE> &ek,
                                TPZFMatrix<STATE> &ef, TPZBndCondT<STATE> &bc) {

    const TPZFMatrix<REAL> &phi = data.phi;
    const TPZFMatrix<REAL> &axes = data.axes;
    int phr = phi.Rows();
    int in, jn;

    STATE v2 = bc.Val2()[0];

    switch (bc.Type()) {
        case 0 : // Dirichlet condition already imposed in the initial solution
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

void TPZNonLinearWellH1::ContributeResidual(const TPZMaterialDataT<STATE> &data, REAL weight, TPZFMatrix<STATE> &ef)
{
    const TPZFMatrix<STATE> &dphix = data.dphix;
    const TPZFMatrix<STATE> &phi = data.phi;
    TPZFNMatrix<1, REAL> Aux(1, 1, 1.);

    REAL psol = data.sol[0][0];
    REAL dpsol = data.dsol[0][0];

    // At this point, how should we compute flow rate? Linear or turbulent?
    // REAL Qsol = - fC * std::pow(std::abs(dpsol), 4./7.); // Non-linear case
    REAL Qsol = - dpsol * fCLin; // Linear case

    bool turbulent = false;
    REAL reynolds = 0.0;

    if (!fIsFirstIteration) { // If it's the first iteration, turbulence = false
      REAL velocity = 4. * Qsol / (M_PI * fDw * fDw);
      reynolds = (fRho * std::abs(velocity) * fDw) / fMu;
      turbulent = reynolds > 1187.38;
    }



    // Check if we are definitely in the turbulent regime
    if (turbulent) {
      Qsol = -fC * std::pow(std::abs(dpsol), 4. / 7.);
      REAL velocity = 4. * Qsol / (M_PI * fDw * fDw);
      reynolds = (fRho * std::abs(velocity) * fDw) / fMu;
      if (reynolds < 1187.38) {
        std::cout << "Warning: Reynolds number oscillating around the critical "
                     "value. Current value: "
                  << reynolds << std::endl;
        DebugStop();
      }
    }
    
    REAL factor = -fCLin * weight * dpsol;
    if (turbulent) {
        factor = -fC * weight * std::pow(std::abs(dpsol), 4./7.);
    }

    ef.AddContribution(0, 0, dphix, 1, Aux, 0, factor);
    ef.AddContribution(0, 0, phi, 0, Aux, 0, fKvw * weight * (fPres - psol)); // Pseudo-resistivity
}

int TPZNonLinearWellH1::VariableIndex(const std::string &name) const {

    if (!strcmp("Solution", name.c_str())) return 1;
    if (!strcmp("Pressure", name.c_str())) return 1;
    if (!strcmp("Derivative", name.c_str())) return 2;
    if (!strcmp("GradPressure", name.c_str())) return 2;
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

int TPZNonLinearWellH1::NSolutionVariables(int var) const {

    if (var == 1) return 1;      // Solution/Pressure
    if (var == 2) return fDim;   // Derivative/GradU
    if (var == 3) return 1;      // KDuDx;
    if (var == 4) return 1;      // KDuDy;
    if (var == 5) return 1;      // KDuDz;
    if (var == 6) return 1;      // NormKDu;
    if (var == 7) return fDim;   // MinusKGradU/Flux;
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

void TPZNonLinearWellH1::Solution(const TPZMaterialDataT<STATE> &data, int var, TPZVec<STATE> &solOut) {

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
            const STATE perm = 1.0;
            solOut[0] = perm * dsoldx(0, 0);
            return;
        }
        case 4: {
            // KDuDy;
            TPZFNMatrix<9, STATE> dsoldx;
            TPZAxesTools<STATE>::Axes2XYZ(data.dsol[0], dsoldx, data.axes);
            const STATE perm = 1.0;
            solOut[0] = perm * dsoldx(1, 0);
            return;
        }
        case 5: {
            // KDuDz;
            TPZFNMatrix<9, STATE> dsoldx;
            TPZAxesTools<STATE>::Axes2XYZ(data.dsol[0], dsoldx, data.axes);
            const STATE perm = 1.0;
            solOut[0] = perm * dsoldx(2, 0);
            return;
        }
        case 6: {
            // NormKDu;
            TPZFNMatrix<9, STATE> dsoldx;
            TPZAxesTools<STATE>::Axes2XYZ(data.dsol[0], dsoldx, data.axes);

            const STATE perm = 1.0;
            STATE res = 0;
            for (int id = 0; id < fDim; id++) {
                res += perm * dsoldx(id, 0) * perm * dsoldx(id, 0);
            }
            solOut[0] = sqrt(res);
            return;
        }
        case 7: {
            // Flux 
            TPZFNMatrix<9, STATE> dsoldx;
            TPZAxesTools<STATE>::Axes2XYZ(data.dsol[0], dsoldx, data.axes);

            REAL signal = (dsoldx(0,0) >= 0.) ? 1. : -1.;

            // Compute the Q and dp/dx in which the flow changes from laminar to turbulent
            REAL Qcrit = (1187.38 * fMu * M_PI * fDw) / (4 * fRho);
            REAL DeltaPcrit = Qcrit/fCLin;

            bool turbulent = false;
            if (!fIsFirstIteration && std::abs(dsoldx(0, 0)) > DeltaPcrit) {
              turbulent = true;
            }

            if (turbulent) {
                solOut[0] = - signal * fC * std::pow(std::abs(dsoldx(0, 0)), 4./7.);
            } else {
                solOut[0] = - fCLin * dsoldx(0, 0);
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
            const STATE perm = 1.0;
            STATE res = 0;
            for (int id = 0; id < fDim; id++) {
                res += - perm * dsoldx(id, 0);
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

void TPZNonLinearWellH1::GetSolDimensions(uint64_t &u_len, uint64_t &du_row, uint64_t &du_col) const {
    u_len=1;
    du_row=fDim;
    du_col=1;
}

void TPZNonLinearWellH1::ErrorNames(TPZVec<std::string> &names) const {
    int nerr = NEvalErrors();
    names.Resize(nerr);
    names[0] = "H1_Norm";
    names[1] = "L2_Norm";
    names[2] = "H1_Semi_Norm";
}

void TPZNonLinearWellH1::Errors(const TPZMaterialDataT<STATE> &data,
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

    // errors[2] - H1 semi-norm: |H1| = K*(grad[u] - grad[u_exact])
    const STATE perm = 1.0;

    TPZManVector<REAL,3> graduDiff(fDim, 0);
    for (int id = 0; id < fDim; id++) {
        graduDiff[id] += fabs(gradu(id) - exact_flux(id, 0));
    }
    diff = 0;
    for (int id = 0; id < fDim; id++) {
        REAL aux = graduDiff[id];
        diff += perm * aux * aux;
    }
    errors[2] = abs(diff);

    // errors[0] - H1 norm
    errors[0] = errors[1] + errors[2];

    // TODO confirm with Phil is the following norms are correct
    // errors[3] - L2 norm of the x-component of the flux
    // errors[4] - L2 norm of the y-component of the flux, if applicable
    // errors[5] - L2 norm of the z-component of the flux, if applicable
    TPZFNMatrix<9, STATE> dsoldx;
    TPZAxesTools<STATE>::Axes2XYZ(dsol, dsoldx, axes);
    TPZManVector<STATE, 3> flux_sol(fDim, 0);
    for (int id = 0; id < fDim; id++) {
        flux_sol[id] = - perm * dsoldx(id, 0);
    }
    return;
    for (int id = 0; id < fDim; id++) {
        diff = fabs(exact_flux[id] - flux_sol[id]);
        errors[3 + id] = diff * diff;
    }
}

int TPZNonLinearWellH1::ClassId() const {
    return Hash("TPZNonLinearWellH1") ^ TBase::ClassId() << 1;
}

TPZMaterial *TPZNonLinearWellH1::NewMaterial() const {
    return new TPZNonLinearWellH1(*this);
}

void TPZNonLinearWellH1::Print(std::ostream &out) const {
    out << "Material Name: " << this->Name() << "\n";
    out << "Material Id: " << this->Id() << "\n";
    out << "Dimension: " << this->Dimension() << "\n\n";
}

void TPZNonLinearWellH1::FillDataRequirements(TPZMaterialData &data) const {
    data.SetAllRequirements(false);
    data.fNeedsSol = true;
}


void TPZNonLinearWellH1::FillBoundaryConditionDataRequirements(int type, TPZMaterialData &data) const {

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
