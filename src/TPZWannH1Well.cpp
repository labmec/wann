/**
 * @file
 * @brief Pressure-based laminar/turbulent well tangent and residual paths.
 * See the matching header for function contracts and docs/ for the architecture.
 */

#include "TPZWannH1Well.h"
#include "pzaxestools.h"

bool TPZWannH1Well::fAssembleRHSOnly = false;
bool TPZWannH1Well::fIsFirstIteration = true;

TPZWannH1Well::TPZWannH1Well() : TPZRegisterClassId(&TPZWannH1Well::ClassId),
                               TBase(), fDim(-1) {}

TPZWannH1Well::TPZWannH1Well(int id, int dim) : TPZRegisterClassId(&TPZWannH1Well::ClassId),
                                              TBase(id), fDim(dim) {
                                              }

[[maybe_unused]] TPZWannH1Well::TPZWannH1Well(int id, REAL Dw, REAL mu, REAL rho, REAL pres, REAL Kvw) : TPZRegisterClassId(&TPZWannH1Well::ClassId),
                                                                                                               TBase(id), fDim(1), fDw(Dw), fMu(mu), fRho(rho), fPres(pres), fKvw(Kvw)
{
    fC = (2.252610888 * pow(fDw, 19. / 7.)) / (pow(fMu, 1. / 7.) * pow(fRho, 3. / 7.));
    fCLin = (M_PI * pow(fDw, 4)) / (128. * fMu);
}

TPZWannH1Well::TPZWannH1Well(const TPZWannH1Well &copy) : TPZMatBase(copy), fDim(copy.fDim)
{
    *this = copy;
}

TPZWannH1Well& TPZWannH1Well::operator=(const TPZWannH1Well &copy) {
    TPZMatBase::operator=(copy);
    fDim = copy.fDim;
    return *this;
}

void TPZWannH1Well::SetDimension(int dim) {
    if (dim > 3 || dim < 1) DebugStop();
    fDim = dim;
}

void TPZWannH1Well::Contribute(const TPZMaterialDataT<STATE> &data, STATE weight, TPZFMatrix<STATE> &ek,
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
        factor = fC * (4./7.) * (1.0/ std::pow(std::abs(dpsol), 3./7.)) * weight;
    }

    ek.AddContribution(0, 0, dphix, 1, dphix, 0, factor);
    ek.AddContribution(0, 0, phi, 0, phi, 1, fKvw * weight); // Pseudo-resistivity
    
    factor = -fCLin * weight * dpsol;
    if (turbulent) {
        factor = -fC * weight * signal * std::pow(std::abs(dpsol), 4./7.);
    }

    ef.AddContribution(0, 0, dphix, 1, Aux, 0, factor);
    ef.AddContribution(0, 0, phi, 0, Aux, 0, fKvw * weight * (fPres - psol)); // Pseudo-resistivity
}

void TPZWannH1Well::ContributeBC(const TPZMaterialDataT<STATE> &data, STATE weight, TPZFMatrix<STATE> &ek,
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

void TPZWannH1Well::ContributeResidual(const TPZMaterialDataT<STATE> &data, REAL weight, TPZFMatrix<STATE> &ef)
{
    const TPZFMatrix<STATE> &dphix = data.dphix;
    const TPZFMatrix<STATE> &phi = data.phi;
    TPZFNMatrix<1, REAL> Aux(1, 1, 1.);

    REAL psol = data.sol[0][0];
    REAL dpsol = data.dsol[0][0];

    const REAL signal = (dpsol >= 0.) ? 1. : -1.;
    REAL Qsol = -dpsol * fCLin; // Laminar estimate for regime selection

    bool turbulent = false;
    REAL reynolds = 0.0;

    if (!fIsFirstIteration) { // If it's the first iteration, turbulence = false
      REAL velocity = 4. * Qsol / (M_PI * fDw * fDw);
      reynolds = (fRho * std::abs(velocity) * fDw) / fMu;
      turbulent = reynolds > 1187.38;
    }



    // Check if we are definitely in the turbulent regime
    if (turbulent) {
      Qsol = -signal * fC * std::pow(std::abs(dpsol), 4. / 7.);
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
        factor = -fC * weight * signal * std::pow(std::abs(dpsol), 4./7.);
    }

    ef.AddContribution(0, 0, dphix, 1, Aux, 0, factor);
    ef.AddContribution(0, 0, phi, 0, Aux, 0, fKvw * weight * (fPres - psol)); // Pseudo-resistivity
}

int TPZWannH1Well::VariableIndex(const std::string &name) const {
    if (name == "Pressure") return 1;
    if (name == "Flux") return 2;
    PZError << __PRETTY_FUNCTION__ << "\nUnsupported field: " << name << "\n";
    DebugStop();
    return -1;
}

int TPZWannH1Well::NSolutionVariables(int var) const {
    if (var == 1 || var == 2) return 1;
    DebugStop();
    return -1;
}

void TPZWannH1Well::Solution(const TPZMaterialDataT<STATE> &data, int var, TPZVec<STATE> &solOut) {
    solOut.Resize(NSolutionVariables(var));
    solOut.Fill(0.);
    if (data.fShapeType == TPZMaterialData::EEmpty) return;

    switch (var) {
        case 1:
            solOut[0] = data.sol[0][0];
            return;
        case 2: {
            // Signed axial flux in the local element direction.
            const REAL dpsol = data.dsol[0](0, 0);

            REAL signal = (dpsol >= 0.) ? 1. : -1.;

            // Compute the Q and dp/ds in which the flow changes from laminar to turbulent
            REAL Qcrit = (1187.38 * fMu * M_PI * fDw) / (4 * fRho);
            REAL DeltaPcrit = Qcrit/fCLin;

            bool turbulent = false;
            if (!fIsFirstIteration && std::abs(dpsol) > DeltaPcrit) {
              turbulent = true;
            }

            if (turbulent) {
                solOut[0] = - signal * fC * std::pow(std::abs(dpsol), 4./7.);
            } else {
                solOut[0] = - fCLin * dpsol;
            }

            return;
        }
        default:
            DebugStop();
    }
}

void TPZWannH1Well::GetSolDimensions(uint64_t &u_len, uint64_t &du_row, uint64_t &du_col) const {
    u_len=1;
    du_row=fDim;
    du_col=1;
}

void TPZWannH1Well::ErrorNames(TPZVec<std::string> &names) const {
    int nerr = NEvalErrors();
    names.Resize(nerr);
    names[0] = "H1_Norm";
    names[1] = "L2_Norm";
    names[2] = "H1_Semi_Norm";
}

void TPZWannH1Well::Errors(const TPZMaterialDataT<STATE> &data,
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

int TPZWannH1Well::ClassId() const {
    return Hash("TPZWannH1Well") ^ TBase::ClassId() << 1;
}

TPZMaterial *TPZWannH1Well::NewMaterial() const {
    return new TPZWannH1Well(*this);
}

void TPZWannH1Well::Print(std::ostream &out) const {
    out << "Material Name: " << this->Name() << "\n";
    out << "Material Id: " << this->Id() << "\n";
    out << "Dimension: " << this->Dimension() << "\n\n";
}

void TPZWannH1Well::FillDataRequirements(TPZMaterialData &data) const {
    data.SetAllRequirements(false);
    data.fNeedsSol = true;
}


void TPZWannH1Well::FillBoundaryConditionDataRequirements(int type, TPZMaterialData &data) const {

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
