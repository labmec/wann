/**
 * @file
 * @brief Flow/pressure well tangent, residual paths, and scalar outputs.
 * See the matching header for function contracts and docs/ for the architecture.
 */

//
// Created by Giovane and Giovanni on 17/09/25
//

#include "TPZWannMixedWell.h"
#include "TPZMaterialDataT.h"
#include "pzlog.h"

#ifdef PZ_LOG
static TPZLogger logger("pz.material.darcy");
#endif

#define USEBLAS

bool TPZWannMixedWell::fAssembleRHSOnly = false;
bool TPZWannMixedWell::fIsFirstIteration = true;

TPZWannMixedWell::TPZWannMixedWell() : TPZRegisterClassId(&TPZWannMixedWell::ClassId),
                                       TBase(), fDim(-1) {}

[[maybe_unused]] TPZWannMixedWell::TPZWannMixedWell(int id, REAL Dw, REAL mu, REAL rho, REAL pres, REAL Kvw) : TPZRegisterClassId(&TPZWannMixedWell::ClassId),
                                                                                                               TBase(id), fDim(1), fDw(Dw), fMu(mu), fRho(rho), fPres(pres), fKvw(Kvw)
{
    REAL c = (2.252610888 * pow(fDw, 19. / 7.)) / (pow(fMu, 1. / 7.) * pow(fRho, 3. / 7.));
    fC = pow(c, -7. / 4.);
    fCLin = 128. * fMu / (M_PI * pow(fDw, 4));
}

/**
         copy constructor
 */
TPZWannMixedWell::TPZWannMixedWell(const TPZWannMixedWell &copy) : TBase(copy), fDim(copy.fDim)
{
    *this = copy;
}
/**
         copy constructor
 */
TPZWannMixedWell &TPZWannMixedWell::operator=(const TPZWannMixedWell &copy)
{
    TBase::operator=(copy);
    fDim = copy.fDim;
    fDw = copy.fDw;
    fMu = copy.fMu;
    fRho = copy.fRho;
    fPres = copy.fPres;
    fKvw = copy.fKvw;
    fC = copy.fC;
    return *this;
}

void TPZWannMixedWell::Contribute(const TPZVec<TPZMaterialDataT<STATE>> &datavec, REAL weight, TPZFMatrix<STATE> &ek,
                                  TPZFMatrix<STATE> &ef)
{
    if (fAssembleRHSOnly)
    {
        ContributeResidual(datavec, weight, ef);
        return;
    }

    TPZFMatrix<REAL> &HDivphiQ = datavec[0].fDeformedDirections;
    TPZFMatrix<REAL> &phiP = datavec[1].phi;
    TPZFMatrix<REAL> &divQ = datavec[0].divphi;
    TPZFNMatrix<1, REAL> Aux(1, 1, 1.);
    const auto &axes = datavec[0].axes;
    TPZFNMatrix<3, REAL> Aux3(HDivphiQ.Rows(), 1, 0.);
    for (int d = 0; d < HDivphiQ.Rows(); d++) {
        Aux3(d, 0) = axes(0, d);
    }

    int nphiQ, nphiP;
    nphiP = phiP.Rows();
    nphiQ = datavec[0].fDeformedDirections.Cols();

    REAL psol = datavec[1].sol[0][0];
    REAL divQsol = datavec[0].divsol[0][0];
    REAL Qsol = 0.;
    for (int d = 0; d < 3; d++) {
        Qsol += datavec[0].sol[0][d] * axes(0, d);
    }

    bool turbulent = false;
    if (!fIsFirstIteration) {
      REAL velocity = 4. * Qsol / (M_PI * fDw * fDw);
      REAL reynolds = (fRho * std::abs(velocity) * fDw) / fMu;
      turbulent = reynolds > 1187.38;
    }
    
    // Tangent matrix
    REAL factor = turbulent ? fC * weight * (pow(std::abs(Qsol), 3. / 4.) + 3. / 4. * pow(std::abs(Qsol), 3. / 4.)) : fCLin * weight;
    ek.AddContribution(0, 0, HDivphiQ, 1, HDivphiQ, 0, factor);      // A
    ek.AddContribution(nphiQ, 0, phiP, 0, divQ, 1, -weight); // B^T
    ek.AddContribution(0, nphiQ, divQ, 0, phiP, 1, -weight); // B
    factor = -fKvw * weight;
    ek.AddContribution(nphiQ, nphiQ, phiP, 0, phiP, 1, factor); // C

    // Residual vector constitutive equation (negative)
    factor = turbulent ? fC * weight * Qsol * pow(std::abs(Qsol), 3. / 4.) : fCLin * Qsol * weight;
    ef.AddContribution(0, 0, HDivphiQ, 1, Aux3, 0, -factor);
    factor = -psol * weight;
    ef.AddContribution(0, 0, divQ, 0, Aux, 1, -factor);

    // Residual vector conservation equation (negative)
    factor = -divQsol * weight;
    ef.AddContribution(nphiQ, 0, phiP, 0, Aux, 0, -factor);
    factor = -fKvw * psol * weight;
    ef.AddContribution(nphiQ, 0, phiP, 0, Aux, 0, -factor);
    factor = fKvw * fPres * weight;
    ef.AddContribution(nphiQ, 0, phiP, 0, Aux, 0, -factor);

    TPZFNMatrix<1, REAL> dualFunctional(1, 1, 0.0);
    if (fForcingFunctionDual) {
      TPZManVector<STATE> res(1);
      fForcingFunctionDual(datavec[1].x, res);
      dualFunctional(0, 0) = res[0];
      ef.AddContribution(0, 0, divQ, 0, dualFunctional, 1, weight);
    }
}

void TPZWannMixedWell::ContributeBC(const TPZVec<TPZMaterialDataT<STATE>> &datavec, REAL weight, TPZFMatrix<STATE> &ek,
                                    TPZFMatrix<STATE> &ef, TPZBndCondT<STATE> &bc)
{

    int dim = Dimension();

    TPZFMatrix<REAL> &phiQ = datavec[0].phi;
    int phrq = phiQ.Rows();

    REAL v2 = bc.Val2()[0];
    REAL v1 = bc.Val1()(0, 0);

    switch (bc.Type())
    {
    case 0: // Dirichlet condition
        for (int iq = 0; iq < phrq; iq++)
        {
            // the contribution of the Dirichlet boundary condition appears in the flow equation
            ef(iq, 0) += (-1.) * v2 * phiQ(iq, 0) * weight;
        }
        break;

    case 1: // Neumann condition
        // for (int iq = 0; iq < phrq; iq++)
        // {
        //     ef(iq, 0) += TPZMaterial::fBigNumber * v2 * phiQ(iq, 0) * weight;
        //     for (int jq = 0; jq < phrq; jq++)
        //     {

        //         ek(iq, jq) += TPZMaterial::fBigNumber * phiQ(iq, 0) * phiQ(jq, 0) * weight;
        //     }
        // }
        break;
    }
}

void TPZWannMixedWell::ContributeResidual(const TPZVec<TPZMaterialDataT<STATE>> &datavec, REAL weight, TPZFMatrix<STATE> &ef)
{
    TPZFMatrix<REAL> &HDivphiQ = datavec[0].fDeformedDirections;
    TPZFMatrix<REAL> &phiP = datavec[1].phi;
    TPZFMatrix<REAL> &divQ = datavec[0].divphi;
    TPZFNMatrix<1, REAL> Aux(1, 1, 1.);
    const auto &axes = datavec[0].axes;
    TPZFNMatrix<3, REAL> Aux3(HDivphiQ.Rows(), 1, 0.);
    for (int d = 0; d < HDivphiQ.Rows(); d++) {
        Aux3(d, 0) = axes(0, d);
    }

    int nphiQ, nphiP;
    nphiP = phiP.Rows();
    nphiQ = datavec[0].fDeformedDirections.Cols();

    REAL psol = datavec[1].sol[0][0];
    REAL divQsol = datavec[0].divsol[0][0];
    REAL Qsol = 0.;
    for (int d = 0; d < 3; d++) {
        Qsol += datavec[0].sol[0][d] * axes(0, d);
    }

    bool turbulent = false;
    if (!fIsFirstIteration) {
        REAL velocity = 4. * Qsol / (M_PI * fDw * fDw);
        REAL reynolds = (fRho * std::abs(velocity) * fDw) / fMu;
        turbulent = reynolds > 1187.38;
    }

    // Residual vector constitutive equation (negative)
    REAL factor = turbulent ? fC * weight * Qsol * pow(std::abs(Qsol), 3. / 4.) : fCLin * Qsol * weight;
    ef.AddContribution(0, 0, HDivphiQ, 1, Aux3, 0, -factor);
    factor = -psol * weight;
    ef.AddContribution(0, 0, divQ, 0, Aux, 1, -factor);

    // Residual vector conservation equation (negative)
    factor = -divQsol * weight;
    ef.AddContribution(nphiQ, 0, phiP, 0, Aux, 0, -factor);
    factor = -fKvw * psol * weight;
    ef.AddContribution(nphiQ, 0, phiP, 0, Aux, 0, -factor);
    factor = fKvw * fPres * weight;
    ef.AddContribution(nphiQ, 0, phiP, 0, Aux, 0, -factor);

    TPZFNMatrix<1, REAL> dualFunctional(1, 1, 0.0);
    if (fForcingFunctionDual) {
      TPZManVector<STATE> res(1);
      fForcingFunctionDual(datavec[1].x, res);
      dualFunctional(0, 0) = res[0];
      ef.AddContribution(0, 0, divQ, 0, dualFunctional, 1, weight);
    }
}

void TPZWannMixedWell::Solution(const TPZVec<TPZMaterialDataT<STATE>> &datavec, int var, TPZVec<STATE> &solOut)
{
    solOut.Resize(this->NSolutionVariables(var));
    solOut.Fill(0.);
    REAL SolP, SolQ;
    const auto &axes = datavec[0].axes;

    SolQ = 0.;
    for (int d = 0; d < 3; d++) {
        SolQ += datavec[0].sol[0][d] * axes(0, d);
    }

    SolP = datavec[1].sol[0][0];

    if (var == 1)
    { // function (state variable Q)
        solOut[0] = SolQ;
        return;
    }

    if (var == 2)
    {
        solOut[0] = SolP;
        return;
    }

    if (var == 3)
    {
        solOut[0] = datavec[0].divsol[0][0];
        return;
    }
}

void TPZWannMixedWell::Errors(const TPZVec<TPZMaterialDataT<STATE>> &data, TPZVec<REAL> &errors)
{

    DebugStop();
}

int TPZWannMixedWell::VariableIndex(const std::string &name) const
{
    if (!strcmp("Flux", name.c_str()))
        return 1;
    if (!strcmp("Pressure", name.c_str()))
        return 2;
    if (!strcmp("DivFlux", name.c_str()))
        return 3;
    if (!strcmp("Divergence", name.c_str()))
        return 3;
    DebugStop();
    return -1;
}

int TPZWannMixedWell::NSolutionVariables(int var) const
{
    if (var == 1)
        return 1;
    if (var == 2)
        return 1;
    if (var == 3)
        return 1;

    DebugStop();
    return -1;
}

int TPZWannMixedWell::ClassId() const
{
    return Hash("TPZWannMixedWell") ^ TBase::ClassId() << 1;
}

TPZMaterial *TPZWannMixedWell::NewMaterial() const
{
    return new TPZWannMixedWell(*this);
}

void TPZWannMixedWell::Print(std::ostream &out) const
{
    out << "Material Name: " << this->Name() << "\n";
    out << "Material Id: " << this->Id() << "\n";
    out << "Dimension: " << this->Dimension() << "\n\n";
    out << "Well Diameter: " << fDw << "\n";
    out << "Viscosity: " << fMu << "\n";
    out << "Density: " << fRho << "\n";
    out << "Reservoir Pressure: " << fPres << "\n";
    out << "Pseudo Resistivity: " << fKvw << "\n";
}

void TPZWannMixedWell::FillDataRequirements(TPZVec<TPZMaterialDataT<STATE>> &datavec) const
{
    int nref = datavec.size();
    for (int i = 0; i < nref; i++)
    {
        datavec[i].SetAllRequirements(false);
        datavec[i].fNeedsSol = true;
    }
}

void TPZWannMixedWell::FillBoundaryConditionDataRequirements(int type, TPZVec<TPZMaterialDataT<STATE>> &datavec) const
{
    // default is no specific data requirements
    int nref = datavec.size();
    for (int iref = 0; iref < nref; iref++)
    {
        datavec[iref].SetAllRequirements(false);
    }
    datavec[0].fNeedsNormal = true;
}
#undef USEBLAS
