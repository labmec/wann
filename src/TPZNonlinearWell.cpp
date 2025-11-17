//
// Created by Giovane and Giovanni on 17/09/25
//

#include "TPZNonlinearWell.h"
#include "TPZMaterialDataT.h"
#include "pzaxestools.h"
#include "pzlog.h"

#ifdef PZ_LOG
static TPZLogger logger("pz.material.darcy");
#endif

#define USEBLAS

TPZNonlinearWell::TPZNonlinearWell() : TPZRegisterClassId(&TPZNonlinearWell::ClassId),
                                       TBase(), fDim(-1) {}

[[maybe_unused]] TPZNonlinearWell::TPZNonlinearWell(int id, REAL Dw, REAL mu, REAL rho, REAL pres, REAL Kvw) : TPZRegisterClassId(&TPZNonlinearWell::ClassId),
                                                                                                               TBase(id), fDim(1), fDw(Dw), fMu(mu), fRho(rho), fPres(pres), fKvw(Kvw)
{
    REAL c = (2.252610888 * pow(fDw, 19. / 7.)) / (pow(fMu, 1. / 7.) * pow(fRho, 3. / 7.));
    fC = pow(c, -7. / 4.);
    fCLin = 128. * fMu / (M_PI * pow(fDw, 4));
}

/**
         copy constructor
 */
TPZNonlinearWell::TPZNonlinearWell(const TPZNonlinearWell &copy) : TBase(copy), fDim(copy.fDim)
{
    *this = copy;
}
/**
         copy constructor
 */
TPZNonlinearWell &TPZNonlinearWell::operator=(const TPZNonlinearWell &copy)
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

void TPZNonlinearWell::Contribute(const TPZVec<TPZMaterialDataT<STATE>> &datavec, REAL weight, TPZFMatrix<STATE> &ek,
                                  TPZFMatrix<STATE> &ef)
{

    TPZFMatrix<REAL> &HDivphiQ = datavec[0].fDeformedDirections;
    TPZFMatrix<REAL> &phiP = datavec[1].phi;
    TPZFMatrix<REAL> &divQ = datavec[0].divphi;
    TPZFNMatrix<1, REAL> Aux(1, 1, 1.);
    TPZFNMatrix<20, REAL> phiQ(HDivphiQ.Cols(), 1);

    int nphiQ, nphiP;
    nphiP = phiP.Rows();
    nphiQ = datavec[0].fDeformedDirections.Cols();
    for (int iq = 0; iq < nphiQ; iq++)
    {
        phiQ(iq, 0) = HDivphiQ(0, iq);
    }

    REAL Qsol = datavec[0].sol[0][0];
    REAL psol = datavec[1].sol[0][0];
    REAL divQsol = datavec[0].divsol[0][0];

    REAL reynolds = (fRho * std::abs(Qsol) * fDw) / fMu;
    bool turbulent = reynolds > 1187.38;

    // Tangent matrix
    REAL factor = turbulent ? fC * weight * (pow(std::abs(Qsol), 3. / 4.) + 3. / 4. * pow(std::abs(Qsol), 3. / 4.)) : fCLin * weight;
    ek.AddContribution(0, 0, phiQ, 0, phiQ, 1, factor);      // A
    ek.AddContribution(nphiQ, 0, phiP, 0, divQ, 1, -weight); // B^T
    ek.AddContribution(0, nphiQ, divQ, 0, phiP, 1, -weight); // B
    factor = -fKvw * weight;
    ek.AddContribution(nphiQ, nphiQ, phiP, 0, phiP, 1, factor); // C

    // Residual vector constitutive equation (negative)
    factor = turbulent ? fC * weight * Qsol * pow(std::abs(Qsol), 3. / 4.) : fCLin * Qsol * weight;
    ef.AddContribution(0, 0, phiQ, 0, Aux, 0, -factor);
    factor = -psol * weight;
    ef.AddContribution(0, 0, divQ, 0, Aux, 1, -factor);

    // Residual vector conservation equation (negative)
    factor = -divQsol * weight;
    ef.AddContribution(nphiQ, 0, phiP, 0, Aux, 0, -factor);
    factor = -fKvw * psol * weight;
    ef.AddContribution(nphiQ, 0, phiP, 0, Aux, 0, -factor);
    factor = fKvw * fPres * weight;
    ef.AddContribution(nphiQ, 0, phiP, 0, Aux, 0, -factor);
}

void TPZNonlinearWell::ContributeBC(const TPZVec<TPZMaterialDataT<STATE>> &datavec, REAL weight, TPZFMatrix<STATE> &ek,
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

void TPZNonlinearWell::Solution(const TPZVec<TPZMaterialDataT<STATE>> &datavec, int var, TPZVec<STATE> &solOut)
{
    solOut.Resize(this->NSolutionVariables(var));
    solOut.Fill(0.);
    REAL SolP, SolQ;

    SolQ = datavec[0].sol[0][0];
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

    if (var == 4)
    {

        if (datavec[1].fShapeType == TPZMaterialData::EEmpty)
            return;
        TPZFNMatrix<9, REAL> dsoldx(3, 1., 0.);
        TPZFNMatrix<9, REAL> dsoldaxes(fDim, 1, 0.);

        dsoldaxes = datavec[1].dsol[0];
        TPZAxesTools<REAL>::Axes2XYZ(dsoldaxes, dsoldx, datavec[1].axes);

        for (int i = 0; i < fDim; i++)
        {
            solOut[i] = dsoldx(i, 0);
        }

        return;
    }
}

void TPZNonlinearWell::Errors(const TPZVec<TPZMaterialDataT<STATE>> &data, TPZVec<REAL> &errors)
{

    DebugStop();
}

int TPZNonlinearWell::VariableIndex(const std::string &name) const
{
    if (!strcmp("Flux", name.c_str()))
        return 1;
    if (!strcmp("Pressure", name.c_str()))
        return 2;
    if (!strcmp("DivFlux", name.c_str()))
        return 3;
    if (!strcmp("GradPressure", name.c_str()))
        return 4;
    DebugStop();
    return -1;
}

int TPZNonlinearWell::NSolutionVariables(int var) const
{
    if (var == 1)
        return 1;
    if (var == 2)
        return 1;
    if (var == 3)
        return 1;
    if (var == 4)
        return 1;

    DebugStop();
    return -1;
}

int TPZNonlinearWell::ClassId() const
{
    return Hash("TPZNonlinearWell") ^ TBase::ClassId() << 1;
}

TPZMaterial *TPZNonlinearWell::NewMaterial() const
{
    return new TPZNonlinearWell(*this);
}

void TPZNonlinearWell::Print(std::ostream &out) const
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

void TPZNonlinearWell::FillDataRequirements(TPZVec<TPZMaterialDataT<STATE>> &datavec) const
{
    int nref = datavec.size();
    for (int i = 0; i < nref; i++)
    {
        datavec[i].SetAllRequirements(false);
        datavec[i].fNeedsSol = true;
    }
}

void TPZNonlinearWell::FillBoundaryConditionDataRequirements(int type, TPZVec<TPZMaterialDataT<STATE>> &datavec) const
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
