#include "TPZMixedDarcyAnisotropic.h"
#include "TPZMaterialDataT.h"
#include "pzaxestools.h"
#include "pzlog.h"

#ifdef PZ_LOG
static TPZLogger logger("pz.material.darcy");
#endif

#define USEBLAS

TPZMixedDarcyAnisotropic::TPZMixedDarcyAnisotropic() : TPZRegisterClassId(&TPZMixedDarcyAnisotropic::ClassId),
                                                       TBase(), fDim(-1) {}

[[maybe_unused]] TPZMixedDarcyAnisotropic::TPZMixedDarcyAnisotropic(int id, int dim) : TPZRegisterClassId(&TPZMixedDarcyAnisotropic::ClassId),
                                                                                      TBase(id), fDim(dim) {}

TPZMixedDarcyAnisotropic::TPZMixedDarcyAnisotropic(const TPZMixedDarcyAnisotropic &copy) : TBase(copy), fDim(copy.fDim)
{
    *this = copy;
}

TPZMixedDarcyAnisotropic &TPZMixedDarcyAnisotropic::operator=(const TPZMixedDarcyAnisotropic &copy)
{
    TBase::operator=(copy);
    fDim = copy.fDim;
    return *this;
}

void TPZMixedDarcyAnisotropic::Contribute(const TPZVec<TPZMaterialDataT<STATE>> &datavec, REAL weight, TPZFMatrix<STATE> &ek,
                                          TPZFMatrix<STATE> &ef)
{

    TPZFMatrix<REAL> &phiQ = datavec[0].fDeformedDirections;
    TPZFMatrix<REAL> &phiP = datavec[1].phi;
    TPZFMatrix<REAL> &divQ = datavec[0].divphi;
    TPZFNMatrix<1, REAL> Aux(1, 1, 1.);

    int nphiQ, nphiP;
    nphiP = phiP.Rows();
    nphiQ = datavec[0].fDeformedDirections.Cols();


    STATE source_term = 0;
    if (this->HasForcingFunction()) {
        TPZManVector<STATE, 1> res(1);
        fForcingFunction(datavec[0].x, res);
        source_term = -res[0];
    }

    // Stiffness matrix
    TPZFMatrix<STATE> InvPerm(3, 3, 0.);
    GetInversePermeability(datavec[0].x, InvPerm);
    TPZFMatrix<REAL> KappaInv_dot_tau;
    InvPerm.Multiply(phiQ, KappaInv_dot_tau);
    ek.AddContribution(0, 0, KappaInv_dot_tau, 1, phiQ, 0, weight); // A
    ek.AddContribution(nphiQ, 0, phiP, 0, divQ, 1, -weight); // B^T
    ek.AddContribution(0, nphiQ, divQ, 0, phiP, 1, -weight); // B

    // Source term
    ef.AddContribution(nphiQ, 0, phiP, 0, Aux, 0, -source_term * weight);
}

void TPZMixedDarcyAnisotropic::ContributeBC(const TPZVec<TPZMaterialDataT<STATE>> &datavec, REAL weight, TPZFMatrix<STATE> &ek,
                                            TPZFMatrix<STATE> &ef, TPZBndCondT<STATE> &bc)
{
    int dim = Dimension();
    TPZFMatrix<REAL> &phiQ = datavec[0].phi;
    int phrq = phiQ.Rows();

    REAL v2 = bc.Val2()[0];
    REAL v1 = bc.Val1()(0, 0);
    REAL u_D = 0;
    REAL normflux = 0.;

    if (bc.HasForcingFunctionBC()) {
        TPZManVector<STATE> res(3);
        TPZFNMatrix<9, STATE> gradu(3, 1);
        bc.ForcingFunctionBC()(datavec[0].x, res, gradu);

        TPZFMatrix<STATE> perm_matrix(3, 3, 0.);
        GetPermeability(datavec[0].x, perm_matrix);

        for (int i = 0; i < dim; i++) {
            for (int j = 0; j < dim; j++) {
                normflux += datavec[0].normal[i] * perm_matrix(i, j) * gradu(j, 0);
            }
        }

        if (bc.Type() == 0 || bc.Type() == 4) {
            v2 = res[0];
            u_D = res[0];
            normflux *= (-1.);
        } else if (bc.Type() == 1 || bc.Type() == 2) {
            v2 = -normflux;
            if (bc.Type() == 2) {
                v2 = -res[0] + v2 / v1;
            }
        } else if (bc.Type() == 5) {
            v2 = res[0];
        } else {
            DebugStop();
        }
    } else {
        v2 = bc.Val2()[0];
    }

    switch (bc.Type()) {
        case 0 :
            for (int iq = 0; iq < phrq; iq++) {
                ef(iq, 0) += (-1.) * v2 * phiQ(iq, 0) * weight;
            }
            break;
        case 1 :
            for (int iq = 0; iq < phrq; iq++) {
                ef(iq, 0) += TPZMaterial::fBigNumber * v2 * phiQ(iq, 0) * weight;
                for (int jq = 0; jq < phrq; jq++) {
                    ek(iq, jq) += TPZMaterial::fBigNumber * phiQ(iq, 0) * phiQ(jq, 0) * weight;
                }
            }
            break;
        case 2 :
            for (int iq = 0; iq < phrq; iq++) {
                ef(iq, 0) += v2 * phiQ(iq, 0) * weight;
                for (int jq = 0; jq < phrq; jq++) {
                    ek(iq, jq) += weight / v1 * phiQ(iq, 0) * phiQ(jq, 0);
                }
            }
            break;
        case 4:
            if (IsZero(bc.Val1()(0, 0))) {

                for (int iq = 0; iq < phrq; iq++) {
                    ef(iq, 0) += TPZMaterial::fBigNumber * normflux * phiQ(iq, 0) * weight;
                    for (int jq = 0; jq < phrq; jq++) {
                        ek(iq, jq) += TPZMaterial::fBigNumber * phiQ(iq, 0) * phiQ(jq, 0) * weight;
                    }
                }

            } else {

                REAL InvKm = 1. / bc.Val1()(0, 0);
                REAL g = normflux;
                for (int in = 0; in < phiQ.Rows(); in++) {
                    ef(in, 0) += (STATE) (InvKm * g - u_D) * phiQ(in, 0) * weight;
                    for (int jn = 0; jn < phiQ.Rows(); jn++) {
                        ek(in, jn) += (STATE) (InvKm * phiQ(in, 0) * phiQ(jn, 0) * weight);
                    }
                }
            }

            break;

        case 5:
            TPZFMatrix<REAL> &phi = datavec[0].fH1.fPhi;
            for (int in = 0; in < phi.Rows(); in++) {
                ef(in, 0) += TPZMaterial::fBigNumber * v2 * phi(in, 0) * weight;
                for (int jn = 0; jn < phi.Rows(); jn++) {
                    ek(in, jn) += TPZMaterial::fBigNumber * phi(in, 0) * phi(jn, 0) * weight;
                }
            }

    }
}

void TPZMixedDarcyAnisotropic::Solution(const TPZVec<TPZMaterialDataT<STATE>> &datavec, int var, TPZVec<STATE> &solOut)
{
    solOut.Resize(this->NSolutionVariables(var));
    solOut.Fill(0.);
    TPZManVector<STATE, 10> SolP, SolQ;
    TPZFMatrix<STATE> perm_matrix(3, 3, 0.);
    GetPermeability(datavec[0].x, perm_matrix);

    SolP = datavec[1].sol[0];
    if(SolP.size() == 0) SolP.Resize(1,0.);

    if (var == 1) {
        for (int i = 0; i < fDim; i++) {
            solOut[i] = datavec[0].sol[0][i];
        }
        return;
    }

    if (var == 2) {
        solOut[0] = SolP[0];
        return;
    }

    if (var == 3) {
        solOut[0] = datavec[0].dsol[0](0, 0);
        solOut[1] = datavec[0].dsol[0](1, 0);
        solOut[2] = datavec[0].dsol[0](2, 0);
        return;
    }

    if (var == 4) {
        solOut[0] = datavec[0].dsol[0](0, 1);
        solOut[1] = datavec[0].dsol[0](1, 1);
        solOut[2] = datavec[0].dsol[0](2, 1);
        return;
    }

    if (var == 5) {
        solOut[0] = datavec[0].divsol[0][0];
        return;
    }

    if (var == 6) {
        TPZVec<STATE> exactSol(1);
        TPZFMatrix<STATE> flux(3, 1);
        if (fExactSol) {
            fExactSol(datavec[0].x, exactSol, flux);
        }
        solOut[0] = exactSol[0];
        return;
    }

    if (var == 7) {
        TPZVec<STATE> exactSol(1);
        TPZFMatrix<STATE> gradu(fDim, 1);

        if (fExactSol) {
            fExactSol(datavec[0].x, exactSol, gradu);
        }

        for (int i = 0; i < fDim; i++) {
            solOut[i] = 0.;
            for (int j = 0; j < fDim; j++) {
                solOut[i] -= perm_matrix(i, j) * gradu(j, 0);
            }
        }

        return;
    }

    if (var == 8) {
        solOut[0] = datavec[1].p;
        return;
    }

    if (var == 9) {

        if(datavec[1].fShapeType == TPZMaterialData::EEmpty) return;
        TPZFNMatrix<9, REAL> dsoldx(3, 1.,0.);
        TPZFNMatrix<9, REAL> dsoldaxes(fDim, 1,0.);

        dsoldaxes = datavec[1].dsol[0];
        TPZAxesTools<REAL>::Axes2XYZ(dsoldaxes, dsoldx, datavec[1].axes);

        for (int i = 0; i < fDim; i++) {
            solOut[i] = dsoldx(i, 0);
        }

        return;
    }

    if (var == 10) {
        solOut[0] = datavec[0].divsol[0][0];
        return;
    }

    if (var == 11) {
        TPZVec<STATE> exactSol(1);
        TPZFMatrix<STATE> flux(3, 1);
        fExactSol(datavec[0].x, exactSol, flux);
        solOut[0] = flux(2, 0);
        return;
    }
    if (var == 12) {
        TPZVec<STATE> exactSol(1);
        TPZFMatrix<STATE> gradu(fDim, 1);
        fExactSol(datavec[0].x, exactSol, gradu);
        for (int i = 0; i < fDim; i++) {
            solOut[i] = 0.;
            for (int j = 0; j < fDim; j++) {
                solOut[i] -= perm_matrix(i, j) * gradu(j, 0);
            }
        }
        return;
    }

    if (var == 13) {
        solOut[0] = datavec[0].divsol[0][0];
        return;
    }

    if (var == 14) {
        solOut[0] = datavec[2].sol[0][0];
        return;
    }
    if (var == 15) {
        solOut[0] = datavec[3].sol[0][0];
        return;
    }
    if (var == 16) {
        STATE infinitesimal = 0.0000000001;
        TPZManVector<REAL, 3> inf = {infinitesimal, infinitesimal, infinitesimal};
        TPZVec<STATE> exactSol(1);
        TPZFMatrix<STATE> gradu(3, 1);
        if (fExactSol) {
            if (datavec[0].x[0] == 0. && datavec[0].x[1] == 0.) {
                fExactSol(inf, exactSol, gradu);
            } else {
                fExactSol(datavec[0].x, exactSol, gradu);
            }
        }
        for (int i = 0; i < 3; i++) {
            solOut[i] = 0.;
            for (int j = 0; j < fDim; j++) {
                solOut[i] -= perm_matrix(i, j) * gradu(j, 0);
            }
        }

        return;
    }

    if (var == 17) {
        TPZVec<STATE> divsigma(1, 0.);
        if (fForcingFunction) {
            fForcingFunction(datavec[0].x, divsigma);
        }
        solOut[0] = divsigma[0];
        return;
    }

    if (var == 18) {
        solOut[0] = datavec[0].sol[0][0];
        return;
    }
}

void TPZMixedDarcyAnisotropic::Errors(const TPZVec<TPZMaterialDataT<STATE>> &data, TPZVec<REAL> &errors)
{
    errors.Resize(NEvalErrors());
    errors.Fill(0.0);

    TPZManVector<STATE, 3> fluxfem(3), pressurefem(1,0);
    fluxfem = data[0].sol[0];
    STATE divsigmafem = data[0].divsol[0][0];

    TPZManVector<STATE,1> divsigma(1,0.);

    TPZManVector<STATE,1> u_exact(1, 0);
    TPZFMatrix<STATE> du_exact(3, 1, 0);
    if (this->fExactSol) {
        this->fExactSol(data[0].x, u_exact, du_exact);
    }
    if (this->fForcingFunction) {
        this->fForcingFunction(data[0].x, divsigma);
    }

    REAL residual = (divsigma[0] - divsigmafem) * (divsigma[0] - divsigmafem);
    if(data[1].sol[0].size())
        pressurefem[0] = data[1].sol[0][0];

    TPZFMatrix<STATE> perm_matrix(3, 3, 0.);
    GetPermeability(data[0].x, perm_matrix);

    TPZManVector<STATE, 3> gradpressurefem(3, 0.);
    this->Solution(data, VariableIndex("GradPressure"), gradpressurefem);

    TPZManVector<STATE, 3> fluxexact(3, 0);
    TPZManVector<STATE, 3> gradpressure(3, 0);
    for (int i = 0; i < 3; i++) {
        gradpressure[i] = du_exact[i];
        fluxexact[i] = 0.;
        for (int j = 0; j < fDim; j++) {
            fluxexact[i] -= perm_matrix(i, j) * gradpressure[j];
        }
    }

    REAL L2flux = 0., L2grad = 0.;
    for (int i = 0; i < 3; i++) {
        L2flux += (fluxfem[i] - fluxexact[i]) * (fluxfem[i] - fluxexact[i]);
        L2grad += (du_exact[i] - gradpressurefem[i]) * (du_exact[i] - gradpressurefem[i]);
    }
    errors[0] = (pressurefem[0] - u_exact[0]) * (pressurefem[0] - u_exact[0]);
    errors[1] = L2flux;
    errors[2] = residual;
    errors[3] = L2grad;
    errors[4] = L2flux + residual;
}

int TPZMixedDarcyAnisotropic::VariableIndex(const std::string &name) const
{
    if (!strcmp("Flux", name.c_str())) return 1;
    if (!strcmp("Pressure", name.c_str())) return 2;
    if (!strcmp("GradFluxX", name.c_str())) return 3;
    if (!strcmp("GradFluxY", name.c_str())) return 4;
    if (!strcmp("DivFlux", name.c_str())) return 5;
    if (!strcmp("ExactPressure", name.c_str())) return 6;
    if (!strcmp("ExactFlux", name.c_str())) return 7;
    if (!strcmp("POrder", name.c_str())) return 8;
    if (!strcmp("GradPressure", name.c_str())) return 9;
    if (!strcmp("Divergence", name.c_str())) return 10;
    if (!strcmp("ExactDiv", name.c_str())) return 11;
    if (!strcmp("Derivative", name.c_str())) return 12;
    if (!strcmp("Permeability", name.c_str())) return 13;
    if (!strcmp("g_average", name.c_str())) return 14;
    if (!strcmp("u_average", name.c_str())) return 15;
    if (!strcmp("ExactFluxShiftedOrigin", name.c_str())) return 16;
    if (!strcmp("EstimatedError", name.c_str())) return 100;
    if (!strcmp("TrueError", name.c_str())) return 101;
    if (!strcmp("EffectivityIndex", name.c_str())) return 102;
    if (!strcmp("ExactDivSigma", name.c_str())) return 17;
    if (!strcmp("normalFlux", name.c_str())) return 18;
    DebugStop();
    return -1;
}

int TPZMixedDarcyAnisotropic::NSolutionVariables(int var) const
{
    if (var == 1) return 3;
    if (var == 2) return 1;
    if (var == 3) return 3;
    if (var == 4) return 3;
    if (var == 5) return 1;
    if (var == 6) return 1;
    if (var == 7) return 3;
    if (var == 8) return 1;
    if (var == 9) return 3;
    if (var == 10 || var == 11) return 1;
    if (var == 12) return 3;
    if (var == 13) return 1;
    if (var == 14) return 1;
    if (var == 15) return 1;
    if (var == 16) return 3;
    if (var == 17) return 1;
    if (var == 18) return 1;
    if (var == 100) return 1;
    if (var == 101) return 1;
    if (var == 102) return 1;
    DebugStop();
    return -1;
}

void TPZMixedDarcyAnisotropic::SetDimension(int dim)
{
    if (dim > 3 || dim < 1) DebugStop();
    fDim = dim;
}

int TPZMixedDarcyAnisotropic::ClassId() const
{
    return Hash("TPZMixedDarcyAnisotropic") ^ TBase::ClassId() << 1;
}

TPZMaterial *TPZMixedDarcyAnisotropic::NewMaterial() const
{
    return new TPZMixedDarcyAnisotropic(*this);
}

void TPZMixedDarcyAnisotropic::Print(std::ostream &out) const
{
    out << "Material Name: " << this->Name() << "\n";
    out << "Material Id: " << this->Id() << "\n";
    out << "Dimension: " << this->Dimension() << "\n\n";
}

void TPZMixedDarcyAnisotropic::FillDataRequirements(TPZVec<TPZMaterialDataT<STATE>> &datavec) const
{
    int nref = datavec.size();
    for (int i = 0; i < nref; i++) {
        datavec[i].SetAllRequirements(false);
        datavec[i].fNeedsNeighborSol = false;
        datavec[i].fNeedsNeighborCenter = false;
        datavec[i].fNeedsNormal = false;
        datavec[i].fNeedsHSize = false;
    }
}

void TPZMixedDarcyAnisotropic::FillBoundaryConditionDataRequirements(int type, TPZVec<TPZMaterialDataT<STATE>> &datavec) const
{
    int nref = datavec.size();
    for (int iref = 0; iref < nref; iref++) {
        datavec[iref].SetAllRequirements(false);
        datavec[iref].fNeedsSol = false;
    }
    datavec[0].fNeedsNormal = true;
    if (type == 50) {
        for (int iref = 0; iref < nref; iref++) {
            datavec[iref].fNeedsSol = false;
        }
    }
}

#undef USEBLAS