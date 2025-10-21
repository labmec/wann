//
// Created by Gustavo Batistela on 5/13/21.
//

#include "TPZWannMixedDarcyNL.h"
#include "TPZMaterialDataT.h"
#include "pzaxestools.h"
#include "pzlog.h"

#ifdef PZ_LOG
static TPZLogger logger("pz.material.darcy");
#endif

#define USEBLAS

TPZWannMixedDarcyNL::TPZWannMixedDarcyNL() : TPZRegisterClassId(&TPZWannMixedDarcyNL::ClassId),
                                         TBase() {}

[[maybe_unused]] TPZWannMixedDarcyNL::TPZWannMixedDarcyNL(int id, int dim) : TPZRegisterClassId(&TPZWannMixedDarcyNL::ClassId),
                                                        TBase(id,dim)
{
}

/**
         copy constructor
 */
TPZWannMixedDarcyNL::TPZWannMixedDarcyNL(const TPZWannMixedDarcyNL &copy) : TBase(copy)
{
    *this = copy;
}
/**
         copy constructor
 */
TPZWannMixedDarcyNL &TPZWannMixedDarcyNL::operator=(const TPZWannMixedDarcyNL &copy)
{
    TBase::operator=(copy);
    return *this;
}

void TPZWannMixedDarcyNL::Contribute(const TPZVec<TPZMaterialDataT<STATE>> &datavec, REAL weight, TPZFMatrix<STATE> &ek,
                                     TPZFMatrix<STATE> &ef)
{

    TPZFMatrix<REAL> &phiQ = datavec[0].fDeformedDirections;
    TPZFMatrix<REAL> &phiP = datavec[1].phi;
    TPZFMatrix<REAL> &divQ = datavec[0].divphi;
    TPZFNMatrix<1,REAL> Aux(1,1,1.);
    
    int nphiQ, nphiP;
    nphiP = phiP.Rows();
    nphiQ = datavec[0].fDeformedDirections.Cols();
    
    TPZVec<STATE> Qsol = datavec[0].sol[0];
    TPZFMatrix<STATE> QsolMat(Qsol.size(), 1);
    for (int i = 0; i < Qsol.size(); ++i) {
        QsolMat(i, 0) = Qsol[i];
    }
    REAL psol = datavec[1].sol[0][0];
    REAL divQsol = datavec[0].divsol[0][0];

    //Tangent matrix
    REAL K = GetPermeability(datavec[0].x);
    REAL factor = weight/K;
    ek.AddContribution(0, 0, phiQ, 1, phiQ, 0, factor);         // A
    ek.AddContribution(nphiQ, 0, phiP, 0, divQ, 1, -weight); // B^T
    ek.AddContribution(0, nphiQ, divQ, 0, phiP, 1, -weight); // B

    // Residual vector constitutive equation (negative)
    ef.AddContribution(0, 0, phiQ, 1, QsolMat, 0, -factor);
    factor = psol * weight;
    ef.AddContribution(0, 0, divQ, 0, Aux, 1, factor);

    // Residual vector conservation equation (negative)
    factor = divQsol * weight;
    ef.AddContribution(nphiQ, 0, phiP, 0, Aux, 0, factor);
}

void TPZWannMixedDarcyNL::ContributeBC(const TPZVec<TPZMaterialDataT<STATE>> &datavec, REAL weight, TPZFMatrix<STATE> &ek,
                                     TPZFMatrix<STATE> &ef, TPZBndCondT<STATE> &bc) {

    int dim = Dimension();

    REAL bigNumber = TPZMaterial::fBigNumber * 1.e-2;

    TPZFMatrix<REAL> &phiQ = datavec[0].phi;
    int phrq = phiQ.Rows();
    TPZManVector<STATE,3> qsol = datavec[0].sol[0];

    REAL v2 = bc.Val2()[0];
    REAL v1 = bc.Val1()(0, 0);
    REAL u_D = 0;
    REAL normflux = 0.;

    if (bc.HasForcingFunctionBC()) {
        TPZManVector<STATE> res(3);
        TPZFNMatrix<9, STATE> gradu(3, 1);
        bc.ForcingFunctionBC()(datavec[0].x, res, gradu);

        const STATE perm = GetPermeability(datavec[0].x);

        for (int i = 0; i < 3; i++) {
            normflux += datavec[0].normal[i] * perm * gradu(i, 0);
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
        case 0 :        // Dirichlet condition
            for (int iq = 0; iq < phrq; iq++) {
                //the contribution of the Dirichlet boundary condition appears in the flow equation
                ef(iq, 0) += (-1.) * v2 * phiQ(iq, 0) * weight;
            }
            break;

        case 1 :            // Neumann condition
            for (int iq = 0; iq < phrq; iq++) {
                REAL qn = qsol[0];
                ef(iq, 0) += bigNumber * (v2-qn) * phiQ(iq, 0) * weight;
                for (int jq = 0; jq < phrq; jq++) {

                    ek(iq, jq) += bigNumber * phiQ(iq, 0) * phiQ(jq, 0) * weight;
                }
            }
            break;
    }
}

void TPZWannMixedDarcyNL::FillDataRequirements(TPZVec<TPZMaterialDataT<STATE>> &datavec) const {
    int nref = datavec.size();
    for (int i = 0; i < nref; i++) {
        datavec[i].SetAllRequirements(false);
        datavec[i].fNeedsSol = true;
    }
}

void
TPZWannMixedDarcyNL::FillBoundaryConditionDataRequirements(int type, TPZVec<TPZMaterialDataT<STATE>> &datavec) const {
    // default is no specific data requirements
    int nref = datavec.size();
    for (int iref = 0; iref < nref; iref++) {
        datavec[iref].SetAllRequirements(false);
        datavec[iref].fNeedsSol = true;
        datavec[iref].fNeedsNormal = true;
    }
}
