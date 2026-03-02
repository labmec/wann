//
// Created by Gi Taraschi 23/02/2026.
//

#include "TPZWannDarcyNL.h"
#include "TPZMaterialDataT.h"
#include "pzaxestools.h"
#include "pzlog.h"

#ifdef PZ_LOG
static TPZLogger logger("pz.material.darcy");
#endif

#define USEBLAS

TPZWannDarcyNL::TPZWannDarcyNL() : TPZRegisterClassId(&TPZWannDarcyNL::ClassId),
                                             TBase() {}

[[maybe_unused]] TPZWannDarcyNL::TPZWannDarcyNL(int id, int dim) : TPZRegisterClassId(&TPZWannDarcyNL::ClassId),
                                                                             TBase(id, dim)
{
}

/**
         copy constructor
 */
TPZWannDarcyNL::TPZWannDarcyNL(const TPZWannDarcyNL &copy) : TBase(copy)
{
    *this = copy;
}
/**
         copy constructor
 */
TPZWannDarcyNL &TPZWannDarcyNL::operator=(const TPZWannDarcyNL &copy)
{
    TBase::operator=(copy);
    return *this;
}

void TPZWannDarcyNL::Contribute(const TPZMaterialDataT<STATE> &data, REAL weight, TPZFMatrix<STATE> &ek,
                                     TPZFMatrix<STATE> &ef)
{

    const TPZFMatrix<REAL> &phi = data.phi;
    const TPZFMatrix<REAL> &dphi = data.dphix;
    const TPZVec<REAL> &x = data.x;
    TPZFNMatrix<1, REAL> Aux(1, 1, 1.);

    TPZFNMatrix dpsol = data.dsol[0];

    STATE source_term = 0;
    if (this->HasForcingFunction()) {
        TPZManVector<STATE, 1> res(1);
        fForcingFunction(x, res);
        source_term = -res[0];
    }

    // Tangent matrix
    REAL K = GetPermeability(x);
    REAL factor = weight * K;
    ek.AddContribution(0, 0, dphi, 1, dphi, 0, factor);

    // Residual vector
    ef.AddContribution(0, 0, dphi, 1, dpsol, 0, -factor);
    ef.AddContribution(0, 0, phi, 0, Aux, 0, -source_term * weight);
}

void TPZWannDarcyNL::ContributeBC(const TPZMaterialDataT<STATE> &data, REAL weight, TPZFMatrix<STATE> &ek,
                                       TPZFMatrix<STATE> &ef, TPZBndCondT<STATE> &bc)
{

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
        const STATE perm = GetPermeability(data.x);
        TPZManVector<STATE, 3> Flux(fDim, 0.);
        for (int id = 0; id < fDim; id++) {
            Flux[id] = - perm * mat_val(id, 0);
        }
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
            // Already done in the initial solution
            break;
        case 1 : // Neumann condition
            for (in = 0; in < phi.Rows(); in++) {
                ef(in, 0) += -v2 * (STATE) (phi(in, 0) * weight);
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

void TPZWannDarcyNL::FillDataRequirements(TPZMaterialData &data) const
{
    data.SetAllRequirements(false);
    data.fNeedsSol = true;
}

void TPZWannDarcyNL::FillBoundaryConditionDataRequirements(int type, TPZMaterialData &data) const
{
    // default is no specific data requirements
    data.SetAllRequirements(false);
    data.fNeedsSol = false;
    data.fNeedsNormal = true;
}

void TPZWannDarcyNL::GetSolDimensions(uint64_t &u_len, uint64_t &du_row, uint64_t &du_col) const {
    u_len=1;
    du_row=fDim;
    du_col=1;
}