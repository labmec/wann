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
    TPZFMatrix<STATE> perm(3, 3, 0.);
    GetPermeability(x, perm);
    TPZFMatrix<STATE> perm_dot_dphi;
    perm.Multiply(dphi, perm_dot_dphi);
    ek.AddContribution(0, 0, perm_dot_dphi, 1, dphi, 0, weight);

    // Residual vector
    ef.AddContribution(0, 0, perm_dot_dphi, 1, dpsol, 0, -weight);
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

    // TODO: Needs refactor
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
            // Nothing to do. Already done in the initial solution
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
                    << "\t 1: Neumann\n";
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

void TPZWannDarcyNL::Solution(const TPZMaterialDataT<STATE> &data, int var, TPZVec<STATE> &solOut) {

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
        case 7: {
            // MinusKGradU/Flux;
            TPZFNMatrix<9, STATE> dsoldx;
            TPZAxesTools<STATE>::Axes2XYZ(data.dsol[0], dsoldx, data.axes);
            TPZFMatrix<STATE> perm(3, 3, 0.);
            GetPermeability(data.x, perm);
            TPZFMatrix<STATE> perm_dot_dsoldx;
            perm.Multiply(dsoldx, perm_dot_dsoldx);
            for (int id = 0; id < 3; id++) {
                solOut[id] = - perm_dot_dsoldx(id, 0);
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

    default: {
            PZError << __PRETTY_FUNCTION__ << "\n Post-processing variable index not implemented!\n";
            DebugStop();
        }
    }
}