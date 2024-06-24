from firedrake import *
import ufl
import numpy as np
import math
import time as tm
import sys

solver_cg_wipc = {"ksp_type": "cg", "pc_type": "bjacobi", "pc_sub_type": "ilu"}

wellbore_id = 111
reservoir_id = 112

reservoir_farfield = 100
wellbore_heel = 101
wellbore_toe = 102
wellbore_cylinder = 103

Pr=10;
Pw=1;

#meshfile = "../geo/mesh3D_rev03.msh"
meshfile = "../geo/mesh2D_rev01.msh"
mesh = Mesh(meshfile)

#CG
V = FunctionSpace(mesh, "CG",2)
Q = VectorFunctionSpace(mesh, "CG",2) # for post procesing
W = V*V #wellbore and reservoir CG spaces

V_DG0_r = FunctionSpace(mesh, "DG", 0)
V_DG0_w = FunctionSpace(mesh, "DG", 0)

# Heaviside step function in reservoir
I_r = Function(V_DG0_r)
par_loop(("{[i] : 0 <= i < f.dofs}", "f[i, 0] = 1.0"),
         dx(reservoir_id),
         {"f": (I_r, WRITE)})
I_cg_r = Function(V)
par_loop(("{[i] : 0 <= i < A.dofs}", "A[i, 0] = fmax(A[i, 0], B[0, 0])"),
         dx,
         {"A": (I_cg_r, RW), "B": (I_r, READ)})

# Heaviside step function in wellbore
I_w = Function(V_DG0_w)
par_loop(("{[i] : 0 <= i < f.dofs}", "f[i, 0] = 1.0"),
         dx(wellbore_id),
         {"f": (I_w, WRITE)})
I_cg_w = Function(V)
par_loop(("{[i] : 0 <= i < A.dofs}", "A[i, 0] = fmax(A[i, 0], B[0, 0])"),
         dx,
         {"A": (I_cg_w, RW), "B": (I_w, READ)})

VTKFile("I_r.pvd").write(I_r)
VTKFile("I_w.pvd").write(I_w)
VTKFile("I_cg_r.pvd").write(I_cg_r)
VTKFile("I_cg_w.pvd").write(I_cg_w)

# We use indicator functions to construct normal unit vector outward to the reservoir domain at the wellbore-reservoir interface:
n_vec = FacetNormal(mesh)
n_int = I_w("+") * n_vec("+") + I_w("-") * n_vec("-")

class MyBC(DirichletBC):
    def __init__(self, V, value, markers):
        # Call superclass init
        # We provide a dummy subdomain id.
        super(MyBC, self).__init__(V, value, 0)
        # Override the "nodes" property which says where the boundary
        # condition is to be applied.
        self.nodes = np.unique(np.where(markers.dat.data_ro_with_halos == 0)[0])


# esse exemplo resolve poço e reservatório tudo junto, mas a condição de contorno no poço também fica no reservatório, mesmo com big number
k = 1 
big = 1.e6
#CG
#u = TrialFunction(V)
#v = TestFunction(V)
#p = Function(V)
#F  = inner(grad(u),grad(v))*dx(reservoir_id) + Constant(0)*v*dx(reservoir_id) + inner(k*grad(u('+')),n_vec('+'))*v('+')*dS(wellbore_cylinder) # reservatório
#F += inner(k*grad(u),grad(v))*dx(wellbore_id) + Constant(0)*v*dx(wellbore_id) + inner(grad(u('-')),n_vec('-'))*v('-')*dS(wellbore_cylinder) #+ big*u('-')*v('-')*dS(wellbore_heel) # poço
#_lhs = lhs(F)
#_rhs = rhs(F)
#p_farfield = DirichletBC(V, Pr, [reservoir_farfield])
#p_well = DirichletBC(V, Pw, [wellbore_heel])
#problem = LinearVariationalProblem(_lhs, _rhs, p, bcs=[p_farfield, p_well]) # condição do poço também fica no reservatório
#problem = LinearVariationalProblem(_lhs, _rhs, p, bcs=[p_farfield])
#solver  = LinearVariationalSolver(problem,solver_parameters=solver_cg_wipc)
#solver.solve()

#vel = Function(W).interpolate(-k*grad(p))
#File("p_all.pvd").write(p)
#File("vel.pvd").write(vel)

#mixed
#sigma, u = TrialFunctions(M)
#tau, v = TestFunctions(M)
#m =  Function(M)
#
#F = inner(sigma, tau)*dx(reservoir_id) + div(tau)*u*dx(reservoir_id) + div(sigma)*v*dx(reservoir_id) + Constant(0)*v*dx(reservoir_id) + inner(k*grad(u('+')),n_vec('+'))*v('+')*dS(wellbore_cylinder)
#F += inner(sigma, tau)*dx(wellbore_id) + div(tau)*u*dx(wellbore_id) + div(sigma)*v*dx(wellbore_id) + Constant(0)*v*dx(wellbore_id) + inner(grad(u('-')),n_vec('-'))*v('-')*dS(wellbore_cylinder)
#F += big*Pr*u*v*dx(reservoir_id) + big*Pr*v*dx(reservoir_id)
#
#_lhs = lhs(F)
#_rhs = rhs(F)
#p_farfield = DirichletBC(M.sub(1), Pr, [reservoir_farfield])
#p_well = DirichletBC(M.sub(1), Pw, [wellbore_heel])
#v_toe = DirichletBC(M.sub(0), as_vector([0,0]), [wellbore_toe])
##v_farfiled = DirichletBC(M.sub(0), as_vector([0,-10]), [wellbore_heel])
#problem = LinearVariationalProblem(_lhs, _rhs, m, bcs=[p_farfield, p_well, v_toe])#, v_farfiled]) # condição do poço também fica no reservatório
##problem = LinearVariationalProblem(_lhs, _rhs, p, bcs=[p_farfield])
#solver  = LinearVariationalSolver(problem,solver_parameters=solver_cg_wipc)
#solver.solve()

#sigma, u = m.subfunctions
#File("p_all.pvd").write(u)
#File("vel.pvd").write(sigma)

#sys.exit("exit")

u_w, u_r = TrialFunction(W)
v_w, v_r = TestFunction(W)
p_wr = Function(W)

k = 100.

vel_w = -k*grad(u_w) # vel inside wellbore -> +
vel_r = -grad(u_r) # vel inside reservoi -> -

F = inner(k*grad(u_w),grad(v_w))*dx(wellbore_id) + Constant(0)*v_w*dx(wellbore_id)
F+= big*(u_w('+')-u_r('+'))*v_w('+')*dS(wellbore_cylinder)
#F+= big*(u_w('+')-u_r('+'))*v_w('+')*dS(wellbore_toe)
F+= inner(grad(u_r),grad(v_r))*dx(reservoir_id) + Constant(0)*v_r*dx(reservoir_id)
F+= big*(u_r('+')-u_w('+'))*v_r('+')*dS(wellbore_cylinder)
#F+= big*(u_r('+')-u_w('+'))*v_r('+')*dS(wellbore_toe)

#F+= (inner(vel_r('-'),n_vec('-'))*v_w('-')*dS(wellbore_cylinder) - inner(vel_w('+'),n_vec('+'))*v_r('+')*dS(wellbore_cylinder))
#F+= (u_w('+')*v_r('+')-u_r('+')*v_w('+'))*dS(wellbore_toe)
#F+= (u_w('+')*v_r('+')-u_r('+')*v_w('+'))*dS(wellbore_cylinder)
#- inner(-grad(u_r('-')),n_vec('-'))*v_w('-')*dS(wellbore_cylinder) + inner(-k*grad(u_w('+')),n_vec('+'))*v_r('+')*dS(wellbore_cylinder)

p_well = DirichletBC(W.sub(0), Pw, [wellbore_heel])
exclude_beyond_wellbore = MyBC(W.sub(0), 0, I_cg_w)

p_farfield = DirichletBC(W.sub(1), Pr, [reservoir_farfield])
exclude_beyond_reservoir = MyBC(W.sub(1), 0, I_cg_r)

_lhs = lhs(F)
_rhs = rhs(F)

problem = LinearVariationalProblem(_lhs, _rhs, p_wr, bcs=[p_well, exclude_beyond_wellbore, p_farfield, exclude_beyond_reservoir])
solver  = LinearVariationalSolver(problem)#, solver_parameters=solver_cg_wipc)
solver.solve()


p_w, p_r = p_wr.subfunctions
VTKFile("p_reservoir.pvd").write(p_r)
VTKFile("p_wellbore.pvd").write(p_w)

pwc = assemble(p_w('-')*dS(wellbore_cylinder))
prc = assemble(p_r('-')*dS(wellbore_cylinder))

vel_w = -k*grad(p_w)
_vel_w = Function(Q).interpolate(vel_w)
flux_w_t = assemble(inner(vel_w('+'),n_vec('+'))*dS(wellbore_toe))
flux_w_h = assemble(inner(vel_w('+'),n_vec('+'))*dS(wellbore_heel))
flux_w_c = assemble(inner(vel_w('+'),n_vec('+'))*dS(wellbore_cylinder))

vel_r = -grad(p_r)
_vel_r = Function(Q).interpolate(vel_r)
flux_r_t = assemble(inner(vel_r('-'),n_vec('-'))*dS(wellbore_toe))
flux_r_h = assemble(inner(vel_r('-'),n_vec('-'))*dS(wellbore_heel))
flux_r_c = assemble(inner(vel_r('-'),n_vec('-'))*dS(wellbore_cylinder))
flux_r_ff = assemble(inner(vel_r,n_vec)*ds)

VTKFile("vel_reservoir.pvd").write(_vel_r)
VTKFile("vel_wellbore.pvd").write(_vel_w)

print('wellbore')
print('pwc='+str(pwc))
print('flux_w (toe)='+str(flux_w_t))
print('flux_w (heel)='+str(flux_w_h))
print('flux_w (cylinder)='+str(flux_w_c))

print('reservoir')
print('prc='+str(prc))
print('flux_r (toe)='+str(flux_r_t))
print('flux_r (heel)='+str(flux_r_h))
print('flux_r (cylinder)='+str(flux_r_c))
print('flux_r (farfiled)='+str(flux_r_ff))
