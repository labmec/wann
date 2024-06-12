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
V = FunctionSpace(mesh, "CG",1)
#W = VectorFunctionSpace(mesh, "CG",2)
W = V*V #wellbore and reservoir CG spaces

#p_order = 1
#V = FunctionSpace(mesh, "DG",p_order-1)
#W = FunctionSpace(mesh, "RTCF",p_order)
#M = W*V

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

File("I_r.pvd").write(I_r)
File("I_w.pvd").write(I_w)
File("I_cg_r.pvd").write(I_cg_r)
File("I_cg_w.pvd").write(I_cg_w)

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

k = 1.

vel_w = -k*grad(u_w) # vel inside wellbore -> +
vel_r = -grad(u_r) # vel inside reservoi -> -

F = inner(k*grad(u_w),grad(v_w))*dx(wellbore_id) + Constant(0)*v_w*dx(wellbore_id)
F+= inner(grad(u_r),grad(v_r))*dx(reservoir_id) + Constant(0)*v_r*dx(reservoir_id)
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
File("p_reservoir.pvd").write(p_r)
File("p_wellbore.pvd").write(p_w)

pwc = assemble(p_w('-')*dS(wellbore_cylinder))
prc = assemble(p_r('-')*dS(wellbore_cylinder))

vel_w = -k*grad(p_w)
#flux_w = assemble(inner(vel_w('+'),n_vec('+'))*dS(wellbore_cylinder))
flux_w = assemble(inner(vel_w('+'),n_vec('+'))*dS(wellbore_toe))

vel_r = -grad(p_r)
#flux_r = assemble(inner(vel_r('-'),n_vec('-'))*dS(wellbore_cylinder))
flux_r = assemble(inner(vel_r('-'),n_vec('-'))*dS(wellbore_toe))

print('wellbore')
print('pwc='+str(pwc))
print('flux_w='+str(flux_w))

print('reservoir')
print('prc='+str(prc))
print('flux_r='+str(flux_r))


sys.exit("exit")

u_w = TrialFunction(V)
u_r = TrialFunction(V)
v_w = TestFunction(V)
v_r = TestFunction(V)
p_w = Function(V)
p_r = Function(V)

max_iter = 1
k = 1
for i in range(max_iter):
    vel_r = -grad(p_r)
    flux_w = inner(vel_r('-'),n_vec('-'))
    #flux_w = Constant(0.)
    F_w = inner(k*grad(u_w),grad(v_w))*dx(wellbore_id) + Constant(0)*v_w*dx(wellbore_id) - flux_w*v_w('-')*dS(wellbore_cylinder) #Constant(0)*v_w('-')*dS(wellbore_cylinder)
    p_well = DirichletBC(V, Pw, [wellbore_heel])
    exclude_beyond_wellbore = MyBC(V, 0, I_cg_w)
    _lhs_w = lhs(F_w)
    _rhs_w = rhs(F_w)
    problem_w = LinearVariationalProblem(_lhs_w, _rhs_w, p_w, bcs=[p_well, exclude_beyond_wellbore])
    #problem_w = LinearVariationalProblem(_lhs_w, _rhs_w, p_w, bcs=[p_well])
    solver_w  = LinearVariationalSolver(problem_w,solver_parameters=solver_cg_wipc)
    solver_w.solve()

    vel_w = -k*grad(p_w)
    #vel_r = -grad(p_r)
    flux_r = inner(vel_w('+'),n_vec('+'))
    #flux_r = Constant(0.)
    F_r = inner(grad(u_r),grad(v_r))*dx(reservoir_id) + Constant(0)*v_r*dx(reservoir_id) + flux_r*v_r('+')*dS(wellbore_cylinder) #Constant(0)*v_r('+')*dS(wellbore_cylinder) 
    p_farfield = DirichletBC(V, Pr, [reservoir_farfield])
    exclude_beyond_reservoir = MyBC(V, 0, I_cg_r)
    _lhs_r = lhs(F_r)
    _rhs_r = rhs(F_r)
    problem_r = LinearVariationalProblem(_lhs_r, _rhs_r, p_r, bcs=[p_farfield, exclude_beyond_reservoir])
    #problem_r = LinearVariationalProblem(_lhs_r, _rhs_r, p_r, bcs=[p_farfield])
    solver_r  = LinearVariationalSolver(problem_r,solver_parameters=solver_cg_wipc)
    solver_r.solve()

    #print("flux_w="+str(assemble(flux_w*dS(wellbore_cylinder))))
    #print("flux_r="+str(assemble(flux_r*dS(wellbore_cylinder))))
    #dflux = assemble(abs(flux_w-flux_r)*dS(wellbore_cylinder))
    #print(dflux)

#vel = Function(W).interpolate(-grad(p))

File("p_reservoir.pvd").write(p_r)
File("p_wellbore.pvd").write(p_w)
#File("vel.pvd").write(vel)
