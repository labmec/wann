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

V = FunctionSpace(mesh, "CG",2)
W = VectorFunctionSpace(mesh, "CG",2)

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

u_w = TrialFunction(V)
u_r = TrialFunction(V)
v_w = TestFunction(V)
v_r = TestFunction(V)
p_w = Function(V)
p_r = Function(V)

#F = inner(grad(u),grad(v))*dx(reservoir_id) + Constant(0)*v*dx(reservoir_id) + inner(grad(u),n_vec)*v('+')*dS(103) + 10000*inner(grad(u),grad(v))*dx(wellbore_id) + Constant(0)*v*dx(wellbore_id)
#F = inner(grad(u),grad(v))*dx(reservoir_id) + Constant(0)*v*dx(reservoir_id) + inner(grad(u('+')),n_vec('+'))*v('+')*dS(103) 
#F += 1*inner(grad(u),grad(v))*dx(wellbore_id) + Constant(0)*v*dx(wellbore_id) + inner(grad(u('-')),n_vec('-'))*v('-')*dS(103)

#F = inner(grad(u),grad(v))*dx(reservoir_id) + Constant(0)*v*dx(reservoir_id) + Constant(50)*v('+')*dS(101) + inner(grad(u),grad(v))*dx(wellbore_id) + Constant(0)*v*dx(wellbore_id)
#F = inner(grad(u),grad(v))*dx(reservoir_id) + Constant(0)*v*dx(reservoir_id)  

F_w = inner(100000*grad(u_w),grad(v_w))*dx + Constant(0)*v_w*dx + Constant(0)*v_w('-')*dS(wellbore_cylinder)
p_well = DirichletBC(V, Pw, [wellbore_heel])
exclude_beyond_wellbore = MyBC(V, 0, I_cg_w)
_lhs_w = lhs(F_w)
_rhs_w = rhs(F_w)
problem_w = LinearVariationalProblem(_lhs_w, _rhs_w, p_w, bcs=[p_well, exclude_beyond_wellbore])
solver_w  = LinearVariationalSolver(problem_w,solver_parameters=solver_cg_wipc)
solver_w.solve()

F_r = inner(grad(u_r),grad(v_r))*dx + Constant(0)*v_r*dx + Constant(0)*v_r('+')*dS(wellbore_cylinder) 
p_farfield = DirichletBC(V, Pr, [reservoir_farfield])
exclude_beyond_reservoir = MyBC(V, 0, I_cg_r)
_lhs_r = lhs(F_r)
_rhs_r = rhs(F_r)
problem_r = LinearVariationalProblem(_lhs_r, _rhs_r, p_r, bcs=[p_farfield, exclude_beyond_reservoir])
solver_r  = LinearVariationalSolver(problem_r,solver_parameters=solver_cg_wipc)
solver_r.solve()

#vel = Function(W).interpolate(-grad(p))

File("p_reservoir.pvd").write(p_r)
File("p_wellbore.pvd").write(p_w)
#File("vel.pvd").write(vel)
