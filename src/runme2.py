from firedrake import *
import ufl
import numpy as np
import math
import time as tm
import sys

solver_cg_wipc = {"ksp_type": "cg", "pc_type": "bjacobi", "pc_sub_type": "ilu"}

reservoir_id = 111
wellbore_id = 112
reservoir_farfield = 100
wellbore_heel = 101
wellbore_toe = 102
wellbore_cylinder = 103

meshfile = "../geo/mesh3D_rev03.msh"
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

#File("I_r.pvd").write(I_r)
#File("I_w.pvd").write(I_w)
#File("I_cg_r.pvd").write(I_cg_r)
#File("I_cg_w.pvd").write(I_cg_w)

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



u = TrialFunction(V)
v = TestFunction(V)
p = Function(V)



#F = inner(grad(u),grad(v))*dx(reservoir_id) + Constant(0)*v*dx(reservoir_id) + inner(grad(u),n_vec)*v('+')*dS(103) + 10000*inner(grad(u),grad(v))*dx(wellbore_id) + Constant(0)*v*dx(wellbore_id)
F = inner(grad(u),grad(v))*dx(reservoir_id) + Constant(0)*v*dx(reservoir_id) + inner(grad(u('+')),n_vec('+'))*v('+')*dS(102) 
F += 1*inner(grad(u),grad(v))*dx(wellbore_id) + Constant(0)*v*dx(wellbore_id) + inner(grad(u('-')),n_vec('-'))*v('-')*dS(102)

Pr=10;
Pw=1;
p_farfield = DirichletBC(V, Pr, [reservoir_farfield])
p_well = DirichletBC(V, Pw, [wellbore_heel])

_lhs = lhs(F)
_rhs = rhs(F)
problem = LinearVariationalProblem(_lhs, _rhs, p, bcs=[p_farfield,p_well])
solver  = LinearVariationalSolver(problem,solver_parameters=solver_cg_wipc)
solver.solve()

vel = Function(W).interpolate(-grad(p))

File("p.pvd").write(p)
File("vel.pvd").write(vel)
