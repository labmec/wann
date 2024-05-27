from firedrake import *
import ufl
import numpy as np
import math
import time as tm
import sys

meshfile = "../geo/mesh3D_rev02.msh"
mesh = Mesh(meshfile)
V = FunctionSpace(mesh, "CG",1)
W = VectorFunctionSpace(mesh, "CG",1)

u = TrialFunction(V)
v = TestFunction(V)
p = Function(V)

F = inner(grad(u),grad(v))*dx + Constant(0)*v*dx

Pr=10;
Pw=1;
p_farfield = DirichletBC(V, Pr, [100])
p_well = DirichletBC(V, Pw, [103])

_lhs = lhs(F)
_rhs = rhs(F)
problem = LinearVariationalProblem(_lhs, _rhs, p, bcs=[p_farfield,p_well])
solver  = LinearVariationalSolver(problem)
solver.solve()

vel = Function(W).interpolate(-grad(p))

File("p.pvd").write(p)
File("vel.pvd").write(vel)
