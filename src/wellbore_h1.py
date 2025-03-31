from firedrake import *
import ufl
import numpy as np
import math
import time as tm
import sys

rho   = 800.
mu    = 0.005
k     = 1e-13 # m2 = 101.3249970000 Milli Darcy
Q_toe = 0
p_heel= 1.177e+7
p_res = 2.206e+7
D     = 0.1
Dr    = 2000
Lw    = 400
Nelem = 100
dpdx  = 0.000001e+7  
p_toe_expected  = 1.219297411e+7 #from Mathematica NDsolve

# K from a vertical wellbore (constant)
K = 2*np.pi*k/(mu*np.log(Dr/D))

mesh = IntervalMesh(Nelem,0,Lw)
V = FunctionSpace(mesh, "CG",1)

u = TrialFunction(V)
v = TestFunction(V)
p = Function(V) #solution
u_n = Function(V) #solution at i-1

#initial condition
#x = SpatialCoordinate(mesh)
#p.interpolate(p_heel+dpdx*x)

x = np.linspace(0,Lw,Nelem+1)
p.dat.data[:] = p_heel + x[:] * dpdx 
u_n.dat.data[:] = p_heel + x[:] * dpdx 

#File("p_wellbore_h1.pvd").write(p)
#sys.exit("exit")


def Re_mask(p_nm1):
    # return a mask based on the Reynolds number
    _dpdx = Dx(p_nm1,0)
    
    # first, compute Reynolds number considering linear flow everywhere
    _Re = (1/32) * rho*( (D**3)/(mu**2) ) * _dpdx 
    
    _mask = ufl.conditional( _Re <= 1187.38, Constant(0), Constant(1))

    return _mask

def C(p_nm1):
    # retrieve the mask that defines linear and non-linear flow
    _mask = Re_mask(p_nm1)

    ## linear case
    _C_linear = np.pi*(D**4)/(128*mu)

    ##non-linear case
    _grad_p_abs = Dx(p_nm1,0)
    _aux = ( 2.252610888*D**(19/7) )/( (mu**(1/7)) * (rho**(3/7)) )
    _C_nonlinear = _aux * _grad_p_abs**(-3/7)
   
    return (1-_mask)*_C_linear + _mask*_C_nonlinear


#F = C(u_n)*inner(grad(u),grad(v))*dx + K*(u-p_res)*v*dx 
F = C(u_n)*inner(grad(u),grad(v))*dx + K*(u_n-p_res)*v*dx 

p_well = DirichletBC(V, p_heel, [1])

_lhs = lhs(F)
_rhs = rhs(F)
problem = LinearVariationalProblem(_lhs, _rhs, p, bcs=[p_well])
solver  = LinearVariationalSolver(problem)

p_file = VTKFile("p_wellbore_h1.pvd")
p_file.write(p,time=0)

num_it=50
ti = tm.time()
for i in range(num_it):
    print("iteration "+str(i+1), end=" ")
    solver.solve()
    error = norm(u_n-p)
    print(" error = "+str(error))
    
    u_n.assign(p)
    p_file.write(p,time=i+1)

tf = tm.time()
print("Exec. time="+str(tf-ti))


print("P heel="+str(p.dat.data[ 0]))
print("P toe ="+str(p.dat.data[-1]))
print("dP toe (Pa)="+str(abs(p.dat.data[-1]-p_toe_expected)))
print("dP toe (%) ="+str(100*abs(p.dat.data[-1]-p_toe_expected)/p_toe_expected))

#solve(F==0,u, bcs=p_well)


#vel = Function(W).interpolate(-grad(p))

#File("p_wellbore_h1.pvd").write(p)
#File("vel.pvd").write(vel)
