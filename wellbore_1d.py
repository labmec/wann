import numpy as np
import matplotlib.pyplot as plt
import sys

# model settings {{{
model = {}

model["fluid_prop"] = {
    "rho": 800., # kg/m3
    "mu": 0.005, # Pa * s
}

model["wellbore_prop"] = {
    "Lw": 400, # m
    "D": 0.1, # m
    #"D": 1., # m
    "Q_toe": 0,
    "p_heel": 1.177e+7 # Pa == 120 kgf/cm2
}

model["reservoir_prop"] = {
    "p_res": 2.206e+7, # Pa == 225 kgf/cm2
    "k": 1e-13, # m2 
    "Dr": 2000, # [m] diameter of the reservoir
}

model["RK_settigns"] = {
    "verbose": 1, 
    "npoints": 21,
    "idp": 0.005, # [Pa] initial delta P
    "eps": 1.e-6,
    "Re_min": 100,
    "rtol": 1e-10,
    "max_n": 100
}
#}}}
def Reynolds(model, Q): #{{{
    rho = model["fluid_prop"]["rho"]
    mu = model["fluid_prop"]["mu"]
    D = model["wellbore_prop"]["D"]
    V = np.sqrt(4*Q/np.pi)
    return rho * V * D / mu
#}}}
def FrictionFactor(model, Q):#{{{
    eps = model["RK_settigns"]["eps"]
    Re_min = model["RK_settigns"]["Re_min"]
    Re = Reynolds(model, Q)
    if Re < eps:
        print("WARNING: Reynolds <= 0! Using mininal Reynolds")
        Re = Re_min        

    if Re<= 1187.38:
        f = 16/Re # laminar flow
    else:
        f= 0.0791/Re**0.25 # turbulent flow
    
    return f
#}}}
def K(model, x, K_x=None, K_v=None): #{{{
    # K comes from a Trained Neural Network 
    #FIXME testing
    if 1:
        k = model["reservoir_prop"]["k"]
        Dr= model["reservoir_prop"]["Dr"]
        D = model["wellbore_prop"]["D"]
        mu= model["fluid_prop"]["mu"] 

        # K from a vertical wellbore (constant)
        K = 2*np.pi*k/(mu*np.log(Dr/D))
    else:
        np.interp(x,K_x,K_v)

    return K
#}}}
def Qfunc(model, x, p): #{{{
    p_res = model["reservoir_prop"]["p_res"]
    
    return -K(model, x) * (p - p_res)
#}}}
def pfunc(model, x, Q): #{{{
    rho = model["fluid_prop"]["rho"]
    D = model["wellbore_prop"]["D"]
    f = FrictionFactor(model, Q)
    # if Q=0, f is obtained with Re=100 (min_Re); 
    # it will return 0 anyway
    return -32*f*rho*Q**2/(np.pi**2 * D**5)
    
    #FIXME testing laminar flow for now
    #mu = model["fluid_prop"]["mu"]
    #return -128*mu*Q/(np.pi*D**4)
#}}}
def read_K(model):#{{{
    _data = np.loadtxt("./pytorch/kpts.txt", delimiter="\t",unpack=False)
    _x = _data[:,0]*model["wellbore_prop"]["Lw"]
    _K = _data[:,1] 

    #plt.plot(x_K,K,linestyle='none',marker='o')
    #plt.show()
    return _x, _K
#}}}
def Q_vertical_wellbore(model):#{{{
    k = model["reservoir_prop"]["k"]
    h = model["wellbore_prop"]["Lw"]
    p_res = model["reservoir_prop"]["p_res"]
    p_heel = model["wellbore_prop"]["p_heel"]
    mu = model["fluid_prop"]["mu"]
    D = model["wellbore_prop"]["D"]
    Dr = model["reservoir_prop"]["Dr"]

    Qw = 2*np.pi*k*h*(p_res-p_heel)/(mu*np.log(Dr/D))
    
    return Qw
#}}}
def test_VerticalWelboreModel(model):#{{{
    model["reservoir_prop"]["k"] = 1e-13 # m2 
    model["wellbore_prop"]["Lw"] = 400 # m 
    model["reservoir_prop"]["p_res"] = 2.206e+7 # Pa 
    model["wellbore_prop"]["p_heel"] = 1.177e+7 # Pa
    model["fluid_prop"]["mu"] = 0.005 # Pa * s
    model["wellbore_prop"]["D"] = 0.1 # m
    model["reservoir_prop"]["Dr"] = 2000 # m 
    
    Q_total = Q_vertical_wellbore(model)
    return Q_total

#}}}

verbose = model["RK_settigns"]["verbose"]
dp = model["RK_settigns"]["idp"] # initial delta p in the toe
npoints = model["RK_settigns"]["npoints"]
Lw = model["wellbore_prop"]["Lw"]
p_res = model["reservoir_prop"]["p_res"]
p_heel= model["wellbore_prop"]["p_heel"]
Q_toe = model["wellbore_prop"]["Q_toe"]
rtol = model["RK_settigns"]["rtol"]
max_n = model["RK_settigns"]["max_n"]
p_toe = p_heel + dp 
#p_toe = p_res # it seems a good initial estimate
rel_error = 1

x = np.linspace(0,Lw,npoints) # x=0, toe; x=Lw, heel
Q = np.zeros(len(x))
p = np.zeros(len(x))

Q[0] = Q_toe
p[0] = p_toe
n = 0

# read K from file
#K_x, K_v = read_K(model)
#K_x = Lw-K_x 
#print(K_x[0])
#plt.plot(K_x,K_v,linestyle='none',marker='o')
#plt.show()

# for the RK solver, x=0 is the toe; x=L is the heel
while rel_error > rtol and n < max_n:
    print("Solving iteration n="+str(n))
    
    if n>0:
        p_toe = p_toe + dp
        p[0] = p_toe
    
    for i in range(0,len(x)-1):
        if verbose>1:
            print("\t solving step i="+str(i))
            print("\tQ[i]="+str(Q[i])+", p[i]="+str(p[i]) )

        h = x[i+1]-x[i]
        
        Q1 = Qfunc(model,  x[i],     p[i] )
        p1 = pfunc(model,  x[i],     Q[i], )
        if verbose>2:
            print("\t\tQ1="+str(Q1)+", p1="+str(p1) )

        Q2 = Qfunc(model,  x[i]+h/2, p[i]+p1*h/2 )
        p2 = pfunc(model,  x[i]+h/2, Q[i]+Q1*h/2 )
        if verbose>2:
            print("\t\tQ2="+str(Q2)+", p2="+str(p2) )
        
        Q3 = Qfunc(model,  x[i]+h/2, p[i]+p2*h/2 )
        p3 = pfunc(model,  x[i]+h/2, Q[i]+Q2*h/2 )
        if verbose>2:
            print("\t\tQ3="+str(Q3)+", p3="+str(p3) )
        
        Q4 = Qfunc(model,  x[i]+h,   p[i]+p3*h )
        p4 = pfunc(model,  x[i]+h,   Q[i]+Q3*h )
        if verbose>2:
            print("\t\tQ4="+str(Q4)+", p4="+str(p4) )

        Q[i+1] = Q[i] + (h/6)*(Q1 + 2*Q2 + 2*Q3 + Q4)
        p[i+1] = p[i] + (h/6)*(p1 + 2*p2 + 2*p3 + p4)

    dp = p_heel-p[-1] 
    rel_error = abs((dp)/p_heel)
    print("Done iteration n="+str(n))
    print("rel_error="+str(rel_error))
    if verbose>0:
        print("BC\t\t\tModel solution")
        print("p_res="+str(p_res)+",\tp[toe]="+str(p[0]))
        print("p_heel="+str(p_heel)+",\tp[heel]="+str(p[-1]))
        print("Q_toe="+str(Q_toe)+",\t\tQ[toe]="+str(Q[0]))
        print("Q_heel=?\t\tQ[heel]="+str(Q[-1]))
    n = n+1

Q_vertical_wellbore = test_VerticalWelboreModel(model)
print("Q_vertical_verbore="+str(Q_vertical_wellbore))

# post-processing
# transforming RK x to FEM x
x_fem = Lw-x  

MPa = 1.e6
fig, axs = plt.subplots(2, 1, layout='constrained')
axs[0].plot(x_fem, p/MPa)
#axs[0].set_xlim(0, 2)
axs[0].set_xlabel('Wellbore distance (m)')
axs[0].set_ylabel('Pressure (MPa)')
axs[0].grid(True)

axs[1].plot(x_fem, Q)
axs[1].set_xlabel('Wellbore distance (m)')
axs[1].set_ylabel('Flow rate (m2/s)')
axs[1].grid(True)

plt.show()
