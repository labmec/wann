import numpy as np
import matplotlib.pyplot as plt
import sys

# model settings {{{
def model_settings():

    model = {}

    model["fluid_prop"] = {
        "rho": 800., # kg/m3
        "mu": 0.005, # Pa * s
    }

    model["wellbore_prop"] = {
        "Lw": 400, # m
        "D": 0.1, # m
        "Q_toe": 0,
        "p_heel": 1.177e+7, # Pa == 120 kgf/cm2
        "f_type": 0, #0: default friction factor (linear or non-linear); 1: linear only; 2: non-linear only  
    }

    model["reservoir_prop"] = {
        "p_res": 2.206e+7, # Pa == 225 kgf/cm2
        "k": 1e-13, # m2 = 101.3249970000 Milli Darcy 
        "Dr": 2000, # [m] diameter of the reservoir
        "K_type": 0, # 0: K=cte from vertical wellbore; 1: K=linear func; 2: parabolic func; 3: K comes from an ANN  
    }

    model["RK_settigns"] = {
        "verbose": 1, 
        "npoints": 21,
        "idp": 0.005e+7, # [Pa] initial delta P
        "eps": 1.e-6,
        "Re_min": 1, #FIXME not sure if 100 is ok
        "rtol": 1e-10,
        "max_n": 100,
        "relax_param": 0.5
    }

    return model
#}}}
# model settings for numerical tests {{{
def model_settings_for_tests():
    # define model parameters to perform numerical tests
    model = model_settings()

    model["fluid_prop"]["rho"] = 800. # kg/m3
    model["fluid_prop"]["mu"]  = 0.005 # Pa * s

    model["wellbore_prop"]["Lw"] = 400 # m
    model["wellbore_prop"]["D"] = 0.1 # m
    model["wellbore_prop"]["Q_toe"] = 0
    model["wellbore_prop"]["p_heel"] = 1.177e+7 # Pa == 120 kgf/cm2
    model["wellbore_prop"]["f_type"] = 0 #0: default friction factor (linear or non-linear); 1: linear only; 2: non-linear only  

    model["reservoir_prop"]["p_res"] = 2.206e+7 # Pa == 225 kgf/cm2
    model["reservoir_prop"]["k"] = 1e-13 # m2 
    model["reservoir_prop"]["Dr"] = 2000 # [m] diameter of the reservoir
    model["reservoir_prop"]["K_type"] = 0 # 0: K=cte from vertical wellbore; 1: K=linear func; 2: parabolic func; 3: K comes from an ANN  

    model["RK_settigns"]["verbose"] = 1 
    model["RK_settigns"]["npoints"] = 21
    model["RK_settigns"]["idp"] = 0.005e+7 # [Pa] initial delta P
    model["RK_settigns"]["eps"] = 1.e-6
    model["RK_settigns"]["Re_min"] = 100
    model["RK_settigns"]["rtol"] = 1e-10
    model["RK_settigns"]["max_n"] = 100
    model["RK_settigns"]["relax_param"] = 0.5

    return model
#}}}
def Reynolds(model, Q): #{{{
    rho = model["fluid_prop"]["rho"]
    mu = model["fluid_prop"]["mu"]
    D = model["wellbore_prop"]["D"]
    V = 4*Q/(np.pi*D**2)
    return rho * V * D / mu
#}}}
def FrictionFactor(model, Re):#{{{

    f_laminar   = 16/Re
    f_turbulent = 0.0791/(Re**0.25) 

    if model["wellbore_prop"]["f_type"]==0: # default
        f = f_laminar if Re<= 1187.38 else f_turbulent 
    else: 
        f = f_laminar if model["wellbore_prop"]["f_type"]==1 else f_turbulent

    #if Re<= 1187.38:
    #    f = 16/Re # laminar flow
    #else:
    #    f= 0.0791/(Re**0.25) # turbulent flow
    #f = 16/Re # laminar flow
    #f= 0.0791/(Re**0.25) # turbulent flow
    return f
#}}}
def K(model, x, K_x=None, K_v=None): #{{{
   
    # x=0, toe; x=Lw, heel
    K_type = model["reservoir_prop"]["K_type"]
    verbose = model["RK_settigns"]["verbose"]  

    if K_type==0:
        # K from a vertical wellbore (constant)
        if verbose>1:
            print("K (constant)")
        k = model["reservoir_prop"]["k"]
        Dr= model["reservoir_prop"]["Dr"]
        D = model["wellbore_prop"]["D"]
        mu= model["fluid_prop"]["mu"] 
        K = 2*np.pi*k/(mu*np.log(Dr/D))
    elif K_type==1:
        # K linear (based on a vertical wellbore) 
        if verbose>1:
            print("K (linear)")
        k = model["reservoir_prop"]["k"]
        Dr= model["reservoir_prop"]["Dr"]
        D = model["wellbore_prop"]["D"]
        Lw= model["wellbore_prop"]["Lw"]
        mu= model["fluid_prop"]["mu"]
        Kvw = 2*np.pi*k/(mu*np.log(Dr/D))
        K = Kvw * x/Lw     
    elif K_type==2:
        # K parabolic (based on a vertical wellbore) 
        if verbose>1:
            print("K (parabolic)")
        k = model["reservoir_prop"]["k"]
        Dr= model["reservoir_prop"]["Dr"]
        D = model["wellbore_prop"]["D"]
        Lw= model["wellbore_prop"]["Lw"]
        mu= model["fluid_prop"]["mu"]
        Kvw = 2*np.pi*k/(mu*np.log(Dr/D))
        K = 4*Kvw * (x/Lw)**2 - 4*Kvw * (x/Lw) + 3*Kvw/2     
    elif K_type==3:
        # K pn  (based on a vertical wellbore) 
        pn=16 # used in the tests
        if verbose>1:
            print("K with p="+str(12) + "(singular function)")
        k = model["reservoir_prop"]["k"]
        Dr= model["reservoir_prop"]["Dr"]
        D = model["wellbore_prop"]["D"]
        Lw= model["wellbore_prop"]["Lw"]
        mu= model["fluid_prop"]["mu"]
        Kvw = 2*np.pi*k/(mu*np.log(Dr/D))
        K = Kvw*(2*x/Lw - 1)**pn + (1/2)*Kvw 
    elif K_type==4:
        sys.exit("K not defined yet")
    else:
        sys.exit("K not defined yet")

    return K
#}}}
def Qfunc(model, x, p): #{{{
    p_res = model["reservoir_prop"]["p_res"]
    
    return -K(model, x) * (p - p_res)
#}}}
def pfunc(model, x, Q): #{{{
    rho = model["fluid_prop"]["rho"]
    D = model["wellbore_prop"]["D"]
    eps = model["RK_settigns"]["eps"]
    Re_min = model["RK_settigns"]["Re_min"]
    verbose = model["RK_settigns"]["verbose"]

    Re = Reynolds(model, Q) 
    if Re < eps:
        if verbose>1:
            print("\tWARNING: Reynolds <= 0! Using mininal Reynolds")
        Re = Re_min        
    f = FrictionFactor(model, Re)
    # if Q=0, f is obtained with Re=100 (min_Re); 
    # it will return 0 anyway
    return -32*f*rho*Q**2/(np.pi**2 * D**5)
    
#}}}
def read_K(model):#{{{
    _data = np.loadtxt("./pytorch/kpts.txt", delimiter="\t",unpack=False)
    _x = _data[:,0]*model["wellbore_prop"]["Lw"]
    _K = _data[:,1] 

    #plt.plot(x_K,K,linestyle='none',marker='o')
    #plt.show()
    return _x, _K
#}}}
def RungeKutta_solver(model):#{{{

    print("Starting Runge-Kutta solver")

    verbose = model["RK_settigns"]["verbose"]
    dp = model["RK_settigns"]["idp"] # initial delta p in the toe
    npoints = model["RK_settigns"]["npoints"]
    Lw = model["wellbore_prop"]["Lw"]
    p_res = model["reservoir_prop"]["p_res"]
    p_heel= model["wellbore_prop"]["p_heel"]
    Q_toe = model["wellbore_prop"]["Q_toe"]
    rtol = model["RK_settigns"]["rtol"]
    max_n = model["RK_settigns"]["max_n"]
    relax_param = model["RK_settigns"]["relax_param"]
    p_toe = p_heel + dp  
    p_rel_error = 1

    # defining x, Q, and p
    x = np.linspace(0,Lw,npoints) # x=0, toe; x=Lw, heel
    Q = np.zeros(len(x))
    p = np.zeros(len(x))

    Q[0] = Q_toe
    p[0] = p_toe
    dp_old = dp
    n = 0

    # for the RK solver, x=0 is the toe; x=L is the heel
    while p_rel_error > rtol and n < max_n:
        print("Solving iteration n="+str(n))
        
        if n>0:
            p_toe += dp 
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

        dp_old = dp
        dp = p_heel-p[-1] # difference between the imposed pressure (BC) p_heel and the estimated by the rK solver (p[-1])
        p_rel_error = abs((dp)/p_heel)
      
        # apply a relaxation and/or a change in dp to avoid error increasing
        dp*=relax_param 
        if abs(dp)>abs(dp_old):
            dp /= 2 # decrease dp (divide it by 2)

        print("Done iteration n="+str(n))
        print("P heel relative error (%) = "+str(p_rel_error*100))
        if verbose>0:
            print("dp="+str(dp)) 
            print("dp old="+str(dp_old))

        if verbose>0:
            print("\tSolution for iteration n="+str(n)+":")
            print("\tBC\t\t\tModel solution")
            print("\tp_res="+str(p_res)+",\tp[toe]="+str(p[0]))
            print("\tp_heel="+str(p_heel)+",\tp[heel]="+str(p[-1]))
            print("\tQ_toe="+str(Q_toe)+",\t\tQ[toe]="+str(Q[0]))
            print("\tQ_heel=?\t\tQ[heel]="+str(Q[-1]))
        
        n = n+1
        
    if n>=max_n:
        sys.exit("ERROR: maximum number of iterations exceeded ("+str(max_n)+")! Check the input parameters.")

    return x, Q, p
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
def test_VerticalWelboreModel():#{{{
    model = model_settings_for_tests()

    model["reservoir_prop"]["k"] = 1e-13 # m2 
    model["reservoir_prop"]["p_res"] = 2.206e+7 # Pa 
    model["reservoir_prop"]["Dr"] = 2000 # m 
    model["reservoir_prop"]["K_type"] = 0
    model["wellbore_prop"]["D"] = 1.0 # m to avoid pressure drop
    model["wellbore_prop"]["Lw"] = 400 # m 
    model["wellbore_prop"]["p_heel"] = 1.177e+7 # Pa
    model["fluid_prop"]["mu"] = 0.005 # Pa * s

    Q_VW = Q_vertical_wellbore(model)
    
    x_RK, Q_RK, p_RK = RungeKutta_solver(model) 
   
    assert np.allclose(Q_VW, Q_RK[-1], rtol= 8.0e-7)
    
    fig, ax = plt.subplots()
    ax.plot(400-x_RK,Q_RK) 
    ax.plot(400-x_RK[-1],Q_VW,"rx",markersize=12)
    ax.set_xlabel("Wellbore distance")
    ax.set_ylabel("Q [m3/s]")
    ax.grid(True)
    plt.show(block=False)

#}}}
def test_FrictionFator():#{{{

    model = model_settings_for_tests()
    model["wellbore_prop"]["f_type"] = 0 # default

    Re_eq = 1187.38
    _eps  = 0.001 

    assert np.allclose(FrictionFactor(model, Re_eq-_eps), 16/(Re_eq-_eps), rtol=1.e-14)
    
    assert np.allclose(FrictionFactor(model, Re_eq+_eps), 0.0791/((Re_eq+_eps)**0.25), rtol=1.e-14)

    assert np.allclose(FrictionFactor(model, Re_eq-_eps), FrictionFactor(model, Re_eq+_eps), rtol=4.e-6)

    Re_min = 100 
    Re_max = 10000
    npoints = 1000

    Re = np.linspace(Re_min,Re_max,npoints)
    f = np.zeros(len(Re))
    g = np.zeros(len(Re))

    for i in range(0,len(Re)):
        f[i] = FrictionFactor(model, Re[i])
        if Re[i] <= Re_eq:
            g[i] = 16/Re[i]
        else:
            g[i] = 0.0791/(Re[i]**0.25)
   
    assert np.allclose(f,g,rtol=1e-12)

    fig, ax = plt.subplots()
    ax.plot(Re,f) 
    plt.xscale("log")
    ax.set_xlabel("Reynolds")
    ax.set_ylabel("Friction factor")
    ax.grid(True)
    plt.show(block=False)

#}}}
def test_Reynolds():#{{{
    D = 0.1 # m
    mu = 0.001 # Pa*s
    rho = 1000 # kg/m3
    V = 0.01 # m/s

    Q = V*np.pi*(D**2)/4
    
    model = model_settings_for_tests()

    model["wellbore_prop"]["D"] = D # m
    model["fluid_prop"]["mu"] = mu # Pa * s
    model["fluid_prop"]["rho"] = rho # kg / m3

    Re_1 = rho*V*D/mu
    Re_2 = Reynolds(model,Q)

    assert np.allclose(Re_1,Re_2,rtol=1e-12)
    
#}}}
def test_RungeKuttaLinearFrictionConstantK():#{{{
    # from Mathematica (using NDSolve)
    # K constant (vertical wellbore), linear friction coeficient
    p_toe_expected  = 1.179124269e+7
    p_heel_expected = 1.177e+7  
    Q_heel_expected = 0.05215536594

    model = model_settings_for_tests()
    model["reservoir_prop"]["K_type"] = 0 # K constant (vertical wellbore)
    model["wellbore_prop"]["f_type"] = 1 # linear only

    # run the RungeKutta model
    x_RK, Q_RK, p_RK = RungeKutta_solver(model)

    #print(p_RK[0])
    #print(p_RK[-1])
    #print(Q_RK[-1])

    assert np.allclose(p_toe_expected,  p_RK[0],  rtol= 1.0e-9)
    assert np.allclose(p_heel_expected, p_RK[-1], rtol= 1.0e-9)
    assert np.allclose(Q_heel_expected, Q_RK[-1], rtol= 2.0e-7, atol=1.0e-12)

#}}}
def test_RungeKuttaLinearFrictionLinearK():#{{{
    # from Mathematica (using NDSolve)
    # K linear, linear friction coeficient
    p_toe_expected  = 1.177708919e+7
    p_heel_expected = 1.177000000e+7  
    Q_heel_expected = 0.02610282486

    model = model_settings_for_tests()
    model["reservoir_prop"]["K_type"] = 1 # K linear 
    model["wellbore_prop"]["f_type"] = 1 # linear only

    # run the RungeKutta model
    x_RK, Q_RK, p_RK = RungeKutta_solver(model)

    print(p_RK[0])
    print(p_RK[-1])
    print(Q_RK[-1])

    assert np.allclose(p_toe_expected,  p_RK[0],  rtol= 2.0e-10)
    assert np.allclose(p_heel_expected, p_RK[-1], rtol= 1.0e-10)
    assert np.allclose(Q_heel_expected, Q_RK[-1], rtol= 7.0e-10, atol=1.0e-12)

#}}}
def test_RungeKuttaLinearFrictionParabolicK():#{{{
    # from Mathematica (using NDSolve)
    # K parabolic, linear friction coeficient
    p_toe_expected  = 1.178770713e+7
    p_heel_expected = 1.177000000e+7  
    Q_heel_expected = 0.04347630488

    model = model_settings_for_tests()
    model["reservoir_prop"]["K_type"] = 2 # K parabolic 
    model["wellbore_prop"]["f_type"] = 1 # linear only

    # run the RungeKutta model
    x_RK, Q_RK, p_RK = RungeKutta_solver(model)

    print(p_RK[0])
    print(p_RK[-1])
    print(Q_RK[-1])

    assert np.allclose(p_toe_expected,  p_RK[0],  rtol= 3.0e-10)
    assert np.allclose(p_heel_expected, p_RK[-1], rtol= 1.0e-10)
    assert np.allclose(Q_heel_expected, Q_RK[-1], rtol= 3.0e-9, atol=1.0e-12)

#}}}
def test_RungeKuttaLinearFrictionSingularK():#{{{
    # from Mathematica (using NDSolve)
    # K singular, linear friction coeficient
    p_toe_expected  = 1.178187991e+7
    p_heel_expected = 1.177000000e+7  
    Q_heel_expected = 0.02916429353

    model = model_settings_for_tests()
    model["reservoir_prop"]["K_type"] = 3 # K singular
    model["wellbore_prop"]["f_type"] = 1 # linear only
    model["RK_settigns"]["npoints"] = 301 # because K is singular here, we need more points

    # run the RungeKutta model
    x_RK, Q_RK, p_RK = RungeKutta_solver(model)

    #print(p_RK[0])
    #print(p_RK[-1])
    #print(Q_RK[-1]-Q_heel_expected)

    assert np.allclose(p_toe_expected,  p_RK[0],  rtol= 3.e-11)
    assert np.allclose(p_heel_expected, p_RK[-1], rtol= 1.e-10)
    assert np.allclose(Q_heel_expected, Q_RK[-1], rtol= 3.0e-9, atol=1.2e-10)

#}}}
def test_RungeKuttaNonLinearFrictionConstantK():#{{{
    # from Mathematica (using NDSolve)
    # K constant (vertical wellbore), non linear friction coeficient
    p_toe_expected  = 1.219297411e+7
    p_heel_expected = 1.176999979e+7  
    Q_heel_expected = 0.05065050899

    model = model_settings_for_tests()
    model["reservoir_prop"]["K_type"] = 0 # K constant (vertical wellbore)
    model["wellbore_prop"]["f_type"] = 2 # non linear only

    # run the RungeKutta model
    x_RK, Q_RK, p_RK = RungeKutta_solver(model)

    #print(p_RK[0])
    #print(p_RK[-1])
    #print(Q_RK[-1]-Q_heel_expected)

    assert np.allclose(p_toe_expected,  p_RK[0],  rtol= 4.0e-8)
    assert np.allclose(p_heel_expected, p_RK[-1], rtol= 2.0e-8)
    assert np.allclose(Q_heel_expected, Q_RK[-1], rtol= 3.0e-8, atol=1.0e-12)

#}}}
def test_RungeKuttaNonLinearFrictionLinearK():#{{{
    # from Mathematica (using NDSolve)
    # K linear, non linear friction coeficient
    p_toe_expected  = 1.185093713e+7
    p_heel_expected = 1.177000000e+7  
    Q_heel_expected = 0.02597136491

    model = model_settings_for_tests()
    model["reservoir_prop"]["K_type"] = 1 # K linear
    model["wellbore_prop"]["f_type"] = 2 # non linear only

    # run the RungeKutta model
    x_RK, Q_RK, p_RK = RungeKutta_solver(model)

    #print(p_RK[0])
    #print(p_RK[-1])
    #print(Q_RK[-1]-Q_heel_expected)

    assert np.allclose(p_toe_expected,  p_RK[0],  rtol= 1.1e-8)
    assert np.allclose(p_heel_expected, p_RK[-1], rtol= 1.0e-10)
    assert np.allclose(Q_heel_expected, Q_RK[-1], rtol= 3.0e-9, atol=1.0e-12)

#}}}
def test_RungeKuttaNonLinearFrictionParabolicK():#{{{
    # from Mathematica (using NDSolve)
    # K parabolic, non linear friction coeficient
    p_toe_expected  = 1.206718661e+7
    p_heel_expected = 1.177000000e+7  
    Q_heel_expected = 0.04267083913

    model = model_settings_for_tests()
    model["reservoir_prop"]["K_type"] = 2 # K parabolic
    model["wellbore_prop"]["f_type"] = 2 # non linear only

    # run the RungeKutta model
    x_RK, Q_RK, p_RK = RungeKutta_solver(model)

    #print(p_RK[0])
    #print(p_RK[-1])
    #print(Q_RK[-1]-Q_heel_expected)

    assert np.allclose(p_toe_expected,  p_RK[0],  rtol= 1.0e-7)
    assert np.allclose(p_heel_expected, p_RK[-1], rtol= 1.0e-10)
    assert np.allclose(Q_heel_expected, Q_RK[-1], rtol= 6.0e-8, atol=1.0e-12)

#}}}
def test_RungeKuttaNonLinearFrictionSingularK():#{{{
    # from Mathematica (using NDSolve)
    # K singular, non linear friction coeficient
    p_toe_expected  = 1.192414383e+7
    p_heel_expected = 1.177000002e+7  
    Q_heel_expected = 0.02887920864

    model = model_settings_for_tests()
    model["reservoir_prop"]["K_type"] = 3 # K singular
    model["wellbore_prop"]["f_type"] = 2 # non linear only
    model["RK_settigns"]["npoints"] = 301 # because K is singular here, we need more points

    # run the RungeKutta model
    x_RK, Q_RK, p_RK = RungeKutta_solver(model)

    #print(p_RK[0])
    #print(p_RK[-1])
    print(Q_RK[-1]-Q_heel_expected)

    assert np.allclose(p_toe_expected,  p_RK[0],  rtol= 2.e-10)
    assert np.allclose(p_heel_expected, p_RK[-1], rtol= 3.e-09)
    assert np.allclose(Q_heel_expected, Q_RK[-1], rtol= 3.0e-9, atol=1.5e-9)

#}}}
def test_h_convergence(): #{{{
    
    model = model_settings_for_tests()
    model["reservoir_prop"]["K_type"] = 3 # 2: parabolic func; 3: singular K 
    model["wellbore_prop"]["f_type"] = 0  # default, linear and non-linear

    Lw    =  model["wellbore_prop"]["Lw"] 
    nele  = 1 # initial number of elements
    ndiv  = 10 # number of simulations with increased mesh resolution

    Q_heel = []
    p_toe  = []
    h = []

    for i in range(ndiv):
        # run the RungeKutta model
        model["RK_settigns"]["npoints"] = nele+1 # nele -> npoints
        x_RK, Q_RK, p_RK = RungeKutta_solver(model)
        
        Q_heel.append( Q_RK[-1] )
        p_toe.append( p_RK[0] )
        h.append( Lw/nele ) 
    
        nele = 2*nele
    
    print("###################################################")
    print("###################################################")
    print("h convergence test - K parabolic function, f linear/nonlinear")
    print("h\t\tQ_heel\t\t\tp_toe")
    for i in range(ndiv-1):
        print(str(h[i]) +"\t"+ str(abs(Q_heel[i]-Q_heel[-1])) + "\t"+ str(abs(p_toe[i]-p_toe[-1])), flush=True )
    print("###################################################")
    print("###################################################")

    e_i = abs(Q_heel[0]-Q_heel[-1])
    e_rk = []
    e_rk.append(e_i)
    for i in range(1,ndiv):
        e_rk.append(e_rk[i-1]/(2**4))

    print(e_rk)

    fig, ax = plt.subplots()
    ax.plot(h[0:-1:],abs(Q_heel[0:-1:]-Q_heel[-1]))
    ax.plot(h[0:-1:],e_rk[0:-1:])
    ax.set_xlabel("h [m]")
    ax.set_ylabel("error")
    ax.set_yscale('log')
    ax.set_xscale('log')
    ax.grid(True)
    plt.show(block=False)
    
#}}}
def test_only_to_show_all_plots():#{{{
    plt.show()
#}}}

def main():
    
    model = model_settings()
    model["reservoir_prop"]["K_type"] = 3# 0: K=cte from vertical wellbore; 1: K=linear func; 2: parabolic func; 3: singular K; 4: K comes from an ANN 
    model["wellbore_prop"]["f_type"] = 0 # 0: default friction factor (linear or non-linear); 1: linear only; 2: non-linear only 
    model["RK_settigns"]["verbose"] = 0 

    model["wellbore_prop"]["D"] = 0.1 # 0.04 to 1 worked fine 
    #model["RK_settigns"]["idp"] = 1.6324689395162981e+7 - 1.177e+7 
    model["RK_settigns"]["idp"] = 0.005e+7 # [Pa] initial delta P
    # 0.50e+7 # [Pa] initial delta P
    model["RK_settigns"]["npoints"] = 301
    #model["RK_settigns"]["rtol"] = 1e-10

    #model["fluid_prop"]["rho"] = 800. # 100 to 5000
    #model["fluid_prop"]["mu"] = 0.050 # 0.0005 to 0.050
    
    model["reservoir_prop"]["k"] = 1e-13 # 0.0001e-13 to 10e-13 

    # calling the RK solver. Solver for Q and p
    x, Q, p = RungeKutta_solver(model)

    # post-processing

    # retrieve the Reynolds number
    Re = Reynolds(model, Q) 

    # transforming RK x to FEM x
    Lw = model["wellbore_prop"]["Lw"]
    x_fem = Lw-x  

    # print final results
    print("Q heel = "+str(Q[-1]))
    print("p heel = "+str(p[-1]))
    print("p toe  = "+str(p[0]))

    # now, plot
    MPa = 1.e6
    fig, axs = plt.subplots(3, 1, layout='constrained')
    axs[0].plot(x_fem, p/MPa)
    #axs[0].set_xlim(0, 2)
    axs[0].set_xlabel('Wellbore distance (m)')
    axs[0].set_ylabel('Pressure (MPa)')
    axs[0].grid(True)

    axs[1].plot(x_fem, Q)
    axs[1].set_xlabel('Wellbore distance (m)')
    axs[1].set_ylabel('Flow rate (m3/s)')
    axs[1].grid(True)
   
    re_lim = np.zeros(Re.shape)
    re_lim[:] = 1187.38
    axs[2].plot(x_fem, Re)
    axs[2].plot(x_fem, re_lim, linestyle='dashed') 
    axs[2].set_xlabel('Wellbore distance (m)')
    axs[2].set_ylabel('Reynolds number []')
    axs[2].grid(True)

    plt.show()

if __name__ == "__main__":
    main()

