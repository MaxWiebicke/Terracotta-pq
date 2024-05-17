# collection of functions for the constitutive models
import numpy
import math
import sys

from scipy.optimize import fsolve
from scipy.integrate import solve_ivp

def eps2pq(eps_1,eps_2):
    """
    This function transforms the principal strains into the triaxial invariants volumetric and shear strain.
    """
    eps_v = eps_1+2*eps_2
    eps_s = 2/3*(eps_1-eps_2)
    return numpy.array([eps_v,eps_s])

# stress and strain increments on different loading paths
def testControl(loadingRate,testType):        
    """
    This function gets the loading rate for strain rate controlled tests.
    It is defined by the type of the loading and the given maximal principle strain rate.
    """
    if testType=="iso":
        deps_v = 3*loadingRate
        deps_s = 0*loadingRate
        
    elif testType=='oed':
        deps_1 = 1*loadingRate
        deps_2 = 0*loadingRate
        deps_v, deps_s = eps2pq(deps_1,deps_2)
        
    elif testType=='txu':
        deps_1 = 1*loadingRate
        deps_2 = -0.5*loadingRate
        deps_v, deps_s = eps2pq(deps_1,deps_2)
        
    elif testType=='relaxation':
        deps_v = 0*loadingRate
        deps_s = 0*loadingRate
    
    return numpy.array([deps_v,deps_s])

def matModel(strainE,phi,modelParam):
    """
    This function determines the instantaneous stiffness matrix defined by the partial derivatives of the elastic stress components with elastic strain.
    """
    epse_v = strainE[0]
    epse_s = strainE[1]
    
    K = modelParam["K"]
    G = modelParam["G"]
    
    M_pp = phi**6*K*epse_v 
    M_pq = phi**6*3*G*epse_s
    M_qp = M_pq
    M_qq = phi**6*3*G*epse_v
        
    M = numpy.array([[M_pp,M_pq],[M_qp,M_qq]])
    
    return M

def initialElasticStrains(stress,phi,modelParam):
    """
    This function calculates the elastic strain components given the elastic stress.
    It is used to set the initial state of the state variables given a stress.
    """
    pe = stress[0]
    qe = stress[1]
    K = modelParam["K"]
    G = modelParam["G"]
    
    x = 1/phi**6 * ( pe/K + math.sqrt( (pe/K)**2 - qe**2/(3*G*K) ) )
    epse_v = x**(1/2)
    epse_s = qe/(3*G*phi**6*epse_v)
        
    return epse_v, epse_s
    
def plastCoefficients(stress,Tg,phi,modelParam):
    """
    This function determines the plastic transport coefficients a, b and c.
    Furthermore, it performs a check on whether Onsager's reciprocity relations are satisfied.
    If that is not the case, it currently throws a warning and exits the calculation.
    """
    eta   = modelParam["eta"]
    Lambda = modelParam["Lambda"]
    p1      = modelParam["p1"]
    M      = modelParam["M"]
    omega  = modelParam["omega"]
    
    alpha, beta, gamma = viscCoefficients(modelParam)
    
    pc = p1*phi**Lambda
    
    a = math.sqrt(eta/alpha)/pc 
    b = -a/M*((stress[1])/(M*(stress[0])))**(1)
    c = math.sqrt(eta/gamma)/(M*omega*pc) + a/M**2
    
    # check Onsager 
    if a*c-b**2 < 0:
        print("Onsager's requirement not satisfied.")
        print("a*c-b**2 > 0 ....")
        # print(r"$\sigma_1 = $", 2/3*stress[1]+stress[0])
        print(stress)
        sys.exit()
    
    return a,b,c
    
def viscCoefficients(modelParam):
    """
    This function fetches the viscous transport coefficients.
    In the current version of the model, this is obsolete as they are only constants and not dependent on other state variables.
    """
    alpha = modelParam["alpha"]
    beta  = modelParam["beta"]
    gamma = modelParam["gamma"]
    
    return alpha, beta, gamma
    
def plasticStrainRate(stress,dStrain,Tm,phi,modelParam):    
    """
    This function determines the plastic strain rate.
    """
    a,b,c = plastCoefficients(stress,Tm,phi,modelParam)
    A = numpy.array([[a,b],[b,c]])
    dStrainP = Tm*numpy.dot(A,stress)
    
    return dStrainP
    
def elasticStrainRate(stress,Tm,phi,dStrain,modelParam):
    """
    This function determines the elastic strain rate.
    """
    dStrainP = plasticStrainRate(stress,dStrain,Tm,phi,modelParam)
    
    dStrainE = dStrain - dStrainP
    
    return dStrainE
    
def TmRate(Tm,dStrain,modelParam):
    """
    This function determines the rate of the meso-related temperature.
    """
    alpha, beta, gamma = viscCoefficients(modelParam)
    dTm = alpha*dStrain[0]**2 + gamma*dStrain[1]**2 - modelParam["eta"]*Tm**2
    
    return dTm

def stateVarRatesComplexLoading(t,stateVar,*data):
    """
    This function calculates the rates for all the variables given the loading path.
    If the loading is not purely strain rate controlled, it calls another function that solves for the strain rate given the respective loading conditions.
    """
    modelParam, testType, dStrainPrev, dt, loadingRate, t_prev = data
    
    alpha, beta, gamma = viscCoefficients(modelParam)
    Gamma = modelParam["Gamma"]
    
    strain  = stateVar[0:2]
    stress  = stateVar[2:4]
    strainE = stateVar[4:6]
    Tm      = stateVar[6]
    pd, qd  = stateVar[7], stateVar[8]
    phi     = stateVar[9]
    
    #get the elastic stiffness matrix
    Me = matModel(strainE,phi,modelParam)
    
    #get the total strain rate
    if testType == "creep":
        args = (stateVar,dStrainPrev,dt,modelParam,testType,loadingRate)
        deps_v, deps_s = fsolve(funcComplexLoading, (1.e-6, 0), args)
        
        dStrain = numpy.array([deps_v, deps_s])
    
    elif testType == "txd":
        args = (stateVar,dStrainPrev,dt,modelParam,testType,loadingRate)
        dp, deps3 = fsolve(funcComplexLoading, (1.e-6, 0), args)
        deps1 = loadingRate
        deps_v = (deps1+2*deps3)
        deps_s = 2/3*(deps1-deps3)
        dStrain = numpy.array([deps_v, deps_s]) 
        
    elif testType == "txuCreep":
        args = (stateVar,dStrainPrev,dt,modelParam,testType,loadingRate)
        dp, deps_s = fsolve(funcComplexLoading, (1.e-6, dStrainPrev[1]), args)
        deps_v = 0
        dStrain = numpy.array([deps_v, deps_s])
        
    elif testType == "pureDev":
        args = (stateVar,dStrainPrev,dt,modelParam,testType,loadingRate)
        deps_v, dq = fsolve(funcComplexLoading, (dStrainPrev[0], dStrainPrev[1]), args)
        deps_s = loadingRate
        dStrain = numpy.array([deps_v, deps_s])
        
    elif testType == "isoStress":
        args = (stateVar,dStrainPrev,dt,modelParam,testType,loadingRate)
        deps_v, dq = fsolve(funcComplexLoading, (dStrainPrev[0], 0), args)
        deps_s = 0
        dStrain = numpy.array([deps_v, deps_s])
        
    else:
        dStrain = testControl(loadingRate,testType)
    
    #determine the rates of the variables
    #elastic strain
    dStrainE = dStrain - plasticStrainRate(stress,dStrain,Tm,phi,modelParam)
    #solid fraction
    dPhi = phi * dStrain[0] 
    #elastic stress
    dpe = Me[0,0]*dStrainE[0] + Me[0,1]*dStrainE[1] + stress[0]*6/phi*dPhi
    dqe = Me[1,0]*dStrainE[0] + Me[1,1]*dStrainE[1] + stress[1]*6/phi*dPhi
    #meso-related temperature
    dTm = TmRate(Tm,dStrain,modelParam)
    
    #Here, we determine the acceleration of the strain (i.e. the rate of the strain rate).
    #We assume that it is linear for each time step.
    if t == 0.:
        ddstrain = numpy.array([0.,0.])
    else:
        ddstrain = (dStrain-dStrainPrev)/dt
    #change of the strain rate within this time step
    dStrainUpdate = dStrainPrev + ddstrain*t 
    #viscous stress
    dpd = 2*alpha/Gamma*( dTm*dStrainUpdate[0] + Tm*ddstrain[0])
    dqd = 2*gamma/Gamma*( dTm*dStrainUpdate[1] + Tm*ddstrain[1])
    
    return numpy.hstack((dStrain[0],dStrain[1],dpe,dqe,dStrainE[0],dStrainE[1],dTm,dpd,dqd,dPhi))
    
def funcComplexLoading(vars, *data):
    """
    This function determines the strain rate for loading paths, that are not purely controlled by the strain rate.
    It solves the constitutive equation dependent on the conditions of the particular loading for a single time step.
    """
    stateVar, dStrainPrev, dt, modelParam, testType, loadingRate = data
    
    strain  = stateVar[0:2]
    stress  = stateVar[2:4]
    strainE = stateVar[4:6]
    Tm      = stateVar[6]
    pd, qd  = stateVar[7], stateVar[8]
    phi     = stateVar[9]
    
    Me = matModel(strainE,phi,modelParam)
    a,b,c = plastCoefficients(stress,Tm,phi,modelParam)
    alpha, beta, gamma = viscCoefficients(modelParam)
    
    Gamma   = modelParam["Gamma"]
    eta    = modelParam["eta"]
    
    if testType=="creep":
        deps_v, deps_s = vars
        dp = 0
        dq = 0
        
    elif testType=="txd":
        deps1 = loadingRate
        
        dp, deps3 = vars
        dq = 3*dp
        
        deps_v = deps1+2*deps3
        deps_s = 2/3*(deps1-deps3)
        
    elif testType=="txuCreep":
        dp, deps_s = vars
        
        deps_v = 0
        dq     = 0
        
    elif testType=="pureDev":
        deps_v, dq = vars
        
        dp = 0
        deps_s = loadingRate
        
    elif testType=="isoStress":
        dp = loadingRate
        deps_s = 0
        
        deps_v, dq = vars
        
    dTm = (alpha*deps_v**2 + gamma*deps_s**2 - eta*Tm**2)/2
    
    dpe = 6*stress[0]*deps_v + Me[0,0]*deps_v + Me[0,1]*deps_s - Tm*(Me[0,0]*(a*stress[0]+b*stress[1])+Me[0,1]*(b*stress[0]+c*stress[1]))
    dpd = 2*alpha/Gamma*(deps_v*dTm + Tm*(deps_v-dStrainPrev[0])/dt)
    dpT = Tm*dTm*2
    
    dqe = 6*stress[1]*deps_v + Me[1,0]*deps_v + Me[1,1]*deps_s - Tm*(Me[1,0]*(a*stress[0]+b*stress[1])+Me[1,1]*(b*stress[0]+c*stress[1]))
    dqd = 2*gamma/Gamma*( deps_s*dTm + Tm*(deps_s-dStrainPrev[1])/dt)
    
    #functions to minimise
    f0 = - dp + dpe + dpd + dpT
    f1 = - dq + dqe + dqd
    
    return [f0, f1]
    
def ratesIntegrationScipyCL(stateVar0,dStrainPrev,modelParam,loadingRate,dt,testType):
    """
    This function integrates the state variables over a time step.
    """
    sol = solve_ivp(stateVarRatesComplexLoading,
                    args=[modelParam, testType, dStrainPrev, dt, loadingRate, numpy.array([0.])],
                    t_span=[0.,dt],
                    t_eval=[dt],
                    y0=stateVar0,
                    vectorized=False,
                    # atol = 1.e-12,
                    # rtol = 1.e-12,
                    method='BDF')
    
    stateVar = sol.y[:,-1]
    
    # computation of the current total strain rate
    dStrainPrev = (stateVar[0:2]-stateVar0[0:2])/dt
    
    return stateVar, dStrainPrev

def integratePath(stateVar,modelParam,testing,loadingRates,targetIncs):
    """
    The function integrates the model along potentially multiple loading paths.
    """
    
    #initialise data
    t = 0
    data = numpy.array([[stateVar[0], stateVar[1], stateVar[2], stateVar[3], stateVar[4], stateVar[5], stateVar[6], stateVar[7], stateVar[8], stateVar[9], t ]])
    
    dStrainPrev = numpy.array([0,0])
    
    for i in range(len(testing)):
        print("---------------------------")
        print("loading: "+testing[i] )
        testType=testing[i]
        # determine step sizes for integration
        if testType=="relaxation" or testType=="creep" or testType=="txuCreep":
            #time increment for the integration step
            dt = loadingRates[i]
            steps = int(targetIncs[i]/dt)
            print("Time increment       [s]   = ", dt)
            print("Time step            [s]   = ", targetIncs[i])
        elif testType=="isoStress":
            #stress increment for the integration step
            stressInc = math.copysign(1.e-1,loadingRates[i])
            dt = stressInc/loadingRates[i]
            steps = int(targetIncs[i]/stressInc)
            print("Stress rate        [kPa/s] = ", loadingRates[i])
            print("time increment       [-]   = ", dt)
            print("Stress increment     [-]   = ", stressInc)
        else:
            #strain increment for the integration step
            epsInc = math.copysign(1.e-4,loadingRates[i])
            dt = epsInc/loadingRates[i]
            steps = int(targetIncs[i]/epsInc)
            print("Strain rate          [1/s] = ", loadingRates[i])
            print("time increment       [-]   = ", dt)
            print("Strain increment     [-]   = ", epsInc)
            
        for step in range(steps):
            print(step, " / ", steps)
            t += dt
            loadingRate = loadingRates[i]
            
            stateVar, dStrainPrev = ratesIntegrationScipyCL(stateVar,dStrainPrev,modelParam,loadingRate,dt,testType)
            
            data = recordData(data,stateVar,t)
            step += 1

    return data

def recordData(data,stateVar,t):
    """
    This function assembles the output of all variables for each time step.
    """
    data = numpy.append(data,[[stateVar[0], stateVar[1], stateVar[2], stateVar[3], stateVar[4], stateVar[5], stateVar[6], stateVar[7], stateVar[8], stateVar[9], t ]], axis=0)
    return data

def printResults(filename,data):
    """
    This function prints the output to a text file.
    """
    outputfile = open(filename+'.txt','w')
    outputfile.write ('#epsp #epsq #p #q #epsev #epses #Tg #pd #qd #e #t \n')
    for i in range(numpy.shape(data)[0]):
        outputfile.write('%1.5e %1.5e %1.5e %1.5e %1.5e %1.5e %1.5e %1.5e %1.5e %1.5e %.8f \n' %(data[i,0], data[i,1], data[i,2], data[i,3], data[i,4], data[i,5], data[i,6], data[i,7], data[i,8], data[i,9], data[i,10]))
    outputfile.close()
    return
