# code for the pure hydrodynamic model for clay 

import numpy, sys
import matplotlib.pyplot as plt

from Terracotta_functions import *

#================================
# Input
#================================

# Material constants for kaolin from Table 1 in [Wiebicke and Einav, 2024]
K       = 13.e6       #[kPa]
G       = 5.e6        #[kPa]
omega   = 0.6         #[-]
M       = 0.73        #[-]
p1      = 0.543*10**6 #[]
Lambda  = 10.0        #[-]

# rheological parameters
eta   = 1.5e2 #[(Ks)^(-1)]
alpha = 6.0e5 #[Ks]
gamma = alpha #[Ks]
beta  = 0     #[Ks]
Gamma = 1     #[K^2/kPa]

modelParam={"K":K, "G":G, "Lambda":Lambda, "p1":p1, "M":M, "eta":eta, "alpha":alpha, "beta":beta, "gamma":gamma, "omega":omega, "Gamma":Gamma}

# initial state
stressE = numpy.array([100,0])        # elastic stress
stressD = numpy.array([0,0])          # viscous stress
Tm      = 0                     # ... mesoscopic temperature
strain  = numpy.array([0,0])          # total strain
phi     = 0.423

stress = stressE + stressD + Tm**2    # total stress

# get the initial elastic strain 
epse_v, epse_s = initialElasticStrains(stressE,phi,modelParam)
stateVar = numpy.array([strain[0],strain[1],stressE[0],stressE[1],epse_v,epse_s,Tm,stressD[0],stressD[1],phi])

print("initial state variables \n", stateVar)

#list for loading conditions
#possible loading paths: txu, iso, oed, txd, relaxation, creep, txuCreep, pureDev, isoStress
#loading rate gives the 
# - strain rate in the case of txu, iso, oed, txd, pureDev
# - stress rate in the case of isoStress
# - time step in the case of relaxation, creep, creepOed, txuCreep
#targetIncs defines the total increment of 
# - strain in the case of txu, iso, oed, txd, pureDev
# - stress in the case of isoStress
# - time in the case of relaxation, creep, creepOed, txuCreep

testing = ["txu","relaxation"]
loadingRates = [2.e-5, 1.e-2]
targetIncs = [0.1, 1.e1]

filenameOutput = 'testOutput'

#================================
# Calculation
#================================

#integration along the given path
data = integratePath(stateVar,modelParam,testing,loadingRates,targetIncs)

#save results to a text file
printResults(filenameOutput,data)

#================================
# Plotting results
#================================

fig, axs = plt.subplots(2, 2, figsize=(10,10))

#plot deviatoric stress-strain response
axs[0,0].plot(data[:,1],data[:,3]+data[:,8], '-', label='Simulation')
axs[0,0].set_xlabel('shear strain $\epsilon_s$ [-]')
axs[0,0].set_ylabel('deviatoric stress $q$ [kPa]')
# axs[0,0].legend()
        
#plot stress path
maxP = numpy.amax(data[:,2]+data[:,6]**2/Gamma+data[:,7])
failure = numpy.arange(0,maxP+10,1)
axs[0,1].plot(failure, failure*M, '', c='r')
axs[0,1].plot(data[:,2]+data[:,6]**2/Gamma+data[:,7], data[:,3]+data[:,8], '')
axs[0,1].set_xlabel('pressure $p$ [kPa]')
axs[0,1].set_ylabel('deviatoric stress $q$ [kPa]')

#plot density evolution
axs[1,0].plot(data[:,1],data[:,9], '-', label='Simulation')
axs[1,0].set_xlabel('shear strain $\epsilon_s$ [-]')
axs[1,0].set_ylabel('solid fraction $\phi$ [-]')
axs[1,0].invert_yaxis()

#plot pressure-density response
pres = numpy.arange(10,500)
phi  = (pres/p1)**(1/Lambda)
phiC = (pres/(p1*omega))**(1/Lambda)

axs[1,1].plot(data[:,2]+data[:,6]**2/Gamma+data[:,7],data[:,9], '-', label='')
axs[1,1].plot(pres, phiC, 'red', label='CSL')
axs[1,1].set_xlabel('pressure $p$ [kPa]')
axs[1,1].set_ylabel('solid fraction $\phi$ [-]')
axs[1,1].invert_yaxis()

plt.tight_layout()
plt.show()



