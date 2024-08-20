# -*- coding: utf-8 -*-

"""
Energy Deposition Rates and Particle Trajectory

"""
import numpy as np

### BEGIN PROGRAM PARAMETERS ###

numParticles = 200                        # Number of particles to simulate
altitude = 1600000.                       # Radial distance to which Titan's atmosphere extends from the surface (m, meters)
T = 250.*11604.                           # average temp of particles in Kelvin; 1eV = 11604 Kelvin
q = 1.6*10**(-19)                         # Charge of one particle (Coulombs)
bulkVel = np.array([1*10**5,0,0],float)   # Bulk (average) velocity; in the positive x-direction (m/s, meters/second)
kb = 1.38*10**(-23)                       # Boltzmann Constant (Joules/Kelvin)

m = 16*1.67*10**(-27)                     # Mass of Oxygen-16 (kilograms)
titanRad=2576000.0                         # Radius of Titan; Meters
delX = 1000.                              # Particle Deposition Step (meters)
dt = .3                                   # timestep; units of seconds
tf = 10000                                  # final time; units of seconds
nstep = int(tf/dt)
a = np.sqrt(kb*T/m)                       # Maxwell scaling factor

# Superparticle Characteristics:
APD = .05                                 # Actual Particle Density. Average CAPS data. Units of O2+/cm^3
OPD = numParticles/((4./3.)*np.pi*((titanRad + 3*altitude)**3 - (titanRad+altitude)**3)*(10**6)) # "Our Particle Density"; units of O2+/cm^3. Density of particles in volume we start them in (a spherical shell around Titan)
SPSF=APD/OPD # "Super Particle Scaling Factor"; the number of particles each simulated particle represents.
PIA = 0 # number of particles that enter the atmosphere of Titan
AtmP = np.zeros((numParticles,7))     # empty array to be used to store information about particles

# Particle Desposition Altitude/Energy Vector
deltaR=20000. # spacing between segments of the atmosphere (meters)
myRvector = np.arange(titanRad,altitude+titanRad,deltaR)
altEn = np.zeros((4,len(myRvector))) # matrix to store cumulative ionization, charge exchange, and energy deposition at each layer of the atmosphere 
altEn[0,:] = myRvector # altEn stores the energy lost to each layer of the atmosphere. layer thickness determined by deltaR


# Titan and Grid Attributes
boxRes = 2*128750.                        # units of meters
rShift = boxRes*np.array([56.,45.,45.])   # shifts position vectors to have Titan centered at r = 0. units of meters.
LB=boxRes*np.array([5,5,5]) - rShift      # lower bound for particle trajectory. units of meters.
UB=boxRes*np.array([106,84,84]) - rShift  # upper bound for particle trajectory. units of meters. 

### END PROGRAM PARAMETERS ###

### polyfitting cross-section data ### 
from titanPolyfit import fitS_N, fitS_E, fitN2, getPolyVal, N2Dens, fitIonAndCharge

### BEGIN FUNCTION DEFINITIONS ###

# Gives Maxwellian distribution of velocities
def MaxwellDist(v):
    return numParticles*4*np.pi*(v**2)*((2*np.pi*kb*T/m)**(-3./2.))*np.exp(-m*(v**2)/(2*kb*T))

def f_ode(values,t):
    pos=values[:3]
    vel=values[3:6]
    E,B = EandB(pos + rShift, bfield, efield) # calculating field values at current position
    a = (q/m)*(E+np.cross(vel,B)) # get acceleration
    return np.concatenate((vel,a))

def magF(r):
    return np.sqrt(r[0]**2 + r[1]**2 + r[2]**2)

def ranNum(): # not strictly necessary, but the code looks a little nicer this way
    return np.random.rand()
    
def randPos(): # calculates random initial position for a particle near Titan
    randRadius = 3*altitude+titanRad+ranNum()*2*altitude # particle starts 3-5 altitudes away from Titan
    randTheta = np.pi*ranNum()
    randPhi = 2*np.pi*ranNum()
    rhat = np.array([np.sin(randTheta)*np.cos(randPhi),np.sin(randTheta)*np.sin(randPhi),np.cos(randTheta)]) # unit vector      
    return randRadius*rhat
    
def randVel(velocityBin,i): # calculates random initial velocity for a particle near Titan.
    # i indicates which velocity bin we're in 
    randTheta = np.pi*ranNum()
    randPhi = 2*np.pi*ranNum()
    vhat = np.array([np.sin(randTheta)*np.cos(randPhi),np.sin(randTheta)*np.sin(randPhi),np.cos(randTheta)]) # unit vector
    return  velocityBin[i]*vhat+bulkVel

def ParticleInfo(r,v,k,n,DFT,altEn):
    x = r[k+1,0,n]
    y = r[k+1,1,n]
    z = r[k+1,2,n]    
        
    ThetaP = np.degrees(np.arccos(z/DFT)) # theta coordinate of the current particle 
        
    # phi coordinate of the current particle:
    if x > 0:
        PhiP = np.degrees(np.arctan(y/x))
    elif x < 0: 
        PhiP = np.degrees(np.arctan(y/x)) + 180
    else:
        PhiP = 90.
    
    # ensuring 0 < PhiP < 360
    if PhiP < 0:
        PhiP += 360
    
    # calculate energies
    energyInit = .5*m*(v[0,0,n]**2 + v[0,1,n]**2 + v[0,2,n]**2)*6.24*10**18 # initial energy in eV
    #print "Initial Energy = ", energyInit
    energyFin = .5*m*(v[k+1,0,n]**2 + v[k+1,1,n]**2 + v[k+1,2,n]**2)*6.24*10**18 # final energy in eV
    #print "Incident Energy = ", energyFin
    
    ### Particle Desposition ###
    
    vhat = v[k+1,:,n]/magF(v[k+1,:,n])
    rhat = -r[k+1,:,n]/magF(r[k+1,:,n]) 
    
    incidentAngle = np.degrees(np.arccos(np.dot(vhat,rhat)))
    numIon = 0.
    numCS = 0.
    numLeft = SPSF
    numLeftCS = SPSF     
    
    # store coordinates of the particle when it hits Titan, as well as the initial energy of the particle and its energy when it enters the atmosphere
    particle_info = np.array([PhiP,ThetaP,energyInit,energyFin,incidentAngle,numIon,numCS])

    energyP = energyFin # storing the incident energy of the particle in a new variable
    p = k+1 # new iterator for the particle
  
    # Energy Deposition
    while energyP > 0:
        particleAlt = magF(r[p,:,n])
        
        # the change in energy is calculated using our cross-sections
        
        deltaE = N2Dens(N2Coef,particleAlt-titanRad)*(getPolyVal(s_nCoef,energyP)+getPolyVal(s_eCoef,energyP))*delX 
        
        energyP -= deltaE # new energy; units of eV
        
        #energyP = np.subtract(energyP,deltaE)
        if energyP < 0:
            break
        # check to see if particle is still in Titan's atmosphere
        
        if particleAlt < titanRad or particleAlt > (titanRad+altitude):
            break
        
        # Net Ionization
        probIonize = N2Dens(N2Coef,particleAlt-titanRad)*np.polyval(ionFit,energyP/(10.**3))*(10.**-21)
        numIon=numLeft*probIonize
        particle_info[5] += numIon
        numLeft -= numIon
        
        # Net Charge Exchange        
        probCS = N2Dens(N2Coef,particleAlt-titanRad)*np.polyval(chExFit,energyP/(10.**3))*(10.**-21)
        numCS = numLeftCS*probCS
        particle_info[6] += numCS
        numLeftCS -= numCS
        
        for gamma in range(len(altEn[0,:])-1):
            
            if particleAlt < altEn[0,gamma+1] and particleAlt > altEn[0,gamma]:
                altEn[1,gamma] += deltaE
                altEn[2,gamma] += numIon
                altEn[3,gamma] += numCS
        
        r[p+1,:,n] = r[p,:,n] + (delX*vhat)
        p += 1 
    return r, particle_info, altEn, p


def particleTrajectory(i,n,PIA,r,v,altEn,AtmP):
# generate initial conditions, store in a single array to be passed to odeint
    init_cond = np.concatenate((randPos(),randVel(i)))
    r[0,:,n] = init_cond[0:3]
    v[0,:,n] = init_cond[3:6]
    
    for k in range(nstep-1):
           
    # use odeint with the function "f_ode" to solve for the trajectory between adjacent timesteps
        
        soln = odeint(f_ode,init_cond,[k*dt, (k+1)*dt])
        
        
        # store the newly calculated trajectory info in our r and v arrays
        r[k+1,:,n] = soln[1,0:3]
        v[k+1,:,n] = soln[1,3:6]
    
        # check to see if the particle is out of bounds
        if r[k+1,0,n] < LB[0] or r[k+1,0,n] > UB[0] or r[k+1,1,n] < LB[1] or r[k+1,1,n] > UB[1] or r[k+1,2,n] < LB[2] or r[k+1,2,n] > UB[2]:
                stopPoint[n] = k+1
                endTimeVector[n] = dt*(k+1)
                break
        
        # check to see if the particle hit the atmosphere of Titan and store information about the particle.
        
        DFT = magF(r[k+1,:,n]) # Distance (of particle From Titan 
        
        if DFT < (altitude + titanRad):
            #print "Particle %s hit the atmosphere!"%(n+1)
            PIA += 1
            # recording the time at which the particle hit Titan's atmosphere
            stopPoint[n] = k+1
            particleHitLogic[n] = 1 
            endTimeVector[n] = dt*(k+1)
            r,AtmP[n,:],altEn,pdStopPoint[n] = ParticleInfo(r,v,k,n,DFT,altEn)
            break
        else:
            init_cond = soln[1,:]
            
    # if the particle didn't hit the atmosphere, we set the initial conditions to our most recent calculation, 
    # so that odeint uses it to get the trajectory for the next timestep
    
    if k+1 == nstep:
        stopPoint[n] = k+1
        endTimeVector[n] = dt*(k+1)
    
    return r,v,AtmP,altEn,PIA

### END FUNCTION DEFINITIONS ###

### BEGIN SETUP CALCULATIONS/VECTOR INITIALIZATION ### 

# calculate polyfit coefficients for each dataset (see titanPolyfit.py)
# s_e = electric cross-section, s_n = nuclear cross section, N2 = nitrogen density profile
s_nCoef = fitS_N()
s_eCoef = fitS_E()
N2Coef = fitN2()
ionFit,chExFit = fitIonAndCharge()

# initialize energy bin, convert energies to velocities
maxE = 4.0*10**3 #eV
minE = 0. # eV
binSize = 10 # eV
energyBin = np.arange(minE,maxE,binSize) # array of energies
energyBin *= 1.602*10**(-19) # convert from eV to Joules
velocityBin = (2*energyBin/m)**(1./2.)

# integrate Maxwellian distribution between the limits of each box
# to find the number of particles in each box.

from scipy.integrate import quad, odeint

N = np.zeros(len(velocityBin),dtype='int')
for i in range(len(N)-1):
    N[i] = np.round(quad(MaxwellDist,velocityBin[i],velocityBin[i+1])[0])

# Initializing Vectors:
n = 0                                  # particle iterator; indicates which particle we're looking at
stopPoint = np.zeros(numParticles)     # a vector that stores the maximum row number for each n value
pdStopPoint = np.zeros(numParticles)    # stopping point for particle deposition in atmosphere
r = np.zeros((nstep,3,numParticles))   # initializing position 3-D array
v = np.zeros((nstep,3,numParticles))   # initializing velocity 3-D array

particleHitLogic = np.zeros((numParticles,numParticles)) # stores boolean data about which particles hit the atmosphere and which did not
endTimeVector = np.zeros(int(np.sum(N)))                 # tells us how much time passed before the particle went out of bounds

# Particle Deposition Altitude/Energy Vector
deltaR = 5000.                                           # spacing between segments of the atmosphere; units of meters
myRvector = np.arange(titanRad,altitude+titanRad,deltaR) # array of positions running from the surface of Titan to the extent of the atmosphere. used only for the calculation of the following energy vector.
altitudeEn=np.zeros((2,len(myRvector)))                  # 2-D array containing positions and corresponding energies
altitudeEn[0,:] = myRvector                              # set all values in the first row equal to the values of myRvector; the second row contains the energies at each point


if 'bfield' not in globals() and 'efield' not in globals():
    from snowden_code_py import build_fields
    print "Obtaining field values . . . "
    bfield, efield = build_fields()

### END SETUP CALCULATIONS/VECTOR INITIALIZATIONS ###

### BEGIN CALCULATION OF PARTICLE TRAJECTORIES ###

# Method 1

from EandBpy import EandB
import sys
from multiprocessing import Queue, Process

if __name__ == '__main__':

    for i in range(len(N)): # pick a velocity bin       
    
        sys.stdout.write("\rCalculating trajectories. %s percent complete." % (100*float(n)/float(np.sum(N))))     
    
        for j in range(N[i]): # iterate through each particle in the current velocity bin
            q = Queue()
            p = Process(target=particleTrajectory,args=())
            p.start()

            p.start()
            p.join()
            n += 1
sys.stdout.write("\r\n %s particles hit the atmosphere of Titan. Godspeed, particles. Godspeed. \n"%PIA)
print("The longest time taken was %s seconds, and the shortest time taken was %s seconds.")%(np.max(endTimeVector),np.min(endTimeVector))
    
 # converting eV to eV/(cm^3*s)
for i in range(len(altEn[0,:])-1):
    altEn[1,i] = SPSF*altEn[1,i]/((4./3.)*np.pi*(altEn[0,i+1]**3 - altEn[0,i]**3)*(10**6))

# measuring distances relative to Titan's surface; converting m to km
altEn[0,:] = (altEn[0,:] - titanRad)/1000

### Graphing Parameters ###

# For each graph, 1 means ON (the graph will be displayed), and 0 means OFF (the graph will not be displayed)

particleDepositionGraph = 1;            # Energy Deposition Graph
incidentAngleHistogram = 1;             # Histogram of incident angles
incidentEnergyHistogram = 1;            # Histogram of incident energies
particleIonGraph = 1;                   # Graph of Ionization Deposition
particleChargeExchange = 1;             # Graph of Charge Exchange Desposition

"""
Plots not currently implemented:

powerFlux3D = 0; # Second 3D without points
particleDepositionTrajectoryGraph = 0;
particleTrajectoryGraph = 0;
powerFlux2D = 0; # 2nd Graph

"""
from titanPlots import plots
plots(particleDepositionGraph,incidentEnergyHistogram,incidentAngleHistogram,particleIonGraph,particleChargeExchange,altEn,AtmP)
            


     
            
            
        
        
        
                             
    
        






