# -*- coding: utf-8 -*-

### Importing Packages and Loading Field Values ### 
from __future__ import division
import numpy as np
from EandBpy import EandB
from scipy.integrate import odeint, quad
import titanPolyfit
import constants as c
import itertools

bulkVel = np.array([1*10**5,0,0],float)   # Bulk (average) velocity; in the positive x-direction (m/s, meters/second)
boxRes = 128750.                          # units of meters
rShift = boxRes*np.array([56.,45.,45.])   # shifts position vectors to have Titan centered at r = 0. units of meters.
LB=boxRes*np.array([5,5,5]) - rShift      # lower bound for particle trajectory. units of meters. 
UB=boxRes*np.array([106,84,84]) - rShift

s_nCoef = titanPolyfit.fitS_N()
s_eCoef = titanPolyfit.fitS_E()
N2Coef = titanPolyfit.fitN2()
ionFit,chExFit = titanPolyfit.fitIonAndCharge()   

# if you don't have bfield.npy or efield.npy, run snowden_code_py.py
bfield = np.load("bfield.npy")
efield = np.load("efield.npy")

class Particle:
    # iterator function that returns particle ID's
    newid = itertools.count().next
    
    def __init__(self,ipos,ivel):
        
        # initialize particle with velocity and position inputs
        self.pos = ipos
        self.vel = ivel
        self.initEn = .5*c.mass*(magF(ivel)**2)*c.JtoEV
        # call newid to give particle a new ID (integer)
        self.id = Particle.newid()
        
        # r and v store the entire history of the particle's trajectory
        self.r = np.zeros((c.nsteps,3))
        self.v = self.r.copy()
        self.r[0,:] = ipos
        self.v[0,:] = ivel
        
        # each particle has its own individual deposition vector
        atmosVector = np.arange(c.titanRad,c.atmosAltitude+c.titanRad,c.deltaR)
        self.depo = np.zeros((4,len(atmosVector)))      
        self.depo[0,:] = atmosVector
        
    def setAltitude(self):
        self.altitude = magF(self.pos)
        
    def setEnergy(self):
        self.energy = .5*c.mass*(magF(self.vel)**2)*c.JtoEV
        
    def setDeltaE(self):
        self.deltaE = calcDeltaE(self.altitude,self.energy)
        
    def setIncidentAngle(self):
        rhat = -self.pos/magF(self.pos)
        vhat = self.vel/magF(self.vel)
        self.incidentAngle = np.degrees(np.arccos(np.dot(rhat,vhat)))
    
    def storeTrajectory(self,pos,vel,tstep):
        self.r[tstep,:] = pos
        self.v[tstep,:] = vel


# createParticles creates a list of particles via the Particle class

def createParticles(numParticles):
    particleList = []
    init_energies = []
    
    # Bin Properties  
    maxE = 4.0*10**3                            # eV
    minE = 0.                                   # eV
    binSize = 10.                               # eV

    energyBin = np.arange(minE,maxE,binSize)
        
    velocityBin = (2*energyBin/(c.JtoEV*c.mass))**(1./2.) # convert energies to velocities  
    
    N = np.zeros(len(velocityBin),dtype='int') 
    
    # this loop integrates the Maxwellian velocity distribution to find 
    # the number of particles that should be in each bin
    for i in range(len(N)-1):
        N[i] = np.round(quad(MaxwellDist,velocityBin[i],velocityBin[i+1],
        args=(numParticles))[0])
    
    # create the particles. We loop through each veloctity
    # bin, and then through the number of particles in each bin.
    for binNum in range(len(N)):
        for i in range(N[binNum]):
            particleList.append(Particle(randPos(),randVel(velocityBin, binNum)))
            init_energies.append(particleList[i].initEn)
            
    # calculate the weights for the particles.
    distE = np.multiply(intensity_drift(velocityBin), energyBin)
    histoE = np.histogram(init_energies,energyBin)[0]
    histoE = np.append(histoE,0)
    
    weights = np.divide(distE*c.dectorArea, histoE)
    counts = np.digitize(init_energies,energyBin)
    
    for w in range(len(weights)):
        if np.isinf(weights[w]) or np.isnan(weights[w]):
           weights[w] = 0
  
    i = 0    
    for p in particleList:
        p.weight = weights[counts[i]]
        i += 1
    return particleList

def ranNum():
    return np.random.rand()

def magF(r): 
    return np.sqrt(r.dot(r))

def intensity_drift(v):
        part_dens = .09*10**-6
        v0 = magF(bulkVel)
        return ((v**2)/c.mass + (2*c.kb*c.T/(c.mass**2)))*part_dens \
        *((2*np.pi*c.kb*c.T/c.mass)**(-3/2))*np.exp(-c.mass*((v-v0)**2)/(2*c.kb*c.T))/(4*np.pi)

# Maxwellian distribution of velocities
def MaxwellDist(v,numParticles):
    return numParticles*4*np.pi*(v**2)*((2*np.pi*c.kb*c.T/c.mass)**(-3./2.))*np.exp(-c.mass*(v**2)/(2*c.kb*c.T))
        
# calculates random initial position for a particle near Titan
# all particles are intially located 3-5 altitudes from Titan
def randPos(): 
    Radius = 3*c.atmosAltitude + c.titanRad 
    randTheta = np.pi*ranNum()
    randPhi = 2*np.pi*ranNum()
    rhat = np.array([np.sin(randTheta)*np.cos(randPhi),np.sin(randTheta)*np.sin(randPhi),np.cos(randTheta)]) # unit vector      
    return Radius*rhat
    
def randVel(velocityBin,i): # calculates random initial velocity for a particle near Titan.
    # i indicates which velocity bin we're in 
    randTheta = np.pi*ranNum()
    randPhi = 2*np.pi*ranNum()
    vhat = np.array([np.sin(randTheta)*np.cos(randPhi),np.sin(randTheta)*np.sin(randPhi),np.cos(randTheta)]) # unit vector
    return  velocityBin[i]*vhat+bulkVel

def outOfBounds(r):
    if r[0] < LB[0] or r[0] > UB[0] or r[1] < LB[1] or r[1] > UB[1] or r[2] < LB[2] or r[2] > UB[2]:
        return True

# lorentz calculates the acceleration of a particle in the magnetosphere
# of Saturn based on the Lorentz force. Note that t is not used in any
# calculations, but it is required as a parameter for odeint.

def lorentz(vals,t,bfield,efield):
    pos = vals[0:3]
    vel = vals[3:6]
    E,B = EandB(pos + rShift, bfield, efield) # calculating field values at current position
    a = (c.q/c.mass)*(E+np.cross(vel,B)) # get acceleration
    return np.concatenate((vel,a))

def runInfo(numParticles,sim):
    while numParticles % sim != 0 and sim < numParticles:
        sim+=1
    runs = int(numParticles/sim)
    return sim,runs
   
def IonAndCE(p,numLeft,numLeftCS):
    # Net Ionization
    probIonize = (titanPolyfit.N2Dens(N2Coef,p.altitude-c.titanRad)*np.polyval(ionFit,p.energy/
    (10.**3))*(10.**-21))
    
    numIon = numLeft*probIonize
    numLeft -= numIon
    
    # Net Charge Exchange        
    probCS = (titanPolyfit.N2Dens(N2Coef,p.altitude-c.titanRad)*np.polyval(chExFit,p.energy/
    (10.**3))*(10.**-21))
    numCS = numLeftCS*probCS
    numLeftCS -= numCS
                
    return numCS,numIon,numLeftCS,numLeft

# calcDeltaE calculates the change in energy of a particle based on 
# cross section data, initial energy, and altitude.
def calcDeltaE(particleAlt,energy):
    return titanPolyfit.N2Dens(N2Coef,particleAlt-c.titanRad)*(titanPolyfit.getPolyVal(s_nCoef,energy)+titanPolyfit.getPolyVal(s_eCoef,energy))*c.delX

# magnetosphere trajectory calculations
def outsideAtmos(p): 
    init_cond = np.concatenate((p.pos,p.vel))
    for tstep in range(c.nsteps-1):
           
    # use odeint with the to solve for the trajectory between adjacent timesteps
        soln = odeint(lorentz, init_cond, [tstep*c.dt, (tstep+1)*c.dt],args=(bfield,efield))
        
        # store the newly calculated trajectory info in our r and v arrays
        p.pos = soln[1,0:3]
        p.vel = soln[1,3:6]
        
        p.storeTrajectory(p.pos,p.vel,tstep)
        
        # check to see if the particle is out of bounds
        if outOfBounds(p.pos) == True: 
            p.hitAtmos = False
            break
        # check to see if the particle hit the atmosphere of Titan
        if magF(p.pos) < c.titanRad + c.atmosAltitude:
            p.hitAtmos = True
            p.hitTime = tstep
            
            p.setEnergy()
            p.incEnergy = p.energy 
            
            p.setIncidentAngle()
            break
        else: 
            init_cond = soln[1]
     
    if tstep+1 == c.nsteps:
        p.hitAtmos = False
    
    return p

# inside Titan's Atmosphere trajectory calculations
def insideAtmos(p):
    # each simulated particle represents a number of actual particles equal to its weight
    numLeft = p.weight
    numLeftCS = p.weight
    
    p.setEnergy()          
    step = p.hitTime + 1          

    while p.energy > 0:
        p.setAltitude()
        p.setDeltaE()
        
        p.energy -= p.deltaE
        vhat = p.vel/magF(p.vel)
        
        if p.energy < 0:
            break
        if p.altitude < c.titanRad or p.altitude > (c.titanRad+c.atmosAltitude):
            break
        
        numCS,numIon,numLeftCS,numLeft = IonAndCE(p,numLeft,numLeftCS)
        
        # store the deposition values for each layer of the atmosphere:
        for i in range(len(p.depo[0,:])-1):
            if p.altitude < p.depo[0,i+1] and p.altitude > p.depo[0,i]:
                
                p.depo[1,i] += p.deltaE
                p.depo[2,i] += numIon
                p.depo[3,i] += numCS 
        
        p.pos += (vhat*c.delX)
        p.storeTrajectory(p.pos,p.vel,step)
        step += 1
        
    return p
