# -*- coding: utf-8 -*-
"""

This file contains all of the polyfit calculations that used to be in 
Energy_Deposition_With_Particle_Trajectory.m. 

We use h5py to import the N2 density profile, which is of format .mat. Note that 
our method of import works only with .mat files of version 7.3 or later.

"""

import numpy as np

def fitS_N():
    # Load values of energy in keV and s_n (nuclear cross-section) in eV*cm^2 #
    energy, s_n = np.loadtxt('stopping_cross_oxygen.dat',skiprows = 2,usecols=(0,1),unpack=True)

    # Calculate polyfit coefficients for s_n; polynomial is of degree 4 #
    s_nCoef = np.polyfit(np.log10(energy*1000), np.log10(s_n/(10**4)), 4)
    return s_nCoef

def fitS_E():
    # Load file with energy and s_e (electric cross-section) data 
    O7gplot = np.loadtxt('O7gplot',float)

    # Add linear data to the end of s_e
    lineEnergyPoints = np.array([O7gplot[0,0]*16.*10**6, O7gplot[2,0]*16.*10**6]) # Use first and third energy values

    lineS_ePoints = np.array([(O7gplot[0,1]*1000*16/(6*10**2))*10**(-15), (O7gplot[2,1]*1000*16/(6*10^2))*10**-15]) # Use first and third s_e values 

    lineEnergyVector = 10.**np.arange(-5,np.log10(O7gplot[2,0]*16.*10**3),.05)*1000 # Energy array for linear data 

    lineS_eCoef=np.polyfit(np.log10(lineEnergyPoints),np.log10(lineS_ePoints),1) # linear coefficents for s_e

    lineS_eVector=10**np.polyval(lineS_eCoef,np.log10(lineEnergyVector)) # s_e array for linear data

    # Obtain s_e in eV*cm^2 and energy in eV data from O7gplot
    s_eEnergy = O7gplot[:,0]*16.*(10**6)
    s_e = (O7gplot[:,1]*1000*16/(6*10**2))*10**(-15)

    # Expanding the imported data with linear data
    s_eEnergy = np.append(lineEnergyVector, s_eEnergy)
    s_e = np.append(lineS_eVector,s_e)

    # Calculate 9th power polyfit for log10(energy [eV]) and log10(Stopping power [eV*m^2])
    s_eCoef = np.polyfit(np.log10(s_eEnergy), np.log10(s_e/(10**4)), 9)
    return s_eCoef

def fitN2():
    # load density profile
    height, n2Density = np.loadtxt('darrels_density_profile.txt',skiprows=1,usecols=(0,1),unpack=True)
    
    #height, n2Density = sio.loadmat('darrels_density_profile_v2.mat')
    n2_coef = np.polyfit(height*1000,np.log10(n2Density*(10**6)),4)
    return n2_coef

def getPolyVal(coef,x):
    return 10**(np.polyval(coef,np.log10(x)))

def N2Dens(coef, x):
    return 10**(np.polyval(coef,x))

def fitIonAndCharge():
    
    ion_energy,ion_area = np.loadtxt('NetIonizationData.txt',skiprows = 1,usecols=(0,1),unpack=True)
    
    charge_energy,charge_area = np.loadtxt('NetChargeExchange.txt', skiprows = 1, usecols = (0,1),unpack=True)
    
    return np.polyfit(ion_energy,ion_area,5),np.polyfit(charge_energy,charge_area,6)
    
    
