# -*- coding: utf-8 -*-
"""

These are values that are used by both titan_util.py and Particle_Trajectory_Multiprocesses.py.
They are stored here to avoid inconsistencies between the files.

"""
from numpy import pi

# Particle Properties and Constants
kb = 1.38*10**(-23) 
temp = 300.
T = temp*11604.
mass = 16*1.67*10**(-27)                  # particle mass (kg)
q = 1.6*10**(-19)                         # Charge of one particle (Coulombs)
JtoEV = 6.24*10**18                       # Joules to electron-volts conversion factor

# Titan Properties
titanRad=2576000.0                        # Radius of Titan; Meters
atmosAltitude = 1600000.                  # Radial distance to which Titan's atmosphere extends from the surface (m, meters)
dectorArea = 4*pi*((titanRad+atmosAltitude)**2)    # surface area of Titan's atmosphere


# Simulation Parameters
tf = 10000
dt = .3
nsteps = int(tf/dt)
delX = 1000.
deltaR=20000. 
