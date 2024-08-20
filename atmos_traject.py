# -*- coding: utf-8 -*-

from __future__ import division
from os import path
if not path.isfile("atmosParticles.npy"):
    print "No particles to work with; run mag_traject.py first."
else:

    import multiprocessing as mp
    from titanUtil import insideAtmos #runInfo
    import constants as c
    import time
    import numpy as np
    import sys
    from titanPlots import plots
    
    
    atmosParticles = np.load("atmosParticles.npy")
    
    allResults = []
    incidentAngles = []
    
    net_depo = np.zeros((4,len(atmosParticles[0].depo[0])))
    
    # make sure the number of runs and the number of simulataneous particles
    # work out (see runInfo in titanUtil.py)
    
    sim = len(atmosParticles)
    runs = int(len(atmosParticles)/sim)
# check to see if there are any particles incident on Titan's atmosphere

    sys.stdout.write("\nBeginning in-atmosphere calculations...")
    
#    if __name__ == "main":
    # create a multiprocessing pool and begin timing calculations
    start = time.time()
    pool = mp.Pool()
    
    # calculate in-atmosphere particle trajectories
    for i in range(runs):
        #sys.stdout.write("\n Taking particles %s to %s for run %s." % (sim*i,
        #sim*(i+1),i+1))
        
        # pass workers to pool
        tempResults = pool.map(insideAtmos,
        atmosParticles[sim*i:sim*(i+1)])
        
        # store results returned from each worker
        allResults.append(tempResults)
    
    # close pool and wait for all calculations to finish
    pool.close()
    pool.join()
    
    allResults = [p for pgroup in allResults for p in pgroup]    
    
    # store layers of atmosphere in first row of net deposition array
    net_depo[0] = atmosParticles[0].depo[0]    
    
    # add up deposition contributions from each particle
    for p in allResults:
        net_depo[1,:] += p.weight*p.depo[1]
        net_depo[2,:] += p.depo[2]
        net_depo[3,:] += p.depo[3]
    
    # convert to stuff/cm^3
    for i in range(len(net_depo[0])-1):
        net_depo[1,i] = net_depo[1,i]/((4/3)*np.pi*(net_depo[0,i+1]**3 - net_depo[0,i]**3)*(10**6))
        net_depo[2,i] = net_depo[2,i]/((4*np.pi*((net_depo[0,i+1] + net_depo[0,i])/2)**2)*(10**6))
        net_depo[3,i] = net_depo[3,i]/((4*np.pi*((net_depo[0,i+1] + net_depo[0,i])/2)**2)*(10**6))
    
    # converting meters to kilometers, shiifting altitude to have 0 = titanRad
    net_depo[0,:] = (net_depo[0,:] - c.titanRad)/1000
    
    # stop timing calculation
    end = time.time()
    
    if end-start > 60:
        print "\n Calculation completed in %s minutes."%((end-start)/60)
    else:
        print "\n Calculation completed in %s seconds."%(end-start)
        
    incidentAngles = [p.incidentAngle for p in atmosParticles]
    incEnergies = [p.incEnergy for p in atmosParticles]
    plots(1,1,1,1,1,net_depo,incidentAngles,incEnergies)
