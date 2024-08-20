# -*- coding: utf-8 -*-
"""

Particle trajectory and deposition rates implementing 
multiprocessing.

"""
from __future__ import division 
from titanUtil import createParticles, runInfo, outsideAtmos
import multiprocessing as mp
import sys
import time
from numpy import save

# Particle Properties 
           
def magTraject(i):
    
    numParticles = 10000         # total number of particles
    sim = 900                    # number of simultaneous trajectory calculations
                                 # sim should always be less than numParticles    
    
    allResults = []         
        
    print "Creating particles..."
    # Create all of our particles (see createParticles in titanUtil)
    particleList = createParticles(numParticles)
    
    print "Created %s particles. Running calculation now...\n"%(len(particleList))
    
    # create a multiprocessing pool and begin timing the calculation
    pool = mp.Pool()
    start = time.time()

    # Get the number of runs and the number of simultaneous particles
    sim, runs = runInfo(len(particleList),sim)    
    print sim, runs
    for i in range(runs):
        sys.stdout.write("Taking particles %s to %s for run %s. \r" % (sim*i,
        sim*(i+1),i+1))
        sys.stdout.flush()
        
        # pass our outer-atmosphere calculations to the pool, using the
        # number of simultaneous particles specified by the user
        tempResults = pool.map(outsideAtmos,
        particleList[sim*i:sim*(i+1)])
        
        # store the results from each parallel run
        allResults.append(tempResults)       
    
    # close the pool, so that no more tasks can be added
    pool.close()
    
    # wait for all processes to finish
    pool.join() 
        
    # stop timing the calculation
    end = time.time()

    # Collect particles that hit the atmosphere and save that data in a file
    atmosParticles = [p for pgroup in allResults for p in pgroup if p.hitAtmos == True]
    save("atmosParticles_run%s",atmosParticles)
    if end-start > 60:
        sys.stdout.write("\nCalculations complete in %s minutes. %s particle(s) hit the atmosphere of Titan!"\
        %((end-start)/60,len(atmosParticles)))
    else:
        sys.stdout.write("\nCalculations complete in %s seconds. %s particle(s) hit the atmosphere of Titan!"\
        %(end-start,len(atmosParticles)))
    
        
    
    
        
    
    
    


    
    
            
            
    
            
            
            
        
        
        
        
        
        






