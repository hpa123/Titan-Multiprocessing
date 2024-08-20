# -*- coding: utf-8 -*-

import numpy as np

def build_fields():    
    bfield_output = np.loadtxt('bfield_output')
    efield_output = np.loadtxt('efield_output')
    
    a = 111
    b = 89
    c = 89
    
    bfield = np.zeros((a,b,c,3),float) 
    efield = bfield.copy()
    
    mm=0
    for i in range(a):
        for j in range(b):
            for k in range(c):
                
                bfield[i,j,k,0] = bfield_output[mm,0]
                bfield[i,j,k,1] = bfield_output[mm,1]
                bfield[i,j,k,2] = bfield_output[mm,2]
                
                efield[i,j,k,0] = efield_output[mm,0]
                efield[i,j,k,1] = efield_output[mm,1]
                efield[i,j,k,2] = efield_output[mm,2]
                
                mm = mm + 1
            
    #efield and bfield are 4D arrays that give the Bx, By, Bz, Ex, Ey, and Ez
    #values in a 3D cube
    
    np.save("bfield",bfield)
    np.save("efield",efield)
