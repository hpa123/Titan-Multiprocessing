# -*- coding: utf-8 -*-


import numpy as np

def EandB(r,bfield,efield):
    boxRes = 128750.*2
    r = r/boxRes
    
    dr = np.zeros((2,3),float)
    dr[0,:] = .5 + (r-np.floor(r))
    dr[1,:] = .5 - (r-np.floor(r))
    Bc = np.zeros((8,3))
    Ec = Bc.copy()
    Bw = Bc.copy()
    Ew = Bc.copy()
    
    for i in range(3):
        if (r[i]-np.floor(r[i])) > .5:
            dr[0,i] = dr[0,i] - 1
            dr[1,i] = dr[1,i] + 1
    
    shift = [.5,.5,.5]
    mm=0
    for i in range(2):
        for j in range(2):
            for k in range(2): 
                corner = np.floor(r+shift) + 1
                corner = corner.astype(int)
                
                Bc[mm, :] = [bfield[corner[0],corner[1],corner[2],0],bfield[corner[0],corner[1],corner[2],1],bfield[corner[0],corner[1],corner[2],2]]              
                Ec[mm, :] = [efield[corner[0],corner[1],corner[2],0],efield[corner[0],corner[1],corner[2],1],efield[corner[0],corner[1],corner[2],2]]
                
                Bw[mm,:] = Bc[mm,:]*(dr[i,0]*dr[j,1]*dr[k,2])               
                Ew[mm,:] = Ec[mm,:]*(dr[i,0]*dr[j,1]*dr[k,2])
                
                mm += 1
                shift[2] = -1*shift[2]
            shift[1] = -1*shift[1]
        shift[0] = -1*shift[0]
  
    Bw=np.squeeze(Bw)
    Ew=np.squeeze(Ew)
    
    B = 0.
    E = 0.
    for i in range(8):
        B += Bw[i,:]
        E += Ew[i,:]
    return E, B
    
