# -*- coding: utf-8 -*-

import matplotlib.pyplot as plt
import numpy as np

def plots(EnDep,incEnHist,incAngHist,altCS,altIon,net_depo,incidentAngles,incidentEnergies):
    
    atmosLayers = net_depo[0,:] 
    enDepResults = net_depo[1,:]   
    ionResults = net_depo[2,:]
    csResults = net_depo[3,:]

    if EnDep == 1:     
        plt.figure(1)
        plt.semilogx(enDepResults,atmosLayers)
        plt.xlabel('Energy Deposition in Atmosphere (eV/cm^3*s)')
        plt.ylabel('Altitude (km)')
        plt.title('Energy Deposition Rate at Altitude Level')
        #plt.axis([np.min(np.log10(altEn[1,:][altEn[1,:] != 0])), np.max(np.log10(altEn[1,:])), 1150, 1600])
        
    if incEnHist == 1: 
        plt.figure(2)
        p2 = plt.hist(incidentEnergies,bins=10)
        plt.xlabel('Energy of O+ Ions Entering Atmosphere of Titan (eV)');
        plt.ylabel('Number of O+ Ions');
        plt.title('Energy of O+ Ions Incident Upon Atmosphere of Titan');
        plt.axis([0,np.max(p2[1]),0,np.max(p2[0])*1.2]);
    
    if incAngHist == 1:
        plt.figure(3)
        p3 = plt.hist(incidentAngles)
        plt.xlabel('Incident Angle of O+ Ions Entering Atmosphere of Titan (deg)');
        plt.ylabel('Number of O+ Ions');
        plt.title('Incident Angle of O+ Ions Entering Atmosphere of Titan');
        plt.axis([0,np.max(p3[1]),0,np.max(p3[0])*1.2]);
    
    if altCS == 1:
        plt.figure(4)
        plt.semilogx(csResults,atmosLayers)
        
        plt.xlabel('Charge Exchange Desposition (1/cm^3*s)')
        plt.ylabel('Altitude (km)')
        #plt.axis([np.min(np.log10(altEn[3,:][altEn[3,:] != 0])), np.max(np.log10(altEn[3,:])), 1150, 1600])
        plt.title('Charge Exchange Deposition in Atmosphere of Titan')
        
    if altIon == 1:
        plt.figure(5)
        plt.semilogx(ionResults,atmosLayers)
        #plt.axis([np.min(np.log10(altEn[4,:][altEn[4,:] != 0])), np.max(np.log10(altEn[4,:])), 1150, 1600])
        plt.xlabel('Ioniziation Desposition (1/cm^3*s)')
        plt.ylabel('Altitude (km)')
        plt.title('Ionization Deposition in Atmosphere of Titan')
    plt.show()
