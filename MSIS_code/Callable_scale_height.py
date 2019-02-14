#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Fri Oct 26 12:44:15 2018

@author: tojo5760\

Wrapping scale height into callable function
"""

from pyglow.pyglow import Point
from datetime import datetime
import matplotlib.pyplot as plt
import numpy as np
import math
import matplotlib

def Scale_Height(lat, lon, dn, altitude):
    m = altitude
    Densfrac = np.zeros(len(m))
    
    Hp= np.zeros(len(m))#Pressure scale heights
    Mbar = np.zeros(len(m))#Mean molecular mass
    Tn = np.zeros(len(m))#Neutral temperature
    He = np.zeros(len(m))
    N2 = np.zeros(len(m))
    O1 = np.zeros(len(m))
    Hp_mbar = np.zeros(len(m))
    Hp_he = np.zeros(len(m))
    Hp_n2 = np.zeros(len(m))
    Hp_o1 = np.zeros(len(m))
    
    k = 1.3806*10**(-23) #Boltzmann's constant
    g0 = 9.807 #m/s^2
    R = 6373 #Radius of earth in km
    Av = 6.022141*10**23 #Avogadro's Constant
    n = 0
    rho = np.zeros(len(m))
    
    
    #Iterate through altitudes in range m
    for x in m:
        pt = Point(dn, lat, lon, x)
        result = pt.run_msis()
        
        #Knudsen number estimate
        n_dens_tot = result.nn['HE']+result.nn['N2']+result.nn['O2'] \
            +result.nn['AR']#+result.nn['H']+result.nn['O']+result.nn['N']# number density per cm^3
    
        g = g0*R**2/(R+x)**2
        m_dens1 = (result.nn['HE']*.004003/Av+result.nn['N2']*.0280134/Av\
                   +result.nn['O2']*.032/Av+result.nn['AR']*.039948/Av)#+result.nn['N']\
                   #*.01401/Av+result.nn['H']*.001/Av+result.nn['O']*.016/Av)#kg/cm^3
        Mbar[n] = m_dens1*Av*1000/(n_dens_tot) #kg/kmol
        Tn[n] = result.Tn_msis
    
        #Follows form of kT/mg for pressure scale height
        meanmass = m_dens1/n_dens_tot #kg/molecule
        Hp_mbar[n] = k*Tn[n]/(meanmass*g)/1000 #Units of km
        Hp_he[n] = k*Tn[n]/(.004003/Av*g)/1000
        Hp_n2[n] = k*Tn[n]/(.0280134/Av*g)/1000
        Hp_o1[n] = k*Tn[n]/(.016/Av*g)/1000
    
        He[n] = result.nn['HE']*.004003/Av #kg/cm^3
        N2[n] = result.nn['N2']*.0280134/Av
        O1[n] = result.nn['O']*.016/Av
    
        n = n+1
        
    F107 = result.f107
    H_tn = np.zeros(len(m))
    H_mass = np.zeros(len(m))
    
    #Finds gradient with three point gradient technique
    for i in range(0,n,1):
        if i==0: #First Point
            coeff1 = (2.0*m[1]-m[2]-m[3])/((m[1]-m[2]*1.0)*(m[1]-m[3]*1.0))
            coeff2 = (2.0*m[1]-m[1]-m[3])/((m[2]-m[1]*1.0)*(m[2]-m[3]*1.0))
            coeff3 = (2.0*m[1]-m[1]-m[2])/((m[3]-m[1]*1.0)*(m[3]-m[2]*1.0))
    
            H_tn[1] = 1.0/Tn[1]*(Tn[1]*coeff1+Tn[2]*coeff2+Tn[3]*coeff3)
            H_mass[1] = -1/Mbar[1]*(Mbar[1]*coeff1+Mbar[2]*coeff2\
                  +Mbar[3]*coeff3)
    
        if i==n-1: #Last point
            coeff1 = (2.0*m[i]-m[i-1]-m[i])/((m[i-2]-m[i-1]*1.0)*(m[i-2]-m[i]*1.0))
            coeff2 = (2.0*m[i]-m[i-2]-m[i])/((m[i-1]-m[i-2]*1.0)*(m[i-1]-m[i]*1.0))
            coeff3 = (2.0*m[i]-m[i-2]-m[i-1])/((m[i]-m[i-2]*1.0)*(m[i]-m[i-1]*1.0))
    
            H_tn[i] = 1.0/Tn[i]*(Tn[i-2]*coeff1+Tn[i-1]*coeff2+Tn[i]*coeff3)
            H_mass[i] = -1.0/Mbar[i]*(Mbar[i-2]*coeff1+Mbar[i-1]*coeff2\
                  +Mbar[i]*coeff3)
    
        else: #Middle Points
            coeff1 = (2*m[i]-m[i]-m[i+1])/((m[i-1]-m[i]*1.0)*(m[i-1]-m[i+1]*1.0))
            coeff2 = (2*m[i]-m[i-1]-m[i+1])/((m[i]-m[i-1]*1.0 )*(m[i]-m[i+1]*1.0))
            coeff3 = (2*m[i]-m[i-1]-m[i])/((m[i+1]-m[i-1]*1.0)*(m[i+1]-m[i]*1.0))
    
            H_tn[i] = 1/Tn[i]*(Tn[i-1]*coeff1+Tn[i]*coeff2+Tn[i+1]*coeff3)
            H_mass[i] = -1/Mbar[i]*(Mbar[i-1]*coeff1+Mbar[i]*coeff2\
                  +Mbar[i+1]*coeff3)
    
    H_temp = 1./H_tn
    H_temp_he = H_temp/.64
    H_mass = 1./H_mass
    
    lnHe = np.log(He)
    lnN2 = np.log(N2)
    lnO1 = np.log(O1)
    H_he_star = np.zeros(len(m))
    H_n2_star = np.zeros(len(m))
    H_o1_star = np.zeros(len(m))
    for i in range(0,n,1):
        if i==0:
            H_he_star[i] = -(lnHe[2]-lnHe[1])/(m[2]-m[1])
            H_n2_star[i] = -(lnN2[2]-lnN2[1])/(m[2]-m[1])
            H_o1_star[i] = -(lnO1[2]-lnO1[1])/(m[2]-m[1])
        if i==len(m)-1:
            H_he_star[i] = -(lnHe[i]-lnHe[i-1])/(m[i]-m[i-1])
            H_n2_star[i] = -(lnN2[i]-lnN2[i-1])/(m[i]-m[i-1])
            H_o1_star[i] = -(lnO1[i]-lnO1[i-1])/(m[i]-m[i-1])
        else:
            H_he_star[i] = -(lnHe[i+1]-lnHe[i-1])/(m[i+1]-m[i-1])
            H_n2_star[i] = -(lnN2[i+1]-lnN2[i-1])/(m[i+1]-m[i-1])
            H_o1_star[i] = -(lnO1[i+1]-lnO1[i-1])/(m[i+1]-m[i-1])
    
    H_he_star = 1./H_he_star
    H_n2_star = 1./H_n2_star
    H_o1_star = 1./H_o1_star
    
    H_he_diff = 1./(1./H_temp_he+1./Hp_he)
    H_n2_diff = 1./(1./H_temp+1./Hp_n2)
    H_o1_diff = 1./(1./H_temp+1./Hp_o1)
    H_tot = 1./(1./H_temp+1./Hp_mbar+1./H_mass)
    
    return (H_tot, H_he_star, H_n2_star, H_o1_star)#Return scale heights in units (km)


def Mean_free_path(lat, lon, dn, altitude, scale):
    m = altitude
    Dens1frac = np.zeros(len(m))
    H1 = np.zeros(len(m))#Pressure scale heights
    M1 = np.zeros(len(m))#Mean molecular mass
    T1 = np.zeros(len(m))#Neutral temperature
    Kn1 = np.zeros(len(m))#Knudsen number for characteristic length
    Kn1At = np.zeros(len(m))#Atmospheric knudsen number
    collision_freq = np.zeros(len(m))
    k = 1.3806*10**(-23) #Boltzmann's constant
    g0 = 9.81 #m/s^2
    R = 6373 #Radius of earth in km
    Rgas = 8.31446 #J/(mol K)
    Av = 6.022141*10**23 #Avogadro's Constant
    n = 0
    rho = np.zeros(len(m))
    Te1=np.zeros(len(m))
    amu_kg = 1.660539*10**(-27) #kg per amu
    
    M_He = 4.0026*amu_kg
    M_O2 = 31.998*amu_kg
    M_N2 = 28.014*amu_kg
    M_AR = 39.948*amu_kg
    D_He = 260*scale
    D_O2 = 346*scale
    D_N2 = 364*scale
    D_AR = 340*scale
    
    pass_value = np.zeros(len(m))
    pass_value2 = np.zeros(len(m))
    
    #Iterate through altitudes in range m
    for x in m:
        pt = Point(dn, lat, lon, x)
        result = pt.run_msis()
        
        #Knudsen number estimate
        n_dens_tot = result.nn['HE']+result.nn['N2']+result.nn['O2'] \
            +result.nn['AR']+result.nn['H']+result.nn['O']+result.nn['N']# number density per cm^3
            
        n_dens_tot_new = result.nn['HE']+result.nn['N2']+result.nn['O2']+result.nn['AR'] #Added Sep. 18 2018
        Te1[n] = n_dens_tot_new
        d_avg = (result.nn['HE']/n_dens_tot*D_He+result.nn['N2']/n_dens_tot*D_N2\
            +result.nn['O2']/n_dens_tot*D_O2+result.nn['AR']/n_dens_tot*D_AR)\
            *10**(-12)+(result.nn['H']+result.nn['O']+result.nn['N'])/n_dens_tot*5.046*10**(-10)# average kinetic diameter of 
            # molecules in atmosphere. Note H and O are ignored due to no data on
            # their kinetic diameters and the breakdown of the model with ions 
            # (hard sphere collisions of molecules)
            
        d_avg = (result.nn['HE']/n_dens_tot_new*D_He+result.nn['N2']/n_dens_tot_new*D_N2\
            +result.nn['O2']/n_dens_tot_new*D_O2+result.nn['AR']/n_dens_tot_new*D_AR)\
            *10**(-12)
            
        
        lam = 1/(math.sqrt(2)*math.pi*d_avg**2*n_dens_tot_new*10**(6))
            # lam = 1/(pi*Nv*D^2*sqrt(2))
        L = 1 #in meters. Chosen characteristic length
        Kn1[n] = lam/L
        
        M_avg = (result.nn['HE']/n_dens_tot_new*M_He+result.nn['N2']/n_dens_tot_new*M_N2\
            +result.nn['O2']/n_dens_tot_new*M_O2+result.nn['AR']/n_dens_tot_new*M_AR) #kg/molecule
       
        Molecular_speed = math.sqrt(8*Rgas*result.Tn_msis/(math.pi*M_avg*Av))
        collision_freq[n] = Molecular_speed/lam
        
        pass_value[n] = Molecular_speed
        pass_value2[n] = d_avg
        
        n = n+1
    
    return (Kn1, collision_freq, pass_value2, pass_value)

