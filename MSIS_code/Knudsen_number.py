#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Wed Feb 13 14:45:56 2019

@author: tojo5760
"""

#==============================================================================
# Plot Knudsen number for a spacecraft with characeristic length of 1m and for
# the atmosphere scale height. Also plot mass percentage of Helium as a 
# function of altitude. All from MSIS model. Definition of msis 
# is in pyglow.py line 290.
# For all plots, "low" and "high" refer to the solar cycle.
#==============================================================================
from pyglow.pyglow import Point
from datetime import datetime
from Collision_Frequency import collision_freq
from Callable_scale_height import Scale_Height
import matplotlib.pyplot as plt
import numpy as np
import math
import matplotlib


dn = datetime(1987, 3, 23, 9, 30)#year,month,day,hour,minute
dn2 = datetime(1981, 3, 23, 9, 30)
lat = 0.
lon = 0.
m = range(100, 600, 5) #Set altitude range and step size (km)
Dens1frac = np.zeros(len(m))
Dens2frac = np.zeros(len(m))
H1 = np.zeros(len(m))#Pressure scale heights
H2 = np.zeros(len(m))
M1 = np.zeros(len(m))#Mean molecular mass
T1 = np.zeros(len(m))#Neutral temperature
T2 = np.zeros(len(m))
M2 = np.zeros(len(m))
Kn1 = np.zeros(len(m))#Knudsen number for characteristic length
Kn2 = np.zeros(len(m))
Kn1At = np.zeros(len(m))#Atmospheric knudsen number
Kn2At = np.zeros(len(m))
k = 1.3806*10**(-23) #Boltzmann's constant
g0 = 9.81 #m/s^2
R = 6373 #Radius of earth in km
Av = 6.022141*10**23 #Avogadro's Constant
n = 0
rho = np.zeros(len(m))

Mass_O2 = 32.00
Mass_He = 4.00
Mass_N2 = 28.01
Mass_Ar = 39.95
Radius_O2 = 346/2.0 #In picometers
Radius_He = 260/2.0
Radius_N2 = 364/2.0
Radius_Ar = 340/2.0
R_gas = 8.3144598 #Universal gas constant

lam = np.zeros(len(m))
lam2 = np.zeros(len(m))
Coll_freq_avg = np.zeros(len(m))
Coll_freq_avg2 = np.zeros(len(m))
Molecular_speed1 = np.zeros(len(m))
Molecular_speed2 = np.zeros(len(m))

#Iterate through altitudes in range m
for x in m:
    pt = Point(dn, lat, lon, x)
    pt2 = Point(dn2, lat, lon, x)
    result = pt.run_msis()
    result2 = pt2.run_msis()
    
    #Knudsen number estimate
    n_dens_tot_new = result.nn['HE']+result.nn['N2']+result.nn['O2']+result.nn['AR']#Number density per cm^3
    
    M1[n] = (result.nn['HE']/n_dens_tot_new*.004003+result.nn['N2']/n_dens_tot_new*.028013\
        +result.nn['O2']/n_dens_tot_new*.032+result.nn['AR']/n_dens_tot_new*.039948) #kg/mol  
        
    Coll_freq_1 = collision_freq(result.nn['HE'], result.nn['N2'], Mass_He, Mass_N2, Radius_He, Radius_N2, result.Tn_msis, result.Tn_msis)
    Coll_freq_2 = collision_freq(result.nn['HE'], result.nn['O2'], Mass_He, Mass_O2, Radius_He, Radius_O2, result.Tn_msis, result.Tn_msis)
    Coll_freq_3 = collision_freq(result.nn['HE'], result.nn['AR'], Mass_He, Mass_Ar, Radius_He, Radius_Ar, result.Tn_msis, result.Tn_msis)
    Coll_freq_4 = collision_freq(result.nn['N2'], result.nn['HE'], Mass_N2, Mass_He, Radius_N2, Radius_He, result.Tn_msis, result.Tn_msis)
    Coll_freq_5 = collision_freq(result.nn['N2'], result.nn['O2'], Mass_N2, Mass_O2, Radius_N2, Radius_O2, result.Tn_msis, result.Tn_msis)
    Coll_freq_6 = collision_freq(result.nn['N2'], result.nn['AR'], Mass_N2, Mass_Ar, Radius_N2, Radius_Ar, result.Tn_msis, result.Tn_msis)
    Coll_freq_7 = collision_freq(result.nn['O2'], result.nn['N2'], Mass_O2, Mass_N2, Radius_O2, Radius_N2, result.Tn_msis, result.Tn_msis)
    Coll_freq_8 = collision_freq(result.nn['O2'], result.nn['AR'], Mass_O2, Mass_Ar, Radius_O2, Radius_Ar, result.Tn_msis, result.Tn_msis)
    Coll_freq_9 = collision_freq(result.nn['O2'], result.nn['HE'], Mass_O2, Mass_He, Radius_O2, Radius_He, result.Tn_msis, result.Tn_msis)
    Coll_freq_10 = collision_freq(result.nn['AR'], result.nn['N2'], Mass_Ar, Mass_N2, Radius_Ar, Radius_N2, result.Tn_msis, result.Tn_msis)
    Coll_freq_11 = collision_freq(result.nn['AR'], result.nn['O2'], Mass_Ar, Mass_O2, Radius_Ar, Radius_O2, result.Tn_msis, result.Tn_msis)
    Coll_freq_12 = collision_freq(result.nn['AR'], result.nn['HE'], Mass_Ar, Mass_He, Radius_Ar, Radius_He, result.Tn_msis, result.Tn_msis)
    Coll_freq_total = Coll_freq_1+Coll_freq_2+Coll_freq_3+Coll_freq_4+Coll_freq_5+Coll_freq_6+Coll_freq_7+Coll_freq_8+Coll_freq_9+Coll_freq_10\
        +Coll_freq_11+Coll_freq_12
    Coll_freq_avg[n] = Coll_freq_total/12
    Molecular_speed1[n] = math.sqrt((8*R_gas*result.Tn_msis)/(math.pi*M1[n]))#m/s
    lam[n] = Molecular_speed1[n]/Coll_freq_avg[n]
    
    
    n_dens_tot2_new = result2.nn['HE']+result2.nn['N2']+result2.nn['O2']+result2.nn['AR'] #Number density per cm^3
    M2[n] = (result2.nn['HE']/n_dens_tot2_new*.004003+result2.nn['N2']/n_dens_tot2_new*.028013\
        +result2.nn['O2']/n_dens_tot2_new*.032+result2.nn['AR']/n_dens_tot2_new*.039948) #kg/mol  
        
    Coll_freq_1 = collision_freq(result2.nn['HE'], result2.nn['N2'], Mass_He, Mass_N2, Radius_He, Radius_N2, result2.Tn_msis, result2.Tn_msis)
    Coll_freq_2 = collision_freq(result2.nn['HE'], result2.nn['O2'], Mass_He, Mass_O2, Radius_He, Radius_O2, result2.Tn_msis, result2.Tn_msis)
    Coll_freq_3 = collision_freq(result2.nn['HE'], result2.nn['AR'], Mass_He, Mass_Ar, Radius_He, Radius_Ar, result2.Tn_msis, result2.Tn_msis)
    Coll_freq_4 = collision_freq(result2.nn['N2'], result2.nn['HE'], Mass_N2, Mass_He, Radius_N2, Radius_He, result2.Tn_msis, result2.Tn_msis)
    Coll_freq_5 = collision_freq(result2.nn['N2'], result2.nn['O2'], Mass_N2, Mass_O2, Radius_N2, Radius_O2, result2.Tn_msis, result2.Tn_msis)
    Coll_freq_6 = collision_freq(result2.nn['N2'], result2.nn['AR'], Mass_N2, Mass_Ar, Radius_N2, Radius_Ar, result2.Tn_msis, result2.Tn_msis)
    Coll_freq_7 = collision_freq(result2.nn['O2'], result2.nn['N2'], Mass_O2, Mass_N2, Radius_O2, Radius_N2, result2.Tn_msis, result2.Tn_msis)
    Coll_freq_8 = collision_freq(result2.nn['O2'], result2.nn['AR'], Mass_O2, Mass_Ar, Radius_O2, Radius_Ar, result2.Tn_msis, result2.Tn_msis)
    Coll_freq_9 = collision_freq(result2.nn['O2'], result2.nn['HE'], Mass_O2, Mass_He, Radius_O2, Radius_He, result2.Tn_msis, result2.Tn_msis)
    Coll_freq_10 = collision_freq(result2.nn['AR'], result2.nn['N2'], Mass_Ar, Mass_N2, Radius_Ar, Radius_N2, result2.Tn_msis, result2.Tn_msis)
    Coll_freq_11 = collision_freq(result2.nn['AR'], result2.nn['O2'], Mass_Ar, Mass_O2, Radius_Ar, Radius_O2, result2.Tn_msis, result2.Tn_msis)
    Coll_freq_12 = collision_freq(result2.nn['AR'], result2.nn['HE'], Mass_Ar, Mass_He, Radius_Ar, Radius_He, result2.Tn_msis, result2.Tn_msis)
    Coll_freq_total2 = Coll_freq_1+Coll_freq_2+Coll_freq_3+Coll_freq_4+Coll_freq_5+Coll_freq_6+Coll_freq_7+Coll_freq_8+Coll_freq_9+Coll_freq_10\
        +Coll_freq_11+Coll_freq_12
    Coll_freq_avg2[n] = Coll_freq_total2/12
    Molecular_speed2[n] = math.sqrt((8*R_gas*result2.Tn_msis)/(math.pi*M2[n]))#m/s
    lam2[n] = Molecular_speed2[n]/Coll_freq_avg2[n]  
    n = n+1
    
F107_1 = result.f107  
F107_2 = result2.f107

[H_tot, H_he_star, H_n2_star, H_o1_star] = Scale_Height(lat, lon, dn, m)
[H_tot2, H_he_star2, H_n2_star2, H_o1_star2] = Scale_Height(lat, lon, dn2, m)
H_tot = H_tot*1000 #Change from km to m
H_tot2 = H_tot2*1000
Kn1At = np.true_divide(lam,H_tot)
Kn2At = np.true_divide(lam2,H_tot2)
#Plot Knudsen number as a function of altitude with charictaristic length 
#of 1 meter
#==============================================================================
fig = plt.figure()
ax = plt.subplot(111)
line1, = ax.plot(lam, m, label='1987 F10.7 = %.1f' %F107_1, color='red')
line2, = ax.plot(lam2, m, label='1981 F10.7 = %.1f' %F107_2, color='blue')
line3, = ax.plot(Kn1At, m, label='Atmospheric, F10.7 = %.1f' %F107_1, color='red', ls='dashed')
line4, = ax.plot(Kn2At, m, label='Atmospheric, F10.7 = %.1f' %F107_2, color='blue', ls='dashed')

plt.xscale('log')
plt.xlabel('Knudsen Number')
plt.ylabel('Altitude (km)')
plt.title('Knudsen Number versus Altitude')
#ax.legend(loc='upper center', bbox_to_anchor=(.5, .2 ), ncol=2)

red_patch = matplotlib.patches.Patch(color='red', label='F10.7 = %.1f' %F107_1)
blue_patch = matplotlib.patches.Patch(color='blue', label='F10.7 = %.1f' %F107_2)
black_line = matplotlib.lines.Line2D([], [], color='black', label='L = 1m')
black_dash = matplotlib.lines.Line2D([], [], ls='dashed', color='black', label='L=Hp scale height')
plt.legend(handles=[red_patch, blue_patch, black_line, black_dash], loc='upper left', bbox_to_anchor=(0, 1))
plt.grid()
plt.show()


fig2 = plt.figure()
ax = plt.subplot(111)
line1, = ax.plot(Coll_freq_avg, m, label='\nu F10.7 = %.1f' %F107_1, color='red')
line2, = ax.plot(Coll_freq_avg2, m, label='\nu F10.7 = %.1f' %F107_2, color='blue')

plt.xscale('log')
plt.xlabel('Collision Frequency')
plt.ylabel('Altitude (km)')
plt.title('Collision Frequency versus Altitude')

red_line = matplotlib.lines.Line2D([], [], color='red', label=r'$\nu$'' F10.7 = %.1f' %F107_1)
blue_line = matplotlib.lines.Line2D([], [], color='blue', label=r'$\nu$'' F10.7 = %.1f' %F107_2)
plt.legend(handles=[red_line, blue_line], loc='upper right', bbox_to_anchor=(1, 1))
plt.grid()
plt.show()

fig3 = plt.figure()
ax = plt.subplot(111)
line1, = ax.plot(H_tot/1000, m, color='red')
line2, = ax.plot(H_tot2/1000, m, color='blue')

plt.xlabel('Scale Height (km)')
plt.ylabel('Altitude (km)')
plt.title('Total Atmospheric Scale Height versus Altitude')

red_line = matplotlib.lines.Line2D([], [], color='red', label='$H_{tot}$ F10.7 = %.1f' %F107_1)
blue_line = matplotlib.lines.Line2D([], [], color='blue', label='$H_{tot}$ F10.7 = %.1f' %F107_2)
plt.legend(handles=[red_line, blue_line], loc='upper right', bbox_to_anchor=(1, .5))
plt.grid()
plt.show()

fig4 = plt.figure()
ax = plt.subplot(111)
line1, = ax.plot(Molecular_speed1, m, color='red')
line2, = ax.plot(Molecular_speed2, m, color='blue')

plt.xlabel('Molecular Speed (m/s)')
plt.ylabel('Altitude (km)')
plt.title('Mean Molecular Speed versus Altitude')

red_line = matplotlib.lines.Line2D([], [], color='red', label='S F10.7 = %.1f' %F107_1)
blue_line = matplotlib.lines.Line2D([], [], color='blue', label='S F10.7 = %.1f' %F107_2)
plt.legend(handles=[red_line, blue_line], loc='upper right', bbox_to_anchor=(1, .5))
plt.grid()
plt.show()