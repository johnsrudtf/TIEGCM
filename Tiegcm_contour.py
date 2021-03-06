#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Thu Dec  7 17:19:23 2017
Plots Contour plots created by Contour_Plots.m
@author: tojo5760
"""
import matplotlib.pyplot as plt
import numpy as np
import math
import matplotlib

#Load data file
data = np.loadtxt("He_Tn_normalized_difference_pdrag400km.txt",delimiter=',')
#data = np.log10(data)

#Load Magnetic Equator
mag_equator = np.loadtxt("Magnetic_equator_lat_lon.txt",delimiter=',')
top = mag_equator[180:360,:]
bottom = mag_equator[0:180,:]
mag_equator = np.concatenate((top,bottom),axis=0)

#Create contour plot mesh of correct size
x = np.linspace(0,143/6, 144)
y = np.linspace(-88.75,88.75,72)
X, Y = np.meshgrid(x, y)
magnetic_x = np.linspace(0,143/6,360)

#Plot contour
plt.figure()
#levels = np.linspace(1.5153e-14,6.8968999999999999e-14,100)
levels = np.linspace(np.amin(data),np.amax(data),100)
myplot = plt.contourf(X, Y, data,levels,cmap='jet')
cont = plt.contour(X,Y,data,10,colors='k')
cbar = plt.colorbar(myplot, format='%.2e')
cbar.ax.set_ylabel('Normalized Difference')
plt.plot(magnetic_x,mag_equator[:,1],'r') #Add magnetic equator line
plt.xticks(np.arange(0.,24.,3.))
plt.title('HE Tn Normalized Difference 400km UT=0.03 with Ion Drag')
plt.xlabel('Local Solar Time [hr]')
plt.ylabel('Latitude [deg]')
plt.show()
ctrmin = np.amin(data)
ctrmax = np.amax(data)
