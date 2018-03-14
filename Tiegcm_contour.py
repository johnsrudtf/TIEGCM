#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Thu Dec  7 17:19:23 2017

@author: tojo5760
"""
import matplotlib.pyplot as plt
import numpy as np
import math
import matplotlib

data = np.loadtxt("He_400km_pdrag_TEST_DELETE.txt",delimiter=',')
#data = np.log10(data)
x = np.linspace(0,143/6, 144)
y = np.linspace(-88.75,88.75,72)
X, Y = np.meshgrid(x, y)

plt.figure()
#levels = np.linspace(1.5153e-14,6.8968999999999999e-14,100)
levels = np.linspace(np.amin(data),np.amax(data),100)
myplot = plt.contourf(X, Y, data,levels,cmap='jet')
cont = plt.contour(X,Y,data,10,colors='k')
cbar = plt.colorbar(myplot, format='%.2e')
cbar.ax.set_ylabel('Log 10 O1/N2 Ratio')
plt.xticks(np.arange(0.,24.,3.))
plt.title('Log 10 O1/N2 Ratio 400km UT=0 with Ion Drag')
plt.xlabel('Local Solar Time [hr]')
plt.ylabel('Latitude [deg]')
plt.show()
ctrmin = np.amin(data)
ctrmax = np.amax(data)