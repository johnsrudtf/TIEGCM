#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Wed Mar  7 13:41:14 2018

@author: tojo5760
"""

import matplotlib.pyplot as plt
import numpy as np
import math
import matplotlib

data = np.loadtxt("He_Tn_Correlation_400km_pdrag.txt",delimiter=',')

x = np.linspace(-143,143,287)
y = np.linspace(-71,71,143)
X, Y = np.meshgrid(x,y)
plt.figure()
levels = np.linspace(np.amin(data),np.amax(data),100)
myplot = plt.contourf(X, Y, data,levels,cmap='jet')
cont = plt.contour(X,Y,data,10,colors='k')
cbar = plt.colorbar(myplot, format='%.2e')
cbar.ax.set_ylabel('Correleation Coefficient')
plt.xticks(np.arange(-143.,143.,26.))
plt.title('Normalized He Tn Correlation plot 400km with Ion Drag')
plt.xlabel('Longitude shift')
plt.ylabel('Latitude shift')
plt.show()