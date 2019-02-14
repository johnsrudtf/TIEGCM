#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Tue Sep 11 10:52:23 2018

Finding collision frequency between two gas species

n1, n2 = number density of species 1 and 2 respectively
m1, m2 = mass of species 1 and 2 respectively
r1, r2 = radii of species 1 and 2 respectively
T1, T2 = temperatures of species 1 and 2 respectively

@author: tojo5760
"""
import numpy as np
import math

def collision_freq(n1,n2,m1,m2,r1,r2,T1,T2):
    
    k = 1.38065*10**(-23)
    n1 = n1*10**6 #Convert from cm^(-3) to m^(-3)
    n2 = n2*10**6
    m1 = m1*1.66054*10**(-27) #Convert from amu to kg
    m2 = m2*1.66054*10**(-27)
    r1 = r1*10**(-12) #Convert from picometers to meters
    r2 = r2*10**(-12)
    
    T_st = (m1*T2+m2*T1)/(m1+m2)
    Mu_st = (m1*m2)/(m1+m2)
    Alpha = math.sqrt((2*k*T_st)/Mu_st)
    Q_st = math.pi*(r1+r2)**2
    Omega_st = (Alpha*Q_st)/math.sqrt(4*math.pi)
    V_st = 16/3.0*(n2*m2)/(m1+m2)*Omega_st
    
    return V_st
