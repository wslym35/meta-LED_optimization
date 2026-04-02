#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Jul 15 12:06:51 2024

@author: wkmills
"""

from makeS4structure import structure1ribbon 
import numpy as np 
import matplotlib.pyplot as plt 

def convergence_test(Fourier_N):  
    sim = structure1ribbon((0.125,0.250), Fourier_N)  
    sim.SetExcitationPlanewave((0,0), 1, 0) 
    N = 100 
    E,H = sim.GetFieldsOnGrid(z=0.200, NumSamples=(N,N), Format='Array')
    I_z = 0
    for i in range(N):
        for j in range(N):
            I_z += np.abs(E[i][j][0])**2 + np.abs(E[i][j][1])**2 + np.abs(E[i][j][2])**2 
    return I_z 

Fourier_N = np.arange(20,200,1) 
results = np.zeros(len(Fourier_N)) 
for i in range(len(Fourier_N)):
    results[i] = convergence_test(Fourier_N[i]) 
    print(i) 
    

plt.plot(Fourier_N, results) 

        