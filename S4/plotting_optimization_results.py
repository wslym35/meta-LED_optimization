#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Sep 12 10:11:59 2024

@author: quadrupole
"""
import matplotlib.pyplot as plt 
import numpy as np 
from os import getcwd 

folder_datename = '2025-11-24' 
file_name = 'data'
data = np.load(folder_datename + '/' + file_name + '.npy', allow_pickle = True) 
FoM = '$D_{s+p}$'#'$D_s + D_p$ '#'- |D_s - D_p|$'#'$D_{s+p}$' #'$(D_s D_p) / (D_s + D_p)$'

def plot_data(data, save = False):
    
    plt.clf() 
    fig, ax = plt.subplots() 
    ti = 1 # Skip the first (label) row  
    while data[ti, 3] == 'Sobol':
        ti += 1
    ax.scatter(np.arange(1,ti), data[1:ti,0], marker='P', label='training') # Scatter plot of directivity, green for Sobol 
    ax.scatter(np.arange(ti,len(data[:,0])), data[ti:,0], marker='o', label='learning') # Scatter plot of directivity, red for BoTorch 
    ax.legend()#, loc = (0.65, 0.75)) 
    ax.set_box_aspect(1) 
    plt.xlabel('Trial #')
    plt.ylabel('FoM') 
    plt.title('Results of ' + folder_datename + ' optimization \nFoM = ' + FoM) 
    if save: 
        plt.savefig(getcwd() + '/' + folder_datename + '/convergence.png', dpi=300, bbox_inches='tight')  
    plt.show() 

plot_data(data, False) 
#plt.plot([i for i in range(len(data))], data[:,0],'o')
#plt.show() 

