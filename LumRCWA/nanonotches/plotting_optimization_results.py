#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Sep 12 10:11:59 2024

@author: quadrupole
"""
import matplotlib.pyplot as plt 
import numpy as np 
import os 
import sys 
import directivity 
from pathlib import Path

current_path = Path(os.getcwd()).resolve()
results_dir = os.path.join(current_path.parent.parent, 'results') #'/home/*/Dropbox/research/computation/Lumerical/RCWA/meta-LED_optimization/results'
date = '2026-04-27' 
filename = 'data.npy'
data = np.load(os.path.join(results_dir, date, filename), allow_pickle = True) 
sys.path.append(os.path.join(results_dir, date)) # For the next line 
from best_params import opt_params 
FoM = '$D_s+D_p$'#'$D_s + D_p$ '#'- |D_s - D_p|$'#'$D_{s+p}$' #'$(D_s D_p) / (D_s + D_p)$'

def plot_optimization(data, save = False):
    
    plt.clf() 
    fig, ax = plt.subplots() 
    ti = 1 # Skip the first (label) row  
    while data[ti, 3] ==  'GenerationStep_0_Sobol':
        ti += 1
    ax.scatter(np.arange(1,ti), data[1:ti,0], marker='P', label='training') # Scatter plot of directivity, green for Sobol 
    ax.scatter(np.arange(ti,len(data[:,0])), data[ti:,0], marker='o', label='learning') # Scatter plot of directivity, red for BoTorch 
    ax.legend()#, loc = (0.65, 0.75)) 
    ax.set_box_aspect(1) 
    plt.xlabel('Trial #')
    plt.ylabel('FoM') 
    plt.title('Results of ' + date + ' optimization \nFoM = ' + FoM) 
    if save: 
        plt.savefig(os.path.join(results_dir, date, 'convergence.png'), dpi=300, bbox_inches='tight')  
    plt.show() 


plot_optimization(data, False) 
#plt.plot([i for i in range(len(data))], data[:,0],'o')
#plt.show() 

opt_params['wavelength_points'] = 1
opt_params['wavelength_FWHM'] = 100e-9
directivity.FoM(opt_params, plot=True) 
