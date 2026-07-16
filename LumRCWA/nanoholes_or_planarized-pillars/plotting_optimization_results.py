#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Sep 12 10:11:59 2024

@author: quadrupole
"""
import matplotlib.pyplot as plt 
import numpy as np 
import os 
import directivity 
from pathlib import Path

current_path = Path(os.getcwd()).resolve()
results_dir = os.path.join(current_path.parent.parent, 'results') #'/home/*/Dropbox/research/computation/Lumerical/RCWA/meta-LED_optimization/results'
date = '2026-06-18' 
data = np.load(os.path.join(results_dir, date, 'data.npy'), allow_pickle = True) 
opt_params = np.load(os.path.join(results_dir, date, 'opt_params.npy'), allow_pickle = True).item() 



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
    plt.title('Results of ' + date + ' optimization \nFoM = ' + opt_params['FoM_definition']) 
    if save: 
        plt.savefig(os.path.join(results_dir, date, 'convergence.png'), dpi=300, bbox_inches='tight')  
    plt.show() 


plot_optimization(data, save=False) 
#plt.plot([i for i in range(len(data))], data[:,0],'o')
#plt.show() 

# Tolerence testing 
# =============================================================================
# opt_params['wavelength_points'] = 5
# opt_params['wavelength_FWHM'] = 100e-9
# =============================================================================
overetch = -100e-9
opt_params['layer_thicknesses'][1] -= overetch
opt_params['layer_thicknesses'][2] += overetch 
QW_fixed = np.array([3.031, 3.041, 3.050]) * 1e-6

directivity.FoM(opt_params, plot=True, QW_fixed=QW_fixed) 
