#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Jul 26 11:08:55 2024

@author: wkmills
"""

import os # for getcwd(), cpu_count(), maybe more 
import gc # for memory management in the for loop 
#import S4
import numpy as np 
import matplotlib.pyplot as plt 
from makeS4structure import struct_1d 

# =============================================================================
# # define a function that takes input angle and returns intensity in the material 
# def I_inplane_angle(sim, X, Y, QWz, s, p, theta_deg):
#     I = 0
#     sim.SetExcitationPlanewave((theta_deg, 0), s, p) 
#     for x in X:
#         for y in Y:
#             if np.real(sim.GetEpsilon(x, y, QWz)) > 1.0: # Only integrate field where the non-air material is 
#                 E, H = sim.GetFields(x, y, QWz) 
#                 I += np.abs(E[0])**2 + np.abs(E[1])**2 # Only integrate in-plane componants  
#                 del(E, H) # save memory 
#                 ' NEED TO ADD SOMETHING FOR IGNORING THE OUTER 20 nm OF THE QW'
#     return I
# =============================================================================

def rms(*args):
    squaresum = 0
    n = len(args) 
    for i in range(n):
        squaresum += args[i]**2
    return np.sqrt(squaresum/n) 
        

def k_emission_1d(sim, params, QW_z): 
    # Takes a sim, parameters, and QW depth (measured from first layer) and returns k_inplane and intensity @ QW depth (both numpy arrays) 
    
    #theta_deg_range = np.linspace(-89, +89, params['reciprocity N']) 
    n_emission = np.sqrt(np.real(rms(params['layer-epsilons'][params['QW layer']][0][0], params['layer-epsilons'][params['QW-layer']][1][1], params['layer-epsilons'][params['QW-layer']][2][2])))
    n_measure = np.sqrt(np.real(rms(params['layer-epsilons'][0][0][0], params['layer-epsilons'][0][1][1], params['layer-epsilons'][0][2][2])))
    k_inplane = np.linspace(-0.99 * n_measure, +0.99 * n_measure, params['reciprocity-N']-1)
    # Insert the exact target angle into k_inplane, since the directivity might be too narrow to easily capture otherwise 
    k_inplane = np.insert(k_inplane, np.searchsorted(k_inplane, np.sin(params['target-angle'])/n_measure), np.sin(params['target-angle'])/n_measure) 
    polar_range = np.rad2deg(np.arcsin(k_inplane / n_measure)) 
    #print('k_inplane min = ' + str(k_inplane_range[0]))
    intensity_k = np.zeros(len(k_inplane))  
    
    if params['polarization'] == 's':
        s = 1
        p = 0
    elif params['polarization'] == 'p':
        s = 0 
        p = 1
    
    # Get x that only includes where the QW is emitting 
    x = np.linspace(0, params['period'], params['reciprocity-N']) 
    if params['ribbon-count'] == 0:
        QW_x = x 
    else: 
        QW_x = np.array([]) 
        for ri in range(params['ribbon-count']):
            QW_x = np.append(QW_x, x[np.where((x >= params['ribbon-centers'][ri] - params['ribbon-widths'][ri]/2 + 0.020) & (x <= params['ribbon-centers'][ri] + params['ribbon-widths'][ri]/2 - 0.020))] ) 
        #print(QW_x) # Confirmed working 2024/09/05 
    
    for i in range(params['reciprocity-N']):
        polar = polar_range[i] 
        sim.SetExcitationPlanewave((polar, 0), s, p) 
        E_abs_inplane = 0
        for x in QW_x:
            E_abs_inplane += np.sqrt(np.sum((np.abs(sim.GetFields(x, 0, QW_z)[0])**2)[0:2])) 
        intensity_k[i] = E_abs_inplane 
# =============================================================================
#     X = np.linspace(-params['period']/2, +params['period']/2, params['reciprocity N'])
#     Y = np.linspace(-params['period']/2, +params['period']/2, params['reciprocity N'])
#     
#     if params['ribbon count'] == 0:
#         for i in range(len(theta_deg_range)):
#             sim.SetExcitationPlanewave((theta_deg_range[i],0), s, p)
#             #E_fields = [[sim.GetFields(x, y, QWz)[0] for x in X] for y in Y]
#             #intensity_k[i] = sum(sum(np.abs(E_fields[0])**2 + np.abs(E_fields[1])**2)) 
#             E_grid = np.abs(sim.GetFieldsOnGrid(QWz,(100,100),'Array')[0])**2 
#             intensity_k[i] = sum(sum(E_grid[...,0] + E_grid[...,1])) # Only add the in-plane componants 
#     else: 
#         print('Need to add feature for integrating only over the QW regions and omitting the outer 20 nm')
#         return 0
#      
# =============================================================================
# =============================================================================
#     # Use S4.SolveInParallel() to parallelize code a bit 
#     # Create a list of simulation clones to be used 
#     clones = list([])
#     for i in range(os.cpu_count() //2): # Use half the available cores
#         clones.append(sim.Clone()) 
#     
#     for theta_chunk = # next len(clones) values in theta_deg_range 
#         for theta in theta_chunk :
#             clones[ ].SetExcitationPlanewave((theta,0),s,p)
#         S4.SolveInParallel(params[''], 
# =============================================================================
# =============================================================================
#     def I_inplane_angle(sim, X, Y, QWz, s, p, theta_deg):
#          I = 0
#          sim.SetExcitationPlanewave((theta_deg, 0), s, p) 
#          for x in X:
#              for y in Y:
#                  if np.real(sim.GetEpsilon(x, y, QWz)) > 1.0: # Only integrate field where the non-air material is 
#                      E, H = sim.GetFields(x, y, QWz) 
#                      I += np.abs(E[0])**2 + np.abs(E[1])**2 # Only integrate in-plane componants  
#                      del(E, H) # save memory 
#                      ' NEED TO ADD SOMETHING FOR IGNORING THE OUTER 20 nm OF THE QW'
#          return I
# =============================================================================
    #print(I_angle(20)) 
    #intensity_k = mp.Pool().map(lambda theta_deg : I_inplane_angle(sim, X, Y, QWz, s, p, theta_deg), theta_deg_range, chunksize=10)
    intensity_k_appod = intensity_k * np.cos(np.deg2rad(polar_range)) 
    #print(params['polarization'] + '-pol Directivity = ' + str(np.max(intensity_k_appod) / np.average(intensity_k_appod))) 
    print(params['polarization'] + '-pol Directivity = ' + str(intensity_k[np.searchsorted(k_inplane, np.sin(params['target-angle'])/n_measure)] / np.average(intensity_k))) 
    return k_inplane, intensity_k

def convergence_test(sim, params, QW_z):
    # Takes sim and params and returns the intensity at the QW depth as excited by a s+p polarized wave at polar_angle = 15 deg 
    sim.SetExcitationPlanewave((15,0),1,1) 
    
    # Get x that only includes where the QW is emitting 
    x = np.linspace(0, params['period'], params['reciprocity-N']) 
    if params['ribbon-count'] == 0:
        QW_x = x 
    else: 
        QW_x = np.array([]) 
        for ri in range(params['ribbon-count']):
            QW_x = np.append(QW_x, x[np.where((x >= params['ribbon-centers'][ri] - params['ribbon-widths'][ri]/2 + 0.020) & (x <= params['ribbon-centers'][ri] + params['ribbon-widths'][ri]/2 - 0.020))] ) 
        #print(QW_x) # Confirmed working 2024/09/05 
    
    E_abs_inplane = 0
    for x in QW_x:
        E_abs_inplane += np.sqrt(np.sum((np.abs(sim.GetFields(x, 0, QW_z)[0])**2)[0:2])) 
    return E_abs_inplane 
    
    
def plot_k_emission(k,Is,Ip, params, save, savepath=None):
    plt.plot(k,Is, label = 's-pol') 
    plt.plot(k,Ip, label = 'p-pol') 
    plt.legend(loc = 'upper center') 
    plt.xlabel('$k_{||}$')
    plt.ylabel('Intensity (a.u.)')
    #plt.ylim([0,7]) 
    plt.title('Emission @ ' + str(params['wavelength']) + ' um into ' + params['layer-names'][0] + ' side')
    if save == True:
        #design_details = '-GaN' + str(params['GaN thickness']) + 'um-ITO' + str(params['ITO thickness']) + 'um-QWdepth' + str(params['QW depth']) + 'um-'
        name = str(params['ribbon-count']) + ' ribbons - ' + str(params['wavelength']) + ' um wavelength - '
        for li in range(params['layer-count']):
            name += str(params['layer-thicknesses'][li]) 
            name += ' um '
            name += params['layer-names'][li] 
            name += ' - ' 
        name += str(params['QW-depth']) + 'um QW depth'
        plt.savefig(savepath + name + '.png', dpi=150)  
    plt.show() 

    
params0 = {'period':1, # um 
          'Fourier N':1,
          'wavelength': 0.475, # Caution: the refactive indicies are hardcoded and will not change when you change this parameter 
          'reciprocity N':100, 
          'polarization': 'p', # 'p' or 's' 
          'QW depth': 0.120 + 0.0646, # um, measured from the first layer in 'layer names' as written below 
          'layer count': 4, 
          'layer names': ['air', 'ITO', 'GaN', 'sapp'], # plane waves are incident from first layer, unless 'orientation' is set to 'back' 
          'layer thicknesses': [0, 0.120, 5.5746, 0], 
          'layer materials': ['vac', 'ITO(475nm)', 'c-GaN(475nm)', 'c-sapp(475nm)'],
          'layer etched': [False, True, True, False], # whether or not to etch through each layer to make the ribbons 
          #'orientation': 'forw', # see note on 'layer names'; switches orientation of the stack so that plane waves are incident from the first layer listed above
          #'emission layer': 'sapp', # which layer to measure momentum-resolved emission in 
          'ribbon count': 0 # number of nanoribbons to etch 
          }

# =============================================================================
# params2 = {'period':0.540,
#             'x1': 0.133,
#             'w1': 0.145,
#             'x2': 0.410,
#             'w2': 0.190, 
#             'Fourier N': 10,
#             'GaN thickness': 5.577, 
#             'ITO thickness': 0.120, 
#             'wavelength': 0.475, 
#             'reciprocity N': 100, 
#             'polarization': 'p',
#             'QW depth': 0.42, 
#             'emission side': 'sapp',
#             'ribbon count': 2
#             }
# =============================================================================

params_Jon = {'period':1,
              'Fourier N': 10,
              'wavelength': 0.700, 
              'reciprocity N': 100,
              'polarization':'p',
              'QW depth': 0.010, 
              'layer count': 3,
              'layer names': ['substrate','perovskite','air'],
              'layer thicknesses': [0, 0.020, 0],
              'layer materials': ['quartz','PTCDA','vac'],
              'orientation': 'forw',
              'emission layer': 'substrate',
              'ribbon count': 0 
              } 

# =============================================================================
# par = params0 
# par['polarization'] = 's'
# k, Is = k_emission_1d(struct_1d(par),par) 
# par['polarization'] = 'p' 
# k, Ip = k_emission(struct_1d(par),par) 
# plot_k_emission_1d(k,Is,Ip, par, save = False, savepath = os.getcwd() + '/uLEDs-Roark/') 
# =============================================================================

params_Larry_case1 = {'period':0.540, # um 
          'Fourier-N':90, # convergence_test(N=90) is 90% of convergence_test(N=900) as of 2024/09/09  
          'wavelength': 0.540, # Caution: the refactive indicies are hardcoded and will not change when you change this parameter 
          'reciprocity-N': 200, 
          'polarization': 'p', # 'p' or 's' 
          #'QW z range': (0 + 0.996-0.070, 0.996-0.070), # um, measured from the first layer in 'layer names' as written below; given twice if QW location is known 
          'layer-count': 3, 
          'layer-names': ['sapp', 'GaN', 'air'], # plane waves are incident from first layer 
          'layer-thicknesses': [0, 0.996, 0], 
          'layer-materials': ['c-sapp(540nm)', 'c-GaN(540nm)', 'vac'],
          'layer-epsilons': [(
                       (3.0726+0.073626j, 0, 0),   # ordinary 
                       (0, 3.0726+0.073626, 0),    # ordinary 
                       (0, 0, 2.9890+0.069160j)    # extraordinary 
                       ),(
                       (2.4194**2, 0, 0),     # ordinary axis
                       (0, 2.4194**2, 0),     # ordinary axis 
                       (0, 0, 2.3121**2)      # extraordinary axis 
                       ), 
                       1.00 + 0.00j 
                       ],
          #'index of emission layer': np.sqrt(3), # Approximate index of sapphire; for k_emission.py 
          'layer-is-etched': [False, True, False], # whether or not to etch through each layer to make the ribbons 
          'QW-layer': 1, # Which layer is the QW in: 0, 1, 2, ... The QW will be optimally placed within [10 nm, 2*lambda/n] of the bottom of this layer 
          'ribbon-count': 2, # number of nanoribbons to etch 
          'ribbon-centers': [0.133, 0.410], # Ribbon centers are assumed to be at y = 0 for 1-d case
          'ribbon-widths': [0.145, 0.190], # Ribbon y-widths are assumed to be period/2 for 1-d case 
          'target-angle': 0 # 1-d simualation, so don't need (0,0) 
          } 

                           
QW_depth_Larry = 0.996 - 0.070  # measured from first layer 
QW_depth_mine = 0.9287
par = params_Larry_case1
par['polarization'] = 's'
k, Is = k_emission_1d(struct_1d(par), par, QW_depth_mine )
par['polarization'] = 'p'
k, Ip = k_emission_1d(struct_1d(par), par, QW_depth_mine ) 
plot_k_emission(k, Is, Ip, par, False)

# =============================================================================
# # Convergence testing results  
# N = [1, 10, 20, 30, 40, 50, 60, 70, 80, 90, 100, 200, 300, 400, 500, 600, 700, 800, 900]
# results = np.zeros(len(N)) 
# for ni in range(len(N)):
#     par['Fourier N'] = N[ni] 
#     print("Fourier N = " + str(par['Fourier N']))
#     intensity = convergence_test(struct_1d(par), par, QW_depth) 
#     print("Intensity = " + str(intensity))
#     results[ni] = intensity 
#  
# plt.plot(N, results) 
# plt.show() 
# =============================================================================
