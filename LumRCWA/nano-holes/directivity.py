#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Mar  3 13:33:33 2026

@author: quadrupole
"""
import sys, os, psutil 

sys.path.append("/opt/lumerical/v241/api/python") 
sys.path.append(os.path.dirname(__file__)) 
import lumapi

import numpy as np
#import h5py
#import matplotlib.pyplot as plt
import gc
import time 
import traceback 


def mem(): # Helper function for tracking where the memory usage spikes 
    return psutil.Process(os.getpid()).memory_info().rss / 1e9

def setup(params, sim, QW_z_limits, wavelengths): 
    
    # Clear anything left over from previous runs 
    sim.switchtolayout() 
    sim.selectall()
    sim.delete() 
    
    # Build the metasurface 
    # z=0 is the bottom of the sapphire, so all z coordinates are > 0 
    sim.addrcwa() # Add rcwa region 
    sim.set("name", "RCWA")
    
    sim.addrcwafieldmonitor() 
    sim.set("name", "monitor") 
    
    configuration = [
        ("RCWA", ( 
            ("x min", 0), ("x max", params['period'][0]),
            ("y min", 0), ("y max", params['period'][1]),
            ("z min", 0), ("z max", sum(params['layer_thicknesses'][:])),
            # set excitation wavelength here (M-point broadband is an option) 
            ("frequency points", len(wavelengths)), ("custom frequency samples", c/wavelengths),
            # Set RCWA interfaces 
            ("interface absolute positions", np.array([sum(params['layer_thicknesses'][0:i]) for i in range(1, params['layer_count'])])), 
            # Set the number of k vectors to use in the RCWA solver 
            ("max number k vectors", params['Fourier_N']), 
            # Record index to find n_sapp 
            ("report index", True),
            )), 
        ("monitor", (
            ("x min", 0), ("x max", params['period'][0]), ("x number points", params['QW_xy_mesh']),
            ("y min", 0), ("y max", params['period'][1]), ("y number points", params['QW_xy_mesh']),
            ("z min", QW_z_limits[0]), ("z max", QW_z_limits[1]), ("z number points", params['QW_z_mesh']),
            # Only need to record Ex and Ey 
            ("output Ex", True), ("output Ey", True),("output Ez", False), 
            ("output Hx", False), ("output Hy", False), ("output Hz", False),
            ))
        ]
    
    for i in range(params['layer_count']):
        
        if not params['layer_is_etched'][i]: 
            sim.addrect()
            sim.set("name", params['layer_names'][i]) 
            configuration.append(
                (params['layer_names'][i], ( 
                    ("x min", 0), ("x max", params['period'][0]),
                    ("y min", 0), ("y max", params['period'][1]),
                    ("z min", sum(params['layer_thicknesses'][0:i])), ("z max", sum(params['layer_thicknesses'][0:i+1])),
                    ("material", params['layer_materials'][i])
                    ))
                )
            
        else: # If the layer is etched 
            # First, add the full solid slab
            sim.addrect()
            sim.set("name", params['layer_names'][i])
            configuration.append(
                (params['layer_names'][i], (
                    ("x min", 0), ("x max", params['period'][0]),
                    ("y min", 0), ("y max", params['period'][1]),
                    ("z min", sum(params['layer_thicknesses'][0:i])),
                    ("z max", sum(params['layer_thicknesses'][0:i+1])),
                    ("material", params['layer_materials'][i])
                ))
            )
            # Then subtract each hole (material = "etch" = air/void)
            for h in range(params['hole_count']):
                sim.addcircle()
                sim.set("name", f"{params['layer_names'][i]}, hole {h}")
                configuration.append(
                    (f"{params['layer_names'][i]}, hole {h}", (
                        ("x", params['hole_centers_x'][h]),
                        ("y", params['hole_centers_y'][h]),
                        ("radius", params['hole_diameters'][h] / 2),
                        ("z min", sum(params['layer_thicknesses'][0:i])),
                        ("z max", sum(params['layer_thicknesses'][0:i+1])),
                        ("material", "etch")
                    ))
                )
    
    for obj, parameters in configuration:
       for name, value in parameters:
           sim.setnamed(obj, name, value)
           
    return 

def QW_xy(params, x, y):
    # Takes params, x, & y as arguments and returns Boolean (whether or not x,y is in QW emitting region) 
    # Checked 2026-03-24 
    nonemitting_thickness = 0.020e-6 
    
    if (params['hole_count'] == 0):
             return True 
         
    for h in range(params['hole_count']):
        dx = x - params['hole_centers_x'][h]
        dy = y - params['hole_centers_y'][h]
        dist = np.sqrt(dx**2 + dy**2)
        # Exclude the hole itself plus the non-emitting skin on its inner wall
        if dist <= params['hole_diameters'][h] / 2 + nonemitting_thickness:
            return False

    return True

def I(params, sim, angles):
    
    # Set table of excitation angles, then get sim results
    #return sim # For troubleshooting
    sim.switchtolayout()
    sim.setnamed("RCWA", "incident angle", "table") 
    sim.setnamed("RCWA", "incident angle table", angles) 
    print(f"Before solve: {mem():.2f} GB")
    sim.save("/tmp/rcwa_session.fsp") # Throwaway file. Saving prevents a lot of issues when running solver loops with lumapi 
    sim.run() 
    print(f"After solve: {mem():.2f} GB")
    
# This is where the memory explosion happens -- Claude suggests calculating I_z_lambda_angle in Lumerical and then just returning I_z_lambda_angle to Python 
    script = """
        Es = getresult("monitor","Es");
        Ep = getresult("monitor","Ep");
        #matlabsave("/home/quadrupole/Dropbox/research/computation/Lumerical/RCWA/meta-LED_optimization/LumRCWA/fields.mat",Es,Ep);
        """
    sim.eval(script) 
    sim.switchtolayout() 
    # Skip matlabsave and instead use sim.getv() to bring variables over to python 
    Es = sim.getv("Es")
    #sim.eval("Es.Ex = []") 
    sim.eval("Es = 0;") # To free up the memory in Lumerical 
    Ep = sim.getv("Ep") 
    sim.close()
    gc.collect() 
    
    # indices are x, y, z, lambda, angle, vector-components 
    # Loop over x, y -> if x,y in QWxy: calculate I_z -> set column of I_ky_ky_z to I_z 
    I_z_lambda_angle = np.zeros((params['QW_z_mesh'], params['wavelength_points'], len(angles))) 
    
    for xi in range(len(Es['x'])):
        for yi in range(len(Es['y'])):
            if QW_xy(params, Es['x'][xi], Es['y'][yi]):
                I_z_lambda_angle += np.abs(Es['Ex'][xi,yi,:,:,:])**2 + np.abs(Es['Ey'][xi,yi,:,:,:])**2
                I_z_lambda_angle += np.abs(Ep['Ex'][xi,yi,:,:,:])**2 + np.abs(Ep['Ey'][xi,yi,:,:,:])**2 
                # Sum intensity of in-plane componants 
    
    return I_z_lambda_angle #, z_range 

def RCWA_sim(params, kx_range, ky_range, QW_z_limits, wavelengths):
    # Initialize 
    I_z_lambda_angle = np.zeros((len(kx_range), len(ky_range), params['QW_z_mesh'])) 
    
    sim = lumapi.FDTD(hide=False) 
    setup(params, sim, QW_z_limits, wavelengths) 
    
    # Find n_sapp, or approximate effective n_sapp, using index result from RCWA object 
# =============================================================================
#     # This works because the Al2O3 model used by Lumerical is real and isotropic
#     sim.switchtolayout() 
#     sim.setnamed("RCWA", "angle phi", 0) 
#     sim.setnamed("RCWA", "angle theta", 0) 
#     sim.run() 
#     n_sapp = np.real(sim.getresult("RCWA", "index")['index_z'][0,0,0,0])  
# =============================================================================
    # But this is even faster    
    n_sapp = np.real(sim.getindex(params['layer_materials'][0], c/wavelengths[0])[0][0])
    

    # Convert kx_range and ky_range to a table of angles to feed the RCWA solver 
    # No need to loop over polarization; s- and p-pol are automatically calculated simultaneously 
    angles = [] 
    index_map = []  # (kxi, kyi) for each angle, used for unfolding data later 
    for kxi in range(len(kx_range)):
        for kyi in range(len(ky_range)):
            kx = kx_range[kxi]
            ky = ky_range[kyi]
            
            if (np.sqrt(kx**2 + ky**2)) <= NA: # inside the measureable circle 
                k0 = np.sqrt(kx**2 + ky**2)
                polar = np.rad2deg(np.arcsin(k0 / n_sapp)) 
                if k0 == 0: 
                    azimuthal = 0
                else: 
                    if ky != 0:
                        azimuthal = np.rad2deg(np.sign(ky) * np.arccos(kx/k0)) 
                    else: 
                        azimuthal = np.rad2deg(1 *  np.arccos(kx/k0))
                angles.append([polar,azimuthal])
                index_map.append((kxi, kyi)) 
    
    I_z_lambda_angle = I(params, sim, np.array(angles))
    sim.close() # Superfluous, but safe 
    gc.collect() 
    
    # Now I just need to unpack I_z_lambda_angle into a kx-by-ky matrix 
    I_kx_ky_z_lambda = np.zeros((len(kx_range),len(ky_range)) + np.shape(I_z_lambda_angle)[:-1])
    for i, (kxi, kyi) in enumerate(index_map):
        I_kx_ky_z_lambda[kxi,kyi,:,:] = I_z_lambda_angle[:,:,i]
        #heatmap[kxi, kyi, :] = E[..., i]
            
    return I_kx_ky_z_lambda

def FoM(params, queue=None):
    
    # Make the array of evenly-spaced wavelengths and the associated weights 
    if params['wavelength_points'] == 1:
        wavelengths = np.array([params['wavelength_center']]) 
    else: 
        wavelengths = np.linspace(params['wavelength_center']-params['wavelength_FWHM']/2, params['wavelength_center']+params['wavelength_FWHM']/2, params['wavelength_points'])
    
    # Make kx and ky range 
    k_inplane = np.linspace(-NA, +NA, params['k_mesh'])
    # Insert the exact target angle into kx, ky since the directivity might be too narrow to easily capture otherwise 
    target_k = params['target_k'] 
    kx_range = np.insert(k_inplane, np.searchsorted(k_inplane, target_k[0]), target_k[0]) 
    ky_range = np.insert(k_inplane, np.searchsorted(k_inplane, target_k[1]), target_k[1]) 
    
    # Stephen says p-GaN is usually 50-300 nm thick 
    QW_z_limits = sum(params['layer_thicknesses'][:-2]) - np.array([0.300e-6, 0.050e-6]) 
    
    I_kx_ky_z_lambda = RCWA_sim(params, kx_range, ky_range, QW_z_limits, wavelengths)
    
    # Apply weights to the different wavelengths calculated 
    def Gaussian_lineshape(x, x0, FWHM):
        return np.array(1 * np.exp(-(x-x0)**2 / (FWHM**2 / (4*np.log(2))))) 
    wavelength_weights = Gaussian_lineshape(wavelengths, params['wavelength_center'], params['wavelength_FWHM']) 
    I_kx_ky_z = np.zeros(np.shape(I_kx_ky_z_lambda)[:-1])
    for i in range(len(wavelengths)):
        I_kx_ky_z += I_kx_ky_z_lambda[:,:,:,i] * wavelength_weights[i] 
    
    # Calculate dual-pol directivity as a function of QW depth 
    I_target_z = I_kx_ky_z[np.where(kx_range == target_k[0])[0][0], np.where(ky_range == target_k[1])[0][0], :]
    I_avg_z = np.trapezoid( np.trapezoid(I_kx_ky_z, x = kx_range, axis = 0) / (kx_range[-1] - kx_range[0]), x = ky_range, axis = 0) / (ky_range[-1] - ky_range[0])
    I_avg_z *= np.size(I_kx_ky_z[:,:,0]) / np.size(np.nonzero(I_kx_ky_z[:,:,0])[0]) # To account for all the I=0 values at k|| > NA; verified 2024-11-05 in S4 version 
    D_z = I_target_z / I_avg_z
    # This ends up the same as z_range = sim.getresult("monitor", "Es")['z']
    z_range = np.linspace(QW_z_limits[0], QW_z_limits[1], params['QW_z_mesh']) 
    
    # mQWs 
    # Find the QW positions that have the most directional LDOS. 
    # Allow for 7-15 nm between QWs (Claude AI), but use uniform spacing for 
    # zi_spacing is the number of indices you can move in z_range without stepping > 15 nm
    max_zi_spacing = int(15e-9 // (z_range[1] - z_range[0])) 
    min_zi_spacing = int(7e-9 // (z_range[1] - z_range[0])) + 1
    
    def find_max_D(D_z, min_zi_spacing, max_zi_spacing):
        best_D = 0.0
        best_depths = None

        for N in range(min_zi_spacing, max_zi_spacing+1):
            if N * (len(params['QW_relative_intensities']) - 1) >= params['QW_z_mesh']:
                continue

            for i in range(params['QW_z_mesh'] - N * (len(params['QW_relative_intensities']) - 1)):
                d = 0 #D_z[i] + D_z[i+N] + D_z[i+2*N] 
                indices = [] 
                indiv_d= []

                for j in range(len(params['QW_relative_intensities'])): 
                    idx = i + j * N
                    d += D_z[idx] * params['QW_relative_intensities'][j] 
                    indices.append(idx) 
                    indiv_d.append(D_z[idx] * params['QW_relative_intensities'][j]) 

                if d > best_D:
                    best_D = d 
                    best_result = {
                        "best D" : float(d), 
                        "indices" : tuple(indices), 
                        "individual directivities" : tuple(indiv_d), 
                        "individual depths" : tuple(round(1e6*float(z_range[i]-params['layer_thicknesses'][0]),3) for i in indices), 
                        "QW spacing" : z_range[N] - z_range[0] 
                            }
        return best_result
    
    best_results = find_max_D(D_z, min_zi_spacing, max_zi_spacing) 
    print(f"Highest directivity is D={best_results['best D']:.3f}, with QWs at depths {' um, '.join(f'{n:.3f}' for n in best_results['individual depths'])} um (from sapphire)")
     
    # For multithreading module 
    if queue != None: queue.put((best_results['best D'], best_results['individual depths']))
    
    return (best_results['best D'], best_results['individual depths']) 

NA = 1.3
#from dOpt import min_mesa_width 
min_wall_width = 50e-9 
#minimum_ribbon_thickness = 0.050e-6 
c = 3e8 # Speed of light 

# =============================================================================
# #These didn't work:  
# params = {'Fourier_N': 50, 
# 	     'wavelength_center': 4.8e-07, 
# 	     'wavelength_FWHM': 2e-08, 
# 	     'wavelength_points': 1, 
# 	     'QW_xy_mesh': 90, 
# 	     'QW_z_mesh': 55, 
# 	     'k_mesh': 24, 
# 	     'layer_count': 5, 
# 	     'layer_names': ['sapp', 'uniform_GaN', 'etched_GaN', 'ITO', 'air'], 
# 	     'layer_thicknesses': [1e-06, 2.738e-06, 8.19e-07, 1.2e-07, 1e-06], 
# 	     'layer_materials': ['Al2O3 - Palik', 'GaN - custom', 'GaN - custom', 'ITO - custom', 'etch'], 
# 	     'layer_is_etched': [False, False, True, True, False], 
# 	     'QW_relative_intensities': [0.45, 0.33, 0.22], 
# 	     'ribbon_count': 3, 
# 	     'notch_count': 3, 
# 	     'target_k': (0, 0), 
# 	     'period': [6.89e-07, 6.63e-07], 
# 	     'ribbon_centers': [1.25e-07, 2.91e-07, 5.73e-07], 
# 	     'ribbon_widths': [1.4e-07, 6.5e-08, 1.44e-07], 
# 	     'notch_centers': [1.07e-07, 3.41e-07, 5.35e-07], 
# 	     'notch_widths': [1.11e-07, 8.5e-08, 1.05e-07]
# 	     }
# =============================================================================

# =============================================================================
# # Test device 
# params = {
#           'Fourier_N' : 50, # N = 50 recommended by Claude after convergence_test, 2026-04-14 
#           'wavelength_center' : 480e-9, 
#           'wavelength_FWHM' : 20e-9, 
#           'wavelength_points' : 1, 
#           'QW_xy_mesh' : 90, # It would seem 90 is the bare minimum, based on a period of 1.5 * 540 nm, a min_mesa_width of 50 nm, and a non-emitting thickness of 20 nm
#           'QW_z_mesh': 55, # N=55 gives at most 20 nm (~wavelength/10 in GaN) between sampled points 
#           'k_mesh': 24, # Should be an even number. Memory contraint is such that k=28 is about as high as you can go when Fourier=50, xy=90, and z=55. But k=24 is much faster and results in roughly the same outcome. 
#           'layer_count' : 5, 
#           'layer_names' : ['sapp', 'uniform_GaN', 'etched_GaN', 'ITO', 'air'], # reciprocity plane waves are incident from first layer 
#           'layer_thicknesses' : [1e-6, 1e-6, 1e-6, 120e-9, 1e-6], # GaN thicknesses can be variable param 
#           'layer_materials' : ["Al2O3 - Palik", "GaN - custom", "GaN - custom", 'ITO - custom', 'etch'],
#           'layer_is_etched' : [False, False, True, True, False], # whether or not to etch through each layer to make the ribbons 
#           'QW_relative_intensities' : [0.45, 0.33, 0.22], # relative intensities of QWs, suggested by Claude AI
#           'ribbon_count' : 3, # number of nanoribbons to etch 
#           'notch_count' : 3, 
#           'target_k' : (0, 0) # (kx, ky) 
#           # The params below will be incorporated into 'var' as fixed or range parameters, then passed to FoM in evaluate() 
#           , 'period' : [0.540e-6, 0.540e-6], # um 
#           'ribbon_centers' : [0.065e-6, 0.215e-6, 0.435e-6], 
#           'ribbon_widths' : [0.050e-6, 0.070e-6, 0.110e-6], 
#           'notch_centers' : [0.065e-6, 0.215e-6, 0.435e-6], 
#           'notch_widths' : [0.050e-6, 0.070e-6, 0.110e-6], 
#           }
# =============================================================================

# =============================================================================
# k_inplane = np.linspace(-NA, +NA, params['k_mesh'])
# # Insert the exact target angle into kx, ky since the directivity might be too narrow to easily capture otherwise 
# target_k = params['target_k'] 
# kx_range = np.insert(k_inplane, np.searchsorted(k_inplane, target_k[0]), target_k[0]) 
# ky_range = np.insert(k_inplane, np.searchsorted(k_inplane, target_k[1]), target_k[1]) 
# 
# tic = time.time()
# try: 
#     D, QWz = FoM(params) 
# except Exception as e:
#     print(type(e))
#     traceback.print_exc()
#     D, QWz = (0,(0,))
# toc = time.time()
# print(f"Took {toc-tic:.0f} seconds with Fourier N = {params['Fourier_N']} and k mesh = {params['k_mesh']}.")
# =============================================================================

# =============================================================================
# # Prasad's unpatterned thin film 
# params['layer_count'] = 3
# params['layer_names'] = ['sapp', 'uniform_GaN', 'air']
# params['layer_materials'] = ["Al2O3 - Palik", 'GaN - custom', 'etch']
# params['layer_thicknesses'] = [1e-6, 1.45e-6, 1e-6] 
# params['layer_epsilons'] = [1.77**2, 2.23**2, 1.00] 
# params['layer_is_etched'] = [False, False, False] 
# params['QW_layer'] = 1 
# params['period'] = [10e-6, 10e-6]
# params['k_mesh'] = 6
# params['Fourier_N'] = 5
# params['ribbon_count'] = 0
# params['notch_count'] = 0 
# params['wavelength_center'] = 530e-9 
# params['QW_z_min_max'] = [1e-6 + 1.35e-6, 1e-6 + 1.35e-6]
# params['QW_z_mesh'] = 1 
# #QW_depth = 1.35 # measured from first layer (sapp.); 0.672 is optimum for thin-film of thickness 0.996
# #directivity.FoM_2d(params_2d)
# =============================================================================
           
# =============================================================================
# # Larry's dual-pol optimized metasurface
# params['wavelength_center'] = 540e-9 
# params['wavelength_points'] = 1 
# params['layer_count'] = 3 
# params['layer_names'] = ['sapp', 'GaN', 'air'] # reciprocity plane waves are incident from first layer 
# params['layer_thicknesses'] = [1e-6, 0.996e-6, 1e-6] # GaN thickness can be variable param 
# params['layer_materials'] = ["Al2O3 - Palik", "GaN - custom", "etch"]
# params['layer_is_etched'] = [False, True, False] # whether or not to etch through each layer to make the ribbons 
# params['ribbon_count'] = 3 # number of nanoribbons to etch 
# params['notch_count'] = 0 
# params['target_k'] = (0, 0) # (kx, ky) 
# # The params below will be incorporated into 'var' as fixed or range parameters, then passed to FoM in evaluate() 
# params['period'] = [0.540e-6, 0.540e-6] # um 
# params['ribbon_centers'] = [0.065e-6, 0.215e-6, 0.435e-6] 
# params['ribbon_widths'] = [0.050e-6, 0.070e-6, 0.110e-6]  
# params['notch_centers'] = [] 
# params['notch_widths'] = [] 
# 
# QW_z_limits = [params['layer_thicknesses'][0] + params['layer_thicknesses'][1] - (0.120e-6 + 0.003e-6),
#                params['layer_thicknesses'][0] + params['layer_thicknesses'][1] - (0.120e-6 + 0.003e-6)]
# params['QW_z_mesh'] = 25 
# params['QW_xy_mesh'] = 50 # params['k_mesh'] = 24 
# #params_2d['notch_depths'] = []
# #QW_depth = 0.996 - (0.120 + 0.003) # measured from first layer (sapp.)
# =============================================================================

        
# =============================================================================
# # Test device
# params['ribbon_count'] = 3
# params['notch_count'] = 3
# params['ribbon_centers'] = [0.250e-6, 0.500e-6, 0.750e-6]
# params['ribbon_widths'] = [0.220e-6, 0.170e-6, 0.110e-6] 
# params['notch_centers'] = [0.250e-6, 0.500e-6, 0.750e-6] 
# params['notch_widths'] = [0.220e-6, 0.170e-6, 0.110e-6] 
# params['notch_depths'] = [0.100e-6, 0.080e-6, 0.060e-6] 
# QW_z_limits = sum(params['layer_thicknesses'][:-2]) - np.array([0.300e-6, 0.050e-6]) 
# =============================================================================


# =============================================================================
# # Convergence testing 
# Fourier_N_range = [50]
# xy_mesh_range = [90] #[110, 125, 140, 155] #[50, 65, 80, 95]
# z_mesh_range = [55] # [25] #[25, 40, 55, 70] # Claude's interpretation can't be trusted because of the way QW depth is parallelized in the FoM 
# k_mesh_range = [16, 20, 22, 24, 26, 28, 30, 32]
# 
# results = np.zeros((len(Fourier_N_range), 
#                     len(xy_mesh_range), 
#                     len(z_mesh_range), 
#                     len(k_mesh_range),
#                     2
#                     ))
# 
# for Fi in range(len(Fourier_N_range)):
#     for xyi in range(len(xy_mesh_range)): 
#         for zi in range(len(z_mesh_range)):
#             for ki in range(len(k_mesh_range)):
#                 
#                 params['Fourier_N'] = Fourier_N_range[Fi] 
#                 params['QW_xy_mesh'] = xy_mesh_range[xyi] 
#                 params['QW_z_mesh'] = z_mesh_range[zi] 
#                 params['k_mesh'] = k_mesh_range[ki] 
#                 
#                 #k_inplane = np.linspace(-0.99 * n_measure, +0.99 * n_measure, params['reciprocity_N'])
#                 k_inplane = np.linspace(-NA, +NA, params['k_mesh'])
#                 # Insert the exact target angle into kx, ky since the directivity might be too narrow to easily capture otherwise 
#                 target_k = params['target_k'] 
#                 kx_range = np.insert(k_inplane, np.searchsorted(k_inplane, target_k[0]), target_k[0]) 
#                 ky_range = np.insert(k_inplane, np.searchsorted(k_inplane, target_k[1]), target_k[1]) 
#                 
#                 
#                 tic = time.time()
#                 try: 
#                     D, QWz = FoM(params) 
#                 except Exception: 
#                     print(Exception) 
#                     D, QWz = (0,0)
#                 toc = time.time()
#                 
#                 
#                 results[Fi, xyi, zi, ki,0] = D
#                 results[Fi, xyi, zi, ki,1] = toc - tic 
#                 with open("results.txt", "w", encoding="utf-8") as f:
#                     f.write(f"Fourier N = {Fourier_N_range[Fi]} \nxy mesh = {xy_mesh_range[xyi]} \nz mesh = {z_mesh_range[zi]} \nk mesh ={k_mesh_range[ki]}")
#                     f.write(f"D = {D} \nTook {toc-tic} seconds to calculate.")
#                     f.write("\n") 
#                 
#                 print(f"Fourier N = {Fourier_N_range[Fi]} \nxy mesh = {xy_mesh_range[xyi]} \nz mesh = {z_mesh_range[zi]} \nk mesh ={k_mesh_range[ki]}")
#                 print(f"Took {toc-tic} seconds to calculate FoM.")
# =============================================================================

# =============================================================================
# plt.imshow((I_kx_ky_z_lambda[:,:,0].T)/np.max(I_kx_ky_z[:,:,0]), extent = [kx_range[0], kx_range[-1], ky_range[0], ky_range[-1]], cmap = 'inferno')
# plt.xlabel('$k_x$/$k_0$')
# plt.ylabel('$k_y$/$k_0$')
# plt.title(f'dual-pol emissions from Larry\'s case 5 \nQW depth {QW_z_min_max[0]*1e6-1:.3f} um ')
# plt.colorbar() 
# plt.show()
# =============================================================================
