#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Mar  3 13:33:33 2026

@author: quadrupole
"""
import sys, os

sys.path.append("/opt/lumerical/v241/api/python") 
sys.path.append(os.path.dirname(__file__)) 
import lumapi

import numpy as np
#import matplotlib.pyplot as plt
import gc

def setup(sim, QW_z_limits, wavelengths): 
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
            ("z min", 0), ("z max", sum(params['layer-thicknesses'][:])),
            # set excitation wavelength here (M-point broadband is an option) 
            ("frequency points", len(wavelengths)), ("custom frequency samples", wavelengths/c),
            # Set RCWA interfaces 
            ("interface absolute positions", np.array([sum(params['layer-thicknesses'][0:i]) for i in range(1, params['layer-count'])])), 
            # Set the number of k vectors to use in the RCWA solver 
            ("max number k vectors", params['Fourier-N']), 
            # Record index to find n_sapp 
            ("report index", True),
            )), 
        ("monitor", (
            ("x min", 0), ("x max", params['period'][0]), ("x number points", params['QW-xy-mesh']),
            ("y min", 0), ("y max", params['period'][1]), ("y number points", params['QW-xy-mesh']),
            ("z min", QW_z_limits[0]), ("z max", QW_z_limits[1]), ("z number points", params['QW-z-mesh']),
            # Only need to record Ex and Ey 
            ("output Ex", True), ("output Ey", True),("output Ez", False), 
            ("output Hx", False), ("output Hy", False), ("output Hz", False),
            ))
        ]
    
    for i in range(params['layer-count']):
        
        if not params['layer-is-etched'][i]: 
            sim.addrect()
            sim.set("name", params['layer-names'][i]) 
            configuration.append(
                (params['layer-names'][i], ( 
                    ("x min", 0), ("x max", params['period'][0]),
                    ("y min", 0), ("y max", params['period'][1]),
                    ("z min", sum(params['layer-thicknesses'][0:i])), ("z max", sum(params['layer-thicknesses'][0:i+1])),
                    ("material", params['layer-materials'][i])
                    ))
                )
            
        else: # If the layer is etched 
            for r in range(params['ribbon-count']):
                sim.addrect()
                sim.set("name", f"{params['layer-names'][i]}, ribbon {r}")
                configuration.append(
                    (f"{params['layer-names'][i]}, ribbon {r}", (
                        ("x", params['ribbon-centers'][r]), ("x span", params['ribbon-widths'][r]),
                        ("y min", 0), ("y max", params['period'][1]),
                        ("z min", sum(params['layer-thicknesses'][0:i])), ("z max", sum(params['layer-thicknesses'][0:i+1])),
                        ("material", params['layer-materials'][i])
                        ))
                    )
                for n in range(params['notch-count']):
                    # notches are carved into both sides of the ribbon, to a depth such that the remaining ribbon width is "min_mesa_width"
                    # checked 2026-03-24 
                    sim.addrect()
                    sim.set("name", f"{params['layer-names'][i]}, ribbon {r}, notch {n}, negative side")
                    configuration.append(
                        (f"{params['layer-names'][i]}, ribbon {r}, notch {n}, negative side", (
                            ("x min", params['ribbon-centers'][r] - params['ribbon-widths'][r]/2), ("x max", params['ribbon-centers'][r] - min_mesa_width/2),
                            ("y", params['notch-centers'][n]), ("y span", params['notch-widths'][n]),
                            ("z min", sum(params['layer-thicknesses'][0:i])), ("z max", sum(params['layer-thicknesses'][0:i+1])),
                            ("material", "etch")
                            ))
                        )
                    
                    sim.addrect()
                    sim.set("name", f"{params['layer-names'][i]}, ribbon {r}, notch {n}, positive side")
                    configuration.append(
                        (f"{params['layer-names'][i]}, ribbon {r}, notch {n}, positive side", (
                            ("x min", params['ribbon-centers'][r] + min_mesa_width/2), ("x max", params['ribbon-centers'][r] + params['ribbon-widths'][r]/2),
                            ("y", params['notch-centers'][n]), ("y span", params['notch-widths'][n]),
                            ("z min", sum(params['layer-thicknesses'][0:i])), ("z max", sum(params['layer-thicknesses'][0:i+1])),
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
    
    if (params['ribbon-count'] == 0) & (params['notch-count'] == 0):
             return True 
         
    def ribbon_x(x):
        result = False 
        for r in range(params['ribbon-count']):
            result = result or (params['ribbon-centers'][r] - params['ribbon-widths'][r]/2 + nonemitting_thickness <= x and x <= params['ribbon-centers'][r] + params['ribbon-widths'][r]/2 - nonemitting_thickness)
            #result = result or (params['ribbon-centers'][r] - params['ribbon-widths'][r]/2 <= x and x <= params['ribbon-centers'][r] + params['ribbon-widths'][r]/2)
        return result 
    
    def notch_y(y):
        result = False 
        for n in range(params['notch-count']):
            result = result or (params['notch-centers'][n] - params['notch-widths'][n]/2 - nonemitting_thickness <= y and y <= params['notch-centers'][n] + params['notch-widths'][n]/2 + nonemitting_thickness)
            #result = result or (params['notch-centers'][n] - params['notch-widths'][n]/2 <= y and y <= params['notch-centers'][n] + params['notch-widths'][n]/2)
        return result 
    
    def notch_x(x):
        result = False 
        for r in range(params['ribbon-count']):
            result = result or (min_mesa_width - 2*nonemitting_thickness <= np.abs(params['ribbon-centers'][r] - x) and np.abs(params['ribbon-centers'][r] - x) <= params['ribbon-widths'][r]/2 - nonemitting_thickness)
        return result 
    
    if ribbon_x(x):
        if notch_y(y):
            if not notch_x(x):
                return True 
        else: 
            return True 
    return False 

def I(sim, angles):
    
    # Set table of excitation angles, then get sim results
    #return sim # For troubleshooting
    sim.switchtolayout()
    sim.setnamed("RCWA", "incident angle", "table") 
    sim.setnamed("RCWA", "incident angle table", angles) 
    sim.run() 
    # indices are x, y, z, lambda, angle, vector-components 
# =============================================================================
#     Es = sim.getresult("monitor", "Es")['E']#[:,:,:,0,0,:]
#     Ep = sim.getresult("monitor", "Ep")['E']#[:,:,:,0,0,:]
# =============================================================================
    Es = np.stack((sim.getresult("monitor", "Es")['Ex'], sim.getresult("monitor", "Es")['Ey'] ), axis=-1)
    Ep = np.stack((sim.getresult("monitor", "Ep")['Ex'], sim.getresult("monitor", "Ep")['Ey'] ), axis=-1)
    x_range = sim.getresult("monitor", "Es")['x']
    y_range = sim.getresult("monitor", "Es")['y']
    #z_range = sim.getresult("monitor", "Es")['z']
    
    # Loop over x, y -> if x,y in QWxy: calculate I_z -> set column of I_ky_ky_z to I_z 
    I_z_lambda_angle = np.zeros((params['QW-z-mesh'], params['wavelength-points'], len(angles))) 
    for xi in range(len(x_range)):
        for yi in range(len(y_range)):
            if QW_xy(params, x_range[xi], y_range[yi]):
                I_z_lambda_angle += np.sum(np.abs(Es[xi,yi,:,:,:,0:2])**2,axis=-1) # Sum intensity of in-plane componants pumped by incident Es
                I_z_lambda_angle += np.sum(np.abs(Ep[xi,yi,:,:,:,0:2])**2,axis=-1) # Sum intensity of in-plane componants pumped by incident Ep
    
    sim.switchtolayout() # Make sure solver gets torn down to avoid license leakage 
    
    return I_z_lambda_angle #, z_range 

def RCWA_sim(params, kx_range, ky_range, QW_z_limits, wavelengths):
    # Initialize 
    I_z_lambda_angle = np.zeros((len(kx_range), len(ky_range), params['QW-z-mesh'])) 
    
    sim = lumapi.FDTD(hide=False) 
    setup(sim, QW_z_limits, wavelengths) 
    
    # Find n_sapp, or approximate effective n_sapp, using index result from RCWA object 
    sim.switchtolayout() 
    sim.setnamed("RCWA", "angle phi", 0) 
    sim.setnamed("RCWA", "angle theta", 0) 
    sim.run() 
    n_sapp = np.real(sim.getresult("RCWA", "index")['index_z'][0,0,0,0]) 
    # This works because the Al2O3 model used by Lumerical is real and isotropic 
    
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
            
    I_z_lambda_angle = I(sim, np.array(angles))
    sim.close() 
    gc.collect() 
    
    # Now I just need to unpack I_z_lambda_angle into a kx-by-ky matrix 
    I_kx_ky_z_lambda = np.zeros((len(kx_range),len(ky_range)) + np.shape(I_z_lambda_angle)[:-1])
    for i, (kxi, kyi) in enumerate(index_map):
        I_kx_ky_z_lambda[kxi,kyi,:,:] = I_z_lambda_angle[:,:,i]
        #heatmap[kxi, kyi, :] = E[..., i]
            
    return I_kx_ky_z_lambda

def FoM(params, queue=None):
    # Make the array of evenly-spaced wavelengths and the associated weights 
    if params['wavelength-points'] == 1:
        wavelengths = np.array([params['wavelength-center']]) 
    else: 
        wavelengths = np.linspace(params['wavelength-center']-params['wavelength-FWHM']/2, params['wavelength-center']+params['wavelength-FWHM']/2, params['wavelength-points'])
    
    k_inplane = np.linspace(-NA, +NA, params['k-mesh'])
    # Insert the exact target angle into kx, ky since the directivity might be too narrow to easily capture otherwise 
    target_k = params['target-k'] 
    kx_range = np.insert(k_inplane, np.searchsorted(k_inplane, target_k[0]), target_k[0]) 
    ky_range = np.insert(k_inplane, np.searchsorted(k_inplane, target_k[1]), target_k[1]) 
    
    # Stephen says p-GaN is usually 50-300 nm thick 
    QW_z_limits = sum(params['layer-thicknesses'][:-2]) - np.array([0.300e-6, 0.050e-6]) 
    
    
    I_kx_ky_z_lambda = RCWA_sim(params, kx_range, ky_range, QW_z_limits, wavelengths)
    
    # Apply weights to the different wavelengths calculated 
    def Gaussian_lineshape(x, x0, FWHM):
        return np.array(1 * np.exp(-(x-x0)**2 / (FWHM**2 / (4*np.log(2))))) 
    wavelength_weights = Gaussian_lineshape(wavelengths, params['wavelength-center'], params['wavelength-FWHM']) 
    I_kx_ky_z = np.zeros(np.shape(I_kx_ky_z_lambda)[:-1])
    for i in range(len(wavelengths)):
        I_kx_ky_z += I_kx_ky_z_lambda[:,:,:,i] * wavelength_weights[i] 
    
    # Calculate dual-pol directivity as a function of QW depth 
    I_target_z = I_kx_ky_z[np.where(kx_range == target_k[0])[0][0], np.where(ky_range == target_k[1])[0][0], :]
    I_avg_z = np.trapezoid( np.trapezoid(I_kx_ky_z, x = kx_range, axis = 0) / (kx_range[-1] - kx_range[0]), x = ky_range, axis = 0) / (ky_range[-1] - ky_range[0])
    I_avg_z *= np.size(I_kx_ky_z[:,:,0]) / np.size(np.nonzero(I_kx_ky_z[:,:,0])[0]) # To account for all the I=0 values at k|| > NA; verified 2024-11-05 in S4 version 
    D_z = I_target_z / I_avg_z
    # This ends up the same as z_range = sim.getresult("monitor", "Es")['z']
    z_range = np.linspace(QW_z_limits[0], QW_z_limits[1], params['QW-z-mesh']) 
    
    # mQWs 
    # Find the QW positions that have the most directional LDOS. Allow for 8-15 nm between QWs 
    
    
    max_z_index = np.argmax(D_z)
    
    D_max = D_z[max_z_index] 
    QW_depth = z_range[max_z_index]
    print(f"Highest directivity is {D_max:.3f} at QW depth of {(QW_depth - params['layer-thicknesses'][0])*1e9:.0f} nm (from sapphire).")
    
    # For multithreading module 
    if queue != None: queue.put((D_max, QW_depth))
    
    return D_max, QW_depth 

NA = 1.3
#from dOpt import min_mesa_width 
min_mesa_width = 50e-9 
#minimum_ribbon_thickness = 0.050e-6 
c = 3e8 # Speed of light 

# Test device 
params = {
          'Fourier-N' : 90, # convergence_test(N=90) is 90% of convergence_test(N=900) as of 2024/09/09  
          'wavelength-center' : 480e-9, 
          'wavelength-FWHM' : 20e-9, 
          'wavelength-points' : 1, 
          'QW-xy-mesh' : 50, # My guess is that 50 is the bare minimum, based on a period of ~500 nm, a min_mesa_width of 50 nm, and a non-emitting thickness of 20 nm
          'QW-z-mesh': 25, # My guess is that 25 is the bare minimum, based on a QW-depth range of 50-300 nm and a QW placement precision of 10 nm
          'k-mesh' : 26, # This is just "how fine of a heatmap do you want?" It should be an even number. 
          'layer-count' : 5, 
          'layer-names' : ['sapp', 'uniform-GaN', 'etched-GaN', 'ITO', 'air'], # reciprocity plane waves are incident from first layer 
          'layer-thicknesses' : [1e-6, 1e-6, 1e-6, 120e-9, 1e-6], # GaN thicknesses can be variable param 
          'layer-materials' : ["Al2O3 - Palik", "GaN - custom", "GaN - custom", 'ITO - custom', 'etch'],
          'layer-is-etched' : [False, False, True, True, False], # whether or not to etch through each layer to make the ribbons 
          'QW-relative-intensities' : [0.45, 0.33, 0.22], # relative intensities of QWs, suggested by Claude AI
          #'QW-spacing' : 15e-9, # distance between QWs 
          'ribbon-count' : 0, # number of nanoribbons to etch 
          'notch-count' : 0, 
          'target-k' : (0, 0) # (kx, ky) 
          # The params below will be incorporated into 'var' as fixed or range parameters, then passed to FoM in evaluate() 
          , 'period' : [1e-6, 1e-6], # um 
          'ribbon-centers' : [], 
          'ribbon-widths' : [], 
          'notch-centers' : [], 
          'notch-widths' : [], 
          }

# =============================================================================
# # Prasad's unpatterned thin film 
# params['layer-count'] = 3
# params['layer-names'] = ['sapp', 'uniform-GaN', 'air']
# params['layer-materials'] = ["Al2O3 - Palik", 'GaN - custom', 'etch']
# params['layer-thicknesses'] = [1e-6, 1.45e-6, 1e-6] 
# params['layer-epsilons'] = [1.77**2, 2.23**2, 1.00] 
# params['layer-is-etched'] = [False, False, False] 
# params['QW-layer'] = 1 
# params['period'] = [10e-6, 10e-6]
# params['k-mesh'] = 6
# params['Fourier-N'] = 5
# params['ribbon-count'] = 0
# params['notch-count'] = 0 
# params['wavelength-center'] = 530e-9 
# params['QW_z-min-max'] = [1e-6 + 1.35e-6, 1e-6 + 1.35e-6]
# params['QW-z-mesh'] = 1 
# #QW_depth = 1.35 # measured from first layer (sapp.); 0.672 is optimum for thin-film of thickness 0.996
# #directivity.FoM_2d(params_2d)
# =============================================================================
           
# =============================================================================
# # Larry's dual-pol optimized metasurface
# params['wavelength-center'] = 540e-9 
# params['wavelength-points'] = 1 
# params['layer-count'] = 3 
# params['layer-names'] = ['sapp', 'GaN', 'air'] # reciprocity plane waves are incident from first layer 
# params['layer-thicknesses'] = [1e-6, 0.996e-6, 1e-6] # GaN thickness can be variable param 
# params['layer-materials'] = ["Al2O3 - Palik", "GaN - custom", "etch"]
# params['layer-is-etched'] = [False, True, False] # whether or not to etch through each layer to make the ribbons 
# params['ribbon-count'] = 3 # number of nanoribbons to etch 
# params['notch-count'] = 0 
# params['target-k'] = (0, 0) # (kx, ky) 
# # The params below will be incorporated into 'var' as fixed or range parameters, then passed to FoM in evaluate() 
# params['period'] = [0.540e-6, 0.540e-6] # um 
# params['ribbon-centers'] = [0.065e-6, 0.215e-6, 0.435e-6] 
# params['ribbon-widths'] = [0.050e-6, 0.070e-6, 0.110e-6]  
# params['notch-centers'] = [] 
# params['notch-widths'] = [] 
# 
# QW_z_limits = [params['layer-thicknesses'][0] + params['layer-thicknesses'][1] - (0.120e-6 + 0.003e-6),
#                params['layer-thicknesses'][0] + params['layer-thicknesses'][1] - (0.120e-6 + 0.003e-6)]
# params['QW-z-mesh'] = 25 
# params['QW-xy-mesh'] = 50 # params['k-mesh'] = 24 
# #params_2d['notch-depths'] = []
# #QW_depth = 0.996 - (0.120 + 0.003) # measured from first layer (sapp.)
# =============================================================================

        
# =============================================================================
# # Test device
# params['ribbon-count'] = 3
# params['notch-count'] = 3
# params['ribbon-centers'] = [0.250e-6, 0.500e-6, 0.750e-6]
# params['ribbon-widths'] = [0.220e-6, 0.170e-6, 0.110e-6] 
# params['notch-centers'] = [0.250e-6, 0.500e-6, 0.750e-6] 
# params['notch-widths'] = [0.220e-6, 0.170e-6, 0.110e-6] 
# params['notch-depths'] = [0.100e-6, 0.080e-6, 0.060e-6] 
# QW_z_limits = sum(params['layer-thicknesses'][:-2]) - np.array([0.300e-6, 0.050e-6]) 
# =============================================================================


#k_inplane = np.linspace(-0.99 * n_measure, +0.99 * n_measure, params['reciprocity-N'])
k_inplane = np.linspace(-NA, +NA, params['k-mesh'])
# Insert the exact target angle into kx, ky since the directivity might be too narrow to easily capture otherwise 
target_k = params['target-k'] 
kx_range = np.insert(k_inplane, np.searchsorted(k_inplane, target_k[0]), target_k[0]) 
ky_range = np.insert(k_inplane, np.searchsorted(k_inplane, target_k[1]), target_k[1]) 

I_kx_ky_z_lambda = FoM(params) 

# =============================================================================
# plt.imshow((I_kx_ky_z_lambda[:,:,0].T)/np.max(I_kx_ky_z[:,:,0]), extent = [kx_range[0], kx_range[-1], ky_range[0], ky_range[-1]], cmap = 'inferno')
# plt.xlabel('$k_x$/$k_0$')
# plt.ylabel('$k_y$/$k_0$')
# plt.title(f'dual-pol emissions from Larry\'s case 5 \nQW depth {QW_z_min_max[0]*1e6-1:.3f} um ')
# plt.colorbar() 
# plt.show()
# =============================================================================
