#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Jun  7 13:14:19 2024

@author: wkmills
"""

#import S4 
import numpy as np 
import makeS4structure as makeS4 
import matplotlib.pyplot as plt
#import multiprocessing as mp 
import gc 

NA = 1.3 # Limit of objective collection 

''' DirectionalityFoM function takes in a simulation object (with already 
defined layers and frequency) and a valid depth range where the QW could be placed
and returns the best directionality and the associated QW depth. It operates 
via reciprocity, i.e., by setting excitation planewaves at various angles and 
measuring the resultant field intensity throughout the material where the QW 
layer could be placed. The directionality at each QW depth is then calculated, 
and the best option (and associated depth) is returned. '''

def FoM(dim, # 1-d or 2-d  
        params, # dict of the chosen and fixed parameters 
        queue = None): # multiprocessing.Queue used to send information back to parent function 
    '''
        params, # Numpy array of physical parameters for making the structure 
        Fourier_N, # Integer number of Fourier coefficients to use in S4 calculations 
        QWrange, # Duple of the minimum and maximum z range to scan for best QW location, measured from the GaN/sapp interface  
        TargetAngle, # Desired angle of emission 
        rec_N): # Integer number of points to use when discretizing theta, phi, x, y, and z 
    '''
        
    match dim:
        case "1-d": # Just use the functions in the 1-d folder (2025-01-08) 
            return 0 
# =============================================================================
#             result = FoM_1d(params)  
#             if queue != None: queue.put(result) 
#             return result 
# =============================================================================
        case "2-d":
            result = FoM_2d(params) 
            if queue != None: queue.put(result)
            return result 
        case _:
            print('Please specify the dimensionality of the simulation as \'1-d\' or \'2-d\' ')
            return 0 

# =============================================================================
#     if mp.get_start_method(True) != 'spawn':
#         mp.set_start_method('spawn')
#         
#     if __name__ == '__main__' or 'directivity':
#         p = mp.Pool(processes = 1) #Just do one process for now, so its easier to handle returned values
#         match dim:
#             case "1-d":
#                 value = p.map(FoM_1d, [params]) 
#                 p.close() 
#                 return value 
#             case "2-d":
#                 value = p.map(FoM_2d, [params])
#                 p.close() 
#                 return value 
#             case _:
#                 print('Please specify the dimensionality of the simulation as \'1-d\' or \'2-d\' ')
#                 p.close() 
#                 return 0 
#     else:
#         print('Current module name is ' + __name__ + ' not \'__main__\' ')
#         return 0 
# =============================================================================

# =============================================================================
# def FoM_1d(params):
#     # Make the sim structure by calling into makeS4structure.py functions 
#     sim = makeS4.struct_1d(params) 
#     # Calculate FoM by calling into Ds_z and Dp_z (modified to take list of tuples defining xy integration range as input)
#     
#     QW_depth = 0
#     for i in range(params['QW-layer'] + 1):
#         QW_depth += params['layer-thicknesses'][i]
#     QW_z_range = np.linspace(QW_depth - 2*params['wavelength']/np.sqrt(params['layer-epsilons'][params['QW-layer']][0][0]), QW_depth-0.010, params['reciprocity-N']) 
#     
#     D_z = {'s':None, 'p':None, 'dual':None}
#     
#     if len(params['polarization']) == 1:
#         pol = params['polarization'][0]
#         D_z[pol] = D_z_1d(sim, params, QW_z_range, pol) 
#     
#     params['polarization'] = 's' 
#     Ds_z = D_z_1d(sim, params, QW_z_range) 
#     params['polarization'] = 'p' 
#     Dp_z = D_z_1d(sim, params, QW_z_range) 
#     print(sim.GetFieldsOnGrid(QW_depth, (10,10), 'Array') )
#     
#     print('s-directivity: ' + str(max(Ds_z)) + '@ QW depth ' + str(QW_z_range[np.argmax(Ds_z)]))
#     print('p-directivity: ' + str(max(Dp_z)) + '@ QW depth ' + str(QW_z_range[np.argmax(Dp_z)]))
# 
#         
#     Dsp_z = Dp_z + Ds_z - np.abs(Dp_z - Ds_z) # combined p-pol and s-pol FoM that Larry used 
#     plt.plot(QW_z_range, Dsp_z) 
#     
#     # sim is passed by reference, so this should work to return memory to Python 
#     del sim 
#     gc.collect() 
#     
#     # Return the best directionality and the associated QW depth 
#     return np.max(Dsp_z), QW_z_range[np.argmax(Dsp_z)]  
# =============================================================================

def FoM_2d(params):
    # Make the sim structure by calling into makeS4structure.py functions 
    sim = makeS4.struct_2d(params) 
    
# =============================================================================
#     'Watch out if you try to re-use this code with the new 5-layer, stop-etched stack' 
#     QW_depth = 0
#     for i in range(params['QW-layer'] + 1):
#         QW_depth += params['layer-thicknesses'][i]
#     QW_z_range = np.linspace(QW_depth - 2*params['wavelength']/np.sqrt(params['layer-epsilons'][params['QW-layer']][0][0]), QW_depth-0.010, params['reciprocity-N']) 
# =============================================================================
    QW_z_range = [1.685] # Set QW at a fixed depth, measured from sapp. 
    
    D_z = {'s' : None, 'p' : None, 'dual' : None}
    
    if len(params['polarization']) == 1: # If only optimizing for one polarization 
        
        pol = params['polarization'][0]
        D_z[pol] = D_z_2d(sim, params, QW_z_range, pol) 
        max_z_index = np.argmax(D_z[pol])
        
        print(pol + '-directivity: ' + str(D_z[pol][max_z_index]) + ' @ QW depth ' + str(QW_z_range[max_z_index]))
        
        # This should free up memory for Python because sim is passed by reference 
        del sim
        gc.collect() 
        
        return D_z[pol][max_z_index], QW_z_range[max_z_index] 
    
    else:
# =============================================================================
#         for pol in params['polarization']:
#             D_z[pol] = D_z_2d(sim, params, QW_z_range, pol) 
#         D_z['dual'] = (D_z['s'] * D_z['p']) / (D_z['s'] + D_z['p']) # combined p-pol and s-pol FoM that Larry used 
# =============================================================================
        D_z['dual'] = D_z_2d(sim, params, QW_z_range) # Combined directivity from both polarizations 
        max_z_index = np.argmax(D_z['dual']) 
    
        #print('s-directivity: ' + str(D_z['s'][max_z_index]) + ' @ QW depth ' + str(QW_z_range[max_z_index]))
        #print('p-directivity: ' + str(D_z['p'][max_z_index]) + ' @ QW depth ' + str(QW_z_range[max_z_index]))
        print('dual-directivity: '+ str(D_z['dual'][max_z_index]) + ' @ QW depth ' + str(QW_z_range[max_z_index]))
        
        # This "should" free up memory for Python because sim is passed by reference 
        # Killing the mp process makes more of an impact though 
        del sim
        gc.collect() 
        
        # Return the best directionality and the associated QW depth 
        return D_z['dual'][max_z_index], QW_z_range[max_z_index]  

def rms(*args):
    squaresum = 0
    n = len(args) 
    for i in range(n):
        squaresum += args[i]**2
    return np.sqrt(squaresum/n) 


# =============================================================================
# def D_z_1d(sim = None, 
#       params = None,  
#       QW_z_range = None 
#       ):
#     
#     # Returns a 1-d array of directivity as a function of QW depth 
#     
#     if (params == None) or (sim== None): 
#         print('Please include the parameters and a simulation as arguments')
#         return 0 
#     
#     n_measure = np.sqrt(np.real(rms(params['layer-epsilons'][0][0][0], params['layer-epsilons'][0][1][1], params['layer-epsilons'][0][2][2])))
#     k_inplane = np.linspace(-0.99 * n_measure, +0.99 * n_measure, params['reciprocity-N']) 
#     # Insert the exact target angle into k_inplane, since the directivity might be too narrow to easily capture otherwise 
#     target_k = np.sin(params['target-angle'])/n_measure 
#     k_inplane = np.insert(k_inplane, np.searchsorted(k_inplane, target_k), target_k) 
#     # Note that now, k_inplane is 1 value longer than reciprocity-N 
#     I_k_z = []
#     for k in k_inplane: 
#         I_k_z.append(emission_k_1d(sim, params, QW_z_range, k)) 
#     I_target_z = I_k_z[np.where(k_inplane == target_k)[0][0]] 
#     I_avg_z = np.trapz(I_k_z, x = k_inplane, axis = 0) / (k_inplane[-1] - k_inplane[0]) # Use trapz because uneven spacing in k_inplane 
#     
#     return I_target_z / I_avg_z 
# =============================================================================

# =============================================================================
# def plot_at_QWz_1d(sim, params, QW_z):
#     
#     # Repeat of D_z_1d, but at fixed QW depth and then making a plot 
#     
#     if (params == None) or (sim== None): 
#         print('Please include the parameters and a simulation as arguments')
#         return 0 
#     
#     n_measure = np.sqrt(np.real(rms(params['layer-epsilons'][0][0][0], params['layer-epsilons'][0][1][1], params['layer-epsilons'][0][2][2])))
#     k_inplane = np.linspace(-0.99 * n_measure, +0.99 * n_measure, params['reciprocity-N']) 
#     # Insert the exact target angle into k_inplane, since the directivity might be too narrow to easily capture otherwise 
#     target_k = np.sin(params['target-angle'])/n_measure 
#     k_inplane = np.insert(k_inplane, np.searchsorted(k_inplane, target_k), target_k) 
#     # Note that now, k_inplane is 1 value longer than reciprocity-N 
#     I_k = []
#     for k in k_inplane: 
#         I_k.append(emission_k_1d(sim, params, [QW_z], k)) 
#         
#     # Normalize 
#     I_k = np.array(I_k)
#     I_k /= np.max(I_k) 
#
#     # Apodization correction 
#      
#     
#     I_target = I_k[np.where(k_inplane == target_k)[0][0]][0] 
#     I_avg = np.trapz(I_k, x = k_inplane, axis = 0)[0] / (k_inplane[-1] - k_inplane[0])
#     D = I_target / I_avg 
#     
#     plt.plot(k_inplane, np.array(I_k), '-o', label = 'D = ' + str(round(D,2)))
#     plt.xlabel('$k_x$')
#     plt.ylabel('Intensity (a.u.)')
#     plt.legend() 
#     plt.show() 
#     
#     return D
# =============================================================================
        

def D_z_2d(sim, params, QW_z_range):
    
    # Returns a 1-d array of directivity as a function of QW depth 
    
    if (params == None) or (sim== None): 
        print('Please include the parameters and a simulation as arguments')
        return 0 
    
    if np.size(params['layer-epsilons'][0]) == 1:
        n_measure = np.sqrt(np.real(params['layer-epsilons'][0])) 
    else: n_measure = np.sqrt(np.real(rms(params['layer-epsilons'][0][0][0], params['layer-epsilons'][0][1][1], params['layer-epsilons'][0][2][2])))
    
    k_inplane = np.linspace(-NA, +NA, params['reciprocity-N'])
    # Insert the exact target angle into kx, ky since the directivity might be too narrow to easily capture otherwise 
    target_k = params['target-k'] 
    kx_range = np.insert(k_inplane, np.searchsorted(k_inplane, target_k[0]), target_k[0]) 
    ky_range = np.insert(k_inplane, np.searchsorted(k_inplane, target_k[1]), target_k[1]) 
    ky, kx = np.meshgrid(ky_range, kx_range) # See silliness-with-numpy-array-indexing.py for logic 
    
    QW_xy = QW_xy_2d(sim, params) 
    #I = {'dual':0}
    D = {}
    
    for pol in params['polarization']: 
        # Use emission_k_2d(kx, ky)
        I_kx_ky_z = []
        for r in range(len(kx_range)): 
            row = [] 
            for c in range(len(ky_range)):
                row.append(emission_k_2d(sim, params, QW_xy, QW_z_range, kx[r,c], ky[r,c], pol)) 
            I_kx_ky_z.append(row) 
        I_kx_ky_z = np.array(I_kx_ky_z)
        
        # Apodization correction 
        apod = []
        for r in range(len(kx_range)):
            row = []
            for c in range(len(ky_range)):
                if (np.sqrt(kx[r,c]**2 + ky[r,c]**2)) > NA: # If outside objective collection region  
                    row.append(1) 
                else:
                    row.append(1/np.cos(np.sqrt(kx[r,c]**2 + ky[r,c]**2) / n_measure * np.pi/2))
            apod.append(row) 
        apod = np.array(apod)     
        I_kx_ky_z *= apod[:,:,None] # Verified visually 2024-09-24 
        
        # Save intensity profile 
        #I['dual'] += I_kx_ky_z  
        I_target_z = I_kx_ky_z[np.where(kx_range == target_k[0])[0][0], np.where(ky_range == target_k[1])[0][0] ]
        I_avg_z = np.trapz( np.trapz(I_kx_ky_z, x = kx_range, axis = 0) / (kx_range[-1] - kx_range[0]), x = ky_range, axis = 0) / (ky_range[-1] - ky_range[0])
        I_avg_z *= np.size(I_kx_ky_z[:,:,0]) / np.size(np.nonzero(I_kx_ky_z[:,:,0])[0]) # To account for all the I=0 values at k|| > NA; verified 2024-11-05
        D[pol] = I_target_z / I_avg_z 
    
# =============================================================================
#     I_kx_ky_z = I['dual'] # Don't need to normalize because it won't affect I_target/I_avg 
#     I_target_z = I_kx_ky_z[np.where(kx_range == target_k[0])[0][0], np.where(ky_range == target_k[1])[0][0] ]
#     I_avg_z = np.trapz( np.trapz(I_kx_ky_z, x = kx_range, axis = 0) / (kx_range[-1] - kx_range[0]), x = ky_range, axis = 0) / (ky_range[-1] - ky_range[0])
#     I_avg_z *= np.size(I_kx_ky_z[:,:,0]) / np.size(np.nonzero(I_kx_ky_z[:,:,0])[0]) # To account for all the I=0 values at k|| > NA; verified 2024-11-05
# =============================================================================
    
    return D['s'] + D['p'] - np.abs(D['s'] - D['p']) # Larry's dual-pol D


def plot_at_QWz_2d(sim, params, QW_z, correct_apod):
    
    # Repeat of D_z_2d, but at a fixed QWz
    # Then make a plot 
    
    if (params == None) or (sim== None): 
        print('Please include the parameters and a simulation as arguments')
        return 0 
    
    if np.size(params['layer-epsilons'][0]) == 1:
        n_measure = np.sqrt(np.real(params['layer-epsilons'][0])) 
    else: n_measure = np.sqrt(np.real(rms(params['layer-epsilons'][0][0][0], params['layer-epsilons'][0][1][1], params['layer-epsilons'][0][2][2])))
    
    #k_inplane = np.linspace(-0.99 * n_measure, +0.99 * n_measure, params['reciprocity-N'])
    k_inplane = np.linspace(-NA, +NA, params['reciprocity-N'])
    # Insert the exact target angle into kx, ky since the directivity might be too narrow to easily capture otherwise 
    target_k = params['target-k'] 
    kx_range = np.insert(k_inplane, np.searchsorted(k_inplane, target_k[0]), target_k[0]) 
    ky_range = np.insert(k_inplane, np.searchsorted(k_inplane, target_k[1]), target_k[1]) 
    ky, kx = np.meshgrid(ky_range, kx_range) # See silliness-with-numpy-array-indexing.py for logic 
    
    QW_xy = QW_xy_2d(sim, params) 
    I = {'dual':0}  
    D = {} 
    
    for pol in params['polarization']:
        # Use emission_k_2d(kx, ky)
        I_kx_ky = []
        for r in range(len(kx_range)): 
            row = [] 
            for c in range(len(ky_range)):
                row.append(emission_k_2d(sim, params, QW_xy, [QW_z], kx[r,c], ky[r,c], pol)[0]) 
            I_kx_ky.append(row) 
    
        # Make an array 
        I_kx_ky = np.array(I_kx_ky)
        
        # Apodization correction 
        if correct_apod:
            apod = []
            for r in range(len(kx_range)):
                row = []
                for c in range(len(ky_range)):
                    if (np.sqrt(kx[r,c]**2 + ky[r,c]**2)) > NA: # If outside objective collection region  
                        row.append(1) 
                    else:
                        row.append(1/np.cos(np.sqrt(kx[r,c]**2 + ky[r,c]**2) / n_measure * np.pi/2))
                apod.append(row) 
            apod = np.array(apod) 
            
            I_kx_ky *= apod 
        
        # Find directivity 
        I_target = I_kx_ky[np.where(kx_range == target_k[0])[0][0], np.where(ky_range == target_k[1])[0][0] ]
        # Note that all the I = 0 values outside of sqrt(kx**2 + ky**2) = 1 are going to pull down this average 
        I_avg = np.trapz( np.trapz(I_kx_ky, x = kx_range, axis = 0) / (kx_range[-1] - kx_range[0]), x = ky_range, axis = 0) / (ky_range[-1] - ky_range[0])
        I_avg *= np.size(I_kx_ky) / np.size(I_kx_ky[np.nonzero(I_kx_ky)]) # To account for all the I=0 values at k|| > NA; verified 2024-11-05
        D[pol] = I_target / I_avg # Confirmed working 2024-10-07 
        
        # Save intensity profile for returning from function 
        #I[pol] = I_kx_ky 
        I['dual'] += I_kx_ky 
        
        # Normalize 
        I_kx_ky /= np.max(I_kx_ky) 
        
        # Transpose it so kx is on the horizontal axis 
        plt.imshow(I_kx_ky.T, extent = [kx_range[0], kx_range[-1], ky_range[0], ky_range[-1]], cmap = 'inferno')
        plt.xlabel('$k_x$/$k_0$')
        plt.ylabel('$k_y$/$k_0$')
        plt.title(pol + '-pol emissions from QW depth ' + str(QW_z) + ' um \n D = ' + str(np.round(D[pol],2)))
        plt.colorbar() 
        plt.show() 
    
    I['dual'] /= np.max(I['dual']) # Normalize 
    # Find directivity 
# =============================================================================
#     I_target = I['dual'][np.where(kx_range == target_k[0])[0][0], np.where(ky_range == target_k[1])[0][0] ]
#     # Note that all the I = 0 values outside of sqrt(kx**2 + ky**2) = 1 are going to pull down this average; we account for that below 
#     I_avg = np.trapz( np.trapz(I['dual'], x = kx_range, axis = 0) / (kx_range[-1] - kx_range[0]), x = ky_range, axis = 0) / (ky_range[-1] - ky_range[0])
#     I_avg *= np.size(I['dual']) / np.size(I['dual'][np.nonzero(I['dual'])]) # To account for all the I=0 values at k|| > NA; verified 2024-11-05
#     D = I_target / I_avg # Confirmed working 2024-10-07 
# =============================================================================
    D = D['s'] + D['p'] - np.abs(D['s'] - D['p']) # Larry's dual-pol D
    
    plt.imshow(I['dual'].T, extent = [kx_range[0], kx_range[-1], ky_range[0], ky_range[-1]], cmap = 'inferno')
    plt.xlabel('$k_x$/$k_0$')
    plt.ylabel('$k_y$/$k_0$')
    plt.title('dual-pol emissions from QW depth ' + str(QW_z) + ' um \n D = ' + str(np.round(D,2)))
    plt.colorbar() 
    plt.show() 
    return I 

# =============================================================================
# def emission_k_1d(sim, params, QW_z_range, kx):
#     # Takes a single kx value and a range of QW depths 
#     # Returns a 1-d array of the emitted intensity at the specified kx as a function of QW depth  
#     
#     
#     if params['polarization'] == 's':
#         s = 1
#         p = 0
#     elif params['polarization'] == 'p':
#         s = 0 
#         p = 1
#         
#     # Get x that only includes where the QW is emitting 
#     x = np.linspace(0, params['period'], params['reciprocity-N']) 
#     if params['ribbon-count'] == 0:
#         QW_x = x 
#     else: 
#         QW_x = np.array([]) 
#         for ri in range(params['ribbon-count']):
#             QW_x = np.append(QW_x, x[np.where((x >= params['ribbon-centers'][ri] - params['ribbon-widths'][ri]/2 + 0.020) & (x <= params['ribbon-centers'][ri] + params['ribbon-widths'][ri]/2 - 0.020))] ) 
#         #print(QW_x) # Confirmed working 2024/09/05 
#     
#     # Send in plane wave and measure total intensity within the emitting region 
#     n_measure = np.sqrt(np.real(rms(params['layer-epsilons'][0][0][0], params['layer-epsilons'][0][1][1], params['layer-epsilons'][0][2][2])))
#     polar = np.rad2deg(np.arcsin(kx / n_measure)) 
#     sim.SetExcitationPlanewave((polar, 0), s, p) 
#     I_inplane_z = np.zeros(len(QW_z_range)) 
#     for zi in range(len(QW_z_range)): 
#         for x in QW_x:
#             I_inplane_z[zi] += np.sqrt(np.sum((np.abs(sim.GetFields(x, 0, QW_z_range[zi])[0])**2)[0:2])) 
#     
#     return I_inplane_z 
# =============================================================================

def QW_xy_2d(sim, params):
    #  Get x, y that only includes where the QW is emitting 
    # Note that ribbons run along y direction 
    # Confirmed that it works 2024-09-16 
    X = np.linspace(0, params['period'][0], 2*params['reciprocity-N'])
    Y = np.linspace(0, params['period'][1], 2*params['reciprocity-N'])
    QW_xy = [] 
    nonemitting_thickness = 0.010 
    
    unpatterned = False 
    if (params['ribbon-count'] == 0) & (params['notch-count'] == 0):
        unpatterned = True 
    
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
            for n in range(params['notch-count']):
                result = result or (params['ribbon-centers'][r] + params['ribbon-widths'][r]/2 - params['notch-depths'][n] - nonemitting_thickness <= x and x <= params['ribbon-centers'][r] + params['ribbon-widths'][r]/2 + nonemitting_thickness)
                #result = result or (params['ribbon-centers'][r] + params['ribbon-widths'][r]/2 - params['notch-depths'][n] <= x and x <= params['ribbon-centers'][r] + params['ribbon-widths'][r]/2)
        return result 
    
    for x in X:
        for y in Y:
            # If (x,y) is in emitting region, append (x,y) to QW_xy 
            emitting = unpatterned 
            if ribbon_x(x):
                if notch_y(y):
                    if not notch_x(x):
                        emitting = True 
                else: 
                    emitting = True 
            if emitting: 
                QW_xy.append((x,y)) 
                
    return QW_xy 

def emission_k_2d(sim, params, QW_xy, QW_z_range, kx, ky, pol): 
    
    # Takes a single kx, ky value and a range of QW depths 
    # Returns a 1-d array of the emitted intensity at the specified kx, ky as a function of QW depth  
    
    if pol == 's':
        s = 1
        p = 0
    elif pol == 'p':
        s = 0 
        p = 1
    
    
# =============================================================================
#     #  Get x, y that only includes where the QW is emitting 
#     # Confirmed that it works 2024-09-16 
#     X = np.linspace(0, params['period'][0], params['reciprocity-N'])
#     Y = np.linspace(0, params['period'][1], params['reciprocity-N'])
#     QW_xy = [] 
#     
#     def ribbon_x(x):
#         result = False 
#         for r in range(params['ribbon-count']):
#             result = result or (params['ribbon-centers'][r] - params['ribbon-widths'][r]/2 + 0.020 <= x and x <= params['ribbon-centers'][r] + params['ribbon-widths'][r]/2 - 0.020)
#         return result 
#     
#     def notch_y(y):
#         result = False 
#         for n in range(params['notch-count']):
#             result = result or (params['notch-centers'][n] - params['notch-widths'][n]/2 - 0.020 <= y and y <= params['notch-centers'][n] + params['notch-widths'][n]/2 + 0.020)
#         return result 
#     
#     def notch_x(x):
#         result = False 
#         for r in range(params['ribbon-count']):
#             for n in range(params['notch-count']):
#                 result = result or (params['ribbon-centers'][r] + params['ribbon-widths'][r]/2 - params['notch-depths'][n] - 0.020 <= x and x <= params['ribbon-centers'][r] + params['ribbon-widths'][r]/2)
#         return result 
#     
#     for x in X:
#         for y in Y:
#             # If (x,y) is in emitting region, append (x,y) to QW_xy 
#             emitting = False 
#             if ribbon_x(x):
#                 if notch_y(y):
#                     if not notch_x(x):
#                         emitting = True 
#                 else: emitting = True 
#             if emitting: 
#                 QW_xy.append((x,y)) 
# =============================================================================
    
    # Send in plane wave and measure total intensity within the emitting region 
    if np.size(params['layer-epsilons'][0]) == 1:
        n_measure = np.sqrt(np.real(params['layer-epsilons'][0])) 
    else: n_measure = np.sqrt(np.real(rms(params['layer-epsilons'][0][0][0], params['layer-epsilons'][0][1][1], params['layer-epsilons'][0][2][2])))
    
    if (np.sqrt(kx**2 + ky**2)) > NA:
        return np.zeros(len(QW_z_range)) # Outside the k|| = 1.3 circle 
    else:
        polar = np.rad2deg(np.arcsin(np.sqrt(kx**2 + ky**2) / n_measure)) 
        if np.sqrt(kx**2 + ky**2) == 0: 
            azimuthal = 0
        else: 
            if ky != 0:
                azimuthal = np.rad2deg(np.sign(ky) * np.arccos(kx/np.sqrt(kx**2 + ky**2))) 
            else: 
                azimuthal = np.rad2deg(1 *  np.arccos(kx/np.sqrt(kx**2 + ky**2)))
        
# =============================================================================
#         azi_z = np.zeros(len(QW_z_range)) 
#         for zi in range(len(QW_z_range)):
#             azi_z[zi] = azimuthal 
#         return azi_z   
# =============================================================================
    sim.SetExcitationPlanewave((np.round(polar,3), np.round(azimuthal, 3)), s, p) 
    I_inplane_z = np.zeros(len(QW_z_range)) 
    #I_outofplane_z = np.zeros(len(QW_z_range)) 
    #I_IP_and_OP_z = np.zeros(len(QW_z_range)) 
    for zi in range(len(QW_z_range)): 
        for xy in QW_xy:
            I_inplane_z[zi] += np.sum((np.abs(sim.GetFields(xy[0], xy[1], QW_z_range[zi])[0])**2)[0:2]) 
            #I_outofplane_z[zi] += np.sum((np.abs(sim.GetFields(xy[0], xy[1], QW_z_range[zi])[0])**2)[2]) 
            #I_IP_and_OP_z[zi] += np.sum((np.abs(sim.GetFields(xy[0], xy[1], QW_z_range[zi])[0])**2)[0:3]) 
    return I_inplane_z 
