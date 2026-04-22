#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Jul 23 13:47:25 2024

@author: wkmills
"""

import S4
import matplotlib.pyplot as plt 
import numpy as np 

# =============================================================================
# # Test apodization factor 
# sim = S4.New(((1,0),(0,1)), 20) 
# sim.SetMaterial('air',1.0) 
# sim.AddLayer('empty', 1.0, 'air')
# sim.SetFrequency(1.0) 
# 
# theta = np.linspace(-90,90,100) 
# =============================================================================

# =============================================================================
# simtest = S4.New(Lattice=((1,0),(0,1)), NumBasis=10) 
# 
# simtest.SetMaterial(Name='vac', Epsilon=1) 
# simtest.SetMaterial(Name = 'sapphire', Epsilon = 1.75**2)
# 
# simtest.AddLayer(Name='top', Thickness=0, Material='vac')
# simtest.AddLayer(Name='slab', Thickness=1, Material='sapphire') 
# simtest.AddLayer('ribbons', 1, 'vac') 
# simtest.SetRegionRectangle(Layer='ribbons', Material='sapphire', Center=(-1.3,0), Angle=0, Halfwidths=(0.2,0.2))
# simtest.AddLayerCopy('bottom',0,'top') 
# simtest.SetFrequency(1/0.550) # Corresponds to free-space wavelength of 550 nm, since S4 has frequency units of inverse length  
# 
# theta = 0
# simtest.SetExcitationPlanewave((theta,0),1,1) 
# 
# f1, b1 = simtest.GetPowerFlux('top') 
# f2, b2 = simtest.GetPowerFlux('bottom') 
# 
# E, H = simtest.GetFieldsOnGrid(1.4, (100,100), 'Array') 
# E_inplane = simtest.GetFields(-.6, 0.2, 1.4 )[0][0:1]
# print('Exyz = ' + str(np.sum(np.abs(E_inplane)**2))) 
# 
# plt.imshow(np.sum(np.abs(E)**2, axis = 2))
# 
# print('Theta = ' + str(theta)) 
# print('Incident from top = ' + str(f1)) 
# print('Transmitted to bottom = ' + str(f2)) 
# print('Incident from bottom = ' + str(b2))
# print('Reflected to top = ' + str(b1)) 
# =============================================================================

# =============================================================================
# simtest = S4.New(Lattice=((0.54,0),(0,0.54)), NumBasis=90) 
# simtest.SetOptions(SubpixelSmoothing = True, DiscretizationResolution = 4) 
# 
# # add the materials 
# simtest.SetMaterial(Name = 'c-sapp(540nm)', Epsilon = (
#                         (3.0726+0.073626j, 0, 0),   # ordinary 
#                         (0, 3.0726+0.073626j, 0),    # ordinary 
#                         (0, 0, 3.0726+0.073626j) #(0, 0, 2.9890+0.069160j)    # extraordinary 
#                         ))     
# simtest.SetMaterial(Name = 'c-GaN(540nm)', Epsilon = (
#                         (2.4194**2, 0, 0),     # ordinary axis
#                         (0, 2.4194**2, 0),     # ordinary axis 
#                         (0, 0, 2.4194**2), #(0, 0, 2.3121**2)      # extraordinary axis 
#                         ))     
# simtest.SetMaterial(Name = 'vac', Epsilon = 1.00 + 0.0j )     
# 
# # add each of the layers 
# simtest.AddLayer('sapp', 0, 'c-sapp(540nm)') 
# simtest.AddLayer('GaN', 0.996, 'c-GaN(540nm)') 
# simtest.AddLayer('air', 0, 'vac') 
# 
# simtest.SetFrequency(1/0.540)  
# angles = [(14, 0), (14, 90), (14, 180), (14, 270)] 
# #angles = [(14,270), (14, 180), (14,90), (14,0)]
# for a in angles: 
#     simtest.SetExcitationPlanewave((a[0],a[1]), 1, 0) 
#     E, H = simtest.GetFieldsOnGrid(0.672, (100,100), 'Array') 
#     field = np.sum((np.abs(E)**2)[:,:,0:2], axis = 2) 
#     plt.imshow(field) 
#     print(sum(sum(field))) 
#     plt.colorbar() 
#     plt.show() 
# =============================================================================

def QW_xy_2d(sim, params):
    #  Get x, y that only includes where the QW is emitting 
    # Note that ribbons run along y direction 
    # Confirmed that it works 2024-09-16 
    X = np.linspace(0, params['period'][0], params['reciprocity-N'])
    Y = np.linspace(0, params['period'][1], params['reciprocity-N'])
    QW_xy = [] 
    
    unpatterned = False 
    if (params['ribbon-count'] == 0) & (params['notch-count'] == 0):
        unpatterned = True 
    
    def ribbon_x(x):
        result = False 
        for r in range(params['ribbon-count']):
            result = result or (params['ribbon-centers'][r] - params['ribbon-widths'][r]/2 + 0.020 <= x and x <= params['ribbon-centers'][r] + params['ribbon-widths'][r]/2 - 0.020)
        return result 
    
    def notch_y(y):
        result = False 
        for n in range(params['notch-count']):
            result = result or (params['notch-centers'][n] - params['notch-widths'][n]/2 - 0.020 <= y and y <= params['notch-centers'][n] + params['notch-widths'][n]/2 + 0.020)
        return result 
    
    def notch_x(x):
        result = False 
        for r in range(params['ribbon-count']):
            for n in range(params['notch-count']):
                result = result or (params['ribbon-centers'][r] + params['ribbon-widths'][r]/2 - params['notch-depths'][n] - 0.020 <= x and x <= params['ribbon-centers'][r] + params['ribbon-widths'][r]/2)
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

def struct_2d(params):
    S = S4.New(Lattice=((params['period'][0],0),(0,params['period'][1])), NumBasis=params['Fourier-N']) 

    for m in range(params['layer-count']):
        S.SetMaterial(Name = params['layer-materials'][m], Epsilon = params['layer-epsilons'][m]) 
    
    S.SetOptions(SubpixelSmoothing = True, DiscretizationResolution = 4)  
    
    # Make an array of the ribbon verticies 
    verts = []
    for r in range(params['ribbon-count']):
        # Begin at the origin 
        x = 0
        y = 0
        ribbon = np.array([x,y]) 
        x += params['ribbon-widths'][r] 
        ribbon = np.vstack([ribbon, [x,y]])
        for n in range(params['notch-count']): # Loop over number of notches 
             # Add the notches that make up the verticies 
             y = params['notch-centers'][n] - params['notch-widths'][n]/2 
             ribbon = np.vstack([ribbon, [x,y]]) 
             x -= params['notch-depths'][n] 
             ribbon = np.vstack([ribbon, [x,y]])
             y += params['notch-widths'][n]
             ribbon = np.vstack([ribbon, [x,y]])
             x += params['notch-depths'][n] 
             ribbon = np.vstack([ribbon, [x,y]]) 
        # Complete the ribbon 
        y = params['period'][1] 
        ribbon = np.vstack([ribbon, [x,y]])
        x = 0
        ribbon = np.vstack([ribbon, [x,y]])
        verts.append(ribbon) 
    
    #print(verts) 
    
    # add each of the layers 
    if params['ribbon-count'] == 0: # If there are no ribbons, add all the layers as specified in 'layer materials'
        for li in range(params['layer-count']):
            S.AddLayer(params['layer-names'][li], params['layer-thicknesses'][li], params['layer-materials'][li]) 
    else: # if there are ribbons 
        for li in range(params['layer-count']):
            # Add etched layers as vacuum, then define the pillars/ribbons using S.SetRegionRectangle 
            if params['layer-is-etched'][li]: # if the layer is etched 
                S.AddLayer(params['layer-names'][li], params['layer-thicknesses'][li], 'vac') 
                for ri in range(params['ribbon-count']): # then add each of the ribbons 
                    S.SetRegionPolygon(params['layer-names'][li], params['layer-materials'][li], (0, 0), 0, tuple(map(tuple, verts[ri])) ) 
            else: # if the layer is not etched, add material as specified in 'layer materials' 
                S.AddLayer(params['layer-names'][li], params['layer-thicknesses'][li], params['layer-materials'][li]) 
    
    S.SetFrequency(1/params['wavelength']) # Since S4 has frequency units of inverse length  
    
    return S 

params = {
          'Fourier-N' : 90, # convergence_test(N=90) is 90% of convergence_test(N=900) as of 2024/09/09  
          #'wavelength' : 0.540, # Caution: the refactive indicies are hardcoded below and will not change when you change this parameter 
          'wavelength' : 0.480, # MUST CHANGE INDICES AS WELL (coded after opt. loop) 
          'reciprocity-N': 26, 
          'polarization' : ['s','p'], # Array containing 's' and/or 'p'  
          'layer-count' : 4, 
          'layer-names' : ['sapp', 'GaN', 'ITO', 'air'], # reciprocity plane waves are incident from first layer 
          'layer-thicknesses' : [0, 0.996, 0.110, 0], # GaN thickness will be variable param 
          'layer-materials' : ['c-sapp(540nm)', 'c-GaN(540nm)', 'ITO(540nm)', 'vac'],
          'layer-epsilons' : [(
                       (3.0726+0.073626j, 0, 0),   # ordinary 
                       (0, 3.0726+0.073626j, 0),    # ordinary 
                       (0, 0, 2.9890+0.069160j)    # extraordinary 
                       ),(
                       (2.4194**2, 0, 0),     # ordinary axis
                       (0, 2.4194**2, 0),     # ordinary axis 
                       (0, 0, 2.3121**2)      # extraordinary axis 
                       ), 
                       3.5097+0.012418j, # ITO index 
                       1.00             # Air index 
                       ],
          'layer-is-etched' : [False, True, True, False], # whether or not to etch through each layer to make the ribbons 
          'QW-layer' : 1, # Which layer is the QW in: 0, 1, 2, ... The QW will be optimally placed within [10 nm, 2*lambda/n] of the bottom of this layer 
          'ribbon-count' : 2, # number of nanoribbons to etch 
          'notch-count' : 1, 
          'target-k' : (0, 0) # (kx, ky) 
          # The params below will be incorporated into 'var' as fixed or range parameters, then passed to FoM in evaluate() 
          , 'period' : [0.377, 0.598], # um 
          'ribbon-centers' : [0.100, 0.250], 
          'ribbon-widths' : [0.055, 0.146],  
          'notch-centers' : [0.136], 
          'notch-widths' : [0.056], 
          'notch-depths' : [0.115]
          }
                           
meta_points = QW_xy_2d(struct_2d(params), params)
for p in meta_points: 
    plt.plot(p[0], p[1], 'o')
plt.xlim([0,params['period'][0]])
plt.ylim([0,params['period'][1]])
plt.show() 

