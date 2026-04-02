#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Jun 28 10:26:21 2024

@author: wkmills

This function takes in a parameter set and returns a S4 structure that can be 
fed as the first argument to directionalityFoM 

"""

import S4 
import numpy as np 

def struct_1d(params):
    S = S4.New(Lattice=((params['period'],0),(0,params['period'])), NumBasis=params['Fourier-N']) 

    for m in range(params['layer-count']):
        S.SetMaterial(Name = params['layer-materials'][m], Epsilon = params['layer-epsilons'][m]) 
    
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
                    S.SetRegionRectangle(params['layer-names'][li], params['layer-materials'][li], (params['ribbon-centers'][ri], 0), 0, (params['ribbon-widths'][ri]/2, params['period']/2)) 
            else: # if the layer is not etched, add material as specified in 'layer materials' 
                S.AddLayer(params['layer-names'][li], params['layer-thicknesses'][li], params['layer-materials'][li]) 
    
    S.SetFrequency(1/params['wavelength']) # Since S4 has frequency units of inverse length  
    
    return S 


def struct_2d(params):
    S = S4.New(Lattice=((params['period'][0],0),(0,params['period'][1])), NumBasis=params['Fourier-N']) 
    #print("HELP!")
    for m in range(params['layer-count']):
        S.SetMaterial(Name = params['layer-materials'][m], Epsilon = params['layer-epsilons'][m]) 
    
    S.SetOptions(SubpixelSmoothing = True, DiscretizationResolution = 4)  
    
    # Make an array of the ribbon verticies 
    verts = []
    for r in range(params['ribbon-count']):
        # Begin at the origin 
        x = 0
        y = 0
        x += params['ribbon-centers'][r] - 1/2 * params['ribbon-widths'][r] 
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
        x -= params['ribbon-widths'][r] 
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
                    #S.SetRegionRectangle(params['layer-names'][li], params['layer-materials'][li], (params['ribbon-centers'][ri], 0), 0, (params['ribbon-widths'][ri]/2, params['period'][0]/2)) 
            else: # if the layer is not etched, add material as specified in 'layer materials' 
                S.AddLayer(params['layer-names'][li], params['layer-thicknesses'][li], params['layer-materials'][li]) 
    
    S.SetFrequency(1/params['wavelength']) # Since S4 has frequency units of inverse length  
    
    return S 

