#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Sep 19 14:05:04 2024

@author: wkmills
"""

import numpy as np 

# Notice the flipped order of meshgrid input and output
x_range = np.array([1, 3, 5])
y_range = np.array([7, 11, 13])
ky, kx = np.meshgrid(y_range, x_range)  

print(kx[:,0]) 
print(ky[0,:]) 

ls = []
for r in range(3):
    row = []
    for c in range(3):
        row.append(kx[r,c] * ky[r,c] * np.array([-1, 0, 1]))
    ls.append(row)

a = np.array(ls) 

target_k = (3, 13) # Given as (x,y) 

# This doesn't have to be silly because we already flipped the meshgrid input and output 
a[np.where(kx[:,0] == target_k[0])[0][0], np.where(ky[0,:] == target_k[1])[0][0] , :] 


