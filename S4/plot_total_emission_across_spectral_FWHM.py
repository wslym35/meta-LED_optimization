#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Jan 26 15:48:07 2026

@author: quadrupole
"""

import matplotlib.pyplot as plt 
import numpy as np 

I = I_total 
kx_range = [-1.3, 1.3] 
ky_range = [-1.3, 1.3] 
QW_z = 1.685 

plt.imshow(I.T, extent = [kx_range[0], kx_range[-1], ky_range[0], ky_range[-1]], cmap = 'inferno')
plt.xlabel('$k_x$/$k_0$')
plt.ylabel('$k_y$/$k_0$')
plt.title('dual-pol emissions from QW depth ' + str(QW_z) + ' um \n Emission across spectral FWHM')
plt.colorbar() 
plt.show() 