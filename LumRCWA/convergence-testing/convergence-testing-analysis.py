#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Apr  9 15:43:04 2026

@author: quadrupole
"""

import numpy as np

seconds_arr     = np.load('convergence-testing/convergence_seconds.npy')
directivity_arr = np.load('convergence-testing/convergence_directivity.npy')

# Index maps
fourier_idx = {v: i for i, v in enumerate([1, 30, 40, 50, 60, 70, 80, 90])}
xy_idx      = {v: i for i, v in enumerate([50, 65, 80, 95])}
z_idx       = {v: i for i, v in enumerate([25, 40, 55, 70])}
k_idx       = {v: i for i, v in enumerate([20, 24, 28])}

# Example lookup
#seconds_arr[fourier_idx[60], xy_idx[65], z_idx[55], k_idx[20]]  # → 2117.78

