#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Jul 22 16:24:04 2024

@author: wkmills
"""

import time
from FoM import directionalityFoM
from makeS4structure import structure1ribbon 

x = (0.20, 0.25) 

tic = time.time() 
D0 = directionalityFoM(sim=structure1ribbon(x, Fourier_N=10), QWregion=(0.120+0.200,0.120+0.300), TargetAngles=(0,0), N=10)
toc= time.time() 

D60 = directionalityFoM(sim=structure1ribbon(x, Fourier_N=10), QWregion=(0.120+0.200,0.120+0.300), TargetAngles=(60,0), N=10)
Dm60 = directionalityFoM(sim=structure1ribbon(x, Fourier_N=10), QWregion=(0.120+0.200,0.120+0.300), TargetAngles=(-60,0), N=10)
Dm180 = directionalityFoM(sim=structure1ribbon(x, Fourier_N=10), QWregion=(0.120+0.200,0.120+0.300), TargetAngles=(60,180), N=10)


ellapsed = toc - tic 

print(ellapsed) 

print("D(0) = "+str(D0)) 
print("D(60) = "+str(D60)) 
print("D(-60) = "+str(Dm60)) 
print("D(60,180) = "+str(Dm180)) 
