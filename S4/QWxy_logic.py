#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Sep 11 15:03:14 2024

@author: wkmills
"""
# =============================================================================
# xy = (7.5, 10) 
# notches = [xy[0]==2, xy[0]==3]
# ribbon = [0<xy[0]<5 or 7<xy[0]<8] 
# result = True 
# =============================================================================
result = True 
x = 2.2 
notches = [x==2, x==3] 
ribbon = [0<x<5 or 7<x<10]
for r in ribbon: 
    for n in notches: 
        result = result and (r and not n) 

print(result) 






result = False 
xy = (2.2, 4.4) 
notches = [x==2, x==3] 
ribbon = [0 < xy[0] < 5 or 7 < xy[0] < 10]
for r in ribbon: 
    for n in notches: 
        result = result or (r and not n) 

print(result) 