#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Sep 12 15:09:01 2024

@author: wkmills
"""

import multiprocessing as mp 
import directivity 
import sys 
    
params_2d = {'period' : [0.540, 0.540], # um 
          'Fourier-N' : 90, # convergence_test(N=90) is 90% of convergence_test(N=900) as of 2024/09/09  
          'wavelength' : 0.540, # Caution: the refactive indicies are hardcoded below and will not change when you change this parameter 
          'reciprocity-N': 6, 
          'polarization' : ['s','p'], # Array containing 's' and/or 'p'  
          'layer-count' : 3, 
          'layer-names' : ['sapp', 'GaN', 'air'], # plane waves are incident from first layer 
          'layer-thicknesses' : [0, 0.996, 0], 
          'layer-materials' : ['c-sapp(540nm)', 'c-GaN(540nm)', 'vac'],
          'layer-epsilons' : [(
                       (3.0726+0.073626j, 0, 0),   # ordinary 
                       (0, 3.0726+0.073626, 0),    # ordinary 
                       (0, 0, 2.9890+0.069160j)    # extraordinary 
                       ),(
                       (2.4194**2, 0, 0),     # ordinary axis
                       (0, 2.4194**2, 0),     # ordinary axis 
                       (0, 0, 2.3121**2)      # extraordinary axis 
                       ), 
                       1.00 + 0.00j 
                       ],
          'layer-is-etched' : [False, True, False], # whether or not to etch through each layer to make the ribbons 
          'QW-layer' : 1, # Which layer is the QW in: 0, 1, 2, ... The QW will be optimally placed within [10 nm, 2*lambda/n] of the bottom of this layer 
          'ribbon-count' : 2, # number of nanoribbons to etch 
          'notch-count' : 2, 
          'target-k' : (0, -0.25) # (kx, ky) 
          # The params below will be incorporated into 'var' as fixed or range parameters during optimization; this is just for testing FoM
          , 'ribbon-centers' : [0.161, 0.396], # Ribbon centers are assumed to be at y = 0 for 1-d case
          'ribbon-widths' : [0.147, 0.163], # Ribbon y-widths are assumed to be period/2 for 1-d case 
          'notch-centers' : [0.215, 0.429], 
          'notch-widths' : [0.049, 0.116], 
          'notch-depths' : [0.008, 0.064]
          }
    

#directivity.FoM('2-d', params_2d) 

try:
    mp.set_start_method('fork') 
except RuntimeError:
    if mp.get_start_method() == 'spawn':
        print('Somehow the start method got set to spawn. Please restart the kernel.')
        sys.exit() 
q = mp.Queue() 
p = mp.Process(target = directivity.FoM, args = ('2-d', params_2d, q,))
p.start() # Starts running the process 
print(p.join()) # Blocks other things from running until p terminates; best practice to use 
print(q.get() ) # Get output of FoM() 
# =============================================================================
# #directivity.FoM('2-d', params_2d) 
# def FoM_child(dim, params):
#     return directivity.FoM(dim, params) 
# 
# if __name__ == '__main__':
#     try: 
#         mp.set_start_method('spawn')
#     except RuntimeError: 
#         if mp.get_start_method() == 'fork':
#             print('Start method was somehow set to fork. Please restart the kernel and try again.') 
#             sys.exit()
#     
#     p = mp.Process(target = FoM_child, args=('2-d',params_2d,))
#     p.start()
#     p.join()
# =============================================================================



# =============================================================================
# def FooM(dic):
#     # All the mp stuff goes here because S4 obj isn't picklable 
#     task_queue = mp.Queue() 
#     done_queue = mp.Queue() 
#     
#     task = [(FooM_child, dic)] 
#     
#     
#     p = mp.Process(target = FooM_child)  
#     returned = p.map(FooM_child, [dic]) 
#     p.close() 
#     print(returned) 
#     
#     return 1 
#      
# 
# def FooM_child(in_q, out_q):
#     
#     string = dic['string']
#     array = dic['array']
#     return string + str(array) 
# 
# 
#         #p = mp.Process(target = ek, args = ([string, 'hello again'], array, dic, sim))
#         #p.start() 
#         #p.join() 
# 
# # See https://docs.python.org/3/library/multiprocessing.html for details 
# 
# # =============================================================================
# # def info(title):
# #     print(title)
# #     print('module name:', __name__)
# #     print('parent process:', os.getppid())
# #     print('process id:', os.getpid())
# # 
# # def f(name):
# #     info('function f')
# #     print('hello', name)
# # 
# dictionary = {'string' : 'hello world ',
#         'array' : np.array([1, 2, 3])
#         }
# 
# FooM(dictionary) 
# 
# # 
# # if __name__ == '__main__':
# #     #mp.set_start_method('spawn') 
# #     print(mp.get_start_method()) 
# #     info('main line')
# #     p = mp.Process(target=f, args=('bob',))
# #     print(p.is_alive())
# #     p.start()
# #     p.join()
# #     print(p.is_alive())
# # =============================================================================
# #if __name__ == '__main__':
# #    FooM(int(1e2)) 
# =============================================================================
