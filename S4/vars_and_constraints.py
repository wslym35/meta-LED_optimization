import numpy as np 

params_2d = {'period' : [0.540, 0.540], # um 
          'Fourier-N' : 9, # convergence_test(N=90) is 90% of convergence_test(N=900) as of 2024/09/09  
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
          #, 'ribbon-centers' : [0.133, 0.410], # Ribbon centers are assumed to be at y = 0 for 1-d case
          #'ribbon-widths' : [0.145, 0.190], # Ribbon y-widths are assumed to be period/2 for 1-d case 
          #'notch-centers' : [0.133, 0.410], 
          #'notch-widths' : [0.100, 0.200], 
          #'notch-depths' : [0.075, 0.085]
          }

min_trench_width = 0.020 
min_mesa_width = 0.050 

var_2d = [ ]
for r in range(params_2d['ribbon-count']):
    var_2d.append({'name' : 'ribbon-centers-' + str(r), 
                   'type' : 'range',
                   'bounds' : [r * params_2d['period'][0]/params_2d['ribbon-count'] + min_mesa_width/2 + min_trench_width/2, (r+1) * params_2d['period'][0]/params_2d['ribbon-count'] - min_mesa_width/2 - min_trench_width/2]
                       })
    var_2d.append({'name' : 'ribbon-widths-' + str(r),
                   'type' : 'range',
                   'bounds' : [min_mesa_width, params_2d['period'][0]/params_2d['ribbon-count'] - min_trench_width/2]
                       })
for n in range(params_2d['notch-count']):
    var_2d.append({'name' : 'notch-centers-' + str(n), 
                   'type' : 'range', 
                   'bounds' : [n * params_2d['period'][1]/params_2d['notch-count'] + min_mesa_width/2 + min_trench_width/2, (n+1) * params_2d['period'][1]/params_2d['notch-count'] - min_mesa_width/2 - min_trench_width/2]
                       })
    var_2d.append({'name' : 'notch-widths-' + str(n),
                   'type' : 'range',
                   'bounds' : [min_trench_width, params_2d['period'][1]/params_2d['notch-count'] - min_mesa_width/2]
                       })
    var_2d.append({'name' : 'notch-depths-' + str(n), 
                   'type' : 'range', 
                   'bounds' : [0.0, params_2d['period'][0]/params_2d['ribbon-count'] - min_mesa_width] 
                       })


constraints = [ ] 
for r in range(params_2d['ribbon-count']):
    constraints.append('ribbon-centers-' + str(r) + ' - 0.5*ribbon-widths-' + str(r) + ' >= ' + str(r * params_2d['period'][0]/params_2d['ribbon-count'] + min_trench_width/2))
    constraints.append('ribbon-centers-' + str(r) + ' + 0.5*ribbon-widths-' + str(r) + ' <= ' + str((r+1) * params_2d['period'][0]/params_2d['ribbon-count'] - min_trench_width/2))
for n in range(params_2d['notch-count']):
    constraints.append('notch-centers-' + str(n) + ' - 0.5*notch-widths-' + str(n) + ' >= ' + str(n * params_2d['period'][1]/params_2d['notch-count'] + min_mesa_width/2))
    constraints.append('notch-centers-' + str(n) + ' + 0.5*notch-widths-' + str(n) + ' <= ' + str((n+1) * params_2d['period'][1]/params_2d['notch-count'] - min_mesa_width/2))
for r in range(params_2d['ribbon-count']):
    for n in range(params_2d['notch-count']):
        # Ensure that each ribbon width minus each notch depth is at least min_mesa_width 
        constraints.append('ribbon-widths-' + str(r) + ' - notch-depths-' + str(n) + ' >= ' + str(min_mesa_width))
