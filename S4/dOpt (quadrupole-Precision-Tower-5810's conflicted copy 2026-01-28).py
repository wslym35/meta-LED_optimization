#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Jun 28 10:15:57 2024

@author: wkmills
"""

import numpy as np 
import pandas as pd 
import multiprocessing as mp 
import sys 
import random # for seeding the SOBOL trial generator 

from ax.service.ax_client import AxClient, ObjectiveProperties
from ax.modelbridge.generation_strategy import GenerationStep, GenerationStrategy
from ax.modelbridge.registry import Models

import time 

import directivity 
import makeS4structure 

dim = "2-d"

# Test my S4 optimization against Larry's results 

params_2d = {
          'Fourier-N' : 90, # convergence_test(N=90) is 90% of convergence_test(N=900) as of 2024/09/09  
          #'wavelength' : 0.540, # Caution: the refactive indicies are hardcoded below and will not change when you change this parameter 
          'wavelength' : 0.480, # MUST CHANGE INDICES AS WELL (coded after opt. loop) 
          'reciprocity-N': 26, 
          'polarization' : ['s','p'], # Array containing 's' and/or 'p'  
          'layer-count' : 4, 
          'layer-names' : ['sapp', 'GaN', 'ITO', 'air'], # reciprocity plane waves are incident from first layer 
          'layer-thicknesses' : [0, 1.740, 0.140, 0], # GaN thickness can be variable param 
          'layer-materials' : ['c-sapp', 'c-GaN', 'ITO', 'vac'],
          'layer-epsilons' : [(
                       (3.0761+0.073668j, 0, 0),   # ordinary 
                       (0, 3.0761+0.073668j, 0),    # ordinary 
                       (0, 0, 2.9890+0.069160j)    # extraordinary 
                       ),(
                       (2.4582**2, 0, 0),     # ordinary axis
                       (0, 2.4582**2, 0),     # ordinary axis 
                       (0, 0, 2.3123**2)      # extraordinary axis 
                       ),
                       3.7430+0.017632j, # ITO index 
                       1.00             # Vacuum index 
                       ],
          'layer-is-etched' : [False, True, True, False], # whether or not to etch through each layer to make the ribbons 
          'QW-layer' : 1, # Which layer is the QW in: 0, 1, 2, ... The QW will be optimally placed within [10 nm, 2*lambda/n] of the bottom of this layer 
          'ribbon-count' : 2, # number of nanoribbons to etch 
          'notch-count' : 2, 
          'target-k' : (0, 0) # (kx, ky) 
          # The params below will be incorporated into 'var' as fixed or range parameters, then passed to FoM in evaluate() 
          , 'period' : [], # um 
          'ribbon-centers' : [], 
          'ribbon-widths' : [],  
          'notch-centers' : [], 
          'notch-widths' : [], 
          'notch-depths' : []
          }
                           
min_trench_width = 0.050  
min_mesa_width = 0.050 
min_period = params_2d['wavelength'] * 0.50 
max_period = params_2d['wavelength'] * 1.50 

var_2d = [ ]
for p in range(2):
    var_2d.append({'name' : 'period-' + str(p), 
                   'type' : 'range',
                   'bounds' : [min_period, max_period]
                       })
for r in range(params_2d['ribbon-count']):
    var_2d.append({'name' : 'ribbon-centers-' + str(r), 
                   'type' : 'range',
                   'bounds' : [r * min_period/params_2d['ribbon-count'] + min_mesa_width/2 + min_trench_width/2, (r+1) * max_period/params_2d['ribbon-count'] - min_mesa_width/2 - min_trench_width/2]
                   #bounds' : [r * params_2d['period'][0]/params_2d['ribbon-count'] + min_mesa_width/2 + min_trench_width/2, (r+1) * params_2d['period'][0]/params_2d['ribbon-count'] - min_mesa_width/2 - min_trench_width/2]
                       })
    var_2d.append({'name' : 'ribbon-widths-' + str(r),
                   'type' : 'range',
                   'bounds' : [min_mesa_width, max_period/params_2d['ribbon-count'] - min_trench_width/2]
                   #'bounds' : [min_mesa_width, params_2d['period'][0]/params_2d['ribbon-count'] - min_trench_width/2]
                       })
for n in range(params_2d['notch-count']):
    var_2d.append({'name' : 'notch-centers-' + str(n), 
                   'type' : 'range', 
                   'bounds' : [n * min_period/params_2d['notch-count'] + min_mesa_width/2 + min_trench_width/2, (n+1) * max_period/params_2d['notch-count'] - min_mesa_width/2 - min_trench_width/2]
                   #'bounds' : [n * params_2d['period'][1]/params_2d['notch-count'] + min_mesa_width/2 + min_trench_width/2, (n+1) * params_2d['period'][1]/params_2d['notch-count'] - min_mesa_width/2 - min_trench_width/2]
                       })
    var_2d.append({'name' : 'notch-widths-' + str(n),
                   'type' : 'range',
                   'bounds' : [min_trench_width, max_period/params_2d['notch-count'] - min_mesa_width/2]
                   #'bounds' : [min_trench_width, params_2d['period'][1]/params_2d['notch-count'] - min_mesa_width/2]
                       })
    var_2d.append({'name' : 'notch-depths-' + str(n), 
                   'type' : 'range', 
                   'bounds' : [0.050, max_period/params_2d['ribbon-count'] - min_mesa_width] 
                   #'bounds' : [0.0, params_2d['period'][0]/params_2d['ribbon-count'] - min_mesa_width] 
                       })
# Allow for optimizable GaN thickness 
#var_2d.append({'name':'GaN-thickness', 'type':'range', 'bounds':[1.0, 2.0]})

constraints = [ ] 
for r in range(params_2d['ribbon-count']):
    constraints.append('ribbon-widths-' + str(r) + ' - ' + str(round(1/params_2d['ribbon-count'],6)) + '*period-0 <= ' + str(-min_trench_width/2))
    constraints.append('ribbon-centers-' + str(r) + ' - 0.5*ribbon-widths-' + str(r) + ' - ' + str(round(r/params_2d['ribbon-count'],6)) + '*period-0 >= ' + str(min_trench_width/2))
    constraints.append('ribbon-centers-' + str(r) + ' + 0.5*ribbon-widths-' + str(r) + ' - ' + str(round((r+1)/params_2d['ribbon-count'],6)) + '*period-0 <= ' + str(-min_trench_width/2))
for n in range(params_2d['notch-count']):
    constraints.append('notch-widths-' + str(n) + ' - ' + str(round(1/params_2d['notch-count'],6)) + '*period-1 <= ' + str(-min_mesa_width/2))
    constraints.append('notch-centers-' + str(n) + ' - 0.5*notch-widths-' + str(n) + ' - ' + str(round(n/params_2d['notch-count'],6)) + '*period-1 >= ' + str(min_mesa_width/2))
    constraints.append('notch-centers-' + str(n) + ' + 0.5*notch-widths-' + str(n) + ' - ' + str(round((n+1)/params_2d['notch-count'],6)) + '*period-1 <= ' + str(-min_mesa_width/2))
    constraints.append('notch-depths-' + str(n) + ' - ' + str(round(1/params_2d['ribbon-count'],6)) + '*period-0 <= ' + str(-min_mesa_width))
for r in range(params_2d['ribbon-count']):
    for n in range(params_2d['notch-count']):
        # Ensure that each ribbon width minus each notch depth is at least min_mesa_width 
        constraints.append('ribbon-widths-' + str(r) + ' - notch-depths-' + str(n) + ' >= ' + str(min_mesa_width))
        # For when notch depth is uniform along a ribbon: 
        #constraints.append('ribbon-widths-' + str(r) + ' - notch-depths-' + str(r) + ' >= ' + str(min_mesa_width))


def D_opt(dim, variables, params, constraints, trials): 
    
    # Set method for multiprocessing. Linux defaults to fork, so this shouldn't ever be an issue 
    try:
        mp.set_start_method('fork') 
    except RuntimeError:
        if mp.get_start_method() == 'spawn':
            print('Somehow the start method got set to spawn. Please restart the kernel.')
            sys.exit() 
    
    # Trials specified as a tuple; unpack here 
    training_count, learning_count = trials 
    
    # Specifies how many training (Sobol model) trials to do, then the rest are Baysian (BoTorch model)  
    gs = GenerationStrategy(
        steps=[
            # 1. Initialization step (does not require pre-existing data and is well-suited for
            # initial sampling of the search space)
            GenerationStep(
                model=Models.SOBOL,
                num_trials=training_count,  # How many trials should be produced from this generation step
                min_trials_observed=training_count,  # How many trials need to be completed to move to next model
                max_parallelism=5,  # Max parallelism for this step
                model_kwargs={"seed": random.randint(1,99), # Random sampling seed, I believe 
                              "fallback_to_sample_polytope": True},  # If all draws fail to pass constraints, fall back on polytope parameter selection model for that training trial 
                model_gen_kwargs={'model_gen_options' : {"max_rs_draws": 1.1e7}},  # Number of draws to attempt 
            ),
            # 2. Bayesian optimization step (requires data obtained from previous phase and learns
            # from all data available at the time of each new candidate generation call)
            GenerationStep(
                model=Models.BOTORCH_MODULAR,
                num_trials=-1,  # No limitation on how many trials should be produced from this step
                max_parallelism=3,  # Parallelism limit for this step, often lower than for Sobol
                # More on parallelism vs. required samples in BayesOpt:
                # https://ax.dev/docs/bayesopt.html#tradeoff-between-parallelism-and-total-number-of-trials
            ),
        ]
    )
    
        
    #'Initialize client' 
    ax_client = AxClient(gs) 

    #'Set up experiment' 
    ax_client.create_experiment(
        name = 'Experiment_name',
        parameters = variables, 
        objectives = {"Directivity": ObjectiveProperties(minimize=False)}, # False indicates to maximize the objective 
        parameter_constraints = constraints,  # Optional.
        #outcome_constraints=["l2norm <= 1.25"],  # Optional.
        )

    #'Define how to evaluate trials' 
    def evaluate(parameterization):
        params['period'] = [] 
        params['ribbon-centers'] = [] 
        params['ribbon-widths'] = []
        params['notch-centers'] = []
        params['notch-widths'] = []
        params['notch-depths'] = []
        for n in [v['name'] for v in variables]:
            if n == 'GaN-thickness':
                params['layer-thicknesses'][1] = np.round(parameterization.get(n), 3)
            # If the last character is a number, use it as an index to set the params
            # Round because we can't fab with greater than 1 nm accuracy anyway; in fact, the limit is more like 10 nm 
            if n[-1].isdigit():
                params[n[:-2]].append(np.round(parameterization.get(n), 3)) 
# =============================================================================
#             else:
#                 print('\n\nEntering the else statement\n\n') 
#                 params[n].append(parameterization.get(n)) 
# =============================================================================
        print("Printing params now: ")
        print(params) # Test that it works 
        # Use child process to free memory after each call to FoM() 
        q = mp.Queue() 
        p = mp.Process(target = directivity.FoM, args = (dim, params, q,))
        p.start() # Starts running the process
        p.join() # Blocks other things from running until p terminates; best practice to use this
        results = q.get() # Get output of FoM(); note, directionality.FoM() returns tuple of best directionality and the associated QW depth 
        return {"Directivity": (results, 0.0)}  # Standard error is 0 (we assume S4 gives an exact solution) 


    results_ls = [['Directivity', 'QW depth', 'Input params']] 

    #'Run optimization loop'
    trial_count = training_count + learning_count 
    for i in range(trial_count):
        parameterization, trial_index = ax_client.get_next_trial() 
        try: 
            D, QW_depth = evaluate(parameterization)['Directivity'][0] 
        except ZeroDivisionError:
            ax_client.abandon_trial(trial_index) 
            print("Trial abandoned. QW_xy was likely empty.")
            continue 
        #D, QW_depth = results["Directivity"][0] 
        #error = results["Directivity"][1] 
        # Add results to my own data list 
        results_ls.append([D, QW_depth, parameterization]) 
        # Local evaluation here can be replaced with deployment to external system.
        # Ax auto-selects an appropriate optimization algorithm based on the search space. For more advance use cases that require a specific optimization algorithm, pass a generation_strategy argument into the AxClient constructor.
        ax_client.complete_trial(trial_index=trial_index, raw_data=D)  
        print('\n')
    


    #Get the output as a dataframe? 
    trials_df = pd.DataFrame(data=results_ls[1:], columns=results_ls[0]) 
    ax_df = ax_client.generation_strategy.trials_as_df 
    trials_df.insert(loc=len(trials_df.columns), column="Generation Model", value=ax_df["Generation Model"])
    print(trials_df) 
    
    best_input, best_output = ax_client.get_best_parameters() 
    print(best_input) 
    means, covariances = best_output 
    print(means) 
    
    # Easier to use np.array than pd.dataframe 
    data = trials_df.to_numpy() 
    p_array = np.reshape(np.repeat(params, trial_count), (trial_count,1))
    data = np.append(data, p_array, 1) 
    data = np.vstack([['\'D\'','\'QW depth\'', '\'variables\'', '\'method\'', '\'fixed parameters\''], data])
    
    return data, ax_client 





# Emission centered at 480 nm instead of 540 nm 
# Wavelength is changed above, so that period bounds can be set accordingly 
params_2d['wavelength'] = 0.480 
params_2d['layer-epsilons'] = [( #sapp, GaN, ITO, air @ 480 nm 
             (3.0761+0.073668j, 0, 0),   # ordinary 
             (0, 3.0761+0.073668j, 0),    # ordinary 
             (0, 0, 2.9890+0.069160j)    # extraordinary 
             ),(
             (2.4582**2, 0, 0),     # ordinary axis
             (0, 2.4582**2, 0),     # ordinary axis 
             (0, 0, 2.3123**2)      # extraordinary axis 
             ), 
             3.7430+0.017632j, # ITO index 
             1.00              # Vacuum index 
             ]

# =============================================================================
# params_2d['wavelength'] = 0.470 
# # Emission centered at 470 nm instead of 540 nm 
# # Wavelength is changed above, so that period bounds can be set accordingly 
# params_2d['layer-epsilons'] = [( #sapp, GaN, ITO, air @ 480 nm 
#              (3.0761+0.073668j, 0, 0),   # ordinary 
#              (0, 3.0761+0.073668j, 0),    # ordinary 
#              (0, 0, 2.9925+0.069200j)    # extraordinary 
#              ),(
#              (2.4668**2, 0, 0),     # ordinary axis
#              (0, 2.4668**2, 0),     # ordinary axis 
#              (0, 0, 2.3123**2)      # extraordinary axis 
#              ), 
#              3.7854+0.019235j, # ITO index 
#              1.00              # Vacuum index 
#              ]
# =============================================================================
                 
# =============================================================================
# # Emission centered at 465 nm instead of 540 nm 
# # Wavelength is changed above, so that period bounds can be set accordingly 
# params_2d['wavelength'] = 0.465 
# params_2d['layer-epsilons'] = [( #sapp, GaN, ITO, air @ 480 nm 
#              (3.0761+0.073668j, 0, 0),   # ordinary 
#              (0, 3.0761+0.073668j, 0),    # ordinary 
#              (0, 0, 2.9925+0.069200j)    # extraordinary 
#              ),(
#              (2.4715**2, 0, 0),     # ordinary axis
#              (0, 2.4715**2, 0),     # ordinary axis 
#              (0, 0, 2.3123**2)      # extraordinary axis 
#              ), 
#              3.8071+0.020153j, # ITO index 
#              1.00              # Vacuum index 
#              ]
# =============================================================================

# =============================================================================
# # Parameters for trial 239 of 2025-10-03 
# # The one with 200 nm of current-spreading, unetched n-GaN  
# params_2d['layer-count'] = 5 
# params_2d['layer-names'] = ['sapp', 'n-GaN-base', 'GaN-ribbons', 'ITO', 'air'] # reciprocity plane waves are incident from first layer 
# params_2d['layer-thicknesses'] = [0, 0.200, 1.570, 0.124, 0] # GaN thickness can be variable param 
# params_2d['layer-materials'] = ['c-sapp', 'c-GaN', 'c-GaN', 'ITO', 'vac']
# params_2d['layer-epsilons'] = [(
#              (3.0761+0.073668j, 0, 0),   # ordinary 
#              (0, 3.0761+0.073668j, 0),    # ordinary 
#              (0, 0, 2.9890+0.069160j)    # extraordinary 
#              ),(
#              (2.4582**2, 0, 0),     # ordinary axis
#              (0, 2.4582**2, 0),     # ordinary axis 
#              (0, 0, 2.3123**2)      # extraordinary axis 
#              ),(
#              (2.4582**2, 0, 0),     # ordinary axis
#              (0, 2.4582**2, 0),     # ordinary axis 
#              (0, 0, 2.3123**2)      # extraordinary axis 
#              ), 
#              3.7430+0.017632j, # ITO index 
#              1.00             # Vacuum index 
#              ]
# params_2d['layer-is-etched'] = [False, False, True, True, False] # whether or not to etch through each layer to make the ribbons 
# params_2d['QW-layer'] = 2 # Which layer is the QW in: 0, 1, 2, ... The QW will be optimally placed within [10 nm, 2*lambda/n] of the bottom of this layer 
# params_2d['period'] = [0.686, 0.451]
# params_2d['ribbon-centers'] = [0.489] 
# params_2d['ribbon-widths'] = [0.343] 
# params_2d['notch-centers'] = [0.155, 0.276] 
# params_2d['notch-widths'] = [0.070, 0.050] 
# params_2d['notch-depths'] = [0.293, 0.000] 
# QW_depth = 1.685 # measured from sapp. 
# =============================================================================
                       
# =============================================================================
# # Parameters for trial 97 of 2025-04-29 
# params_2d['layer-thicknesses'] = [0, 1.77, .110, 0] 
# params_2d['period'] = [0.338, 0.624] 
# params_2d['ribbon-centers'] = [0.200] 
# params_2d['ribbon-widths'] = [0.145] # Test with 0.15 for stability 
# params_2d['notch-centers'] = [0.160] 
# params_2d['notch-widths'] = [0.175] 
# params_2d['notch-depths'] = [0.00] 
# QW_depth = 1.46 # measured from sapp. 
# =============================================================================

# =============================================================================
# # Parameters for highly sensitive mode (see notes from 2025-04-25) 
# params_2d['layer-thicknesses'][1] = 1.9999999999993134
# params_2d['period'] = [0.5557414710969972, 0.5525891006157623]
# params_2d['ribbon-centers'] = [0.3460707981203343] 
# params_2d['ribbon-widths'] = [0.3693413459534806]
# params_2d['notch-centers'] = [0.3983340637735282] 
# params_2d['notch-widths'] = [0.2585100736847604] 
# params_2d['notch-depths'] = [0.31934134595321495] 
# QW_depth = 1.6094703441535727 # measured from sapp. 
# =============================================================================

# =============================================================================
# # Results of 2025-02-13 optimization 
# #params_2d['wavelength'] = 0.46 # Check for stability 
# params_2d['ribbon-count'] = 1
# params_2d['notch-count'] = 2 
# params_2d['period'] = [0.424 * 48/54, 0.673 * 48/54]
# params_2d['ribbon-centers'] = [0.281 * 48/54]
# params_2d['ribbon-widths'] = [0.159 * 48/54 + 0.005] # finalized 
# params_2d['notch-centers'] = [0.153 * 48/54, 0.448 * 48/54] 
# params_2d['notch-widths'] = [0.063 * 48/54, 0.193 * 48/54] # finalized 
# params_2d['notch-depths'] = [0.129 * 48/54, 0.0 * 48/54]# finalized 
# QW_depth = 1.685 # measured from sapp. 
# # 0.784 for GaN thickness of 0.996*48/54 and ITO thickness of 0.140*48/54 
# # 1.685 for GaN thickness of 0.996*48/54*2 and ITO thickness of 0.140*48/54 
# params_2d['layer-thicknesses'] = [0, 0.996 * 48/54 * 2, 0.140 * 48/54, 0]
# params_2d['target-k'] = (0.0,0.0)
# =============================================================================
                 
# =============================================================================
# # Optimization routine 
# training_count = 10 * len(var_2d) 
# learning_count = 3 * training_count 
# #params_2d['reciprocity_N'] = 26 # When NA=1.0 instead of 1.3, decrease reciprocity_N from 26 to 20 to keep the same grittiness 
# data, client = D_opt(dim, var_2d, params_2d, constraints, (training_count, learning_count)) 
# =============================================================================

# =============================================================================
# # Results of 2025-06-24 optimization #2
# params_2d['period'] = [0.720, 0.692]
# params_2d['ribbon-centers'] = [0.242, 0.552]
# params_2d['ribbon-widths'] = [0.054, 0.134] 
# params_2d['notch-centers'] = [0.199] 
# params_2d['notch-widths'] = [0.050]
# params_2d['notch-depths'] = [0.0]
# QW_depth = 1.685 # measured from sapp. 
# #params_2d['reciprocity-N'] = 52
# =============================================================================
# =============================================================================
# # Results of 2025-06-24 optimization #3
# params_2d['period'] = [0.609, 0.679]
# params_2d['ribbon-centers'] = [0.208, 0.435]
# params_2d['ribbon-widths'] = [0.053, 0.211] 
# params_2d['notch-centers'] = [0.178, 0.549] 
# params_2d['notch-widths'] = [0.174, 0.148]
# params_2d['notch-depths'] = [0.003, 0.0]
# QW_depth = 1.685 # measured from sapp. 
# #params_2d['reciprocity-N'] = 52
# =============================================================================
# =============================================================================
# # Results of 2025-06-24 optimization #5
# params_2d['period'] = [0.72, 0.72]
# params_2d['ribbon-centers'] = [0.209, 0.582]
# params_2d['ribbon-widths'] = [0.077, 0.226] 
# params_2d['notch-centers'] = [0.253] 
# params_2d['notch-widths'] = [0.128]
# params_2d['notch-depths'] = [0.027]
# QW_depth = 1.685 # measured from sapp. 
# #params_2d['reciprocity-N'] = 52
# =============================================================================
# =============================================================================
# # Results of 2025-06-24 optimization #6
# params_2d['period'] = [0.526, 0.485]
# params_2d['ribbon-centers'] = [0.190, 0.410]
# params_2d['ribbon-widths'] = [0.080, 0.076] 
# params_2d['notch-centers'] = [0.341] 
# params_2d['notch-widths'] = [0.100]
# params_2d['notch-depths'] = [0.026]
# QW_depth = 1.685 # measured from sapp. 
# #params_2d['reciprocity-N'] = 52
# =============================================================================
# =============================================================================
# # Results of 2025-06-24 optimization #7
# params_2d['ribbon-count'] = 2 # number of nanoribbons to etch 
# params_2d['notch-count'] = 1 
# params_2d['target-k'] = (0.3, 0) # (kx, ky) 
# params_2d['period'] = [0.617, 0.528]
# params_2d['ribbon-centers'] = [0.160, 0.530]
# params_2d['ribbon-widths'] = [0.050, 0.124] 
# params_2d['notch-centers'] = [0.233] 
# params_2d['notch-widths'] = [0.184]
# params_2d['notch-depths'] = [0.0]
# QW_depth = 1.685 # measured from sapp. 
# #params_2d['reciprocity-N'] = 52
# =============================================================================

# =============================================================================
# # Results of 2025-02-13 optimization 
# params_2d['period'] = [0.424, 0.673]
# params_2d['ribbon-centers'] = [0.281]
# params_2d['ribbon-widths'] = [0.159] 
# params_2d['notch-centers'] = [0.153, 0.448] 
# params_2d['notch-widths'] = [0.063, 0.193]
# params_2d['notch-depths'] = [0.129, 0.0]
# QW_depth = 0.864 # measured from sapp. 
# #params_2d['reciprocity-N'] = 52
# =============================================================================

# =============================================================================
# # Results of 2025-01-13 optimization 
# params_2d['period'] = [0.499, 0.599]
# params_2d['ribbon-centers'] = [0.073, 0.260, 0.413]
# params_2d['ribbon-widths'] = [0.116, 0.066, 0.049] 
# params_2d['notch-centers'] = [0.123, 0.310, 0.550] 
# params_2d['notch-widths'] = [0.120, 0.148, 0.044]
# params_2d['notch-depths'] = [0, 0.015, 0.019]
# QW_depth = 0.986 # measured from sapp. 
# #params_2d['reciprocity-N'] = 52
# =============================================================================

# =============================================================================
# # Results of 2025-02-03 optimization 
# params_2d['period'] = [0.698, 0.543] #[0.6976447484269739, 0.5430398130230606]
# params_2d['ribbon-centers'] = [0.194, 0.538] #[0.19218935960903763, 0.5325503608491272]
# params_2d['ribbon-widths'] = [0.246, 0.236] #[0.2438067206740379, 0.23792181830853223] 
# params_2d['notch-centers'] = [0.466] #[0.4661413100268692] 
# params_2d['notch-widths'] = [0.040] #[0.04051819493528455]
# params_2d['notch-depths'] = [0.003] #[0.0027578413719311357]
# QW_depth = 0.550 # measured from sapp. 
# =============================================================================

# =============================================================================
# # Results of 2025-01-27 optimization (1 ribbon, 1 notch) 
# params_2d['period'] = [0.469, 0.499]
# params_2d['ribbon-centers'] = [0.285] 
# params_2d['ribbon-widths'] = [0.338] 
# params_2d['notch-centers'] = [0.276]
# params_2d['notch-widths'] = [0.355]
# params_2d['notch-depths'] = [0.308] 
# QW_depth = 0.619 
# =============================================================================

# =============================================================================
# params_2d['layer-epsilons'] = [(
#              (3.0726+0.073626j, 0, 0),   # ordinary 
#              (0, 3.0726+0.073626j, 0),    # ordinary 
#              (0, 0, 3.0726+0.073626j)    # extraordinary 
#              ),(
#              (2.4194**2, 0, 0),     # ordinary axis
#              (0, 2.4194**2, 0),     # ordinary axis 
#              (0, 0, 2.4194**2)      # extraordinary axis 
#              ), 
#              1.00 + 0.00j 
#              ] 
# =============================================================================

# =============================================================================
# # Unpatterned thin film with optimized QW depth 
# params_2d['layer-count'] = 3
# params_2d['period'] = [10, 10]
# params_2d['reciprocity-N'] = 52
# params_2d['Fourier-N'] = 5
# params_2d['ribbon-count'] = 0
# params_2d['notch-count'] = 0 
# params_2d['layer-thicknesses'] = [0, 0.160 + 0.002 + 0.854, 0] 
# params_2d['wavelength'] = 0.560
# QW_depth = .854 # measured from first layer (sapp.); 0.672 is optimum for thin-film of thickness 0.996
# #directivity.FoM_2d(params_2d)
# =============================================================================

# =============================================================================
# # The first epi Stephen grew 
# params_2d['layer-count'] = 3
# params_2d['period'] = [10, 10]
# params_2d['reciprocity-N'] = 52
# params_2d['Fourier-N'] = 5
# params_2d['ribbon-count'] = 0
# params_2d['notch-count'] = 0 
# params_2d['layer-thicknesses'] = [0, 0.160 + 0.002 + 0.854, 0] 
# params_2d['wavelength'] = 0.510
# QW_depth = .854 # measured from first layer (sapp.); 0.672 is optimum for thin-film of thickness 0.996
# #directivity.FoM_2d(params_2d)
# =============================================================================


# =============================================================================
# # Prasad's unpatterned thin film 
# params_2d['layer-count'] = 3
# params_2d['period'] = [10, 10]
# params_2d['reciprocity-N'] = 52
# params_2d['Fourier-N'] = 5
# params_2d['ribbon-count'] = 0
# params_2d['notch-count'] = 0 
# params_2d['layer-thicknesses'] = [0, 1.45, 0] 
# params_2d['wavelength'] = 0.530 
# params_2d['layer-epsilons'] = [1.77**2, 2.23**2, 1.00] 
# QW_depth = 1.35 # measured from first layer (sapp.); 0.672 is optimum for thin-film of thickness 0.996
# #directivity.FoM_2d(params_2d)
# =============================================================================

# =============================================================================
# # Thin-film optimization 
# params_2d['layer-count'] = 3
# params_2d['period'] = [10, 10]
# params_2d['reciprocity-N'] = 52
# params_2d['Fourier-N'] = 5
# params_2d['ribbon-count'] = 0
# params_2d['notch-count'] = 0 
# params_2d['layer-thicknesses'] = [0, 0.996, 0] 
# params_2d['wavelength'] = 0.530 
# params_2d['layer-epsilons'] = [(
#              (3.0726+0.073626j, 0, 0),   # ordinary 
#              (0, 3.0726+0.073626j, 0),    # ordinary 
#              (0, 0, 3.0726+0.073626j)    # extraordinary 
#              ),(
#              (2.4194**2, 0, 0),     # ordinary axis
#              (0, 2.4194**2, 0),     # ordinary axis 
#              (0, 0, 2.4194**2)      # extraordinary axis 
#              ), 
#              1.00 + 0.00j 
#              ] 
# QW_depth = 0.558 # measured from first layer (sapp.); 0.558 is optimum for thin-film of thickness 0.996
# #directivity.FoM_2d(params_2d) 
# =============================================================================

# Larry's dual-pol optimized metasurface
params_2d['wavelength'] = 0.540 # MUST CHANGE INDICES AS WELL (coded after opt. loop) 
params_2d['layer-count'] = 3 
params_2d['layer-names'] = ['sapp', 'GaN', 'air'] # reciprocity plane waves are incident from first layer 
params_2d['layer-thicknesses'] = [0, 0.996, 0] # GaN thickness can be variable param 
params_2d['layer-materials'] = ['c-sapp', 'c-GaN', 'vac']
params_2d['layer-epsilons'] = [(
             (3.0726+0.073626j, 0, 0),   # ordinary 
             (0, 3.0726+0.073626j, 0),    # ordinary 
             (0, 0, 2.9890+0.069160j)    # extraordinary 
             ),(
             (2.4194**2, 0, 0),     # ordinary axis
             (0, 2.4194**2, 0),     # ordinary axis 
             (0, 0, 2.3121**2)      # extraordinary axis 
             ),
             1.00             # Vacuum index 
             ]
params_2d['layer-is-etched'] = [False, True, False] # whether or not to etch through each layer to make the ribbons 
params_2d['QW-layer'] = 1 # Which layer is the QW in: 0, 1, 2, ... The QW will be optimally placed within [10 nm, 2*lambda/n] of the bottom of this layer 
params_2d['ribbon-count'] = 3 # number of nanoribbons to etch 
params_2d['notch-count'] = 0 
params_2d['target-k'] = (0, 0) # (kx, ky) 
# The params below will be incorporated into 'var' as fixed or range parameters, then passed to FoM in evaluate() 
params_2d['period'] = [0.540, 0.540] # um 
params_2d['ribbon-centers'] = [0.065, 0.215, 0.435] 
params_2d['ribbon-widths'] = [0.050, 0.070, 0.110]  
params_2d['notch-centers'] = [] 
params_2d['notch-widths'] = [] 
params_2d['notch-depths'] = []
QW_depth = 0.996 - (0.120 + 0.003) # measured from first layer (sapp.)

# =============================================================================
# params_2d = {'Fourier-N': 90, 'wavelength': 0.54, 'reciprocity-N': 26, 'polarization': ['s', 'p'], 'layer-count': 3, 'layer-names': ['sapp', 'GaN', 'air'], 'layer-thicknesses': [0, 0.996, 0], 'layer-materials': ['c-sapp(540nm)', 'c-GaN(540nm)', 'vac'], 'layer-epsilons': [(((3.0726+0.073626j), 0, 0), (0, (3.0726+0.073626j), 0), (0, 0, (2.989+0.06916j))), ((5.85349636, 0, 0), (0, 5.85349636, 0), (0, 0, 5.34580641)), (1+0j)], 'layer-is-etched': [False, True, False], 'QW-layer': 1, 'ribbon-count': 3, 'notch-count': 3, 'target-k': (0, 0), 'period': [0.5207941891518663, 0.6438343045404686], 'ribbon-centers': [0.07411327637296769, 0.3037281582148236, 0.47516329593492007], 'ribbon-widths': [0.035201465216677386, 0.049492082704456625, 0.03520146521676472], 'notch-centers': [0.15621531344420592, 0.3341998392656123, 0.47503417351713345], 'notch-widths': [0.0867918135828172, 0.07609373344123979, 0.0571007672417993], 'notch-depths': [3.982266383994796e-16, 0.0052014652168501, 0.0052014652168062715]}
# tic = time.time() 
# try: 
#     D, QW_depth = directivity.FoM(dim, params_2d) 
# except ZeroDivisionError: 
#     print("Error encountered")
# toc = time.time() 
# print('Total time = ' + str(toc-tic) + ' seconds')
# =============================================================================

# =============================================================================
# # 2025-10-31 optimization, trial 87 
# params_2d['period'] = [0.720, 0.456]
# params_2d['ribbon-centers'] = [0.185, 0.514]
# params_2d['ribbon-widths'] = [0.092, 0.093]
# params_2d['notch-centers'] = [0.121, 0.310]
# params_2d['notch-widths'] = [0.163, 0.114]
# params_2d['notch-depths'] = [0.042, 0.042] 
# QW_depth = 1.685 
# =============================================================================

# =============================================================================
# #2025-11-04 optimization, trial 120 
# params_2d['period'] = [0.691, 0.683]
# params_2d['ribbon-centers'] = [0.162, 0.506]
# params_2d['ribbon-widths'] = [0.234, 0.080]
# params_2d['notch-centers'] = [0.122, 0.505]
# params_2d['notch-widths'] = [0.106, 0.165]
# params_2d['notch-depths'] = [0.024, 0.019] 
# QW_depth = 1.685 
# =============================================================================

# =============================================================================
# #2025-11-04 optimization, trial 478 
# params_2d['period'] = [0.720, 0.719]
# params_2d['ribbon-centers'] = [0.133, 0.446]
# params_2d['ribbon-widths'] = [0.216, 0.109]
# params_2d['notch-centers'] = [0.082, 0.483]
# params_2d['notch-widths'] = [0.050, 0.173]
# params_2d['notch-depths'] = [0.059, 0.031] 
# QW_depth = 1.685 
# =============================================================================

# =============================================================================
# #2025-11-04 optimization, trial 406 
# params_2d['wavelength'] = 0.470 #What was actually fabbed 
# params_2d['layer-thicknesses'] = [0, 1.740, 0.140, 0] # What was actually fabbed
# params_2d['period'] = [0.720, 0.720]
# params_2d['ribbon-centers'] = [0.145, 0.475]
# params_2d['ribbon-widths'] = [0.215, 0.108]
# params_2d['notch-centers'] = [0.094, 0.554]
# params_2d['notch-widths'] = [0.060, 0.197]
# params_2d['notch-depths'] = [0.058, 0.037] 
# QW_depth = 1.685 
# =============================================================================

# =============================================================================
# #2025-11-10 optimization, trial 275 (Ds=22, Dp=0.4) 
# # (all "good" results from this run are unbalanced like this; 
# # I'm changing the figure of merit to Larry's: Ds + Dp - |Ds - Dp|)
# params_2d['target-k'] = (0.3, 0) 
# params_2d['period'] = [0.599, 0.666]
# params_2d['ribbon-centers'] = [0.164, 0.427]
# params_2d['ribbon-widths'] = [0.221, 0.204]
# params_2d['notch-centers'] = [0.181, 0.401]
# params_2d['notch-widths'] = [0.082, 0.086]
# params_2d['notch-depths'] = [0.010, 0.008] 
# QW_depth = 1.685 
# =============================================================================

# =============================================================================
# #2025-11-14, trial 269 (sent to Roark)
# params_2d['wavelength'] = 0.470 
# params_2d['target-k'] = (0, 0) 
# params_2d['ribbon-count'] = 1
# params_2d['notch-count'] = 1
# params_2d['period'] = [0.698, 0.691]
# params_2d['ribbon-centers'] = [0.188]
# params_2d['ribbon-widths'] = [0.220] #[0.220 - 0.010]
# params_2d['notch-centers'] = [0.331 - 0.005] #[0.331] # 
# params_2d['notch-widths'] = [0.444 ]#+ 0.010] #[0.444] was optimized for 20 nm non-emitting perimeter 
# params_2d['notch-depths'] = [0.170] #[0.170 - 0.010]  
# QW_depth = 1.685 
# =============================================================================

# =============================================================================
# params_2d['wavelength'] = 0.478
# =============================================================================

# =============================================================================
# #2025-11-19, trial 414 (kx=0.3)
# params_2d['target-k'] = (0.3, 0) 
# params_2d['ribbon-count'] = 2
# params_2d['notch-count'] = 2
# params_2d['period'] = [0.664, 0.705]
# params_2d['ribbon-centers'] = [0.195, 0.524]
# params_2d['ribbon-widths'] = [0.143, 0.143]
# params_2d['notch-centers'] = [0.186, 0.534]
# params_2d['notch-widths'] = [0.219, 0.197]
# params_2d['notch-depths'] = [0.093, 0.090] 
# QW_depth = 1.685 
# =============================================================================

# =============================================================================
# #2025-11-19, trial 407 (kx=0.7)
# params_2d['target-k'] = (0.7, 0) 
# params_2d['ribbon-count'] = 2
# params_2d['notch-count'] = 2
# params_2d['period'] = [0.510, 0.651]
# params_2d['ribbon-centers'] = [0.152, 0.407]
# params_2d['ribbon-widths'] = [0.155, 0.154]
# params_2d['notch-centers'] = [0.202, 0.504]
# params_2d['notch-widths'] = [0.198, 0.097]
# params_2d['notch-depths'] = [0.104, 0.087] 
# QW_depth = 1.685 
# =============================================================================

#D, QW_depth = directivity.FoM(dim, params_2d) 

#params_2d['wavelength'] = 0.530 

# =============================================================================
# # Plot at QW depth 
# tic = time.time() 
# I = directivity.plot_at_QWz_2d(makeS4structure.struct_2d(params_2d), params_2d, QW_depth, correct_apod=True)
# toc = time.time() 
# print('Total time: ' + str(toc - tic))
# =============================================================================

# Plot multiple emission wavelengths
# =============================================================================
# weights = [                              [0.468, 0.52], [0.469, 0.59], [0.470, 0.65],
#            [0.471, 0.72], [0.472, 0.79], [0.473, 0.85], [0.474, 0.91], [0.475, 0.96],
#            [0.476, 0.99], [0.477, 1.00], [0.478, 1.00], [0.479, 0.98], [0.480, 0.95], 
#            [0.481, 0.91], [0.482, 0.87], [0.483, 0.82], [0.484, 0.78], [0.485, 0.73], 
#            [0.486, 0.68], [0.487, 0.63], [0.488, 0.59], [0.489, 0.55], [0.490, 0.51] 
#            ] # Generated from 2026-01-07 data 
# =============================================================================
spectral_FWHM = 0.020 
wavelength_range = np.arange(params_2d['wavelength'] - spectral_FWHM/2, params_2d['wavelength'] + spectral_FWHM/2, 0.001)
weights = np.stack((wavelength_range, np.exp(-(wavelength_range - 0.540)**2 / (spectral_FWHM**2/(4*np.log(2))))), axis=1)
run = True
if run: 
    input('You sure you want to recalculate everything?') 
    I_total = 0 
    for w in weights:
        print('Calculating wavelength' + str(w[0]))
        params_2d['wavelength'] = w[0] 
        I_total += w[1] * directivity.plot_at_QWz_2d(makeS4structure.struct_2d(params_2d), params_2d, QW_depth, correct_apod=True)['dual']

# =============================================================================
# # Check sensitivity to ITO thickness
# for ITO in [0.140]:
#     params_2d['layer-thicknesses'][2] = ITO
#     print("ITO thickness = " + str(ITO) + " um")
#     directivity.plot_at_QWz_2d(makeS4structure.struct_2d(params_2d), params_2d, QW_depth, correct_apod = True)
# =============================================================================

# =============================================================================
# tic = time.time() 
# Ds = directivity.D_z_2d(makeS4structure.struct_2d(params_2d), params_2d, [0.550], 's')
# Dp = directivity.D_z_2d(makeS4structure.struct_2d(params_2d), params_2d, [0.550], 'p')
# Dsp = Ds + Dp - np.abs(Ds - Dp)
# toc = time.time() 
# print('Fourier N: ' + str(params_2d['Fourier-N']))
# print('Reciprocity N: ' + str(params_2d['reciprocity-N']))
# print('Dual directivity = ' + str(Dsp)) 
# print('Total time = ' + str(toc-tic) + ' seconds')
# =============================================================================

