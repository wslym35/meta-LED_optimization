#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Jun 28 10:15:57 2024

@author: wkmills
"""

import numpy as np 
import pandas as pd 
#import multiprocessing as mp 
import sys, os
from pathlib import Path
import random # for seeding the SOBOL trial generator 
from datetime import date 

from ax.service.ax_client import AxClient, ObjectiveProperties
#from ax.modelbridge.generation_strategy import GenerationStep, GenerationStrategy
from ax.generation_strategy.generation_strategy import (GenerationStrategy, GenerationStep)
from ax.adapter.registry import Generators

import time 

import warnings
warnings.filterwarnings("ignore", category=UserWarning, message=".*CUDA initialization.*")
# The GPU I'm using is too old to be useful, no need to rub it in :'( 

from directivity import FoM, min_wall_width
min_hole_diameter = 0.050e-6 

dim = "2-d"

params_2d = {
          'Fourier_N' : 50, # N = 50 recommended by Claude after convergence_test, 2026-04-14. Note, this doesn't seem to affect memory bottleneck like the other mesh sizes do. 
          'wavelength_center' : 480e-9, 
          'wavelength_FWHM' : 20e-9, 
          'wavelength_points' : 1, 
          'QW_xy_mesh' : 90, # It would seem 90 is the bare minimum, based on a period of 1.5 * 540 nm, a min_mesa_width of 50 nm, and a non-emitting thickness of 20 nm
          'QW_z_mesh': 55, # N=55 gives at most 20 nm (~wavelength/10 in GaN) between sampled points 
          'k_mesh': 24, # Memory contraint is such that k=28 is about as high as you can go, but k=24 gives roughly the same D result 
          'layer_count' : 5, 
          'layer_names' : ['sapp', 'uniform_GaN', 'etched_GaN', 'ITO', 'air'], # reciprocity plane waves are incident from first layer 
          'layer_thicknesses' : [1e-6, None, None, 0.120e-6, 1e-6], # GaN thicknesses can be variable param 
          'layer_materials' : ["Al2O3 - Palik", "GaN - custom", "GaN - custom", 'ITO - custom', 'etch'],
          'QW_relative_intensities' : [0.45, 0.33, 0.22], # relative intensities of QWs 
          'layer_is_etched' : [False, False, True, True, False], # whether or not to etch through each layer to make the ribbons 
          #'QW_layer' : 2, # Which layer is the QW in: 0, 1, 2, ... The QW will be optimally placed within 50-300 nm of the bottom of this layer (assuming orientation is layer 0 on top)
          'geometry' : 'pillars', # Either 'holes' or 'pillars' 
          'hole_count' : 1, # For now, should be a perfect square (1, 4, 9, ...)
          'target_k' : (0, 0) # (kx, ky) 
          # The params below will be incorporated into 'var' as fixed or range parameters, then passed to FoM in evaluate() 
          , 'period' : [], # um 
          'hole_centers_x' : [],     # x position of each hole center
          'hole_centers_y' : [],     # y position of each hole center
          'hole_diameters' : [],     # every hole is allowed a different diameter 
          'FoM_definition' : ['$D_s+D_p$', '$D_s D_p / (D_s + D_p)$', '$D_s + D_p - |D_s - D_p|$'][0] 
          }
        

min_period = params_2d['wavelength_center'] * 0.50 
max_period = params_2d['wavelength_center'] * 1.50 

var_2d = [ ]
scale = 1e9 
for p in range(2):
    var_2d.append({'name' : 'period_' + str(p), 
                   'value_type' : 'int', # Values are multiplied by `scale` so that they're large enough for Ax to work with 
                   'type' : 'range',
                   'bounds' : [round(scale * el) for el in [min_period, max_period]]
                       })
for h in range(params_2d['hole_count']):
    # Slot index in a sqrt(N) x sqrt(N) grid
    # For hole_count=4, this gives a 2x2 grid of slots
    # This portion of the code would need updating if you want hole_counts other than 1, 4, 9, 16...
    grid_n = int(np.round(np.sqrt(params_2d['hole_count'])))  # e.g. 2 for 4 holes
    slot_x = h % grid_n        # 0, 1, 0, 1
    slot_y = h // grid_n       # 0, 0, 1, 1

    var_2d.append({'name': 'hole_centers_x_' + str(h),
                   'value_type': 'int',
                   'type': 'range',
                   'bounds': [round(scale * el) for el in [
                       slot_x / grid_n * min_period + min_hole_diameter/2 + min_wall_width,
                       (slot_x + 1) / grid_n * max_period - min_hole_diameter/2 - min_wall_width
                   ]]})
    var_2d.append({'name': 'hole_centers_y_' + str(h),
                   'value_type': 'int',
                   'type': 'range',
                   'bounds': [round(scale * el) for el in [
                       slot_y / grid_n * min_period + min_hole_diameter/2 + min_wall_width,
                       (slot_y + 1) / grid_n * max_period - min_hole_diameter/2 - min_wall_width
                   ]]})
    var_2d.append({'name': 'hole_diameters_' + str(h),
                   'value_type': 'int',
                   'type': 'range',
                   'bounds': [round(scale * el) for el in [
                       min_hole_diameter,
                       max_period / grid_n - min_wall_width  # max diameter is slot width minus wall
                   ]]})

# Allow for optimizable GaN thickness 
var_2d.append({'name':'uniform_GaN_thickness', 'value_type' : 'int', 'type':'range', 'bounds':[round(scale * el) for el in [2.0e-6, 4.0e-6]]}) # Dro suggests UID thickness of >=2 um to recover morphology. Double that is a convinient max.  
var_2d.append({'name':'etched_GaN_thickness', 'value_type' : 'int', 'type':'range', 'bounds':[round(scale * el) for el in [0.500e-6, 1.100e-6]]}) # This is probably close to the max you can etch and still get reasonable verticality 

constraints = [ ] 
for h in range(params_2d['hole_count']):
    grid_n = int(np.round(np.sqrt(params_2d['hole_count'])))
    slot_x = h % grid_n
    slot_y = h // grid_n

    # Hole must fit within the unit cell (accounts for variable period)
    constraints.append(f'hole_centers_x_{h} + 0.5*hole_diameters_{h} - period_0 <= {round(scale * -min_wall_width)}')
    constraints.append(f'hole_centers_y_{h} + 0.5*hole_diameters_{h} - period_1 <= {round(scale * -min_wall_width)}')
    constraints.append(f'hole_centers_x_{h} - 0.5*hole_diameters_{h} >= {round(scale * min_wall_width)}')
    constraints.append(f'hole_centers_y_{h} - 0.5*hole_diameters_{h} >= {round(scale * min_wall_width)}')

    # Diameter can't exceed slot width (keeps Sobol trials sane)
    constraints.append(f'hole_diameters_{h} - {round(1/grid_n, 9)}*period_0 <= {round(scale * -min_wall_width)}')
    
def D_opt(dim, variables, params, constraints, trials): 
    
# =============================================================================
#     # Set method for multiprocessing. Linux defaults to fork, so this shouldn't ever be an issue 
#     try:
#         mp.set_start_method('fork') 
#     except RuntimeError:
#         if mp.get_start_method() == 'spawn':
#             print('Somehow the start method got set to spawn. Please restart the kernel.')
#             sys.exit() 
# =============================================================================
    
    # Trials specified as a tuple; unpack here 
    training_count, learning_count = trials 
    
    # Specifies how many training (Sobol model) trials to do, then the rest are Baysian (BoTorch model)  
    gs = GenerationStrategy(
        steps=[
            # 1. Initialization step (does not require pre-existing data and is well-suited for
            # initial sampling of the search space)
            GenerationStep(
                generator=Generators.SOBOL,
                num_trials=training_count,  # How many trials should be produced from this generation step
                min_trials_observed=training_count,  # How many trials need to be completed to move to next model
                max_parallelism=1,  # Max parallelism for this step
                generator_kwargs={"seed": random.randint(1,99), # Random sampling seed, I believe 
                              "fallback_to_sample_polytope": True},  # If all draws fail to pass constraints, fall back on polytope parameter selection model for that training trial 
                generator_gen_kwargs={'model_gen_options' : {"max_rs_draws": 1.1e7}},  # Number of draws to attempt 
            ),
            # 2. Bayesian optimization step (requires data obtained from previous phase and learns
            # from all data available at the time of each new candidate generation call)
            GenerationStep(
                generator=Generators.BOTORCH_MODULAR,
                num_trials=-1,  # No limitation on how many trials should be produced from this step
                max_parallelism=1,  # Parallelism limit for this step, often lower than for Sobol
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
    
    def holes_overlap(params):
        for i in range(params['hole_count']):
            for j in range(i + 1, params['hole_count']):
                dx = params['hole_centers_x'][i] - params['hole_centers_x'][j]
                dy = params['hole_centers_y'][i] - params['hole_centers_y'][j]
                dist = np.sqrt(dx**2 + dy**2)
                min_dist = (params['hole_diameters'][i] + params['hole_diameters'][j]) / 2 + min_wall_width
                if dist < min_dist:
                    return True
        return False

    #'Define how to evaluate trials' 
    def evaluate(parameterization):
        params['period'] = [] 
        params['hole_centers_x'] = []
        params['hole_centers_y'] = []
        params['hole_diameters'] = []
        #params['notch_depths'] = []
        for n in [v['name'] for v in variables]:
            if n[-9:] == 'thickness': 
                if params['layer_count'] == 5:
                    if n[:3] == 'uni': 
                        params['layer_thicknesses'][1] = round(parameterization.get(n) / scale, 9) 
                    elif n[:3] == 'etc': 
                        params['layer_thicknesses'][2] = round(parameterization.get(n) / scale, 9) 
                    elif n[:3] == 'ITO': 
                        params['layer_thicknesses'][3] = round(parameterization.get(n) / scale, 9) 
                    else:
                        raise RuntimeWarning("Uh oh, {n} is not a valid thickness parameter name!")
                elif params['layer_count'] == 4:
                    if n[:3] == 'tot': 
                        params['layer_thicknesses'][1] = round(parameterization.get(n) / scale, 9) 
                    elif n[:3] == 'ITO': 
                        params['layer_thicknesses'][2] = round(parameterization.get(n) / scale, 9) 
                    else:
                        raise RuntimeWarning("Uh oh, {n} is not a valid thickness parameter name!")
# =============================================================================
#             if n == 'GaN_thickness':
#                 raise ValueError("This case (only one GaN layer) needs to be checked and probably modified. 2026/04/20")
# =============================================================================
                #params['layer_thicknesses'][1] = np.round(parameterization.get(n), 3)
            # If the last character is a number, use it as an index to set the params
            # Round because we can't fab with greater than 1 nm accuracy anyway; in fact, the limit is more like 10 nm 
            elif n[-1].isdigit():
                params[n[:-2]].append(round(parameterization.get(n) / scale, 9)) 
            else:
                raise RuntimeWarning("Uh oh, {n} is not a valid parameter name!")
                
        # Check nonlinear overlap constraint 
        if holes_overlap(params):
            raise RuntimeError("Holes overlap.")    

        print("\n Printing params now: ")
        print(params) # Test that it works, and provide backup output channel 
        
# =============================================================================
#         # Use child process to free memory after each call to FoM() 
#         q = mp.Queue() 
#         p = mp.Process(target = FoM, args = (params, q,))
#         p.start() # Starts running the process
#         p.join() # Blocks other things from running until p terminates; best practice to use this
#         results = q.get() # Get output of FoM(); note, directionality.FoM() returns tuple of best directionality and the associated QW depths 
# =============================================================================
        return {"Directivity": (FoM(params), 0.0)}  # Standard error is 0 (we assume S4 gives an exact solution) 


    def parameterization_to_sim_params(parameterization):
        """Convert a raw Ax parameterization dict to the physical/grouped param dict
        that was actually passed to FoM() — mirrors the logic in evaluate()."""
        sim = {k: v for k, v in params.items()}  # start from a copy of fixed params
        sim['period'] = []
        sim['hole_centers_x'] = []
        sim['hole_centers_y'] = []
        sim['hole_diameters'] = []
        for n in [v['name'] for v in variables]:
            if n[-9:] == 'thickness':
                if params['layer_count'] == 5:
                    if n[:3] == 'uni':
                        sim['layer_thicknesses'][1] = round(parameterization.get(n) / scale, 9)
                    elif n[:3] == 'etc':
                        sim['layer_thicknesses'][2] = round(parameterization.get(n) / scale, 9)
                    elif n[:3] == 'ITO':
                        sim['layer_thicknesses'][3] = round(parameterization.get(n) / scale, 9)
                elif params['layer_count'] == 4:
                    if n[:3] == 'tot':
                        sim['layer_thicknesses'][1] = round(parameterization.get(n) / scale, 9)
                    elif n[:3] == 'ITO':
                        sim['layer_thicknesses'][2] = round(parameterization.get(n) / scale, 9)
            elif n[-1].isdigit():
                sim[n[:-2]].append(round(parameterization.get(n) / scale, 9))
        return sim

    results_ls = [['Directivity', 'QW depth', 'Simulated params']]

    #'Run optimization loop'
    trial_count = training_count + learning_count 
    for i in range(trial_count):
        parameterization, trial_index = ax_client.get_next_trial() 
        try: 
            D, QW_depths = evaluate(parameterization)['Directivity'][0] 
        except (ZeroDivisionError, ValueError) as e:
            ax_client.abandon_trial(trial_index) 
            print(f"Trial abandoned. Error: {e} (division by zero might mean QWxy was empty)")
            continue 
        # Store the physical/grouped sim params for this trial (varies row-to-row)
        sim_params = parameterization_to_sim_params(parameterization)
        results_ls.append([D, QW_depths, sim_params]) 
        ax_client.complete_trial(trial_index=trial_index, raw_data=D)  
        print('\n')
    

    #Get the output as a dataframe? 
    trials_df = pd.DataFrame(data=results_ls[1:], columns=results_ls[0]) 
    ax_df = ax_client.get_trials_data_frame()
    trials_df.insert(loc=len(trials_df.columns), column="Generation Model", value=ax_df["generation_node"])
    print(trials_df) 
    
    best_ax_params, best_output = ax_client.get_best_parameters() 
    print(best_ax_params) 
    means, covariances = best_output 
    print(means) 

    # Convert best Ax params (period_0, period_1, ...) to grouped physical params
    # so that opt_params keys match what FoM() actually received (period, hole_centers_x, ...)
    best_sim_params = parameterization_to_sim_params(best_ax_params)
    
    # Easier to use np.array than pd.dataframe 
    data = trials_df.to_numpy() 
    data = np.vstack([['\'D\'', '\'QW depths\'', '\'simulated params\'', '\'method\''], data])
    
    # Save results 
    current_path = Path(os.getcwd()).resolve()
    results_dir = os.path.join(current_path.parent.parent, 'results')
    def make_unique_dir(base_path):
        if not os.path.exists(base_path):
            os.makedirs(base_path)
            return base_path
        counter = 1
        while True:
            new_path = f"{base_path}({counter})"
            if not os.path.exists(new_path):
                os.makedirs(new_path)
                return new_path
            counter += 1
    
    date_dir = make_unique_dir(os.path.join(results_dir, str(date.today()) )) 
    np.save(os.path.join(date_dir, 'opt_params.npy'), best_sim_params)
    np.save(os.path.join(date_dir, 'data.npy'), data) 
    
    return data, ax_client 

#D, QW_depth = FoM(dim, params_2d) 

# Optimization routine 
training_count = 10 * len(var_2d) 
learning_count = 3 * training_count 
#params_2d['reciprocity_N'] = 26 # When NA=1.0 instead of 1.3, decrease reciprocity_N from 26 to 20 to keep the same grittiness 
data, client = D_opt(dim, var_2d, params_2d, constraints, (training_count, learning_count)) 



# After optimization, I should check sensitivity to ITO thickness and wavelength