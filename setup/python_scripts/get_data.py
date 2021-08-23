from sys import flags
import numpy as np
import pandas as pd
import os
import subprocess
import fileinput
import shutil

from biopandas.pdb import PandasPdb

"""
Collection of functions to read necessary input provided and calculate information required about the specific structure.
"""

def read_input(control_file):
    """Reads control.txt and outputs a dictionary with control parameters."""
    control_dict = {} # Dictionary to store parameters
    functions = [] # List to store functional relationship of chan
    with open(control_file, 'r') as cf:
        lines = cf.readlines()
        for line in lines:
            control = line.split('=')
            # Create control dictionary
            parameter_name = control[0].strip()
            if parameter_name == 'windows1': # Get the sequence of windows
                sequence = control[1].strip().partition('seq ')[2].partition('))')[0]
                sequence_parms = sequence.split()
                control_dict['windows'] = np.arange(float(sequence_parms[0]), float(sequence_parms[2]), float(sequence_parms[1]))
                control_dict['windows'] =  np.append(control_dict['windows'], float(sequence_parms[2]))
            elif parameter_name == 'windows1+':
                additional_sequence = control[1].strip().partition('(')[2].partition(')')[0]
                if additional_sequence:
                    additional_windows = [float(x) for x in additional_sequence.split()]
                    additional_windows = np.array(additional_windows)
                    control_dict['windows'] = np.concatenate([control_dict['windows'], additional_windows])
                    control_dict['windows'] = np.sort(control_dict['windows'])
            elif parameter_name == 'res_pos':
                control_dict['residue_position'] = int(control[1].strip())
            elif parameter_name == 'res_mut':
                control_dict['residue_mutant'] = control[1].strip()
            elif parameter_name == 'chains':
                chains = control[1].strip().split(',')
                if all(isinstance(x, str) for x in chains):
                    control_dict['chains'] = chains
                elif chains[0] == 'tripeptide':
                    control_dict['chains'] = chains
                elif chains[0] == 'all':
                    control_dict['chains'] = chains
                else:
                    raise AttributeError("\'chains\' can only be \'tripeptide\', \'all\' or list of")
            elif parameter_name == 'function_GB':
                control_dict['function_GB'] = control[1].strip()
            elif parameter_name == 'function_ele':
                control_dict['function_ele'] = control[1].strip()
            elif parameter_name == 'function_Rlj':
                control_dict['function_Rlj'] = control[1].strip()
            elif parameter_name == 'function_epsilonlj':
                control_dict['function_epsilonlj'] = control[1].strip()
    parameters = ['windows', 'residue_position', 'residue_mutant', 'chains', 'function_GB', 'function_ele', 'function_Rlj', 'function_epsilonlj']
    if not all(x in control_dict.keys() for x in parameters):
        diff_list = [x for x in parameters if x not in control_dict.keys()]
        raise Exception(f'Missing parameters {diff_list}')

    return control_dict

def create_mutant(wt_structure, chains):
    if chains == 'tripeptide':
        fjdk
    elif chains == 'all':
        gjgj
    else:
        jgi