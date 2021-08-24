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
                residue_position = control[1].strip().split(',')
                if len(residue_position) > 1:
                    control_dict['residue_position'] = [int(x) for x in residue_position]
                elif len(residue_position) == 1 and residue_position[0]:
                    control_dict['residue_position'] = int(residue_position[0])
            elif parameter_name == 'res_mut':
                residue_mutant = control[1].strip().split(',')
                if len(residue_position) > 1:
                    control_dict['residue_mutant'] = residue_mutant
                elif len(residue_mutant) == 1 and residue_mutant[0]:
                    control_dict['residue_mutant'] = residue_mutant[0]
            elif parameter_name == 'chains':
                chains = control[1].strip().split(',')
                if len(chains) > 1:
                    control_dict['chains'] = chains
                elif len(chains) == 1 and chains[0]:
                    control_dict['chains'] = chains[0]
            # Later add checks for the different options of functions
            elif parameter_name == 'function_GB':
                control_dict['function_GB'] = control[1].strip()
            elif parameter_name == 'function_ele':
                control_dict['function_ele'] = control[1].strip()
            elif parameter_name == 'function_Rlj':
                control_dict['function_Rlj'] = control[1].strip()
            elif parameter_name == 'function_epsilonlj':
                control_dict['function_epsilonlj'] = control[1].strip()
    
    # Check for unexpected number of values in the control parameters
    if isinstance(control_dict['chains', list]):
        number_chains = len(control_dict['chains'])
    else:
        number_chains = 1
    if isinstance(control_dict['residue_mutant'], list):
        number_mutant = len(control_dict['residue_mutant'])
    else:
        number_mutant = 1
    if isinstance(control_dict['residue_position'], list):
        number_position = len(control_dict['residue_position'])
    else:
        number_position = 1

    if number_mutant > 1 and number_mutant != number_chains:
        raise Exception("\'res_mut\' parameter must be either a single amino acid or same number of amino acids as number of chains")
    if number_position > 1 and number_position != number_chains:
        raise Exception("\'res_pos\' parameter must be either a single residue number or same number of residue numbers as number of chains")
    parameters = ['windows', 'residue_position', 'residue_mutant', 'chains', 'function_GB', 'function_ele', 'function_Rlj', 'function_epsilonlj']
    if not all(x in control_dict.keys() for x in parameters):
        diff_list = [x for x in parameters if x not in control_dict.keys()]
        raise Exception(f'Missing parameters {diff_list}')

    return control_dict

# Add more functions if more input is necessary