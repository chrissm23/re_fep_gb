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
            elif parameter_name == 'windows1+': # Get custom windows
                additional_sequence = control[1].strip().partition('(')[2].partition(')')[0]
                if additional_sequence:
                    additional_windows = [float(x) for x in additional_sequence.split()]
                    additional_windows = np.array(additional_windows)
                    control_dict['windows'] = np.concatenate([control_dict['windows'], additional_windows])
                    control_dict['windows'] = np.sort(control_dict['windows'])
            elif parameter_name == 'res_pos': # Get position of the mutation(s)
                residue_position = control[1].strip().split(',')
                if len(residue_position) > 1:
                    control_dict['residue_position'] = [int(x) for x in residue_position]
                elif len(residue_position) == 1 and residue_position[0]:
                    control_dict['residue_position'] = int(residue_position[0])
            elif parameter_name == 'res_mut': # Get amino acid(s) to mutate to
                residue_mutant = control[1].strip().split(',')
                if len(residue_position) > 1:
                    control_dict['residue_mutant'] = residue_mutant
                elif len(residue_mutant) == 1 and residue_mutant[0]:
                    control_dict['residue_mutant'] = residue_mutant[0]
            elif parameter_name == 'chains': # Get chain(s) to mutate
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
    if isinstance(control_dict['chains'], list):
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

def make_executable(path):
    """Make file executable"""
    mode = os.stat(path).st_mode
    mode |= (mode & 0o444) >> 2 # Copy R bits to X
    os.chmod(path, mode)

def replace_in_file(path, replace_dict):
    """Replace in file changes[0] to changes[1]"""
    with fileinput.FileInput(path, inplace=True, backup='.bak') as f:
        for line in f:
            new_line = line
            for change in replace_dict:
                new_line = new_line.replace(change, replace_dict[change])
            print(new_line, end='')

def details_str_to_pd(string):
    """Convert table from printLJType turned into a string to pandas dataframe"""
    lines = string.strip().split('\n')
    lines = lines[2:]
    line_0 = lines[0].split()
    line_0.remove('GB')
    line_0.remove('GB')
    line_0[-2] = 'GB ' + line_0[-2]
    line_0[-1] = 'GB ' + line_0[-1]
    line_0.remove('LJ')
    line_0.remove('LJ')
    line_0[6] = 'LJ ' + line_0[6]
    line_0[7] = 'LJ ' + line_0[7]
    table_headers = line_0
    table_list = []
    for line in lines[1:]:
        table_list.append(line.split())

    table_np = np.array(table_list)
    table_pd = pd.DataFrame(table_np, columns=table_headers)
    table_pd['ATOM'] = pd.to_numeric(table_pd['ATOM'], downcast="integer")
    table_pd['RES'] = pd.to_numeric(table_pd['RES'], downcast="integer")
    table_pd['At.#'] = pd.to_numeric(table_pd['At.#'], downcast="integer")
    table_pd['LJ Radius'] = pd.to_numeric(table_pd['LJ Radius'], downcast="float")
    table_pd['LJ Depth'] = pd.to_numeric(table_pd['LJ Depth'], downcast="float")
    table_pd['Mass'] = pd.to_numeric(table_pd['Mass'], downcast="float")
    table_pd['Charge'] = pd.to_numeric(table_pd['Charge'], downcast="float")
    table_pd['GB Radius'] = pd.to_numeric(table_pd['GB Radius'], downcast="float")
    table_pd['GB Screen'] = pd.to_numeric(table_pd['GB Screen'], downcast="float")
    return table_pd

def ljtypes_str_to_pd(string):
    """Convert table from printLJTypes turned into a string to pandas dataframe"""
    lines = string.strip().split('\n')
    lines = lines[2:]
    table_list = []
    for line in lines:
        row_list = line.split()
        table_list.append([int(row_list[1]), int(row_list[-1])])
    
    table_np = np.array(table_list)
    table_pd = pd.DataFrame(table_np, columns=['ATOM', 'LJ Type'])
    table_pd['ATOM'] = pd.to_numeric(table_pd['ATOM'], downcast="integer")
    table_pd['LJ Type'] = pd.to_numeric(table_pd['LJ Type'], downcast="integer")
    return table_pd

def ljmatrix_str_to_pd(string):
    """Convert table from printLJMatrix turned into a string to pandas dataframe"""
    lines = string.strip().split('\n')
    lines = lines[2:]
    table_list = []
    for line in lines:
        row_list = line.split()
        table_list.append([row_list[1].strip().partition('[')[2].partition(']')[0], row_list[3].strip().partition('[')[2].partition(']')[0], 
            row_list[-2], row_list[-1]])

    table_np = np.array(table_list)
    table_pd = pd.DataFrame(table_np, columns=['Atom Type 1', 'Atom Type 2', 'R_ij', 'Eps_ij'])
    table_pd['Atom Type 1'] = pd.to_numeric(table_pd['Atom Type 1'], downcast="integer")
    table_pd['Atom Type 2'] = pd.to_numeric(table_pd['Atom Type 2'], downcast="integer")
    table_pd['R_ij'] = pd.to_numeric(table_pd['R_ij'], downcast="float")
    table_pd['Eps_ij'] = pd.to_numeric(table_pd['Eps_ij'], downcast="float")
    return table_pd