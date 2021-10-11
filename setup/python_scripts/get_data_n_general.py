import numpy as np
import pandas as pd
import os
import fileinput

"""
Collection of functions to read necessary input provided and calculate information required about the specific structure.
"""

def read_input(control_file):
    """Reads control.txt and outputs a dictionary with control parameters."""
    control_dict = {} # Dictionary to store parameters
    possible_functions = ['constant', 'linear', 'quadratic', 'sqrt', 'root6']
    possible_intermediates = ['GLY', 'ALA']
    with open(control_file, 'r') as cf:
        lines = cf.readlines()
        for line in lines:
            control = line.split('=')
            # Create control dictionary
            parameter_name = control[0].strip()
            if parameter_name == 'n_windows_per_topology': # Get the sequence of windows
                control_dict['n_windows'] = int(control[1].strip())
            elif parameter_name == 'res_pos': # Get position of the mutation(s)
                residue_position = control[1].strip().split(',')
                if len(residue_position) > 1:
                    control_dict['residue_position'] = [int(x) for x in residue_position]
                elif len(residue_position) == 1 and residue_position[0]:
                    control_dict['residue_position'] = int(residue_position[0])
            elif parameter_name == 'intermediate':
                intermediate = control[1].strip()
                if intermediate in possible_intermediates:
                    control_dict['intermediate'] = intermediate
                else:
                    raise Exception('Intermediate state not recognized')
            elif parameter_name == 'include_mut':
                control_dict['include_mut'] = bool(int(control[1].strip()))
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
            elif parameter_name == 'function_GB':
                functional = control[1].strip()
                if functional in possible_functions:
                    control_dict['function_GB'] = functional
                else:
                    raise Exception('Function for function_GB not recognized.')
            elif parameter_name == 'function_ele':
                functional = control[1].strip()
                if functional in possible_functions:
                    control_dict['function_ele'] = functional
                else:
                    raise Exception('Function for function_ele not recognized.')
            elif parameter_name == 'function_Rlj':
                functional = control[1].strip()
                if functional in possible_functions:
                    control_dict['function_Rlj'] = functional
                else:
                    raise Exception('Function for function_Rlj not recognized.')
            elif parameter_name == 'function_epsilonlj':
                functional = control[1].strip()
                if functional in possible_functions:
                    control_dict['function_epsilonlj'] = functional
                else:
                    raise Exception('Function for function_epsilonlj not recognized.')
    
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
    parameters = ['n_windows', 'residue_position', 'residue_mutant', 'chains', 'function_GB', 'function_ele', 'function_Rlj', 'function_epsilonlj']
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
        if len(row_list) == 8:
            table_list.append([row_list[1].strip().partition('[')[2].partition(']')[0], row_list[3].strip().partition('[')[2].partition(']')[0], 
                row_list[4], row_list[5], row_list[6], row_list[7]])

    table_np = np.array(table_list)
    table_pd = pd.DataFrame(table_np, columns=['Atom Type 1', 'Atom Type 2', 'A_ij', 'B_ij', 'R_ij', 'Eps_ij'])
    table_pd['Atom Type 1'] = pd.to_numeric(table_pd['Atom Type 1'], downcast="integer")
    table_pd['Atom Type 2'] = pd.to_numeric(table_pd['Atom Type 2'], downcast="integer")
    table_pd['A_ij'] = pd.to_numeric(table_pd['A_ij'], downcast="float")
    table_pd['B_ij'] = pd.to_numeric(table_pd['B_ij'], downcast="float")
    table_pd['R_ij'] = pd.to_numeric(table_pd['R_ij'], downcast="float")
    table_pd['Eps_ij'] = pd.to_numeric(table_pd['Eps_ij'], downcast="float")
    return table_pd

def get_multiplier(window, functional, truncate=False, ele_or_GB=None):
    """Calculates multiplier of the parameter according to window from 0 to 1 and functional. If truncate=True multiplier goes to 0 faster than window"""
    y_1 = 1
    x_1 = 1
    a = 1
    if truncate == True:
        x_0 = 0.2
        a = 1/(x_0^2 - 2*x_0 + 1)
        b = -2*a*x_0
        c = 1-(a+b)
    else:
        x_0 = 0
        b = 0
        c = 0
    if window < x_0:
        multiplier = 0
    else:
        if functional == 'constant':
            multiplier = 1
        if functional == 'linear':
            multiplier = y_1/(x_1-x_0)*(window - x_0)
        if functional == 'quadratic':
            multiplier = window*(a*window + b) + c
        if functional == 'sqrt':
            multiplier = np.sqrt((window - x_0)/(1 - x_0))
        if functional == 'root6':
            if truncate == True:
                raise Exception('I have not yet implemented the truncate function with the 6th root scaling.')
            else:
                multiplier = np.power(window, 1/6)

    return multiplier