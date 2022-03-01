import numpy as np
import pandas as pd
import fileinput
import os

import parmed as pmd
import pytraj as pt
from biopandas.pdb import PandasPdb

import setup.python_scripts.get_data_n_general as get_data_n_general
import setup.python_scripts.create_parm as create_parm
import setup.python_scripts.MutationReGbFe as MutationReGbFe

def read_gb_modifiers():
    """Read recalculate_gb.txt to obtain modifications to original GB radii"""


if __name__ == '__main__':
    # List input files and determine if everything needed is there
    files = [x for x in os.listdir() if os.path.isfile(x)] # List input files
    pdb_files = [x for x in files if x[-4:]==".pdb"]
    control_file = 'control.txt'
    if len(pdb_files) == 0:
        raise FileNotFoundError('No PDB file found.')
    if len(pdb_files) > 1:
        raise Exception('More than one PDB file found. Only one is acceptable.')
    if not control_file in files:
        raise FileNotFoundError('No control file \'control.txt\' found.')
    pdb_file = pdb_files[0]

    control_dict = get_data_n_general.read_input(control_file) # Get required parameters from control file
    control_dict['Rgb_modifiers'] = read_gb_modifiers()

    mutation = MutationReGbFe.MutationReGbFe(control_dict, pdb_file) # Get information about structure
    create_parm.modify_og_GBRadius(mutation.gb_modifiers, mutation.include_mut) # Modify original GB radius
    print("Modifying intermediate parameter files...")
    # Update GB radii in intermediate parameter files
    create_parm.create_intermediate_parms(mutation.functions, mutation.windows, mutation.leap_residue_position, 
        mutation.intermediate, mutation.include_mut, recalculation=True)
    print("Intermediate parameters finished.")

    # Recalculate using updated parameter files and already simulated REMD trajectories