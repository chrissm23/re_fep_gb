import numpy as np
import pandas as pd
import os

import setup.python_scripts.get_data_n_general as get_data_n_general
import setup.python_scripts.MutationReGbFe as MutationReGbFe

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

# Create required directories to store parameter files, bash scripts and pdbs
parms_n_pdbs_dir = 'setup/parms_n_pdbs/'
parms_dir = 'setup/parms_n_pdbs/parms/'
pdbs_dir = 'setup/parms_n_pdbs/pdbs/'
parms_windows_dir = 'setup/parms_n_pdbs/parms/parms_windows'
rst_windows_dir = 'setup/parms_n_pdbs/parms/rst_windows'
for x in [parms_n_pdbs_dir, parms_dir, pdbs_dir, parms_windows_dir, rst_windows_dir]:
    if not os.path.exists(x):
        os.makedirs(x)

mutation = MutationReGbFe.MutationReGbFe(control_dict, pdb_file)
mutation.create_mutant() # Create mutant PDB and WT PDB if tripeptide was selected
mutation.create_parms() # Create unmodified (end points) parameter files
mutation.create_FE_dir() # Create directory with all that is required to run on cow
# Restart output files from heating must match those of group file