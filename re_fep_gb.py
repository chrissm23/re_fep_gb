from genericpath import isfile
import numpy as np
import pandas as pd
import os
import subprocess
import fileinput
import shutil

from biopandas.pdb import PandasPdb

import setup.python_scripts.get_data as get_data
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

control_dict = get_data.read_input(control_file) # Get required parameters from control file

wt_structure = PandasPdb().read_pdb(pdb_file) # Read PDB file into pandas dataframe
get_data.create_mutant(wt_structure, control_dict['chains'])
mt_structure = PandasPdb().read_pdb('./setup/leap/mutant.pdb')
wt_structure = PandasPdb().read_pdb(pdb_file) # Read PDB file into pandas dataframe
mutation = MutationReGbFe(control_dict, wt_structure, mt_structure)
mutation.create_og_parms() # Create unmodified (end points) parameter files
mutation.create_inter_parms() # Create intermediate parameter files according to functions in control file
mutation.create_RE_files() # Create files required for replica exchage
mutation.create_free_energy_dir() # Create directory with all that is required to run on cow