import numpy as np
import pandas as pd
import fileinput
import shutil
import os
import subprocess
import parmed

from biopandas.pdb import PandasPdb
from parmed.tools.actions import parm

import setup.python_scripts.get_data_n_general as get_data_n_general

def create_og_parms(path_wt_pdb, path_mt_pdb):
    """Moves, fills out template and executes resulting bash file to create parameter files of the WT and mutant structures given as imput"""
    # Destination paths
    wt_bash_path = 'setup/parms_n_pdbs/parms/gen_parm_wt.sh'
    mt_bash_path = 'setup/parms_n_pdbs/parms/gen_parm_mt.sh'
    tmpl_path = 'setup/tmpls/leap_tmpls/gen_parm_pdb.tmpl'

    # Output paths
    wt_parm_path = 'setup/parms_n_pdbs/parms/parms_windows/wt_0.parm7'
    wt_rst_path = 'setup/parms_n_pdbs/parms/rst_windows/wt_0.rst7'
    mt_parm_path = 'setup/parms_n_pdbs/parms/parms_windows/mt_0.parm7'
    mt_rst_path = 'setup/parms_n_pdbs/parms/rst_windows/mt_0.rst7'

    # Dictionaries of fields to replace in template
    replace_dict_wt = {
        "%pdb_path%": path_wt_pdb,
        "%name_structure_path_parm%": wt_parm_path,
        "%name_structure_path_rst%": wt_rst_path
    }
    replace_dict_mt = {
        "%pdb_path%": path_mt_pdb,
        "%name_structure_path_parm%": mt_parm_path,
        "%name_structure_path_rst%": mt_rst_path
    }
    
    # Copy templates
    shutil.copyfile(tmpl_path, wt_bash_path)
    shutil.copyfile(tmpl_path, mt_bash_path)
    # Fill out templates
    get_data_n_general.replace_in_file(wt_bash_path, replace_dict_wt)
    get_data_n_general.replace_in_file(mt_bash_path, replace_dict_mt)
    # Make bash files executable and run
    print("Creating WT and mutant parameter files...")
    get_data_n_general.make_executable(wt_bash_path)
    get_data_n_general.make_executable(mt_bash_path)
    subprocess.call(wt_bash_path)
    subprocess.call(mt_bash_path)
    print("Parameter files created.")

def get_new_dataframes(parmed_object, residue_mask):
    """Creates new LJ atom types and outputs new dataframes"""
    residue_details = get_data_n_general.details_str_to_pd(str(parmed.tools.printDetails(parmed_object, f':{residue_mask}')))
    atom_numbers = residue_details['ATOM'].tolist() # Get atom numbers of mutated residues
    structure_types = get_data_n_general.ljtypes_str_to_pd(str(parmed.tools.printLJTypes(parmed_object)))
    structure_types_groupby = structure_types.groupby('LJ Type')
    old_atoms_types = structure_types['LJ Type'].unique().tolist()
    for atom_type in old_atoms_types:
        atom_numbers_type = structure_types_groupby.get_group(atom_type)['ATOM'] # Get list of atom numbers in atom_numbers with atom_type
        # Change atom_type for atom_numbers_type

def create_intermediate_parms(functions, residue_position):
    """Creates intermediate parameter files using scaling of function on residues given by residue_position"""
    # Leap generated parameter files
    wt_parm_path = 'setup/parms_n_pdbs/parms/parms_windows/wt_0.parm7'
    wt_rst_path = 'setup/parms_n_pdbs/parms/rst_windows/wt_0.rst7'
    mt_parm_path = 'setup/parms_n_pdbs/parms/parms_windows/mt_0.parm7'
    mt_rst_path = 'setup/parms_n_pdbs/parms/rst_windows/mt_0.rst7'

    # Import to ParmEd
    wt_parmed = parmed.amber.AmberParm(wt_parm_path, wt_rst_path)
    mt_parmed = parmed.amber.AmberParm(mt_parm_path, mt_rst_path)

    #Create mutated residues masks
    if not isinstance(residue_position, list):
        residue_mask = f'{residue_position}'
    else:
        residue_mask_list = [str(x) for x in residue_position]
        residue_mask = '.'.join(residue_mask_list)

    # Get atom numbers, charge and GB Radius
    if not isinstance(residue_position, list):
        residue_details_wt = get_data_n_general.details_str_to_pd(str(parmed.tools.printDetails(wt_parmed, f':{residue_position}')))
        structure_types_wt = get_data_n_general.ljtypes_str_to_pd(str(parmed.tools.printLJTypes(wt_parmed)))
        residue_LJmatrix_wt = get_data_n_general.ljmatrix_str_to_pd(str(parmed.tools.printLJMatrix(wt_parmed, f':{residue_position}')))
        