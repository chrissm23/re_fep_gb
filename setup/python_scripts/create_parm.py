import numpy as np
import pandas as pd
import fileinput
import shutil
import os
import subprocess

from pandas.core.groupby.groupby import get_groupby
import parmed

from biopandas.pdb import PandasPdb
from parmed import amber
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

def get_new_LJParms(parmed_object, residue_mask, functions, windows):
    """Creates new LJ atom types and outputs new AmberParm bojects with modified LJ  matrix"""
    residue_details = get_data_n_general.details_str_to_pd(str(parmed.tools.printDetails(parmed_object, f':{residue_mask}')))
    atom_numbers = residue_details['ATOM'].tolist() # Get atom numbers of mutated residues

    # Get residue numbers grouped by LJ type
    structure_types_old = get_data_n_general.ljtypes_str_to_pd(str(parmed.tools.printLJTypes(parmed_object)))
    mutated_types_old = structure_types_old[structure_types_old['ATOM'].isin(atom_numbers)]
    atom_types_old = mutated_types_old['LJ Type'].unique().tolist()
    atomtype_old_group = mutated_types_old.groupby('LJ Type')
    masks_by_atomtype_int = []
    for atom_type_old in atom_types_old:
        masks_by_atomtype_int.append(atomtype_old_group.get_group(atom_type_old)['ATOM'].tolist())
    masks_by_atomtype_str = [[str(masks_by_atomtype_int[i][j]) for j in range(len(masks_by_atomtype_int[i]))] for i in range(len(masks_by_atomtype_int))]
    masks_by_atomtype_str = [','.join(x) for x in masks_by_atomtype_str]

    # Create new LJ atom types
    for mask in masks_by_atomtype_str:
        parmed.tools.addLJType(parmed_object, f'@{mask}').execute()
    
    # Group not mutating atoms by LJ atom type
    nonmutating_types = structure_types_old[~structure_types_old['ATOM'].isin(atom_numbers)]
    atom_types_old_nonmut = nonmutating_types['LJ Type'].unique().tolist()
    atomtype_group = nonmutating_types.groupby('LJ Type')
    masks_by_nonmutated_atomtype_int = []
    for atom_type_old in atom_types_old_nonmut:
        masks_by_nonmutated_atomtype_int.append(atomtype_group.get_group(atom_type_old)['ATOM'].tolist())
    masks_by_nonmutated_atomtype_str = [[str(masks_by_nonmutated_atomtype_int[i][j]) for j in range(len(masks_by_nonmutated_atomtype_int[i]))] 
        for i in range(len(masks_by_nonmutated_atomtype_int))]
    masks_by_nonmutated_atomtype_str = [','.join(x) for x in masks_by_nonmutated_atomtype_str]

    # Deep copy parmed object and modify LJ matrix elements
    new_parms = [parmed.amber.AmberParm.from_structure(parmed_object, copy=True) for i in range(len(windows))]
    new_ljmatrix = get_data_n_general.ljmatrix_str_to_pd(str(parmed.tools.printLJMatrix(parmed_object, f':{residue_mask}')))
    #print(str(parmed.tools.printLJMatrix(parmed_object, f':{residue_mask}')))
    for i in range(len(windows)):
        multiplier_R = get_data_n_general.get_multiplier(windows[i], functions[-2], truncate=False)
        multiplier_eps = get_data_n_general.get_multiplier(windows[i], functions[-1], truncate=False)
        for mask_mutant in masks_by_atomtype_str:
            atom_type_mut = get_data_n_general.ljtypes_str_to_pd(str(parmed.tools.printLJTypes(parmed_object, f'@{mask_mutant}')))['LJ Type'][0]
            for mask_nonmutant in masks_by_nonmutated_atomtype_str:
                atom_type_nonmut = get_data_n_general.ljtypes_str_to_pd(str(parmed.tools.printLJTypes(parmed_object, f'@{mask_nonmutant}')))['LJ Type'][0]

                # Get LJ Radius for atom type pair and calculate new LJ Radius
                R_ij = new_ljmatrix[((new_ljmatrix['Atom Type 1'] == atom_type_mut) & (new_ljmatrix['Atom Type 2'] == atom_type_nonmut)) | 
                    ((new_ljmatrix['Atom Type 2'] == atom_type_mut) & (new_ljmatrix['Atom Type 1'] == atom_type_nonmut))]['R_ij'].iloc[0]
                new_R_ij = multiplier_R*R_ij
                # Get LJ epsilon for atom type pair and calculate new LJ epsilon
                Eps_ij = new_ljmatrix[((new_ljmatrix['Atom Type 1'] == atom_type_mut) & (new_ljmatrix['Atom Type 2'] == atom_type_nonmut)) | 
                    ((new_ljmatrix['Atom Type 2'] == atom_type_mut) & (new_ljmatrix['Atom Type 1'] == atom_type_nonmut))]['Eps_ij'].iloc[0]
                new_Eps_ij = multiplier_eps*Eps_ij

                # Change LJ parameters for atom type pair
                parmed.tools.changeLJPair(new_parms[i], f'@{mask_mutant}', f'@{mask_nonmutant}', f'{new_R_ij}', f'{new_Eps_ij}').execute()
        #print(str(parmed.tools.printLJMatrix(new_parms[i], f':{residue_mask}')))

    return new_parms

def get_new_Parms(parms_list, residue_mask, propty, functional, windows, truncate):
    """Changes charge of mutating residues according to functional for the different windows in each parameter file"""
    propty_pd_to_parmed = {
        'Charge': 'CHARGE',
        'GB Radius': 'RADII',
        'GB Screen': 'SCREEN'
    }
    GBRadius_minimum = 0.2
    for i in range(len(windows)):
        residue_details = get_data_n_general.details_str_to_pd(str(parmed.tools.printDetails(parms_list[i], f':{residue_mask}')))
        atom_numbers = residue_details['ATOM'].tolist() # Get atom numbers of mutated residues
        multiplier = get_data_n_general.get_multiplier(windows[i], functional, truncate)
        for atom in atom_numbers:
            value = residue_details[residue_details['ATOM'] == atom][propty].iloc[0]
            # Take into account the lowest limit of 0.1 Angstrom for GB Radius
            if propty != 'GB Radius':
                new_value = multiplier*value
            elif propty == 'GB Radius':
                new_value = multiplier*(value-GBRadius_minimum) + GBRadius_minimum
            parmed.tools.change(parms_list[i], propty_pd_to_parmed[propty], f'@{atom}', f'{new_value}').execute()
        #print(str(parmed.tools.printDetails(parms_list[i], f':{residue_mask}')))
    return parms_list

def create_intermediate_parms(functions, windows, residue_position):
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
        residue_mask = ','.join(residue_mask_list)
    residue_mask = residue_mask + '&!@CA,C,O,N'

    # Create parms with modified LJ matrix according to windows and functions
    wt_parms_LJ = get_new_LJParms(wt_parmed, residue_mask, functions[-2:], windows[:-1])
    mt_parms_LJ = get_new_LJParms(mt_parmed, residue_mask, functions[-2:], windows[:-1])

    # Change charge of mutating residues according to windows and functions
    wt_parms_ele = get_new_Parms(wt_parms_LJ, residue_mask, 'Charge', functions[1], windows[:-1], truncate=True)
    mt_parms_ele = get_new_Parms(mt_parms_LJ, residue_mask, 'Charge', functions[1], windows[:-1], truncate=True)

    # Change GB Radius of mutating residues according to windws and functions
    wt_parms_GB = get_new_Parms(wt_parms_ele, residue_mask, 'GB Radius', functions[1], windows[:-1], truncate=False)
    mt_parms_GB = get_new_Parms(mt_parms_ele, residue_mask, 'GB Radius', functions[1], windows[:-1], truncate=False)

    for i in range(len(wt_parms_GB)):
        #print(parmed.tools.printDetails(wt_parms_GB[i], f':{residue_mask}'))
        parmed.tools.outparm(wt_parms_GB[i], f'setup/parms_n_pdbs/parms/parms_windows/wt_{i+1}.parm7').execute()
        parmed.tools.outparm(mt_parms_GB[i], f'setup/parms_n_pdbs/parms/parms_windows/mt_{i+1}.parm7').execute()
    #print(parmed.tools.printDetails(wt_parmed, f':{residue_mask}'))