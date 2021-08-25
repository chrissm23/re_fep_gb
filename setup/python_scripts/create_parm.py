import numpy as np
import pandas as pd
import fileinput
import shutil
import os
import subprocess

from biopandas.pdb import PandasPdb

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