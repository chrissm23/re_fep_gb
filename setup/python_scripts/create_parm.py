import numpy as np
import pandas as pd
import fileinput
import shutil
import os
import subprocess

from biopandas.pdb import PandasPdb

def create_og_parms(path_wt_pdb, path_mt_pdb):
    """Moves, fills out template and executes resulting bash file to create parameter files of the WT and mutant structures given as imput"""
    wt_bash_path = 'setup/leap/parms/gen_parm_wt.sh'
    mt_bash_path = 'setup/leap/parms/gen_parm_mt.sh'
    tmpl_path = 'setup/tmpls/leap_tmpls/gen_parm_pdb.tmpl'