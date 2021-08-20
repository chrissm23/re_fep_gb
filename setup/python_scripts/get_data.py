from biopandas.pdb import PandasPdb
import pandas as pd
import fileinput
from shutil import copyfile
import os
import subprocess

"""
Collection of functions to read necessary input provided and calculate information required about the specific structure
"""

def read_input(self):
    """Reads control.txt and number name of PDB from outside setup directory"""