import numpy as np
import pandas as pd
import os
import subprocess
import fileinput
import shutil

from biopandas.pdb import PandasPdb

"""
Collection of functions to read necessary input provided and calculate information required about the specific structure
"""

def read_input(self):
    """Reads control.txt and number name of PDB from outside setup directory"""