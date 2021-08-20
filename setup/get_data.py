from biopandas.pdb import PandasPdb
import pandas as pd
import fileinput
from shutil import copyfile
import os
import subprocess

def read_input(self):
    """Reads control.txt and gets name of PDB from outside setup directory"""