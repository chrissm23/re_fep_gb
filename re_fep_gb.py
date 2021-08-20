import numpy as np
import pandas as pd
import os
import subprocess
import fileinput
import shutil

from biopandas.pdb import PandasPdb

import setup.python_scripts.get_data as get_data
import setup.python_scripts.get_data as MutationReGbFe

residue_position, residue_mutation, chains, functions, windows = get_data.read_input()