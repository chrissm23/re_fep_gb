from biopandas.pdb import PandasPdb
import numpy as np
import pandas as pd
import fileinput
from shutil import copyfile
import os
import subprocess

class MutantReGbFe:
    """Class with required information from structure and mutation as well as methods to calculate data to be filled into the templates for RE-FEP-GB"""
    chain_to_number = {
        1: "A",
        2: "B",
        3: "C",
        4: "D",
        5: "E",
        6: "F",
        7: "G",
        8: "H",
        9: "I",
        10: "J",
        11: "K",
        12: "L",
        13: "M",
        14: "N",
        15: "O",
        16: "P",
        17: "Q",
        18: "R",
        19: "S",
        20: "T",
        21: "U",
        22: "V",
        23: "W",
        24: "X",
        25: "Y",
        26: "Z"
    }
    
    def __init__(self) -> None:
        self.residue_position = None
        self.residue_mutation = None
        self.chains = None
        self.functions = None
        self.windows = None