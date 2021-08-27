import numpy as np
import pandas as pd
import fileinput
import shutil
import os
import subprocess

from biopandas.pdb import PandasPdb

import setup.python_scripts.create_parm as create_parm

class MutationReGbFe:
    """Class with required information from structure and mutation as well as methods to calculate data to be filled into the templates for RE-FEP-GB"""
    # Dictionary to change from alphabetic chain ID to numeric one
    chain_to_number = {
        "A": 0,
        "B": 1,
        "C": 2,
        "D": 3,
        "E": 4,
        "F": 5,
        "G": 6,
        "H": 7,
        "I": 8,
        "J": 9,
        "K": 10,
        "L": 11,
        "M": 12,
        "N": 13,
        "O": 14,
        "P": 15,
        "Q": 16,
        "R": 17,
        "S": 18,
        "T": 19,
        "U": 20,
        "V": 21,
        "W": 22,
        "X": 23,
        "Y": 24,
        "Z": 25
    }
    
    def __init__(self, control_dict, wt_structure_path) -> None:
        # Control parameters
        self.residue_position = control_dict['residue_position']
        self.residue_mutant = control_dict['residue_mutant']
        self.chains_to_mutate_str = control_dict['chains']
        self.chains_to_mutate_int = None
        self.functions = [control_dict['function_GB'], control_dict['function_ele'], control_dict['function_Rlj'], control_dict['function_epsilonlj']]
        self.windows = control_dict['windows']
        # Read PDBs into pandas dataframe and delete hydrogens
        self.wt_pdb_path = wt_structure_path
        self.mt_pdb_path = None
        self.wt_structure = PandasPdb().read_pdb(wt_structure_path) # Read PDB file into pandas dataframe
        self.mt_structure = PandasPdb().read_pdb(wt_structure_path) # Read PDB file into pandas dataframe which then will be change to mutant
        self.wt_structure.df['ATOM'] = self.wt_structure.df['ATOM'][self.wt_structure.df['ATOM']['element_symbol']!="H"]
        self.mt_structure.df['ATOM'] = self.mt_structure.df['ATOM'][self.mt_structure.df['ATOM']['element_symbol']!="H"]
        # Information about PDBs
        self.chains_pdb = None
        self.number_chains = None
        self.first_residues = None
        self.last_residues = None
        self.leap_first = None
        self.leap_last = None
        self.leap_residue_position = None

        self.get_structure_info()

    def get_structure_info(self):
        """Reads relevant structural information from PDB"""
        self.chains_pdb = self.wt_structure.df['ATOM']['chain_id'].unique().tolist()
        self.number_chains = self.wt_structure.df['ATOM']['chain_id'].nunique()
        self.first_residues = []
        self.last_residues = []
        # Get fist and last residue number of each chain
        for chain in self.chains_pdb:
            res_min = self.wt_structure.df['ATOM']['residue_number'].loc[self.wt_structure.df['ATOM']['chain_id'] == chain].min()
            self.first_residues.append(res_min)
            res_max = self.wt_structure.df['ATOM']['residue_number'].loc[self.wt_structure.df['ATOM']['chain_id'] == chain].max()
            self.last_residues.append(res_max)
        
        # Chainge from alphabetic chain ID's to numeric ones
        if self.chains_to_mutate_str != 'all' and self.chains_to_mutate_str != 'tripeptide':
            try:
                for x in self.chains_to_mutate_str:
                    self.chain_to_number[x]
            except:
                raise Exception('Chains must be selected using chain IDs from the PDB (e.g. A,B,C). If there are more that 26 chains please add chain IDs '
                'to dictionary \"chain_to_number\" in setup/python_scripts/MutationReGbFe.py')
            else:
                self.chains_to_mutate_int = [self.chain_to_number[x] for x in self.chains_to_mutate_str]

        # Get new mutation positions, initial and final residues for each chain after running leap
        if self.chains_to_mutate_str != 'tripeptide':
            self.leap_first = [0]*len(self.first_residues)
            self.leap_last = [0]*len(self.last_residues)
            for i in range(len(self.first_residues)):
                if i == 0:
                    self.leap_first[i] = 1
                    self.leap_last[i] = self.last_residues[i] - (self.first_residues[i] - 1)
                else:
                    self.leap_first[i] = self.leap_last[i-1] + 1
                    self.leap_last[i] = self.last_residues[i] - self.first_residues[i] + self.leap_first[i]
            if not isinstance(self.residue_position, list):
                self.leap_residue_position = [self.residue_position - self.first_residues[i] + self.leap_first[i] for i in range(len(self.leap_first))
                if i in self.chains_to_mutate_int]
            else:
                self.leap_residue_position = [self.residue_position[i] - self.first_residues[self.chains_to_mutate_int[i]] + 
                    self.leap_first[self.chains_to_mutate_int[i]] for i in range(len(self.residue_position))]
        else:
            self.leap_first = 1
            if self.residue_position == self.last_residues[0]:
                self.leap_last = 2
                self.leap_residue_position = 2
            elif self.residue_position == self.first_residues[0]:
                self.leap_last = 2
                self.leap_residue_position = 1
            else:
                self.leap_last = 3
                self.leap_residue_position = 2

    def change_residue_name(self, chain, residue_position, residue_mutant):
        """Changes the residue name in the PDB to the desired mutation"""
        self.mt_structure.df['ATOM']['residue_name'].loc[
            (self.mt_structure.df['ATOM']['chain_id'] == chain) & 
            (self.mt_structure.df['ATOM']['residue_number'] == residue_position)] = residue_mutant

    def delete_side_chain(self, chain, residue_position):
        """Deletes side chain of mutated residue for leap to autocomplete"""
        self.mt_structure.df['ATOM'] = self.mt_structure.df['ATOM'].drop(
            self.mt_structure.df['ATOM'][(self.mt_structure.df['ATOM']['chain_id'] == chain) & 
            (self.mt_structure.df['ATOM']['residue_number'] == residue_position) &
            (self.mt_structure.df['ATOM']['atom_name']!="CA") & 
            (self.mt_structure.df['ATOM']['atom_name']!="C") & (self.mt_structure.df['ATOM']['atom_name']!="O") &
            (self.mt_structure.df['ATOM']['atom_name']!="N")].index)
    
    def create_mutant(self):
        """Creates mutant PDB and, in case of tripeptide, also WT PDB"""
        # If tripeptide is selected, then only use one amino acid before and after the residue position of first chain
        if self.chains_to_mutate_str == 'tripeptide':
            # Change residue name
            self.change_residue_name("A", self.residue_position, self.residue_mutant)

            # Delete side chain
            self.delete_side_chain("A", self.residue_position)

            # Select only the three peptides from WT and mutant structures
            self.wt_structure.df['ATOM'] = self.wt_structure.df['ATOM'][
                (self.wt_structure.df['ATOM']['residue_number'] >= self.residue_position - 1) & 
                (self.wt_structure.df['ATOM']['residue_number'] <= self.residue_position + 1)]
            self.mt_structure.df['ATOM'] = self.mt_structure.df['ATOM'][
                (self.mt_structure.df['ATOM']['residue_number'] >= self.residue_position - 1) & 
                (self.mt_structure.df['ATOM']['residue_number'] <= self.residue_position + 1)]

            # Write WT tripeptide PDB
            self.wt_pdb_path = 'setup/parms_n_pdbs/pdbs/wt_tripeptide.pdb'
            self.wt_structure.to_pdb(path=self.wt_pdb_path, records=None, gz=False, append_newline=True)

        # If all is selected, mutate all chains
        elif self.chains_to_mutate_str == 'all':
            
            for chain in self.chains_pdb:
                self.change_residue_name(chain, self.residue_position, self.residue_mutant) # Change residue name
                self.delete_side_chain(chain, self.residue_position) # Delete side chain

        # If specific chains are given
        else:
            if isinstance(self.residue_position, list) and isinstance(self.residue_mutant, list): # If there are matching lists
                for chain, position, aminoacid in zip(self.chains_to_mutate_str, self.residue_position, self.residue_mutant):
                    self.change_residue_name(chain, position, aminoacid)
                    self.delete_side_chain(chain, position)
            elif isinstance(self.residue_position, list): # If there is only a matching list for position
                for chain, position in zip(self.chains_to_mutate_str, self.residue_position):
                    self.change_residue_name(chain, position, self.residue_mutant)
                    self.delete_side_chain(chain, position)
            elif isinstance(self.residue_mutant, list): # If there is only a matching list for mutant aminoacid
                for chain, aminoacid in zip(self.chains_to_mutate_str, self.residue_mutant):
                    self.change_residue_name(chain, self.residue_position, aminoacid)
                    self.delete_side_chain(chain, self.residue_position)
            else: # If there are single values for position and mutant aminoacid
                for chain in self.chains_to_mutate_str:
                    self.change_residue_name(chain, self.residue_position, self.residue_mutant)
                    self.delete_side_chain(chain, self.residue_position)
        
        self.mt_pdb_path = 'setup/parms_n_pdbs/pdbs/mt_tripeptide.pdb'
        self.mt_structure.to_pdb(path=self.mt_pdb_path, records=None, gz=False, append_newline=True)

    def create_parms(self):
        """Creates parameter files of WT and mutant"""
        create_parm.create_og_parms(self.wt_pdb_path, self.mt_pdb_path) # Create parameter files from original WT and mutant structures

        create_parm.create_intermediate_parms(self.functions, self.leap_residue_position) # Add other information required to create modified parameter files