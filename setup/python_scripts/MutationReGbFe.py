from genericpath import exists
import numpy as np
import pandas as pd
import fileinput
import shutil
import os
import subprocess

from biopandas.pdb import PandasPdb

import setup.python_scripts.get_data_n_general as get_data_n_general
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
        self.windows = control_dict['windows'].tolist()
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

        print("Creating intermediate parameter files...")
        create_parm.create_intermediate_parms(self.functions, self.windows, self.leap_residue_position) # Create intermediate parameter files
        print("Intermediate parameters finished.")

    def create_RE_n_equil_files(self):
        """Creates necessary files for hamiltonian replica exchange using the intermediate parameter files"""
        # Get names of parameter files and order them in decreasing order for WT and concatenate increasing order of MT
        def sortParmPaths_numerically(parm_path):
            return int(parm_path[3:-6])
        
        parm_files = [x for x in os.listdir('setup/parms_n_pdbs/parms/parms_windows') if os.path.isfile(f'setup/parms_n_pdbs/parms/parms_windows/{x}')]
        parm_files_wt = [x for x in parm_files if x[:2]=='wt']
        parm_files_mt = [x for x in parm_files if x[:2]=='mt']
        parm_files_wt.sort(key=sortParmPaths_numerically, reverse=True)
        parm_files_wt = [parm_files_wt[-1]] + parm_files_wt[:-1]
        parm_files_mt.sort(key=sortParmPaths_numerically)
        parm_files_mt = parm_files_mt[1:] + [parm_files_mt[0]]
        parm_files = parm_files_wt + parm_files_mt
        parm_files_str = '\n'.join(parm_files)

        equilibration_dir = 'FE/equilibration'
        re_dir = 'FE/RE'

        # Create directories for equilibration and copy templates to them
        for i in range(len(parm_files)):
            if not os.path.exists(f'{equilibration_dir}/{i}'):
                os.makedirs(f'{equilibration_dir}/{i}')
            shutil.copyfile('setup/tmpls/equilibration_tmpl/heat.in', f'{equilibration_dir}/{i}/heat.in')
            shutil.copyfile('setup/tmpls/equilibration_tmpl/equilibration.in', f'{equilibration_dir}/{i}/equilibration.in')
            shutil.copyfile(f'setup/parms_n_pdbs/parms/parms_windows/{parm_files[i]}', f'{equilibration_dir}/{i}/topology.parm7')
            if i == 0:
                shutil.copyfile('FE/minimization/minimization.rst7', f'{equilibration_dir}/{i}/minimization.rst7')

        # Copy files to RE directory
        shutil.copyfile('setup/tmpls/re_tmpls/groupfile.ref', f'{re_dir}/groupfile.ref')
        shutil.copyfile('setup/tmpls/re_tmpls/mdin.ref', f'{re_dir}/mdin.ref')
        shutil.copyfile('setup/tmpls/re_tmpls/hamiltonians.tmpl', f'{re_dir}/hamiltonians.dat')
        replace_dict = {
            "%files%": parm_files_str
        }
        get_data_n_general.replace_in_file(f'{re_dir}/hamiltonians.dat', replace_dict)
        shutil.copyfile('setup/tmpls/re_tmpls/generate_remd_inputs.sh', f'{re_dir}/generate_remd_inputs.sh')
        get_data_n_general.make_executable(f'{re_dir}/generate_remd_inputs.sh')
        subprocess.call(f'{re_dir}/generate_remd_inputs.sh')
        for x in os.listdir('setup/parms_n_pdbs/parms/parms_windows'):
            shutil.copyfile(f'setup/parms_n_pdbs/parms/parms_windows/{x}', f'{re_dir}/{x}')

    def create_FE_dir(self):
        """Creates directory with all required files and scripts to setup and run equilibration and replica exchange"""
        print("Generating output files and directories...")
        # Create directories
        fe_dir = 'FE'
        if not os.path.exists(fe_dir):
            os.makedirs(fe_dir)
        
        equilibration_dir = 'FE/equilibration'
        re_dir = 'FE/RE'
        minimization_dir = 'FE/minimization'

        for x in [equilibration_dir, re_dir, minimization_dir]:
            if not os.path.exists(x):
                os.makedirs(x)

        # Copy files for minimization of WT
        shutil.copyfile('setup/tmpls/minimize_tmpl/minimization.in', f'{minimization_dir}/minimization.in')
        shutil.copyfile('setup/parms_n_pdbs/parms/parms_windows/wt_0.parm7', f'{minimization_dir}/topology.parm7')
        shutil.copyfile('setup/parms_n_pdbs/parms/rst_windows/wt_0.rst7', f'{minimization_dir}/coordinates.rst7')
        shutil.copyfile('setup/tmpls/minimize_tmpl/minimize.sh', f'{minimization_dir}/minimize.sh')
        # Minimize WT
        get_data_n_general.make_executable(f'{minimization_dir}/minimize.sh')
        subprocess.call(f'{minimization_dir}/minimize.sh')

        # Copy slurm template for equilibration
        shutil.copyfile('setup/tmpls/free_energy_tmpls/submit_equilibration.tmpl', 'FE/submit_equilibration.slurm')
        replace_dict_equil = {
            '%residue_number%': self.residue_position,
            '%aminoacid_mutant%': self.residue_mutant
        }
        get_data_n_general.replace_in_file('FE/submit_equilibration.slurm', replace_dict_equil)
        
        # Copy slurm template for remd
        n_replicas = 2*self.windows - 1
        n_nodes = n_replicas//8 + (n_replicas%8 > 0)
        n_replica_by_nodes = 8
        shutil.copyfile('setup/tmpls/free_energy_tmpls/submit_remd.tmpl', 'FE/submit_remd.slurm')
        replace_dict_remd = {
            '%residue_number%': self.residue_position,
            '%aminoacid_mutant%': self.residue_mutant,
            '%n_nodes%': n_nodes,
            '%n_replicas%': n_replicas,
            '%n_replica_by_nodes%': n_replica_by_nodes
        }
        get_data_n_general.replace_in_file('FE/submit_remd.slurm', replace_dict_remd)
        
        self.create_RE_n_equil_files()
        print("Done.")