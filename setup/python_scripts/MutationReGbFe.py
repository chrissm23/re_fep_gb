import numpy as np
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
        self.n_windows = control_dict['n_windows']
        self.windows = np.linspace(0, 1, self.n_windows).tolist()
        self.windows.reverse()
        self.intermediate = control_dict['intermediate']
        self.include_mut = control_dict['include_mut']
        if 'Rgb_modifiers' in control_dict.keys():
            self.gb_modifiers = control_dict['Rgb_modifiers']
        else:
            self.gb_modifiers = None
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
        self.missing_residues = None
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
        self.missing_residues = []
        # Get fist and last residue number of each chain and also missing residues
        for chain in self.chains_pdb:
            res_min = self.wt_structure.df['ATOM']['residue_number'].loc[self.wt_structure.df['ATOM']['chain_id'] == chain].min()
            self.first_residues.append(res_min)
            res_max = self.wt_structure.df['ATOM']['residue_number'].loc[self.wt_structure.df['ATOM']['chain_id'] == chain].max()
            self.last_residues.append(res_max)
            full_chain_set = set(range(res_min, res_max))
            chain_residues = self.wt_structure.df['ATOM']['residue_number'].loc[self.wt_structure.df['ATOM']['chain_id'] == chain].unique().tolist()
            self.missing_residues.append(full_chain_set - set(chain_residues))
        
        # Change from alphabetic chain ID's to numeric ones
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
            if self.chains_to_mutate_str == 'all':
                self.chains_to_mutate_int = [self.chain_to_number[x] for x in self.wt_structure.df['ATOM']['chain_id'].unique()]
            self.leap_first = [0]*len(self.first_residues)
            self.leap_last = [0]*len(self.last_residues)
            for i in range(len(self.first_residues)):
                if i == 0:
                    self.leap_first[i] = 1
                    self.leap_last[i] = self.last_residues[i] - (self.first_residues[i] - 1) - len(self.missing_residues[i])
                else:
                    self.leap_first[i] = self.leap_last[i-1] + 1
                    self.leap_last[i] = self.last_residues[i] - self.first_residues[i] + self.leap_first[i] - len(self.missing_residues[i])
            if not isinstance(self.residue_position, list):
                for i in self.chains_to_mutate_int:
                    if self.residue_position in self.missing_residues[i]:
                        raise Exception(f"No residue {self.residue_position} "
                        f"found in chain {list(self.chain_to_number.keys())[list(self.chain_to_number.values()).index(i)]} to mutate.")
                self.leap_residue_position = [self.residue_position - self.first_residues[i] + self.leap_first[i] - 
                    sum(x < self.residue_position for x in self.missing_residues[i]) for i in range(len(self.leap_first))
                    if i in self.chains_to_mutate_int]
            else:
                self.leap_residue_position = [self.residue_position[i] - self.first_residues[self.chains_to_mutate_int[i]] + 
                    self.leap_first[self.chains_to_mutate_int[i]] - 
                    sum(x < self.residue_position[i] for x in self.missing_residues[self.chains_to_mutate_int[i]]) 
                    for i in range(len(self.residue_position))]
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
        self.mt_structure.df['ATOM'].loc[
            (self.mt_structure.df['ATOM'].chain_id == chain) & 
            (self.mt_structure.df['ATOM'].residue_number == residue_position), "residue_name"] = residue_mutant

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

            # Write WT PDB
            self.wt_pdb_path = 'setup/parms_n_pdbs/pdbs/wt.pdb'
            self.wt_structure.to_pdb(path=self.wt_pdb_path, records=None, gz=False, append_newline=True)

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

            # Write WT PDB
            self.wt_pdb_path = 'setup/parms_n_pdbs/pdbs/wt.pdb'
            self.wt_structure.to_pdb(path=self.wt_pdb_path, records=None, gz=False, append_newline=True)
        
        self.mt_pdb_path = 'setup/parms_n_pdbs/pdbs/mt.pdb'
        self.mt_structure.to_pdb(path=self.mt_pdb_path, records=None, gz=False, append_newline=True)

    def create_parms(self):
        """Creates parameter files of WT and mutant"""
        create_parm.create_og_parms(self.wt_pdb_path, self.mt_pdb_path) # Create parameter files from original WT and mutant structures
        if self.gb_modifiers is not None:
            create_parm.modify_og_GBRadius(self.gb_modifiers, self.include_mut) # Modify original GB radius
            shutil.copyfile('setup/recalculate/wt_0.parm7', 'setup/parms_n_pdbs/parms/parms_windows/wt_0.parm7')
            shutil.copyfile('setup/recalculate/mt_0.parm7', 'setup/parms_n_pdbs/parms/parms_windows/mt_0.parm7')
        else:
            shutil.copyfile('setup/parms_n_pdbs/parms/parms_windows/wt_0_og.parm7', 'setup/parms_n_pdbs/parms/parms_windows/wt_0.parm7')
            shutil.copyfile('setup/parms_n_pdbs/parms/parms_windows/mt_0_og.parm7', 'setup/parms_n_pdbs/parms/parms_windows/mt_0.parm7')

        print("Creating intermediate parameter files...")
        # Create intermediate parameter files
        create_parm.create_intermediate_parms(self.functions, self.windows, self.leap_residue_position, self.intermediate, self.include_mut)
        print("Intermediate parameters finished.")

    def create_RE_n_equil_files(self):
        """Creates necessary files for hamiltonian replica exchange using the intermediate parameter files"""
        # Get names of parameter files and order them in decreasing order for WT and concatenate increasing order of MT
        def sortParmPaths_numerically(parm_path):
            return int(parm_path[3:-6])
        
        parm_files = [x for x in os.listdir('setup/parms_n_pdbs/parms/parms_windows') if os.path.isfile(f'setup/parms_n_pdbs/parms/parms_windows/{x}')]
        parm_files_wt = [x for x in parm_files if x[:2]=='wt']
        parm_files_mt = [x for x in parm_files if x[:2]=='mt']
        parm_files_wt.sort(key=sortParmPaths_numerically)
        parm_files_mt.sort(key=sortParmPaths_numerically)
        parm_files_wt_str = '\n'.join(parm_files_wt)
        parm_files_mt_str = '\n'.join(parm_files_mt)

        equilibration_dir = 'FE/equilibration'
        re_dir = 'FE/RE'
        sasa_dir = 'FE/SASA'

        def copy_re_equil_files(wt_or_mt):
            # Check which files are being copied
            if wt_or_mt == 'WT':
                equil_dir = equilibration_dir + '/WT'
                rep_dir = re_dir + '/WT'
                surf_dir = sasa_dir + '/WT'
                parm_files_topology = parm_files_wt
                parm_files_topology_str = parm_files_wt_str
            elif wt_or_mt == 'MT':
                equil_dir = equilibration_dir + '/MT'
                rep_dir = re_dir + '/MT'
                surf_dir = sasa_dir + '/MT'
                parm_files_topology = parm_files_mt
                parm_files_topology_str = parm_files_mt_str
            # Create directories for equilibration and copy templates to them
            if not os.path.exists(f'{equil_dir}'):
                os.makedirs(f'{equil_dir}')
            shutil.copyfile('setup/tmpls/equilibration_tmpls/heat.in', f'{equil_dir}/heat.in')
            shutil.copyfile('setup/tmpls/equilibration_tmpls/equilibration.in', f'{equil_dir}/equilibration.in')
            shutil.copyfile(f'setup/parms_n_pdbs/parms/parms_windows/{parm_files_topology[0]}', f'{equil_dir}/topology.parm7')
            shutil.copyfile(f'FE/minimization/{wt_or_mt}/minimization.rst7', f'{equil_dir}/minimization.rst7')
            if wt_or_mt != 'MT' or self.include_mut:
                shutil.copyfile('setup/tmpls/cpptraj_tmpls/create_snapshots.sh', f'{equil_dir}/create_snapshots.sh')
                get_data_n_general.make_executable(f'{equil_dir}/create_snapshots.sh')
                shutil.copyfile('setup/tmpls/equilibration_tmpls/sample.in', f'{equil_dir}/sample.in')
            # Copy files to RE directory
            shutil.copyfile('setup/tmpls/re_tmpls/groupfile.ref', f'{rep_dir}/groupfile.ref')
            shutil.copyfile('setup/tmpls/re_tmpls/mdin.ref', f'{rep_dir}/mdin.ref')
            shutil.copyfile('setup/tmpls/re_tmpls/hamiltonians.tmpl', f'{rep_dir}/hamiltonians.dat')
            replace_dict = {
                "%files%": parm_files_topology_str
            }
            get_data_n_general.replace_in_file(f'{rep_dir}/hamiltonians.dat', replace_dict)
            shutil.copyfile('setup/tmpls/re_tmpls/generate_remd_inputs.sh', f'{rep_dir}/generate_remd_inputs.sh')
            replace_dict_genremd = {
                '%wt_or_mt%': wt_or_mt
            }
            get_data_n_general.replace_in_file(f'{rep_dir}/generate_remd_inputs.sh', replace_dict_genremd)
            get_data_n_general.make_executable(f'{rep_dir}/generate_remd_inputs.sh')
            if wt_or_mt != 'MT' or self.include_mut:
                subprocess.call(f'{rep_dir}/generate_remd_inputs.sh')
            for x in parm_files_topology:
                shutil.copyfile(f'setup/parms_n_pdbs/parms/parms_windows/{x}', f'{rep_dir}/{x}')
            # Copy files to SASA directory
            shutil.copyfile('setup/tmpls/sasa_tmpls/sasa.in', f'{surf_dir}/sasa.in')
            shutil.copyfile(f'setup/parms_n_pdbs/parms/parms_windows/{wt_or_mt.lower()}_0.parm7', f'{surf_dir}/topology.parm7')
        
        for x in ['WT', 'MT']:
            copy_re_equil_files(x)
                 

    def create_FE_dir(self):
        """Creates directory with all required files and scripts to setup and run equilibration and replica exchange"""
        print("Generating output files and directories...")
        # Create directories
        fe_dir = 'FE'
        if not os.path.exists(fe_dir):
            os.makedirs(fe_dir)
        
        # Path of directories
        equilibration_dir = 'FE/equilibration'
        equil_wt_dir = 'FE/equilibration/WT'
        equil_mt_dir = 'FE/equilibration/MT'

        re_dir = 'FE/RE'
        re_wt_dir = 'FE/RE/WT'
        re_mt_dir = 'FE/RE/MT'

        minimization_dir = 'FE/minimization'
        min_wt_dir = 'FE/minimization/WT'
        min_mt_dir = 'FE/minimization/MT'

        sasa_dir = 'FE/SASA'
        sasa_wt_dir = 'FE/SASA/WT'
        sasa_mt_dir = 'FE/SASA/MT'

        for x in [equilibration_dir, equil_wt_dir, equil_mt_dir, re_dir, re_wt_dir, re_mt_dir, 
            minimization_dir, min_wt_dir, min_mt_dir, sasa_dir, sasa_wt_dir, sasa_mt_dir]:
            if not os.path.exists(x):
                os.makedirs(x)

        # Copy files for minimization of WT and MT
        for x,y in zip([min_wt_dir, min_mt_dir], ['wt', 'mt']):
            shutil.copyfile(f'setup/parms_n_pdbs/parms/parms_windows/{y}_0.parm7', f'{x}/topology.parm7')
            shutil.copyfile(f'setup/parms_n_pdbs/parms/rst_windows/{y}_0.rst7', f'{x}/coordinates.rst7')
            shutil.copyfile('setup/tmpls/minimize_tmpl/minimization.in', f'{x}/minimization.in')
            shutil.copyfile('setup/tmpls/minimize_tmpl/minimize.sh', f'{x}/minimize.sh')
            replace_dict_min = {
                '%wt_or_mt%': y.upper()
            }
            get_data_n_general.replace_in_file(f'{x}/minimize.sh', replace_dict_min)
            # Minimize WT and WT
            get_data_n_general.make_executable(f'{x}/minimize.sh')
            subprocess.call(f'{x}/minimize.sh')

        # Copy slurm template for equilibration
        for wt_or_mt in ['WT', 'MT']:
            shutil.copyfile('setup/tmpls/free_energy_tmpls/submit_equilibration.tmpl', f'FE/equilibration/submit_equilibration_{wt_or_mt}.slurm')
            replace_dict_equil = {
                '%residue_number%': str(self.residue_position),
                '%aminoacid_mutant%': str(self.residue_mutant),
                '%wt_or_mt%': wt_or_mt
            }
            get_data_n_general.replace_in_file(f'FE/equilibration/submit_equilibration_{wt_or_mt}.slurm', replace_dict_equil)
        
        # Copy slurm template for remd
        n_replicas = self.n_windows
        n_nodes = n_replicas//8 + (n_replicas%8 > 0)
        n_replica_by_nodes = 8
        for wt_or_mt in ['WT', 'MT']:
            shutil.copyfile('setup/tmpls/free_energy_tmpls/submit_remd.tmpl', f'FE/RE/submit_remd_{wt_or_mt}.slurm')
            replace_dict_remd = {
                '%residue_number%': str(self.residue_position),
                '%aminoacid_mutant%': str(self.residue_mutant),
                '%n_nodes%': str(n_nodes),
                '%n_replicas%': str(n_replicas),
                '%n_replica_by_nodes%': str(n_replica_by_nodes),
                '%wt_or_mt%': wt_or_mt
            }
            get_data_n_general.replace_in_file(f'FE/RE/submit_remd_{wt_or_mt}.slurm', replace_dict_remd)
        
        # Copy bash file to run SASA
        shutil.copyfile('setup/tmpls/free_energy_tmpls/surface_area.sh', 'FE/SASA/surface_area.sh')
        get_data_n_general.make_executable('FE/SASA/surface_area.sh')
        # Copy python scripts to calculate free energy differences
        shutil.copyfile('setup/python_scripts/FE_diff.py', 'FE/FE_diff.py')
        shutil.copyfile('setup/python_scripts/fe_recalculate.py', 'FE/fe_recalculate.py')
        
        self.create_RE_n_equil_files()

        print("Done.")