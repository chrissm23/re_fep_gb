import numpy as np
import pandas as pd
import fileinput
import os
import shutil
import csv

import parmed as pmd
import pytraj as pt
from biopandas.pdb import PandasPdb

import setup.python_scripts.get_data_n_general as get_data_n_general
import setup.python_scripts.create_parm as create_parm
import setup.python_scripts.MutationReGbFe as MutationReGbFe

def read_gb_modifiers():
    """Read recalculate_gb.txt to obtain modifications to original GB radii"""
    control_file = 'recalculate_gb.txt'
    if not os.path.isfile(f'./{control_file}'):
        raise FileNotFoundError('File with GB modifications should be called recalculate_gb.txt')
    with open(control_file, 'r') as cf:
        lines = cf.readlines()
        for line in lines:
            control = line.split('=')
            parameter_name = control[0].strip()
            if parameter_name == 'gb_radius_mod':
                modifiers = control[1].strip().split(';')
                if len(modifiers) > 1:
                    gb_modification = modifiers
                if len(modifiers) == 1 and modifiers[0]:
                    gb_modification = modifiers[0]
            else:
                raise Exception('Only allowed control key is \"gb_radius_mod = <mask1>+0.1;...\"')
    
    return gb_modification

def get_force_constant(wt_or_mt):
    """Read mdin.ref file and obtain force constant used for restraints"""
    if not os.path.isfile(f'./FE/RE/{wt_or_mt}/mdin.ref'):
        raise FileNotFoundError(f'File ./FE/RE/{wt_or_mt}/mdin.ref not found')
    with open(f'./FE/RE/{wt_or_mt}/mdin.ref', 'r') as mdin:
        lines = mdin.readlines()
        for line in lines:
            elements = line.split(',')
            for element in elements:
                if element.strip()[0:12] == 'restraint_wt':
                    force_constant = float(element.split('=')[1].strip())
    return force_constant


def print_remlog(remlog_np, snapshot_dir):
    """Writes the recalculated energies in the same format as the remlog file to use with FE_diff.py"""
    print("Printing recalculated remlog file...")
    pairs_columns = [remlog_np[:,i:i+2] for i in range(remlog_np.shape[1]) if i%2==0]
    replica_column_1 = np.arange(1,remlog_np.shape[0]+1).T
    # Read temperature and number of exchanges from mdin.ref
    with open(f'./FE/RE/{wt_or_mt}/mdin.ref', 'r') as mdin:
        lines = mdin.readlines()
        for line in lines:
            elements = line.split(',')
            for element in elements:
                if element.strip()[0:5] == 'temp0':
                    temperature = float(element.split('=')[1].strip())
                if element.strip()[0:8] == 'numexchg':
                    n_exchanges = float(element.split('=')[1].strip())
    shutil.copyfile('setup/tmpls/free_energy_tmpls/remlog_recalculated.tmpl', f'{snapshot_dir}/remlog_recalculated.out')
    replace_dict = {
        '%n_exchanges%': str(n_exchanges)
    }
    get_data_n_general.replace_in_file(f'{snapshot_dir}/remlog_recalculated.out', replace_dict)
    # Create pandas dataframes to print
    for i in range(len(pairs_columns)):
        exchange_separator = f'# exchange        {i}'
        if i%2 == 0:
            replica_column_2 = np.array([replica_column_1[i] + 1 if i%2==0 else replica_column_1[i] - 1 for i in range(len(replica_column_1))]).T
        else:
            replica_column_2 = np.array([36] + 
                [replica_column_1[i] + 1 if i%2==1 else replica_column_1[i] - 1 for i in range(1,len(replica_column_1[1:]))] + [1]).T
        temperature_column = np.full_like(replica_column_1, temperature)
        left_fe = np.zeros_like(replica_column_1)
        right_fe = np.zeros_like(replica_column_1)
        success_column = np.full_like(replica_column_1, False)
        rate_column = np.zeros_like(replica_column_1)
        exchange_np = np.stack(replica_column_1, replica_column_2, pairs_columns, temperature_column, 
            left_fe, right_fe, success_column, rate_column, axis=1)
        exchange_df = pd.DataFrame(data=exchange_np, 
            columns=['Replica', 'Neighbor', 'Temperature', 'PotE(x_1)', 'PotE(x_2)', 'left_fe', 'right_fe', 'Success', 'Success_rate'])
        # Append to file
        with open(f'{snapshot_dir}/remlog_recalculated.out', 'a') as remlog:
            remlog.write(exchange_separator)
        exchange_df.to_csv(f'{snapshot_dir}/remlog_recalculated.out', mode='a', index=False, header=False, sep="\t", quoting=csv.QUOTE_NONE)


if __name__ == '__main__':
    # List input files and determine if everything needed is there
    files = [x for x in os.listdir() if os.path.isfile(x)] # List input files
    pdb_files = [x for x in files if x[-4:]==".pdb"]
    control_file = 'control.txt'
    if len(pdb_files) == 0:
        raise FileNotFoundError('No PDB file found.')
    if len(pdb_files) > 1:
        raise Exception('More than one PDB file found. Only one is acceptable.')
    if not control_file in files:
        raise FileNotFoundError('No control file \'control.txt\' found.')
    pdb_file = pdb_files[0]

    control_dict = get_data_n_general.read_input(control_file) # Get required parameters from control file
    control_dict['Rgb_modifiers'] = read_gb_modifiers()

    mutation = MutationReGbFe.MutationReGbFe(control_dict, pdb_file) # Get information about structure
    create_parm.modify_og_GBRadius(mutation.gb_modifiers, mutation.include_mut) # Modify original GB radius
    print("Modifying intermediate parameter files...")
    # Update GB radii in intermediate parameter files
    create_parm.create_intermediate_parms(mutation.functions, mutation.windows, mutation.leap_residue_position, 
        mutation.intermediate, mutation.include_mut, recalculation=True)
    print("Intermediate parameters finished.")

    # Recalculate using updated parameter files and already simulated REMD trajectories
    re_dir = './FE/RE/'
    wt_dir = 'WT/'
    mt_dir = 'MT/'
    path_to_parms = './setup/recalculate/'

    if mutation.include_mut == True:
        loop_dirs = [wt_dir, mt_dir]
    else:
        loop_dirs = [wt_dir]
    
    for wt_or_mt in loop_dirs:
        dir_to_file = {
            'WT/': 'wt',
            'MT/': 'mt'
        }
        print(f"Recalculating {dir_to_file[wt_or_mt]}...")
        # Get pdb structure
        if wt_or_mt == 'WT/':
            pdb_structure = mutation.wt_structure
        elif wt_or_mt == 'MT/':
            pdb_structure = mutation.mt_structure
        # Get directories of REMD simulations
        remd_dir = re_dir + wt_or_mt
        snapshot_dirs = [x for x in os.listdir(remd_dir) if os.path.isdir(remd_dir + x)]
        snapshot_dirs.sort(key=int)
        # Get parameter file paths
        parm_files = [x for x in os.listdir(path_to_parms) if os.path.isfile(path_to_parms + x) and x[:2] == dir_to_file[wt_or_mt]]
        parm_files.sort(key=lambda x: int(x[x.find(f'{dir_to_file[wt_or_mt]}_')+len(f'{dir_to_file[wt_or_mt]}_'):x.rfind('.parm7')]))
        parm_paths = [path_to_parms + x for x in parm_files]

        for snapshot in snapshot_dirs:
            print(f"Recalculating snapshot: {snapshot}")
            # Get trajectory file paths
            snapshot_dir = f'{remd_dir + snapshot}/'
            traj_files = [x for x in os.listdir(snapshot_dir) if os.path.isfile(snapshot_dir + x) and x[-3:] == '.nc']
            traj_files.sort(key=lambda x: int(x[-6:-3]))
            traj_paths = [snapshot_dir + x for x in traj_files]
            # Read reference structure for restraints
            ref_structure = pt.iterload(f'{snapshot_dir}/restart.rst7', parm_paths[0])

            # Recalculate E_1(x_1) for all exchanges
            # Import trajectories with corresponding parameters
            trajs_straight = [pt.iterload(traj_paths[i], parm_paths[i]) for i in range(len(traj_paths))]
            print("Recalculating E_1(x_1)...")
            e1_x1 = [pt.esander(x, igb=5) for x in trajs_straight] # Calculate energies of trajectory at corresponding Hamiltonian
            n_residues = pdb_structure.df['ATOM']['residue_number'].nunique() # Number of residues
            rmsd_ref = [pt.rmsd(x, '@CA', ref=ref_structure, nofit=True) for x in trajs_straight]
            force_constant = get_force_constant(wt_or_mt[:-1])
            e1_x1_restraints = [force_constant*n_residues*x**2 for x in rmsd_ref] # Calculate restraint energy
            # Calculate total energy for E_1(x_1) column
            e1_x1_total = [e1_x1_restraints[i] + e1_x1[i]['vdw'] + e1_x1[i]['elec'] + e1_x1[i]['gb'] + e1_x1[i]['bond'] + 
                e1_x1[i]['angle'] + e1_x1[i]['dihedral'] + e1_x1[i]['vdw_14'] + e1_x1[i]['elec_14'] for i in range(len(e1_x1))]
            e1_x1_total = [x.round(decimals=2) for x in e1_x1_total]
            e1_x1_np = np.stack(e1_x1_total, axis=1).T
            print("Recalculation finished...")

            # Recalculate E_1(x_2) for even exchanges
            # Import trajectories with parameters shifted to the right
            trajs_even_right = [pt.iterload(traj_paths[i], parm_paths[i+1]) for i in range(len(traj_paths)-1)]
            trajs_even_right.append(pt.iterload(traj_paths[-1], parm_paths[0]))
            print("Recalculating E_1(x_2) for even changes (1/2)...")
            e1_x2_even1 = [pt.esander(x[::2], igb=5) for x in trajs_even_right]
            e1_x2_even1_total = [e1_x1_restraints[i][::2] + e1_x2_even1[i]['vdw'] + e1_x2_even1[i]['elec'] + e1_x2_even1[i]['gb'] + 
                e1_x2_even1[i]['bond'] + e1_x2_even1[i]['angle'] + e1_x2_even1[i]['dihedral'] + e1_x2_even1[i]['vdw_14'] + e1_x2_even1[i]['elec_14'] 
                for i in range(len(e1_x2_even1)) if i%2 == 0]
            e1_x2_even1_total = [x.round(decimals=2) for x in e1_x2_even1_total]
            e1_x2_even1_np = np.stack(e1_x2_even1_total, axis=1).T
            print("Recalculation finished...")
            # Import trajectories with parameters shifted to the left
            trajs_even_left = [pt.iterload(traj_paths[i+1], parm_paths[i]) for i in range(len(traj_paths)-1)]
            trajs_even_left = [pt.iterload(traj_paths[0], parm_paths[-1])] + trajs_even_left
            print("Recalculating E_1(x_2) for even changes (2/2)...")
            e1_x2_even2 = [pt.esander(x[::2], igb=5) for x in trajs_even_left]
            e1_x2_even2_total = [e1_x1_restraints[i][::2] + e1_x2_even2[i]['vdw'] + e1_x2_even2[i]['elec'] + e1_x2_even2[i]['gb'] + 
                e1_x2_even2[i]['bond'] + e1_x2_even2[i]['angle'] + e1_x2_even2[i]['dihedral'] + e1_x2_even2[i]['vdw_14'] + e1_x2_even2[i]['elec_14'] 
                for i in range(len(e1_x2_even2)) if i%2 != 0]
            e1_x2_even2_total = [x.round(decimals=2) for x in e1_x2_even2_total]
            e1_x2_even2_np = np.stack(e1_x2_even2_total, axis=1).T
            print("Recalculation finished...")
            # Intercalate array rows to form columns of even exchanges
            e1_x2_even_np = np.empty((e1_x2_even1_np.shape[0] + e1_x2_even2_np.shape[0], e1_x2_even1_np.shape[1]), dtype=float)
            e1_x2_even_np[1::2,:] = e1_x2_even1_np
            e1_x2_even_np[::2,:] = e1_x2_even2_np

            # Recalculate E_1(x_2) for odd exchanges
            # Import trajectories with parameters shifted to the right
            trajs_odd_right = [pt.iterload(traj_paths[i], parm_paths[i-1]) for i in range(1,len(traj_paths))]
            trajs_odd_right.append(pt.iterload(traj_paths[0], parm_paths[-1]))
            print("Recalculating E_1(x_2) for odd changes (1/2)...")
            e1_x2_odd1 = [pt.esander(x[1::2], igb=5) for x in trajs_odd_right]
            e1_x2_odd1_total = [e1_x1_restraints[i+1][1::2] + e1_x2_odd1[i]['vdw'] + e1_x2_odd1[i]['elec'] + e1_x2_odd1[i]['gb'] + 
                e1_x2_odd1[i]['bond'] + e1_x2_odd1[i]['angle'] + e1_x2_odd1[i]['dihedral'] + e1_x2_odd1[i]['vdw_14'] + e1_x2_odd1[i]['elec_14'] 
                for i in range(len(e1_x2_odd1)-1) if i%2 != 0] + [e1_x1_restraints[0][1::2] + e1_x2_odd1[-1]['vdw'] + e1_x2_odd1[-1]['elec'] + 
                e1_x2_odd1[-1]['gb'] + e1_x2_odd1[-1]['bond'] + e1_x2_odd1[-1]['angle'] + e1_x2_odd1[-1]['dihedral'] + e1_x2_odd1[-1]['vdw_14'] + 
                e1_x2_odd1[-1]['elec_14']]
            e1_x2_odd1_total = [x.round(decimals=2) for x in e1_x2_odd1_total]
            e1_x2_odd1_np = np.stack(e1_x2_odd1_total, axis=1).T
            print("Recalculation finished...")
            # Import trajectories with parameters shifted to the left
            trajs_odd_left = [pt.iterload(traj_paths[i-1], parm_paths[i]) for i in range(len(traj_paths))]
            print("Recalculating E_1(x_2) for odd changes (2/2)...")
            e1_x2_odd2 = [pt.esander(x[1::2], igb=5) for x in trajs_odd_left]
            e1_x2_odd2_total = [e1_x1_restraints[i-1][1::2] + e1_x2_odd2[i]['vdw'] + e1_x2_odd2[i]['elec'] + e1_x2_odd2[i]['gb'] + 
                e1_x2_odd2[i]['bond'] + e1_x2_odd2[i]['angle'] + e1_x2_odd2[i]['dihedral'] + e1_x2_odd2[i]['vdw_14'] + e1_x2_odd2[i]['elec_14'] 
                for i in range(len(e1_x2_odd2)) if i%2 == 0]
            e1_x2_odd2_total = [x.round(decimals=2) for x in e1_x2_odd2_total]
            e1_x2_odd2_np = np.stack(e1_x2_odd2_total, axis=1).T
            print("Recalculation finished...")
            # Intercalate array rows to form columns of odd exchanges
            e1_x2_odd_np = np.empty((e1_x2_odd1_np.shape[0] + e1_x2_odd2_np.shape[0], e1_x2_odd1_np.shape[1]), dtype=float)
            e1_x2_odd_np[1::2,:] = e1_x2_odd1_np
            e1_x2_odd_np[::2,:] = e1_x2_odd2_np

            # Intercalate array columns to form E_1(x_2)
            e1_x2_np = np.empty((e1_x2_odd_np.shape[0], e1_x2_odd_np.shape[1] + e1_x2_even_np.shape[1]), dtype=float)
            e1_x2_np[:,::2] = e1_x2_even_np
            e1_x2_np[:,1::2] = e1_x2_odd_np
            # Intercalate array columns to form whole remlog file
            remlog_np = np.empty((e1_x1_np.shape[0], e1_x1_np.shape[1] + e1_x2_np.shape[1]), dtype=float)
            remlog_np[:,::2] = e1_x1_np
            remlog_np[:,1::2] = e1_x2_np

            print_remlog(remlog_np, snapshot_dir)
    
    print("Done.")