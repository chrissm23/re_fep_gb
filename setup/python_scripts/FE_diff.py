import numpy as np
from scipy.constants import k as k_B
from scipy.constants import N_A
import pandas as pd
import os

import pymbar

def correct_splitting_errors(line):
    """Correct line splitting errors when two  columns merge"""
    for i in range(len(line)):
        if line[i].count('.') > 1 or (line[i].count('.') == 1 and line[i][-10:] == '**********'):
            new_element = line[i][line[i].find('.') + 3:]
            line.insert(i+1, new_element)
            line[i] = line[i][:line[i].find('.') + 3]
        if line[i] == '********************' or (line[i].count('.') == 1 and line[i][:10] == '**********'):
            new_element = line[i][10:]
            line.insert(i+1, new_element)
            line[i] = line[i][:10]
    return line

def read_remlog(path):
    """Read remlog file in path and extract energy"""
    print("Parsing remlog file...")
    lines = []
    with open(path, 'r') as remlog:
        lines = remlog.readlines()
    
    numexchg = int(lines[1][13:].strip()) # Number of exchanges
    lines_exchg = [[] for _ in range(numexchg)] # List of lines for each exchange

    # Split lines into exchanges
    counter_exchg = -1
    exchg_flag = False
    for line in lines:
        if line[:10] == "# exchange":
            exchg_flag = False
        if exchg_flag == True:
            new_line = line.strip().split()
            if len(new_line) < 9:
                new_line = correct_splitting_errors(new_line)
            lines_exchg[counter_exchg].append(new_line)
        if line[:10] == "# exchange":
            exchg_flag = True
            counter_exchg += 1

    # Convert into pandas dataframe
    headers = ['Replica', 'Neighbor', 'Temperature', 'PotE(x_1)', 'PotE(x_2)', 'left_fe', 'right_fe', 'Success', 'Success_rate']
    replicas_pd = [pd.DataFrame(np.array(line_exchg), columns=headers) for line_exchg in lines_exchg if line_exchg]
    n_exchgs = len(replicas_pd[0].index)
    # Discard incomplete exchanges
    if len(replicas_pd[-1].index) != n_exchgs:
        if len(replicas_pd)%2 == 0:
            replicas_pd = replicas_pd[:-2]
        else:
            replicas_pd = replicas_pd[:-1]
    # Transform entries into numeric
    for replica in replicas_pd:
        replica['Replica'] = pd.to_numeric(replica['Replica'])
        replica['Neighbor'] = pd.to_numeric(replica['Neighbor'])
        replica['Temperature'] = pd.to_numeric(replica['Temperature'])
        replica['PotE(x_1)'] = pd.to_numeric(replica['PotE(x_1)'], errors='coerce')
        replica['PotE(x_2)'] = pd.to_numeric(replica['PotE(x_2)'], errors='coerce')
        replica['left_fe'] = pd.to_numeric(replica['left_fe'], errors='coerce')
        replica['right_fe'] = pd.to_numeric(replica['right_fe'], errors='coerce')
        replica['Success_rate'] = pd.to_numeric(replica['Success_rate'])
    
    return replicas_pd

def DeltaG_FEP(replicas_pd):
    """Calculate DeltaG using free energy perturbation"""
    print("Calculating free energy differences using FEP...")
    n_replicas = len(replicas_pd[0].index)
    exp_sum_backward_1 = np.zeros(int(n_replicas/2)) # Array for sum of exponentials for backward free energy difference
    exp_sum_backward_2 = np.zeros(int(n_replicas/2)-1) # Array for sum of exponentials for backward free energy difference
    exp_sum_forward_1 = np.zeros(int(n_replicas/2)) # Array for sum of exponentials for forward free energy difference
    exp_sum_forward_2 = np.zeros(int(n_replicas/2)-1) # Array for sum of exponentials for forward free energy difference
    Temp = replicas_pd[0]['Temperature'].iloc[0] # Temperature
    beta = 4184/(k_B*Temp*N_A)
    n_exchanges = len(replicas_pd)
    nan_counter_f1 = 0
    nan_counter_f2 = 0
    nan_counter_b1 = 0
    nan_counter_b2 = 0
    # Array of Boltzmann factors for distributions
    exp_deltaE_forward_1 = np.zeros((n_exchanges, int(n_replicas/2)))
    exp_deltaE_forward_2 = np.zeros((n_exchanges, int(n_replicas/2)-1))
    exp_deltaE_backward_1 = np.zeros((n_exchanges, int(n_replicas/2)))
    exp_deltaE_backward_2 = np.zeros((n_exchanges, int(n_replicas/2)-1))

    for j in range(n_exchanges):
        if j%2 == 0: # Exchanges starting with 1 to last replica
            # Half of forward exponential sums
            E_ip1_forward_2 = replicas_pd[j][replicas_pd[j].index % 2 == 0]['PotE(x_2)'].to_numpy()
            E_i_forward_2 = replicas_pd[j][replicas_pd[j].index % 2 == 1]['PotE(x_1)'].to_numpy()
            if not (np.isnan(np.sum(E_ip1_forward_2)) or np.isnan(np.sum(E_i_forward_2))):
                DeltaE_forward_2 = -(E_ip1_forward_2[1:] - E_i_forward_2[:-1])*beta
                exp_forward_2 = np.exp(DeltaE_forward_2)
                exp_deltaE_forward_2[j,:] = exp_forward_2
                exp_sum_forward_2 += exp_forward_2
            else:
                nan_counter_f2 += 1

            # Half of backward exponential sums
            E_im1_backward_2 = replicas_pd[j][replicas_pd[j].index % 2 == 1]['PotE(x_2)'].to_numpy()
            E_i_backward_2 = replicas_pd[j][replicas_pd[j].index % 2 == 0]['PotE(x_1)'].to_numpy()
            if not (np.isnan(np.sum(E_im1_backward_2)) or np.isnan(np.sum(E_i_backward_2))):
                DeltaE_backward_2 = -(E_im1_backward_2[:-1] - E_i_backward_2[1:])*beta
                exp_backward_2 = np.exp(DeltaE_backward_2)
                exp_deltaE_backward_2[j,:] = exp_backward_2
                exp_sum_backward_2 += exp_backward_2
            else:
                nan_counter_b2 += 1

        elif j%2 == 1: # Exchanges starting with 1 to 2
            # Half of forward exponential sums
            E_ip1_forward_1 = replicas_pd[j][replicas_pd[j].index % 2 == 1]['PotE(x_2)'].to_numpy()
            E_i_forward_1 = replicas_pd[j][replicas_pd[j].index % 2 == 0]['PotE(x_1)'].to_numpy()
            if not (np.isnan(np.sum(E_ip1_forward_1)) or np.isnan(np.sum(E_i_forward_1))):
                DeltaE_forward_1 = -(E_ip1_forward_1 - E_i_forward_1)*beta
                exp_forward_1 = np.exp(DeltaE_forward_1)
                exp_deltaE_forward_1[j,:] = exp_forward_1
                exp_sum_forward_1 += np.exp(DeltaE_forward_1)
            else:
                nan_counter_f1 += 1

            # Half of backward exponential sums
            E_im1_backward_1 = replicas_pd[j][replicas_pd[j].index % 2 == 0]['PotE(x_2)'].to_numpy()
            E_i_backward_1 = replicas_pd[j][replicas_pd[j].index % 2 == 1]['PotE(x_1)'].to_numpy()
            if not (np.isnan(np.sum(E_im1_backward_1)) or np.isnan(np.sum(E_i_backward_1))):
                DeltaE_backward_1 = -(E_im1_backward_1 - E_i_backward_1)*beta
                exp_backward_1 = np.exp(DeltaE_backward_1)
                exp_deltaE_backward_1[j,:] = exp_backward_1
                exp_sum_backward_1 += exp_backward_1
            else:
                nan_counter_b1 += 1

    # Exponential averages
    n_data = n_exchanges/2
    av_exp_backward_1 = exp_sum_backward_1/(n_data - nan_counter_b1)
    av_exp_backward_2 = exp_sum_backward_2/(n_data - nan_counter_b2)
    av_exp_forward_1 = exp_sum_forward_1/(n_data - nan_counter_f1)
    av_exp_forward_2 = exp_sum_forward_2/(n_data - nan_counter_f2)
    # Free energy differences
    DeltaGs_backward_1 = -np.log(av_exp_backward_1)/beta
    DeltaGs_backward_2 = -np.log(av_exp_backward_2)/beta
    DeltaGs_forward_1 = -np.log(av_exp_forward_1)/beta
    DeltaGs_forward_2 = -np.log(av_exp_forward_2)/beta
    # Total free energy difference
    DeltaG_backward = np.sum(np.concatenate([DeltaGs_backward_1, DeltaGs_backward_2]))
    DeltaG_forward = np.sum(np.concatenate([DeltaGs_forward_1, DeltaGs_forward_2]))

    # Delete remaning zeros and intercalate arrays with Boltzmann factors
    # Forward
    exp_deltaE_forward_1 = exp_deltaE_forward_1[~np.all(exp_deltaE_forward_1 == 0, axis=1)]
    exp_deltaE_forward_2 = exp_deltaE_forward_2[~np.all(exp_deltaE_forward_2 == 0, axis=1)]
    row1, col1 = np.shape(exp_deltaE_forward_1)
    row2, col2 = np.shape(exp_deltaE_forward_2)
    exp_deltaE_fwd = np.zeros((max(row1,row2), col1+col2))
    exp_deltaE_fwd[:row1,::2] = exp_deltaE_forward_1
    exp_deltaE_fwd[:row2,1::2] = exp_deltaE_forward_2
    # Backward
    exp_deltaE_backward_1 = exp_deltaE_backward_1[~np.all(exp_deltaE_backward_1 == 0, axis=1)]
    exp_deltaE_backward_2 = exp_deltaE_backward_2[~np.all(exp_deltaE_backward_2 == 0, axis=1)]
    row1, col1 = np.shape(exp_deltaE_backward_1)
    row2, col2 = np.shape(exp_deltaE_backward_2)
    exp_deltaE_bwd = np.zeros((max(row1,row2), col1+col2))
    exp_deltaE_bwd[:row1,::2] = exp_deltaE_backward_1
    exp_deltaE_bwd[:row2,1::2] = exp_deltaE_backward_2

    return [DeltaG_backward, DeltaG_forward, exp_deltaE_fwd, exp_deltaE_bwd]

def DeltaG_BAR(replicas_pd):
    """Calculate DeltaG using BAR"""
    print("Computing free energy differences using BAR...")
    n_replicas = len(replicas_pd[0].index)
    n_exchg = len(replicas_pd)
    Temp = replicas_pd[0]['Temperature'].iloc[0] # Temperature
    beta = 4184/(k_B*Temp*N_A)

    # Rearrange potential energy into numpy arrays
    u_11_a = np.stack([replicas_pd[j][replicas_pd[j].index % 2 == 1]['PotE(x_1)'].iloc[:-1].to_numpy() for j in range(n_exchg) if j%2 == 0],
        axis=1) # First half of u_1(x_1)
    u_11_b = np.stack([replicas_pd[j][replicas_pd[j].index % 2 == 0]['PotE(x_1)'].to_numpy() for j in range(n_exchg) if j%2 == 1], 
        axis=1) # Second half of u_1(x_1)
    u_21_a = np.stack([replicas_pd[j][replicas_pd[j].index % 2 == 0]['PotE(x_2)'].iloc[1:].to_numpy() for j in range(n_exchg) if j%2 == 0], 
        axis=1) # First half of u_2(x_1)
    u_21_b = np.stack([replicas_pd[j][replicas_pd[j].index % 2 == 1]['PotE(x_2)'].to_numpy() for j in range(n_exchg) if j%2 == 1], 
        axis=1) # Second half of u_2(x_1)
    u_12_a = np.stack([replicas_pd[j][replicas_pd[j].index % 2 == 1]['PotE(x_2)'].iloc[:-1].to_numpy() for j in range(n_exchg) if j%2 == 0], 
        axis=1) # First half of u_1(x_2)
    u_12_b = np.stack([replicas_pd[j][replicas_pd[j].index % 2 == 0]['PotE(x_2)'].to_numpy() for j in range(n_exchg) if j%2 == 1],
        axis=1) # Second half of u_1(x_2)
    u_22_a = np.stack([replicas_pd[j][replicas_pd[j].index % 2 == 0]['PotE(x_1)'].iloc[1:].to_numpy() for j in range(n_exchg) if j%2 == 0],
        axis=1) # First half of u_2(x_2)
    u_22_b = np.stack([replicas_pd[j][replicas_pd[j].index % 2 == 1]['PotE(x_1)'].to_numpy() for j in range(n_exchg) if j%2 == 1],
        axis=1) # Second half of u_2(x_2)
    # Concatenate arrays u_i(x_i) and u_i(x_j) or u_j(x_i) and u_j(x_j)
    u_1_a = np.concatenate([u_11_a, u_12_a], axis=1)
    u_1_b = np.concatenate([u_11_b, u_12_b], axis=1)
    u_2_a = np.concatenate([u_21_a, u_22_a], axis=1)
    u_2_b = np.concatenate([u_21_b, u_22_b], axis=1)

    # Concatenate pairs of rows for BAR
    u_ns_a = [np.stack([u_1_a[i,:], u_2_a[i,:]], axis=0)*beta for i in range(int(n_replicas/2 - 1))]
    u_ns_b = [np.stack([u_1_b[i,:], u_2_b[i,:]], axis=0)*beta for i in range(int(n_replicas/2))]
    # Filter out NaN
    n_exchgs_a = [np.array([n_exchg/2, n_exchg/2]) for _ in range(int(n_replicas/2 - 1))]
    n_exchgs_b = [np.array([n_exchg/2, n_exchg/2]) for _ in range(int(n_replicas/2))]
    for u_ns,n_exchanges in zip([u_ns_a, u_ns_b],[n_exchgs_a, n_exchgs_b]):
        nans_idx = [np.argwhere(np.isnan(x)) for x in u_ns]
        nans_columns = [list(set(x[:,1].tolist())) for x in nans_idx]
        count_nans = [np.array([sum(map(lambda x : x < n_exchg/2, y)),sum(map(lambda x : x >= n_exchg/2, y))]) for y in nans_columns]
        new_exchgs = [x-y for x,y in zip(n_exchanges, count_nans)]
        if u_ns == u_ns_a:
            n_exchgs_a = new_exchgs
            u_ns_a = [np.delete(u_ns_a[i], nans_columns[i], axis=1) for i in range(len(u_ns_a))]
        elif u_ns == u_ns_b:
            n_exchgs_b = new_exchgs
            u_ns_b = [np.delete(u_ns_b[i], nans_columns[i], axis=1) for i in range(len(u_ns_b))]
    # Calculate DeltaG with BAR
    mbars_a = [pymbar.MBAR(u_ns_a[i], n_exchgs_a[i], maximum_iterations=100000).getFreeEnergyDifferences(return_dict=True) for i in range(len(u_ns_a))]
    mbars_b = [pymbar.MBAR(u_ns_b[i], n_exchgs_b[i], maximum_iterations=100000).getFreeEnergyDifferences(return_dict=True) for i in range(len(u_ns_b))]
    DeltaGs_a = [x['Delta_f'][0,1] for x in mbars_a]
    DeltaGs_b = [x['Delta_f'][0,1] for x in mbars_b]
    DeltaG = sum(DeltaGs_a) + sum(DeltaGs_b)

    return DeltaG/beta

def get_SASA(path):
    """Read remlog file in path and extract energy"""
    print("Parsing output file of GBSA calculation and averaging SASA contribution...")
    lines = []
    with open(path, 'r') as sasa_output:
        lines = sasa_output.readlines()
    
    # Filter lines with ESURF only
    esurf_lines = []
    for line in lines:
        if 'A V E R A G E S   O V E R' in line:
            break
        if 'ESURF=' in line:
            esurf_lines.append(float(line.split('=')[1].strip()))
    E_surf = np.average(esurf_lines)

    return E_surf

def get_averageE(replicas_pd, replica1, replica2):
    """Calculate average potential energy between replica1 and replica2"""
    print(f'Computing potential energy difference between replica {replica1} and replica {replica2}')
    epot_sum1 = 0
    epot_sum2 = 0
    
    # Add all PotE(x_1)
    for replica in replicas_pd:
        replica1_data = replica[replica['Replica'] == replica1]
        epot_sum1 += replica1_data['PotE(x_1)'].iloc[0]
        replica2_data = replica[replica['Replica'] == replica2]
        epot_sum2 += replica2_data['PotE(x_1)'].iloc[0]
    # Calculate average
    n_exchgs = len(replicas_pd)
    epot1_av = epot_sum1/n_exchgs
    epot2_av = epot_sum2/n_exchgs

    return epot2_av - epot1_av

if __name__ == '__main__':
    FE_dir = './'
    snapshot_windows = [name for name in os.listdir('./RE/WT') if os.path.isdir(os.path.join('./RE/WT', name))]
    snapshot_windows.sort()

    fep_avs_f = []
    fep_avs_b = []
    bar_avs = []
    sasas = []
    
    # Get DeltaG for WT and MT
    if os.path.exists('./RE/MT/0/rem.log'):
        wt_mt_loop = ['WT', 'MT']
    else:
        wt_mt_loop = ['WT']
    for wt_or_mt in wt_mt_loop:
        fep_energies = []
        bar_energies = []

        counter_fep_errors = 0
        counter_bar_errors = 0

        for snapshot in snapshot_windows:
            # Parse remlog file
            print(f'Calculating DeltaG of {wt_or_mt}, snapshot {snapshot}, from H-REMD...\n')
            replicas_pd = read_remlog(FE_dir + f'/RE/{wt_or_mt}/{snapshot}/rem.log')
            if snapshot == snapshot_windows[0]:
                Temp = replicas_pd[0]['Temperature'].iloc[0] # Temperature
                beta = 4184/(k_B*Temp*N_A)
                n_replicas = len(replicas_pd[0].index)
                exp_deltaEs_forward = [[] for _ in range(n_replicas-1)]
                exp_deltaEs_backward = [[] for _ in range(n_replicas-1)]

            # Calculate DeltaG using FEP
            try:
                [DeltaG_backward, DeltaG_forward, exponentials_fwd, exponentials_bwd] = DeltaG_FEP(replicas_pd)
            except:
                counter_fep_errors += 1
            else:
                fep_energies.append([DeltaG_forward, DeltaG_backward])
                print(f'Forward DeltaG = {round(DeltaG_forward, 2)}, Backward DeltaG = {round(DeltaG_backward, 2)}\n')
                for i in range(len(exp_deltaEs_forward)):
                    exp_deltaEs_forward[i].extend(exponentials_fwd[exponentials_fwd[:,i] != 0, i].tolist())
                    exp_deltaEs_backward[i].extend(exponentials_bwd[exponentials_bwd[:,i] != 0, i].tolist())
            
            # Calculate DeltaG using BAR
            try:
                DeltaG = DeltaG_BAR(replicas_pd)
            except:
                counter_bar_errors += 0
            else:
                bar_energies.append(DeltaG)
                print(f'Forward DeltaG = {round(DeltaG, 2)}\n')

        # Get FEP deltaG over all samples
        fep_forward_allavg = [np.average(exp_deltaEs_forward[i]) for i in range(len(exp_deltaEs_forward))]
        fep_forward_allstd = [np.std(exp_deltaEs_forward[i]) for i in range(len(exp_deltaEs_forward))]
        filtered_exp_deltaEs_forward = [list(filter(lambda b_factor: 
            (b_factor <= fep_forward_allavg[i]+5*fep_forward_allstd[i]) and (b_factor >= fep_forward_allavg[i]-5*fep_forward_allstd[i]), 
            exp_deltaEs_forward[i])) 
            for i in range(len(exp_deltaEs_forward))]
        fep_forward = sum(-np.log(fep_forward_allavg)/beta)
        fep_backward_allavg = [np.average(exp_deltaEs_backward[i]) for i in range(len(exp_deltaEs_backward))]
        fep_backward = sum(-np.log(fep_backward_allavg)/beta)
        fep_backward_allstd = [np.std(exp_deltaEs_backward[i]) for i in range(len(exp_deltaEs_backward))]

        # Get averages of DeltaG and print warnings
        bar_avs.append(np.average(bar_energies))
        print(f'BAR {wt_or_mt}: Average DeltaG = {round(bar_avs[-1], 2)}\n')

        fep_avs_f.append(np.average([fep_energies[i][0] for i in range(len(fep_energies))]))
        fep_avs_b.append(np.average([fep_energies[i][1] for i in range(len(fep_energies))]))
        print(f'FEP {wt_or_mt}: Average Forward DeltaG = {round(fep_avs_f[-1], 2)}, Average Backward DeltaG = {round(fep_avs_b[-1], 2)}')
        print(f'FEP {wt_or_mt}: Forward DeltaG = {round(fep_forward, 2)}, Backward DeltaG = {round(fep_backward, 2)}\n')

        if abs(fep_avs_f[-1]) < abs(fep_avs_b[-1]) - 1 or abs(fep_avs_f[-1]) > abs(fep_avs_b[-1]) + 1:
            print(f'Convergence error\n')

        if counter_fep_errors > 0 or counter_bar_errors > 0:
            print(f'Warning: {counter_fep_errors}/{len(snapshot_windows)} explosions\n')

    for wt_or_mt in ['WT', 'MT']:
        E_surf = get_SASA(FE_dir + f'/SASA/{wt_or_mt}/sasa.out')
        sasas.append(E_surf)
        print(f'{wt_or_mt} E_surf = {round(E_surf, 2)}\n')

    if len(bar_avs) > 1:
        G_diff_fep_f = fep_avs_f[0] - fep_avs_f[1]
        G_diff_fep_b = fep_avs_b[0] - fep_avs_b[1]
        G_diff_bar = bar_avs[0] - bar_avs[1]
    else:
        G_diff_fep_f = fep_avs_f[0]
        G_diff_fep_b = fep_avs_b[0]
        G_diff_bar = bar_avs[0]
    sasa_diff = sasas[1] - sasas[0]

    print(f'DeltaG: {round(G_diff_bar, 2)}')
    print(f'DeltaE_surf: {round(sasa_diff, 2)}')
    print(f'DeltaG_total: {round(G_diff_bar + sasa_diff, 2)}')