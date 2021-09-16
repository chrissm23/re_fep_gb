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
    exp_sum_forward_1 = np.zeros(int(n_replicas/2)) # Array for sum of exponentials for forward free energy difference
    exp_sum_forward_2 = np.zeros(int(n_replicas/2)-1) # Array for sum of exponentials for forward free energy difference
    exp_sum_backward_1 = np.zeros(int(n_replicas/2)) # Array for sum of exponentials for backward free energy difference
    exp_sum_backward_2 = np.zeros(int(n_replicas/2)-1) # Array for sum of exponentials for backward free energy difference
    Temp = replicas_pd[0]['Temperature'].iloc[0] # Temperature
    beta = 4184/(k_B*Temp*N_A)
    n_exchanges = len(replicas_pd)

    for j in range(n_exchanges):
        if j%2 == 0: # Exchanges starting with 1 to last replica
            # Half of forward exponential sums
            E_ip1_forward_2 = replicas_pd[j][replicas_pd[j].index % 2 == 0]['PotE(x_2)'].to_numpy()
            E_i_forward_2 = replicas_pd[j][replicas_pd[j].index % 2 == 1]['PotE(x_1)'].to_numpy()
            DeltaE_forward_2 = -(E_ip1_forward_2[1:] - E_i_forward_2[:-1])*beta
            exp_sum_forward_2 += np.exp(DeltaE_forward_2)
            # Half of backward exponential sums
            E_im1_backward_2 = replicas_pd[j][replicas_pd[j].index % 2 == 1]['PotE(x_2)'].to_numpy()
            E_i_backward_2 = replicas_pd[j][replicas_pd[j].index % 2 == 0]['PotE(x_1)'].to_numpy()
            DeltaE_backward_2 = -(E_im1_backward_2[:-1] - E_i_backward_2[1:])*beta
            exp_sum_backward_2 += np.exp(DeltaE_backward_2)

        elif j%2 == 1: # Exchanges starting with 1 to 2
            # Half of forward exponential sums
            E_ip1_forward_1 = replicas_pd[j][replicas_pd[j].index % 2 == 1]['PotE(x_2)'].to_numpy()
            E_i_forward_1 = replicas_pd[j][replicas_pd[j].index % 2 == 0]['PotE(x_1)'].to_numpy()
            DeltaE_forward_1 = -(E_ip1_forward_1 - E_i_forward_1)*beta
            exp_sum_forward_1 += np.exp(DeltaE_forward_1)
            # Half of backward exponential sums
            E_im1_backward_1 = replicas_pd[j][replicas_pd[j].index % 2 == 0]['PotE(x_2)'].to_numpy()
            E_i_backward_1 = replicas_pd[j][replicas_pd[j].index % 2 == 1]['PotE(x_1)'].to_numpy()
            DeltaE_backward_1 = -(E_im1_backward_1 - E_i_backward_1)*beta
            exp_sum_backward_1 += np.exp(DeltaE_backward_1)

    # Exponential averages
    n_data = n_exchanges/2
    av_exp_forward_1 = exp_sum_forward_1/n_data
    av_exp_forward_2 = exp_sum_forward_2/n_data
    av_exp_backward_1 = exp_sum_backward_1/n_data
    av_exp_backward_2 = exp_sum_backward_2/n_data
    # Free energy differences
    DeltaGs_forward_1 = -np.log(av_exp_forward_1)/beta
    DeltaGs_forward_2 = -np.log(av_exp_forward_2)/beta
    DeltaGs_backward_1 = -np.log(av_exp_backward_1)/beta
    DeltaGs_backward_2 = -np.log(av_exp_backward_2)/beta
    # Total free energy difference
    DeltaG_forward = np.sum(np.concatenate([DeltaGs_forward_1, DeltaGs_forward_2]))
    DeltaG_backward = np.sum(np.concatenate([DeltaGs_backward_1, DeltaGs_backward_2]))

    return DeltaG_forward, DeltaG_backward

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
    # Calculate DeltaG with BAR
    mbars_a = [pymbar.MBAR(x, 2*[n_exchg/2], maximum_iterations=100000).getFreeEnergyDifferences(return_dict=True) for x in u_ns_a]
    mbars_b = [pymbar.MBAR(x, 2*[n_exchg/2], maximum_iterations=100000).getFreeEnergyDifferences(return_dict=True) for x in u_ns_b]
    DeltaGs_a = [x['Delta_f'][1,0] for x in mbars_a]
    DeltaGs_b = [x['Delta_f'][1,0] for x in mbars_b]
    DeltaG = sum(DeltaGs_a) + sum(DeltaGs_b)

    return DeltaG

if __name__ == '__main__':
    replicas_pd = read_remlog('/home/christiansustay/Desktop/rem_mt_vacuum.log')
    [DeltaG_forward, DeltaG_backward] = DeltaG_FEP(replicas_pd)
    print(DeltaG_forward, DeltaG_backward)
    DeltaG = DeltaG_BAR(replicas_pd)
    print(DeltaG)