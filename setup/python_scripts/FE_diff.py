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
    replicas_pd = [pd.DataFrame(np.array(line_exchg), columns=headers) for line_exchg in lines_exchg]
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
    exp_sum_forward = np.zeros(len(replicas_pd[0].index)) # Array for sum of exponentials for forward free energy difference
    exp_sum_backward = np.zeros(len(replicas_pd[0].index)) # Array for sum of exponentials for backward free energy difference
    Temp = replicas_pd[0]['Temperature'].iloc[0] # Temperature
    beta = 4184/(k_B*Temp*N_A)

    for j in range(len(replicas_pd)):
        if j%2 == 0: # Exchanges starting with 1 to last replica
            
            E_ip1_forward_2 = replicas_pd[j][replicas_pd[j].index % 2 == 0]['PotE(x_2)'].to_numpy()
            E_i_forward_2 = replicas_pd[j][replicas_pd[j].index % 2 == 1]['PotE(x_1)'].to_numpy()

        elif j%2 == 1: # Exchanges starting with 1 to 2
            
            E_ip1_forward_1 = replicas_pd[j][replicas_pd[j].index % 2 == 1]['PotE(x_2)'].to_numpy()
            E_i_forward_1 = replicas_pd[j][replicas_pd[j].index % 2 == 0]['PotE(x_1)'].to_numpy()

def DeltaG_BAR(replicas_pd):
    """Calculate DeltaG using BAR"""
    sdfsfd

if __name__ == '__main__':
    replicas_pd = read_remlog('/home/christiansustay/Desktop/rem_wt.log')