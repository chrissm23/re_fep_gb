import numpy as np
import pandas as pd
import os

# Obtain directories with both termodynamic pathways
only_dirs = [d for d in os.listdir() if os.path.isdir(os.path.join('./', d))]
#only_dirs.sort()
headers = ['DeltaG', 'DeltaE_surf', 'DeltaG_total']
results_paths = []

for current_dir in only_dirs:
    print(current_dir)
    # Obtain directories of mutations
    mutants_dir = os.path.join('./', current_dir, 'Mutations')
    mutations = [m for m in os.listdir(mutants_dir) if os.path.isdir(os.path.join(mutants_dir, m))]
    mutations.sort()
    results = []

    for mutant in mutations:
        #print(mutant)
        # Read results files and store as dataframe
        fe_file = os.path.join(mutants_dir, mutant, 'FE/FE_diff.out')
        with open(fe_file, 'r') as f:
            lines = [line.strip() for line in f.readlines()]
        results.append([float(line[line.find(':')+1:].strip()) for line in lines[-3:]])
    results_paths.append(pd.DataFrame(np.array(results), columns=headers))
    pd_mutations = pd.DataFrame({'Mutations': mutations})

# Calculate DeltaDeltaG
DeltaDeltaG = results_paths[0]['DeltaG_total'] - results_paths[1]['DeltaG_total']
print(pd_mutations.join(DeltaDeltaG))