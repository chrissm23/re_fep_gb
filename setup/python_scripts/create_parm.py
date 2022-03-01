import shutil
import subprocess

import parmed

import setup.python_scripts.get_data_n_general as get_data_n_general

def create_og_parms(path_wt_pdb, path_mt_pdb):
    """Moves, fills out template and executes resulting bash file to create parameter files of the WT and mutant structures given as imput"""
    # Destination paths
    wt_bash_path = 'setup/parms_n_pdbs/parms/gen_parm_wt.sh'
    mt_bash_path = 'setup/parms_n_pdbs/parms/gen_parm_mt.sh'
    tmpl_path = 'setup/tmpls/leap_tmpls/gen_parm_pdb.tmpl'

    # Output paths
    wt_parm_path = 'setup/parms_n_pdbs/parms/parms_windows/wt_0_og.parm7'
    wt_rst_path = 'setup/parms_n_pdbs/parms/rst_windows/wt_0.rst7'
    mt_parm_path = 'setup/parms_n_pdbs/parms/parms_windows/mt_0_og.parm7'
    mt_rst_path = 'setup/parms_n_pdbs/parms/rst_windows/mt_0.rst7'

    # Dictionaries of fields to replace in template
    replace_dict_wt = {
        "%pdb_path%": path_wt_pdb,
        "%name_structure_path_parm%": wt_parm_path,
        "%name_structure_path_rst%": wt_rst_path
    }
    replace_dict_mt = {
        "%pdb_path%": path_mt_pdb,
        "%name_structure_path_parm%": mt_parm_path,
        "%name_structure_path_rst%": mt_rst_path
    }
    
    # Copy templates
    shutil.copyfile(tmpl_path, wt_bash_path)
    shutil.copyfile(tmpl_path, mt_bash_path)
    # Fill out templates
    get_data_n_general.replace_in_file(wt_bash_path, replace_dict_wt)
    get_data_n_general.replace_in_file(mt_bash_path, replace_dict_mt)
    # Make bash files executable and run
    print("Creating WT and mutant parameter files...")
    get_data_n_general.make_executable(wt_bash_path)
    get_data_n_general.make_executable(mt_bash_path)
    subprocess.call(wt_bash_path)
    subprocess.call(mt_bash_path)
    print("Parameter files created.")

def modify_og_GBRadius(modifiers, include_mut):
    """Modifies GB radius of atoms with type atom_types from parmed_object through either a multiplier or a set value"""
    # Leap generated parameter files
    wt_parm_path = 'setup/parms_n_pdbs/parms/parms_windows/wt_0_og.parm7'
    wt_rst_path = 'setup/parms_n_pdbs/parms/rst_windows/wt_0.rst7'
    mt_parm_path = 'setup/parms_n_pdbs/parms/parms_windows/mt_0_og.parm7'
    mt_rst_path = 'setup/parms_n_pdbs/parms/rst_windows/mt_0.rst7'

    # Import to ParmEd
    wt_parmed = parmed.amber.AmberParm(wt_parm_path, wt_rst_path)
    if include_mut:
        mt_parmed = parmed.amber.AmberParm(mt_parm_path, mt_rst_path)

    def Rgb_modify(mask, proportion):
        """Change initial GB radius of mask by proportion"""
        if include_mut:
            topologies = [wt_parmed, mt_parmed]
        else:
            topologies = [wt_parmed]
        for wt_or_mt in topologies:
            if mask == 'all':
                residues_details = get_data_n_general.details_str_to_pd(str(parmed.tools.printDetails(wt_or_mt, ':*')))
            else:
                residues_details = get_data_n_general.details_str_to_pd(str(parmed.tools.printDetails(wt_or_mt, f'{mask}')))
            atom_numbers = residues_details['ATOM'].tolist() # Get atom numbers of mutated residues

            for atom in atom_numbers:
                value = residues_details[residues_details['ATOM'] == atom]['GB Radius'].iloc[0]
                new_value = value + value*proportion
                parmed.tools.change(wt_or_mt, 'RADII', f'@{atom}', f'{new_value}').execute()
    
    print("Modifying original GB Radius...")
    if not isinstance(modifiers, list):
        if '+' in modifiers:
            [mask, proportion] = modifiers.split('+')
            proportion = float(proportion)
            Rgb_modify(mask, proportion)
        elif '-' in modifiers:
            [mask, proportion] = modifiers.split('-')
            proportion = -float(proportion)
            Rgb_modify(mask, proportion)
    else:
        for modifier in modifiers:
            if '+' in modifier:
                [mask, proportion] = modifier.split('+')
                proportion = float(proportion)
                Rgb_modify(mask, proportion)
            elif '-' in modifier:
                [mask, proportion] = modifier.split('-')
                proportion = -float(proportion)
                Rgb_modify(mask, proportion)

    wt_parm_path_new = 'setup/recalculate/wt_0.parm7'
    mt_parm_path_new = 'setup/recalculate/mt_0.parm7'
    parmed.tools.outparm(wt_parmed, wt_parm_path_new).execute()
    if include_mut:
        parmed.tools.outparm(mt_parmed, mt_parm_path_new).execute()

def get_new_LJParms(parmed_object, residue_mask, functions, windows):
    """Creates new LJ atom types and outputs new AmberParm bojects with modified LJ matrix"""
    residue_details = get_data_n_general.details_str_to_pd(str(parmed.tools.printDetails(parmed_object, f':{residue_mask}')))
    atom_numbers = residue_details['ATOM'].tolist() # Get atom numbers of mutated residues

    # Get residue numbers grouped by LJ type
    structure_types_old = get_data_n_general.ljtypes_str_to_pd(str(parmed.tools.printLJTypes(parmed_object)))
    mutated_types_old = structure_types_old[structure_types_old['ATOM'].isin(atom_numbers)]
    atom_types_old = mutated_types_old['LJ Type'].unique().tolist()
    atomtype_old_group = mutated_types_old.groupby('LJ Type')
    masks_by_atomtype_int = []
    for atom_type_old in atom_types_old:
        masks_by_atomtype_int.append(atomtype_old_group.get_group(atom_type_old)['ATOM'].tolist())
    masks_by_atomtype_str = [[str(masks_by_atomtype_int[i][j]) for j in range(len(masks_by_atomtype_int[i]))] for i in range(len(masks_by_atomtype_int))]
    masks_by_atomtype_str = [','.join(x) for x in masks_by_atomtype_str]

    # Create new LJ atom types
    for mask in masks_by_atomtype_str:
        parmed.tools.addLJType(parmed_object, f'@{mask}').execute()
    
    # Group not mutating atoms by LJ atom type
    nonmutating_types = structure_types_old[~structure_types_old['ATOM'].isin(atom_numbers)]
    atom_types_old_nonmut = nonmutating_types['LJ Type'].unique().tolist()
    atomtype_group = nonmutating_types.groupby('LJ Type')
    masks_by_nonmutated_atomtype_int = []
    for atom_type_old in atom_types_old_nonmut:
        masks_by_nonmutated_atomtype_int.append(atomtype_group.get_group(atom_type_old)['ATOM'].tolist())
    masks_by_nonmutated_atomtype_str = [[str(masks_by_nonmutated_atomtype_int[i][j]) for j in range(len(masks_by_nonmutated_atomtype_int[i]))] 
        for i in range(len(masks_by_nonmutated_atomtype_int))]
    masks_by_nonmutated_atomtype_str = [','.join(x) for x in masks_by_nonmutated_atomtype_str]

    # Deep copy parmed object and modify LJ matrix elements
    new_parms = [parmed.amber.AmberParm.from_structure(parmed_object, copy=True) for i in range(len(windows))]
    #print(str(parmed.tools.printLJMatrix(parmed_object, f':{residue_mask}')))
    new_ljmatrix = get_data_n_general.ljmatrix_str_to_pd(str(parmed.tools.printLJMatrix(parmed_object, f':{residue_mask}')))
    for i in range(len(windows)):
        multiplier_R = get_data_n_general.get_multiplier(windows[i], functions[-2], truncate=False)
        multiplier_eps = get_data_n_general.get_multiplier(windows[i], functions[-1], truncate=False)
        for mask_mutant in masks_by_atomtype_str:
            atom_type_mut = get_data_n_general.ljtypes_str_to_pd(str(parmed.tools.printLJTypes(parmed_object, f'@{mask_mutant}')))['LJ Type'][0]
            for mask_nonmutant in masks_by_nonmutated_atomtype_str:
                atom_type_nonmut = get_data_n_general.ljtypes_str_to_pd(str(parmed.tools.printLJTypes(parmed_object, f'@{mask_nonmutant}')))['LJ Type'][0]

                LJRadius_minimum = 0
                LJEps_minimum = 0
                # Get LJ Radius for atom type pair and calculate new LJ Radius
                R_ij = new_ljmatrix[((new_ljmatrix['Atom Type 1'] == atom_type_mut) & (new_ljmatrix['Atom Type 2'] == atom_type_nonmut)) | 
                    ((new_ljmatrix['Atom Type 2'] == atom_type_mut) & (new_ljmatrix['Atom Type 1'] == atom_type_nonmut))]['R_ij'].iloc[0]
                if R_ij <= LJRadius_minimum:
                    new_R_ij = R_ij
                else:
                    new_R_ij = multiplier_R*(R_ij - LJRadius_minimum) + LJRadius_minimum
                # Get LJ epsilon for atom type pair and calculate new LJ epsilon
                Eps_ij = new_ljmatrix[((new_ljmatrix['Atom Type 1'] == atom_type_mut) & (new_ljmatrix['Atom Type 2'] == atom_type_nonmut)) | 
                    ((new_ljmatrix['Atom Type 2'] == atom_type_mut) & (new_ljmatrix['Atom Type 1'] == atom_type_nonmut))]['Eps_ij'].iloc[0]
                if Eps_ij <= LJEps_minimum:
                    new_Eps_ij = Eps_ij
                else:
                    new_Eps_ij = multiplier_eps*(Eps_ij - LJEps_minimum) + LJEps_minimum

                # Change LJ parameters for atom type pair
                parmed.tools.changeLJPair(new_parms[i], f'@{mask_mutant}', f'@{mask_nonmutant}', f'{new_R_ij}', f'{new_Eps_ij}').execute()
        #print(str(parmed.tools.printLJMatrix(new_parms[i], f':{residue_mask}')))

    return new_parms

def get_new_Parms(parms_list, residue_mask, propty, functional, windows, truncate):
    """Changes charge or GB radius of mutating residues according to functional for the different windows in each parameter file"""
    propty_pd_to_parmed = {
        'Charge': 'CHARGE',
        'GB Radius': 'RADII',
        'GB Screen': 'SCREEN'
    }
    GBRadius_minimum = 0.1
    for i in range(len(windows)):
        residue_details = get_data_n_general.details_str_to_pd(str(parmed.tools.printDetails(parms_list[i], f':{residue_mask}')))
        atom_numbers = residue_details['ATOM'].tolist() # Get atom numbers of mutated residues
        multiplier = get_data_n_general.get_multiplier(windows[i], functional, truncate)
        for atom in atom_numbers:
            value = residue_details[residue_details['ATOM'] == atom][propty].iloc[0]
            # Take into account the lowest limit of 0.1 Angstrom for GB Radius
            if propty != 'GB Radius':
                new_value = multiplier*value
            elif propty == 'GB Radius':
                new_value = multiplier*(value - GBRadius_minimum) + GBRadius_minimum
                #print(multiplier)
            parmed.tools.change(parms_list[i], propty_pd_to_parmed[propty], f'@{atom}', f'{new_value}').execute()
        #print(str(parmed.tools.printDetails(parms_list[i], f':{residue_mask}')))
    return parms_list

def change_atom_charge(residue, parm_inter, window, functional, new_charge, atom, extra_charge=0):
    old_charge = parmed.tools.netCharge(parm_inter, f':{residue}@{atom}').execute()
    old_charge = round(old_charge, 4)
    diff_charge = new_charge - old_charge + extra_charge
    multiplier = get_data_n_general.get_multiplier(window, functional, truncate=True)
    inter_charge = old_charge + (1 - multiplier)*diff_charge
    parmed.tools.change(parm_inter, 'Charge', f':{residue}@{atom}', f'{inter_charge}').execute()

def get_CA_Parms(parms_list, residue_position, functional, windows):
    """Changes charge of CA carbon to keep residue charge constant"""
    def compensate_residue(residue, parm_inter, window, end=None):
        """Compensate the charge of each mutating residue to GLY"""
        # Charge of GLY hydrogens bound to CA
        if end == None:
            charge_GLY_N = -0.4157
            charge_GLY_H = 0.2719
            charge_GLY_CA = -0.0252
            charge_GLY_HAs = 0.1396
            charge_GLY_C = 0.5973
            charge_GLY_O = -0.5679
        elif end == 'C_terminus':
            charge_GLY_N = -0.3821
            charge_GLY_H = 0.2681
            charge_GLY_CA = -0.2493
            charge_GLY_HAs = 0.2112
            charge_GLY_C = 0.7231
            charge_GLY_O = -0.7855
            charge_GLY_OXT = -0.7855
        elif end == 'N_terminus':
            charge_GLY_N = 0.2943
            charge_GLY_H1 = 0.1642
            charge_GLY_H2 = 0.1642
            charge_GLY_H3 = 0.1642
            charge_GLY_CA = -0.0100
            charge_GLY_HAs = 0.1790
            charge_GLY_C = 0.6163
            charge_GLY_O = -0.5722

        charges = {
            'N': charge_GLY_N,
            'CA': charge_GLY_CA,
            'C': charge_GLY_C,
            'O': charge_GLY_O
        }

        if end == None:
            charges['H'] = charge_GLY_H
        elif end == 'C_terminus':
            charges['H'] = charge_GLY_H
            charges['OXT'] = charge_GLY_OXT
        elif end == 'N_terminus':
            charges['H1'] = charge_GLY_H1
            charges['H2'] = charge_GLY_H2
            charges['H3'] = charge_GLY_H3

        for i in charges.keys():
            if i != 'CA':
                change_atom_charge(residue, parm_inter, window, functional, charges[i], i)
            else:
                change_atom_charge(residue, parm_inter, window, functional, charges[i], i, charge_GLY_HAs)

    for i in range(len(parms_list)):
        if not isinstance(residue_position, list):
            no_oxt = int(str(parmed.tools.printDetails(parms_list[i], f':{residue_position}@OXT')).strip().split('\n')[0].split()[-2])
            no_h1 = int(str(parmed.tools.printDetails(parms_list[i], f':{residue_position}@H1')).strip().split('\n')[0].split()[-2])
            no_h2 = int(str(parmed.tools.printDetails(parms_list[i], f':{residue_position}@H2')).strip().split('\n')[0].split()[-2])
            no_h3 = int(str(parmed.tools.printDetails(parms_list[i], f':{residue_position}@H3')).strip().split('\n')[0].split()[-2])
            if no_oxt == 1:
                compensate_residue(residue_position, parms_list[i], windows[i], end='C_terminus')
            elif no_h1 == 1 and no_h2 == 1 and no_h3 == 1:
                compensate_residue(residue_position, parms_list[i], windows[i], end='N_terminus')
            else:
                compensate_residue(residue_position, parms_list[i], windows[i], end=None)
        else:
            for residue in residue_position:
                no_oxt = int(str(parmed.tools.printDetails(parms_list[i], f':{residue}@OXT')).strip().split('\n')[0].split()[-2])
                no_h1 = int(str(parmed.tools.printDetails(parms_list[i], f':{residue}@H1')).strip().split('\n')[0].split()[-2])
                no_h2 = int(str(parmed.tools.printDetails(parms_list[i], f':{residue}@H2')).strip().split('\n')[0].split()[-2])
                no_h3 = int(str(parmed.tools.printDetails(parms_list[i], f':{residue}@H3')).strip().split('\n')[0].split()[-2])
                if no_oxt == 1:
                    compensate_residue(residue, parms_list[i], windows[i], end='C_terminus')
                elif no_h1 == 1 and no_h2 == 1 and no_h3 == 1:
                    compensate_residue(residue, parms_list[i], windows[i], end='N_terminus')
                else:
                    compensate_residue(residue, parms_list[i], windows[i], end=None)
    return parms_list

def get_CB_Parms(parms_list, residue_position, functional, windows):
    """Changes charge of CB carbon to keep residue charge constant"""  
    def compensate_residue(residue, parm_inter, window, end=None):
        """Compensate the charge of each mutating residue to ALA"""
        # Charge of ALA hydrogens bound to CB
        if end == None:
            charge_ALA_N = -0.4157
            charge_ALA_H = 0.2719
            charge_ALA_CA = 0.0337
            charge_ALA_HA = 0.0823
            charge_ALA_CB = -0.1825
            charge_ALA_HBs = 0.1809
            charge_ALA_C = 0.5973
            charge_ALA_O = -0.5679
        elif end == 'C_terminus':
            charge_ALA_N = -0.3821
            charge_ALA_H = 0.2681
            charge_ALA_CA = -0.1747
            charge_ALA_HA = 0.1067
            charge_ALA_CB = -0.2093
            charge_ALA_HBs = 0.2292
            charge_ALA_C = 0.7731
            charge_ALA_O = -0.8055
            charge_ALA_OXT = -0.8055
        elif end == 'N_terminus':
            charge_ALA_N = 0.1414
            charge_ALA_H1 = 0.1997
            charge_ALA_H2 = 0.1997
            charge_ALA_H3 = 0.1997
            charge_ALA_CA = 0.0962
            charge_ALA_HA = 0.0889
            charge_ALA_CB = -0.0597
            charge_ALA_HBs = 0.0900
            charge_ALA_C = 0.6163
            charge_ALA_O = -0.5722
        
        charges = {
            'N': charge_ALA_N,
            'CA': charge_ALA_CA,
            'HA': charge_ALA_HA,
            'CB': charge_ALA_CB,
            'C': charge_ALA_C,
            'O': charge_ALA_O
        }

        if end == None:
            charges['H'] = charge_ALA_H
        elif end == 'C_terminus':
            charges['H'] = charge_ALA_H
            charges['OXT'] = charge_ALA_OXT
        elif end == 'N_terminus':
            charges['H1'] = charge_ALA_H1
            charges['H2'] = charge_ALA_H2
            charges['H3'] = charge_ALA_H3

        for i in charges.keys():
            if i != 'CB':
                change_atom_charge(residue, parm_inter, window, functional, charges[i], i)
            else:
                change_atom_charge(residue, parm_inter, window, functional, charges[i], i, charge_ALA_HBs)

    for i in range(len(parms_list)):
        if not isinstance(residue_position, list):
            no_oxt = int(str(parmed.tools.printDetails(parms_list[i], f':{residue_position}@OXT')).strip().split('\n')[0].split()[-2])
            no_h1 = int(str(parmed.tools.printDetails(parms_list[i], f':{residue_position}@H1')).strip().split('\n')[0].split()[-2])
            no_h2 = int(str(parmed.tools.printDetails(parms_list[i], f':{residue_position}@H2')).strip().split('\n')[0].split()[-2])
            no_h3 = int(str(parmed.tools.printDetails(parms_list[i], f':{residue_position}@H3')).strip().split('\n')[0].split()[-2])
            if no_oxt == 1:
                compensate_residue(residue_position, parms_list[i], windows[i], end='C_terminus')
            elif no_h1 == 1 and no_h2 == 1 and no_h3 == 1:
                compensate_residue(residue_position, parms_list[i], windows[i], end='N_terminus')
            else:
                compensate_residue(residue_position, parms_list[i], windows[i], end=None)
        else:
            for residue in residue_position:
                no_oxt = int(str(parmed.tools.printDetails(parms_list[i], f':{residue}@OXT')).strip().split('\n')[0].split()[-2])
                no_h1 = int(str(parmed.tools.printDetails(parms_list[i], f':{residue}@H1')).strip().split('\n')[0].split()[-2])
                no_h2 = int(str(parmed.tools.printDetails(parms_list[i], f':{residue}@H2')).strip().split('\n')[0].split()[-2])
                no_h3 = int(str(parmed.tools.printDetails(parms_list[i], f':{residue}@H3')).strip().split('\n')[0].split()[-2])
                if no_oxt == 1:
                    compensate_residue(residue, parms_list[i], windows[i], end='C_terminus')
                elif no_h1 == 1 and no_h2 == 1 and no_h3 == 1:
                    compensate_residue(residue, parms_list[i], windows[i], end='N_terminus')
                else:
                    compensate_residue(residue, parms_list[i], windows[i], end=None)
    return parms_list

def create_intermediate_parms(functions, windows, residue_position, intermediate, include_mut, recalculation=False):
    """Creates intermediate parameter files using scaling of function on residues given by residue_position"""
    # Leap generated parameter files
    wt_parm_path = 'setup/parms_n_pdbs/parms/parms_windows/wt_0.parm7'
    wt_rst_path = 'setup/parms_n_pdbs/parms/rst_windows/wt_0.rst7'
    mt_parm_path = 'setup/parms_n_pdbs/parms/parms_windows/mt_0.parm7'
    mt_rst_path = 'setup/parms_n_pdbs/parms/rst_windows/mt_0.rst7'

    # Import to ParmEd
    wt_parmed = parmed.amber.AmberParm(wt_parm_path, wt_rst_path)
    if include_mut:
        mt_parmed = parmed.amber.AmberParm(mt_parm_path, mt_rst_path)

    #Create mutated residues masks
    if not isinstance(residue_position, list):
        residue_mask = f'{residue_position}'
    else:
        residue_mask_list = [str(x) for x in residue_position]
        residue_mask = ','.join(residue_mask_list)
    # Select mask of intermediate
    if intermediate == 'GLY':
        residue_mask_nobackbone = residue_mask + '&!@CA,C,O,N,H,H1,H2,H3,OXT'
    elif intermediate == 'ALA':
        residue_mask_nobackbone = residue_mask + '&!@CA,C,O,N,H,HA,CB,H1,H2,H3,OXT'

    # Create parms with modified LJ matrix according to windows and functions
    wt_parms_LJ = get_new_LJParms(wt_parmed, residue_mask_nobackbone, functions[-2:], windows[1:])
    if include_mut:
        mt_parms_LJ = get_new_LJParms(mt_parmed, residue_mask_nobackbone, functions[-2:], windows[1:])

    # Change charge of mutating residues according to windows and functions
    wt_parms_ele = get_new_Parms(wt_parms_LJ, residue_mask_nobackbone, 'Charge', functions[1], windows[1:], truncate=True)
    if include_mut:
        mt_parms_ele = get_new_Parms(mt_parms_LJ, residue_mask_nobackbone, 'Charge', functions[1], windows[1:], truncate=True)

    # Change GB Radius of mutating residues according to windws and functions
    wt_parms_GB = get_new_Parms(wt_parms_ele, residue_mask_nobackbone, 'GB Radius', functions[0], windows[1:], truncate=False)
    if include_mut:
        mt_parms_GB = get_new_Parms(mt_parms_ele, residue_mask_nobackbone, 'GB Radius', functions[0], windows[1:], truncate=False)

    if intermediate == 'GLY':
        wt_parms_CA = get_CA_Parms(wt_parms_GB, residue_position, functions[1], windows[1:])
        if include_mut:
            mt_parms_CA = get_CA_Parms(mt_parms_GB, residue_position, functions[1], windows[1:])
    elif intermediate == 'ALA':
        wt_parms_CA = get_CB_Parms(wt_parms_GB, residue_position, functions[1], windows[1:])
        if include_mut:
            mt_parms_CA = get_CB_Parms(mt_parms_GB, residue_position, functions[1], windows[1:])

    #print(parmed.tools.printDetails(wt_parmed, f':{residue_mask}&!@C,O,N,H'))
    for i in range(len(wt_parms_CA)):
        #print(parmed.tools.printDetails(wt_parms_CA[i], f':{residue_mask}&!@C,O,N,H'))
        parmed.tools.HMassRepartition(wt_parms_CA[i]).execute()
        if recalculation:
            parmed.tools.outparm(wt_parms_CA[i], f'setup/recalculate/wt_{i+1}.parm7').execute()
        else:
            parmed.tools.outparm(wt_parms_CA[i], f'setup/parms_n_pdbs/parms/parms_windows/wt_{i+1}.parm7').execute()
        if include_mut:
            parmed.tools.HMassRepartition(mt_parms_CA[i]).execute()
            if recalculation:
                parmed.tools.outparm(mt_parms_CA[i], f'setup/recalculate/mt_{i+1}.parm7').execute()
            else:
                parmed.tools.outparm(mt_parms_CA[i], f'setup/parms_n_pdbs/parms/parms_windows/mt_{i+1}.parm7').execute()