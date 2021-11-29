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
    wt_parm_path = 'setup/parms_n_pdbs/parms/parms_windows/wt_0.parm7'
    wt_rst_path = 'setup/parms_n_pdbs/parms/rst_windows/wt_0.rst7'
    mt_parm_path = 'setup/parms_n_pdbs/parms/parms_windows/mt_0.parm7'
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
    wt_parm_path = 'setup/parms_n_pdbs/parms/parms_windows/wt_0.parm7'
    wt_rst_path = 'setup/parms_n_pdbs/parms/rst_windows/wt_0.rst7'
    mt_parm_path = 'setup/parms_n_pdbs/parms/parms_windows/mt_0.parm7'
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
                residues_details = get_data_n_general.details_str_to_pd(str(parmed.tools.printDetails(wt_or_mt, f':{mask}')))
            atom_numbers = residues_details['ATOM'].tolist() # Get atom numbers of mutated residues

            for atom in atom_numbers:
                value = residues_details[residues_details['ATOM'] == atom]['GB Radius'].iloc[0]
                if proportion > 0:
                    new_value = value + value*proportion
                else:
                    new_value = value - value*proportion
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

    parmed.tools.outparm(wt_parmed, wt_parm_path).execute()
    if include_mut:
        parmed.tools.outparm(mt_parmed, mt_parm_path).execute()

def get_new_LJParms(parmed_object, residue_mask, functions, windows):
    """Creates new LJ atom types and outputs new AmberParm bojects with modified LJ  matrix"""
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

def get_CA_Parms(parms_list, residue_position, functional, windows):
    """Changes charge of CA carbon to keep residue charge constant"""
    def compensate_residue(residue, parm_inter, window):
        """Compensate the charge of each mutating residue to GLY"""
        # Charge of GLY hydrogens bound to CA
        charge_GLY_hydrogens = 0.1396
        charge_GLY_CA = -0.0252
        charge_CA = parmed.tools.netCharge(parm_inter, f':{residue}@CA').execute()
        charge_CA = round(charge_CA, 4)
        diff_charge_CA = charge_GLY_CA - charge_CA
        diff_charge = charge_GLY_hydrogens + diff_charge_CA
        multiplier = get_data_n_general.get_multiplier(window, functional, truncate=True)
        new_charge_CA = charge_CA + (1 - multiplier)*diff_charge
        # Change CA charge
        parmed.tools.change(parm_inter, 'Charge', f':{residue}@CA', f'{new_charge_CA}').execute()

    for i in range(len(parms_list)):
        if not isinstance(residue_position, list):
            compensate_residue(residue_position, parms_list[i], windows[i])
        else:
            for residue in residue_position:
                compensate_residue(residue, parms_list[i], windows[i])
    return parms_list

def get_CB_Parms(parms_list, residue_position, functional, windows):
    """Changes charge of CB carbon to keep residue charge constant"""
    def compensate_residue(residue, parm_inter, window):
        """Compensate the charge of each mutating residue to ALA"""
        # Charge of ALA hydrogens bound to CB
        charge_ALA_hydrogens = 0.1809
        charge_ALA_CB = -0.1825
        charge_CB = parmed.tools.netCharge(parm_inter, f':{residue}@CB').execute()
        charge_CB = round(charge_CB, 4)
        diff_charge_CB = charge_ALA_CB - charge_CB
        diff_charge = charge_ALA_hydrogens + diff_charge_CB
        multiplier = get_data_n_general.get_multiplier(window, functional, truncate=True)
        new_charge_CB = charge_CB + (1 - multiplier)*diff_charge
        # Change CB charge
        parmed.tools.change(parm_inter, 'Charge', f':{residue}@CB', f'{new_charge_CB}').execute()
        # Charge of CA to get that of ALA
        charge_ALA_CA = 0.0337
        charge_CA = parmed.tools.netCharge(parm_inter, f':{residue}@CA').execute()
        charge_CA = round(charge_CA, 4)
        diff_charge_CA = charge_ALA_CA - charge_CA
        new_charge_CA = charge_CA + (1 - multiplier)*diff_charge_CA
        parmed.tools.change(parm_inter, 'Charge', f':{residue}@CA', f'{new_charge_CA}').execute()
        # Charge of HA to get that of ALA
        charge_ALA_HA = 0.0823
        charge_HA = parmed.tools.netCharge(parm_inter, f':{residue}@HA').execute()
        charge_HA = round(charge_HA, 4)
        diff_charge_HA = charge_ALA_HA - charge_HA
        new_charge_HA = charge_HA + (1 - multiplier)*diff_charge_HA
        parmed.tools.change(parm_inter, 'Charge', f':{residue}@HA', f'{new_charge_HA}').execute()

    for i in range(len(parms_list)):
        if not isinstance(residue_position, list):
            compensate_residue(residue_position, parms_list[i], windows[i])
        else:
            for residue in residue_position:
                compensate_residue(residue, parms_list[i], windows[i])
    return parms_list

def create_intermediate_parms(functions, windows, residue_position, intermediate, include_mut):
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
        residue_mask_nobackbone = residue_mask + '&!@CA,C,O,N,H'
    elif intermediate == 'ALA':
        residue_mask_nobackbone = residue_mask + '&!@CA,C,O,N,H,HA,CB'

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
        parmed.tools.outparm(wt_parms_CA[i], f'setup/parms_n_pdbs/parms/parms_windows/wt_{i+1}.parm7').execute()
        if include_mut:
            parmed.tools.HMassRepartition(mt_parms_CA[i]).execute()
            parmed.tools.outparm(mt_parms_CA[i], f'setup/parms_n_pdbs/parms/parms_windows/mt_{i+1}.parm7').execute()