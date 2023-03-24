#!/usr/bin/python
#### Functions for DISCO by Liliana C. Gallegos

import numpy as np

periodictable = ["","H","He","Li","Be","B","C","N","O","F","Ne","Na","Mg","Al","Si","P","S","Cl","Ar","K","Ca","Sc","Ti","V","Cr","Mn","Fe","Co","Ni","Cu","Zn","Ga","Ge","As","Se","Br","Kr","Rb","Sr","Y","Zr",
"Nb","Mo","Tc","Ru","Rh","Pd","Ag","Cd","In","Sn","Sb","Te","I","Xe","Cs","Ba","La","Ce","Pr","Nd","Pm","Sm","Eu","Gd","Tb","Dy","Ho","Er","Tm","Yb","Lu","Hf","Ta","W","Re","Os","Ir","Pt","Au","Hg","Tl",
"Pb","Bi","Po","At","Rn","Fr","Ra","Ac","Th","Pa","U","Np","Pu","Am","Cm","Bk","Cf","Es","Fm","Md","No","Lr","Rf","Db","Sg","Bh","Hs","Mt","Ds","Rg","Uub","Uut","Uuq","Uup","Uuh","Uus","Uuo"]

## Bondi VDW radii in Angstrom
bondi = {"Bq": 0.00, "H": 1.09,"He": 1.40,
         "Li":1.81,"Be":1.53,"B":1.92,"C":1.70,"N":1.55,"O":1.52,"F":1.47,"Ne":1.54,
         "Na":2.27,"Mg":1.73,"Al":1.84,"Si":2.10,"P":1.80,"S":1.80,"Cl":1.75,"Ar":1.88,
         "K":2.75,"Ca":2.31,"Ni": 1.63,"Cu":1.40,"Zn":1.39,"Ga":1.87,"Ge":2.11,"As":1.85,"Se":1.90,"Br":1.83,"Kr":2.02,
         "Rb":3.03,"Sr":2.49,"Pd": 1.63,"Ag":1.72,"Cd":1.58,"In":1.93,"Sn":2.17,"Sb":2.06,"Te":2.06,"I":1.98,"Xe":2.16,
         "Cs":3.43,"Ba":2.68,"Pt":1.72,"Au":1.66,"Hg":1.55,"Tl":1.96,"Pb":2.02,"Bi":2.07,"Po":1.97,"At":2.02,"Rn":2.20,
         "Fr":3.48,"Ra":2.83, "U":1.86 }

## covalent radii in Angstrom (taken from Pyykko and Atsumi, Chem. Eur. J. 15, 2009, 188-197)
rcov = {"H": 0.32,"He": 0.46,
        "Li":1.33,"Be":1.02,"B":0.85,"C":0.75,"N":0.71,"O":0.63,"F":0.64,"Ne":0.67,
        "Na":1.55,"Mg":1.39,"Al":1.26, "Si":1.16,"P":1.11,"S":1.03,"Cl":0.99, "Ar":0.96,
        "K":1.96,"Ca":1.71,"Sc": 1.48, "Ti": 1.36, "V": 1.34, "Cr": 1.22, "Mn":1.19, "Fe":1.16, "Co":1.11, "Ni":1.10,"Zn":1.18, "Ga":1.24, "Ge":1.21, "As":1.21, "Se":1.16, "Br":1.14, "Kr":1.17,
        "Rb":2.10, "Sr":1.85,"Y":1.63, "Zr":1.54, "Nb":1.47, "Mo":1.38, "Tc":1.28, "Ru":1.25,"Rh":1.25,"Pd":1.20,"Ag":1.28,"Cd":1.36, "In":1.42,"Sn":1.40,"Sb":1.40,"Te":1.36,"I":1.33,"Xe":1.31,
        "Bi":1.46} # Bismuth added manually

def read_xyz_file(filename): #, look_for_charge=True
    """
    From xyz2mol.py
    """

    atomic_symbols = []
    xyz_coordinates = []
    #charge = 0
    title = ""

    with open(filename, "r") as file:
        for line_number, line in enumerate(file):
            if line_number == 0: num_atoms = int(line)
            elif line_number == 1:
                title = line
                if "charge=" in line:
                    charge = int(line.split("=")[1])
            else:
                atomic_symbol, x, y, z = line.split()
                atomic_symbols.append(atomic_symbol)
                xyz_coordinates.append([float(x), float(y), float(z)])
    return atomic_symbols, xyz_coordinates


def log2xyz(log_file):
    ''' Converts log files into xyz files. '''
    name, ext = os.path.splitext(log_file)

    #### Create xyz file
    xyz_file = f'{name}.xyz'
    gv_command = ['python', '-m', 'goodvibes', f'{log_file}', '--xyz']
    rename_command = ['mv', 'Goodvibes_output.xyz', f'{xyz_file}']
    subprocess.call(gv_command)
    subprocess.call(rename_command)

    #### Grab charge
    charge = 0
    with open(log_file, "r") as file:
        for line_number, line in enumerate(file):
            if "Charge =" in line:
                charge = int(line.split(' ')[4])

    if not os.path.exists(xyz_file): xyz_file = 'None'
    return xyz_file, charge

def par_xyz_data(file):
    """
    Opens a Gaussian NBO or NMR calculation file, seperates the coordinates and the atoms and recasts the coordinates as floats.
    Inputs: Gaussian (filename)_nbo.log or (filename)_NMR.log file.
        NOTE: Case sensitive.
    Returns: atom and coordinates
    """
    # Reads ALL data --> for xyz coords
    outfile = open(file, 'r')
    data = outfile.readlines()
    outfile.close()
    if file.find("nbo".upper()) > -1 or file.find("nmr".upper())> -1 or file.find("NBO".lower()) > -1 or file.find("NMR".lower())> -1 :
        # Get xyz data
        for i in range(0,len(data)):
            if data[i].find("Symbolic Z-matrix") > -1:
                start = i+2
        for j in range(start,len(data)):
            if data[j].find("Input orientation") > -1 or data[j].find("Symmetry") > -1 or data[j].find(" Stoichiometry") > -1:
                end = j-1
                break
    else:
        iaxes_list = []
        #### For optimization calculations
        for i in range(0,len(data)):
            if data[i].find("Standard orientation") > -1:
                start = i+5
                iaxes_list.append(start)
        #print('start:', iaxes_list[-1])
        for j in range(iaxes_list[-1],len(data)):
            if data[j].find("Rotational constants") > -1:
                end = j-2
                break
    xyz_data = data[start:end]

    atoms = []  # element
    coordinates = []  # xyz coords
    for atom in xyz_data:
        atom_data = atom.split()
        if len(atom_data) > 4:
            atoms.append(atom_data[1])
            x, y, z = atom_data[3], atom_data[4], atom_data[5]
            coordinates.append( [float(x), float(y), float(z)] )
        else:
            atoms.append(atom_data[0])
            x, y, z = atom_data[1], atom_data[2], atom_data[3]
            coordinates.append( [float(x), float(y), float(z)] )

    if atoms[0].isnumeric() == True: atoms = [ periodictable[int(a)] for a in atoms ]
    return atoms, coordinates


def gen_connectivity(coords, mol_atoms, vdw_frac, cov_frac):

    ## Use VDW radii to infer atom connectivity
    bonded_atoms, bond_distance = [], []
    for i, elem_i in enumerate(mol_atoms):
        for j, elem_j in enumerate(mol_atoms):
            if j > i:
                vdw_ij = bondi[elem_i] + bondi[elem_j]
                rcov_ij = rcov[elem_i] + rcov[elem_j]
                dist_ij = np.linalg.norm(np.array(coords[i])-np.array(coords[j]))

                if dist_ij / vdw_ij < vdw_frac or dist_ij / rcov_ij < cov_frac:
                    bonded_atoms.append([str(i+1)+elem_i,str(j+1)+elem_j])
                    bond_distance.append(round(dist_ij,5))
                    #print(dist_ij, vdw_ij, rcov_ij)
                else: pass
    return bonded_atoms, bond_distance


def calculate_distance(atom1_coord, atom2_coord):
    """
    Calculates the distance between two points in 3D space (xyz coordinates).
    Inputs: coordinates of two atoms
    Return: distance between the atoms
    """
    x_distance = atom1_coord[0] - atom2_coord[0]
    y_distance = atom1_coord[1] - atom2_coord[1]
    z_distance = atom1_coord[2] - atom2_coord[2]
    atom_distance = np.sqrt(x_distance**2 + y_distance**2 + z_distance**2)
    return atom_distance

def bond_check(distance, min_value=0, max_value=1.55):
    """
    Defaults for minimum values = 0, and maximum values = 1.55 for C-C bond reference.
    """
    if distance > min_value and distance <= max_value:
        return True
    else:
        return False

def nmr_shielding(file):
    """
    Opens NMR calculation file and collects nmr shielding tensors.
    Inputs: Gaussian (filename)_NMR.log file.
        NOTE: Case sensitive.
    Returns: Atoms, atom number, and NMR shielding tensors.
    """
    outfile = open(file, 'r')
    data = outfile.readlines()
    outfile.close()

    if file.find("nmr".upper()) > -1 or file.find("NMR".lower())> -1:
        # Get NMR shielding tensor data
        for i in range(0,len(data)):
            if data[i].find("shielding tensors") > -1:
                start = i+2

        nmr_values = []
        for j in range(start,len(data),5):
            nmr = data[j]
            nmr_values.append(nmr)
            if data[j].find("End of Minotr F.D.") > -1:
                break
        nmr_values = nmr_values[:-1]

        atoms = []
        numbers = []
        nmr_shielding = []
        for line in nmr_values:
            nmr_line = line.split()
            number = nmr_line[0]
            numbers.append(number)
            atom = nmr_line[1]
            atoms.append(atom)
            tensors = float(nmr_line[4])
            nmr_shielding.append(tensors)
        return atoms, numbers, nmr_shielding
    else:
        pass

def scale_chem_shift(slope, intercept, nmr_isotropic):
    """
    Utlizing input scaled factors (slope and intercept) to compute and return scaled chemical shifts.
    """
    new_chem_shift = (intercept-nmr_isotropic)/(-slope)
    #new_chem_shift = intercept+(slope*nmr_isotropic)
    return new_chem_shift


def atom_charge(file):
    """
    Opens NBO calculation file and collects atom charges.
    Inputs: Gaussian (filename)_nbo.log file.
        NOTE: Case sensitive.
    Returns: Atoms, atom number, and natural charges.
    """
    outfile = open(file, 'r')
    data = outfile.readlines()
    outfile.close()

    if file.find("nbo".upper()) > -1 or file.find("NBO".lower()) > -1:
        # Get natural charge data
        for i in range(0,len(data)):
            if data[i].find("Natural Population Analysis") > -1:
                start = i+6

        for j in range(start,len(data)):
            if data[j].find("===") > -1:
                end = j
                break
        charge_data = data[start:end]

        atoms = []
        numbers = []
        natural_charges = []
        for atom in charge_data:
            atom_data = atom.split()
            atom = atom_data[0]
            atoms.append(atom)
            number = atom_data[1]
            numbers.append(number)
            charge = float(atom_data[2])
            natural_charges.append( charge)
        return atoms, numbers, natural_charges
    else:
        pass

def homo_lumo_values(file):
    """
    Opens NBO calculation file and collects HOMO/LUMO values.
    Inputs: Gaussian (filename)_nbo.log file.
        NOTE: Case sensitive.
    Returns: HOMO and LUMO values.
    """
    outfile = open(file, 'r')
    data = outfile.readlines()
    outfile.close()

    if file.find("nbo".upper()) > -1 or file.find("NBO".lower()) > -1:
        # Get LUMO energy data
        for j in range(0,len(data)):
            if data[j].find("Alpha virt. eigenvalues") > -1:
                end = j
                break
        hl_data = data[end-1:end+1]  # Grabs both HOMO/LUMO lines

        homo = hl_data[0]   # For homo, grab first line
        homo_values = homo.split()[-1]  # HOMO last value from line
        homo_values = float(homo_values)

        lumo = hl_data[1]   # For lumo, grab second line
        lumo_values = lumo.split()[4]  # LUMO first value from line
        lumo_values = float(lumo_values)
        return homo_values, lumo_values
    else: pass

def write_tab_csv(dictionary, file_name=False):
    """
    Takes dictionary, forms a pandas table, sorts by 'Name' column, and writes out to csv file.
    Input: {'Name': [list], 'var1': [list1], ...} - needs same length.
    Returns: DataFrame. Also, prints table and writes csv file.
    """
    df_unsort = pd.DataFrame(dictionary)
    df = df_unsort.sort_values(by=['Name'])
    if file_name != False:
        print(df)
        #write to csv file on working directory
        df.to_csv(file_name+'.csv', index=False)
    else:
        return df
