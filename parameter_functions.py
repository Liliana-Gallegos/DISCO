#### Author: Liliana C. Gallegos ####
import numpy as np
import pandas as pd

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

    if file.find("nbo") > -1 or file.find("NMR")> -1:
        # Get xyz data
        for i in range(0,len(data)):
            if data[i].find("Symbolic Z-matrix") > -1:
                start = i+2
        for j in range(start,len(data)):
            if data[j].find("Input orientation") > -1 or data[j].find(" Stoichiometry") > -1:
                end = j-1
                break
        xyz_data = data[start:end]

        atoms = []  # element
        coordinates = []  # xyz coords
        for atom in xyz_data:
            atom_data = atom.split()
            atom_n = atom_data[0]
            atoms.append(atom_n)
            x, y, z = atom_data[1], atom_data[2], atom_data[3]
            coordinates.append( [float(x), float(y), float(z)] )
        return atoms, coordinates

    else:
        print('File not parsed: ', file)
        pass

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

    if file.find("NMR") > -1:
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

    if file.find("nbo") > -1:
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

    if file.find("nbo") > -1:
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
