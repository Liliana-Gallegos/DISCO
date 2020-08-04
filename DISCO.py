#!/usr/bin/python
#### Author: Liliana C. Gallegos ####
import os, sys
import re
import pandas as pd
import numpy as np
import argparse

#import parameter_functions as pf

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

## Main code:
def main():
    parser = argparse.ArgumentParser(description='The script analysis ...')
    parser.add_argument('file', help="Analysis NBO and NMR files.", nargs='*')
    parser.add_argument('--distance', default="False", help="Gives bond length of given X atom in Angstrom.")
    parser.add_argument('--nmr', default="False", help="Gives NMR shielding tensor values of given X or (number)X atom. ")
    parser.add_argument('--charge', default="False", help="Gives atom charge value of given X or (number)X atom.")
    parser.add_argument('--fmo', default="False", help="Gives the frontier molecular orbitals for molecule - options: HOMO, LUMO, or both.")
    parser.add_argument('--csv', default="False", help="Prints out csv file of collected data - options: charge, nmr, fmo, nbo, or all.")
    args = parser.parse_args()

    files = args.file
    # else:
    #     print("\nNo files were found.\n")

    s_nmr_name = []
    s_fmo_name = []
    s_chrg_name = []
    num_atom_fmo = []
    num_atom_chrg = []
    num_atom_nmr = []
    s_distance = []
    s_charge = []
    s_nmr = []
    s_homo = []
    s_lumo = []

    for file in files:

        # Get molecule name
        file_name = os.path.basename(file)
        split_filename = file_name.split('.')
        name = split_filename[0]
        print('\n',name,':\n')
        molecule_name = name.split('_')
        molecule = molecule_name[0]

        # Frontier Molecular Orbitals: HOMO/LUMO values
        if homo_lumo_values(file) == None: pass
        else:
            if args.fmo != "False":
                s_fmo_name.append(molecule)
                homo, lumo = homo_lumo_values(file)
                if args.fmo.upper() == "HOMO":
                    homo_len = round(homo, 4)
                    s_homo.append(homo_len)
                    print(F'HOMO-value: {homo:.4f}')
                elif args.fmo.upper() == "LUMO":
                    lumo_len = round(lumo, 4)
                    s_lumo.append(lumo_len)
                    print(F'LUMO-value: {lumo:.4f}')
                elif args.fmo.lower() == "both":
                    homo_len = round(homo, 4)
                    s_homo.append(homo_len)
                    lumo_len = round(lumo, 4)
                    s_lumo.append(lumo_len)
                    print(F'LUMO-value: {lumo:.4f}')
                    print(F'HOMO-value: {homo:.4f}')

        # Atom charges
        if atom_charge(file) == None: pass
        else:
            if args.charge != "False":
                element = args.charge
                el = re.split('(\d+)', element)
                if len(el) == 1:
                    atoms, numbers, natural_charges = atom_charge(file)
                    # for i in range(len(atoms)):
                        # print('one',atoms[i])
                        # if atoms[i].find(element) > -1:
                        #     print(i)
                    atom_i = [i for i in range(len(atoms)) if element in atoms[i]]

                    for j in atom_i:
                        s_chrg_name.append(molecule)
                        chrg_len = round(natural_charges[j], 5)
                        s_charge.append(chrg_len)
                        n_a = numbers[j]+atoms[j]
                        num_atom_chrg.append(n_a)
                        print(F'{numbers[j]}{atoms[j]}-charge: {natural_charges[j]:.5f}')
                else:
                    s_chrg_name.append(molecule)
                    n, a = el[1], el[2]  # number of atom, atom
                    atoms, numbers, natural_charges = atom_charge(file)
                    j = numbers.index(n) # find by atom number
                    chrg_len = round(natural_charges[j], 5)
                    s_charge.append(chrg_len)
                    n_a = numbers[j]+atoms[j]
                    num_atom_chrg.append(n_a)
                    print(F'{numbers[j]}{atoms[j]}-charge: {natural_charges[j]:.5f}')

        # NMR shielding
        if nmr_shielding(file) == None: pass
        else:
            if args.nmr != "False":
                element = args.nmr
                el = re.split('(\d+)', element)
                if len(el) == 1:
                    atoms, numbers, shielding_tensor = nmr_shielding(file)
                    atom_i = [i for i in range(len(atoms)) if element in atoms[i]]     #atoms.index(element)
                    for j in atom_i:
                        s_nmr_name.append(molecule)
                        nmr_len = round(shielding_tensor[j], 4)
                        s_nmr.append(nmr_len)
                        n_a = numbers[j]+atoms[j]
                        num_atom_nmr.append(n_a)
                        print(F'{numbers[j]}{atoms[j]}-shielding: {shielding_tensor[j]:.4f}')
                else:
                    n, a = el[1], el[2]  # number of atom, atom
                    s_nmr_name.append(molecule)
                    atoms, numbers, shielding_tensor = nmr_shielding(file)
                    j = numbers.index(n)
                    nmr_len = round(shielding_tensor[j], 4)
                    s_nmr.append(nmr_len)
                    n_a = numbers[j]+atoms[j]
                    num_atom_nmr.append(n_a)
                    print(F'{numbers[j]}{atoms[j]}-shielding: {shielding_tensor[j]:.4f}')
            # else: pass

        # Find bond lengths
        if par_xyz_data(file) == None: pass
        else:
            atoms, coords = par_xyz_data(file)
            num_atoms = len(atoms)

            for num1 in range(0,num_atoms):
                for num2 in range(0,num_atoms):
                    if num1<num2:
                        # Bond distance
                        if args.distance != "False":
                            element = args.distance
                            atom_distance = calculate_distance(coords[num1], coords[num2])
                            if atoms[num1].find(element) > -1 or atoms[num2].find(element) > -1:
                                if element == "Pd":
                                    if bond_check(atom_distance, min_value=0, max_value=2.5) is True:
                                        s_len = round(atom_distance, 3)
                                        s_distance.append(s_len)
                                        print(F'{num1+1}{atoms[num1]} to {num2+1}{atoms[num2]} : {atom_distance:.3f}')
                                if element == "N":
                                    if bond_check(atom_distance, min_value=0, max_value=1.5) is True:
                                        s_len = round(atom_distance, 3)
                                        s_distance.append(s_len)
                                        print(F'{num1+1}{atoms[num1]} to {num2+1}{atoms[num2]} : {atom_distance:.3f}')
                                if element == "Cl":
                                    if bond_check(atom_distance, min_value=1.5, max_value=1.8) is True:
                                        s_len = round(atom_distance, 3)
                                        s_distance.append(s_len)
                                        print(F'{num1+1}{atoms[num1]} to {num2+1}{atoms[num2]} : {atom_distance:.3f}')
                                if element == "H":
                                    if bond_check(atom_distance, min_value=0, max_value=1.5) is True:
                                        s_len = round(atom_distance, 3)
                                        s_distance.append(s_len)
                                        print(F'{num1+1}{atoms[num1]} to {num2+1}{atoms[num2]} : {atom_distance:.3f}')
                                if element == "C":
                                    if bond_check(atom_distance, min_value=0, max_value=1.8) is True:
                                        s_len = round(atom_distance, 3)
                                        s_distance.append(s_len)
                                        print(F'{num1+1}{atoms[num1]} to {num2+1}{atoms[num2]} : {atom_distance:.3f}')

    if args.csv != False:
        charge_data = {'Name': s_chrg_name, 'atom_number': num_atom_chrg, 'charge': s_charge}
        fmo_data = {'Name': s_fmo_name, 'HOMO': s_homo, 'LUMO': s_lumo}
        nmr_data = {'Name': s_nmr_name, 'atom_number': num_atom_nmr, 'NMR': s_nmr}
        if args.csv.lower() == "charge":
            write_tab_csv(charge_data, 'Charge_data')
        if args.csv.lower() == "fmo":
            if len(fmo_data['HOMO']) == 0:
                del fmo_data['HOMO']
            if len(fmo_data['LUMO']) == 0:
                del fmo_data['LUMO']
            write_tab_csv(fmo_data, 'FMO_data')
        if args.csv.lower() == "nmr":
            write_tab_csv(nmr_data, 'NMR_data')

        if args.csv.lower() == "nbo":
            if len(fmo_data['HOMO']) == 0:
                del fmo_data['HOMO']
            if len(fmo_data['LUMO']) == 0:
                del fmo_data['LUMO']
            df_fmo = pd.DataFrame(fmo_data)
            df_chrg = pd.DataFrame(charge_data)
            df_nbo = df_chrg.join(df_fmo.set_index('Name'), on='Name').sort_values(by=['Name'])
            print(df_nbo)
            # write to csv
            df_nbo.to_csv('NBO_data.csv', index=False)

        if args.csv.lower() == "all":
            if len(fmo_data['HOMO']) == 0:
                del fmo_data['HOMO']
            if len(fmo_data['LUMO']) == 0:
                del fmo_data['LUMO']
            df_fmo = pd.DataFrame(fmo_data)
            df_chrg = pd.DataFrame(charge_data)
            df_nmr = pd.DataFrame(nmr_data)
            df_merged = pd.merge(df_nmr, df_chrg, on=['Name', 'atom_number'])
            df_all = df_merged.join(df_fmo.set_index('Name'), on='Name')
            print(df_all)
            # write to csv
            df_all.to_csv('All_data.csv', index=False)

if __name__ == "__main__":
    main()

#TEST FINAL PRINTS
# print('nmr file',s_nmr_name)
# print('chrg file',s_chrg_name)
# print('fmo file',s_fmo_name)
# print('nmr values',s_nmr)
# print('distance',s_distance)
# print('charge',s_charge)
# print('homo',s_homo)
# print('lumo',s_lumo)
# print('num. and atom-chrg', num_atom_chrg)
# print('num. and atom-nmr', num_atom_nmr)
# print('num. and atom-fmo', num_atom_fmo)
