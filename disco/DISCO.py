#!/usr/bin/python
#### Author: Liliana C. Gallegos ####
import os, sys
import re
import pandas as pd
import numpy as np
import subprocess
from glob import glob
from optparse import OptionParser
from disc_functions import *

class disco:
    def __init__(self, *args, **kwargs):
        self.file = args[0]

        if 'options' in kwargs:
            self.options = kwargs['options']
        else:
            self.options = set_options(kwargs)

        file = self.file
        options = self.options

        #### Atomic charges
        charge_list, num_element_list = [], []
        if options.charge != False:
            atoms, numbers, natural_charges = atom_charge(file)

            elements = options.charge
            el = re.split('(\d+)', elements)
            #### parse by element
            if len(el) == 1:
                atom_i = [ i for i, a in enumerate(atoms) if a == elements ]
                for j in atom_i:
                    charge_list.append(round(natural_charges[j], 5))
                    num_element_list.append(str(numbers[j]+atoms[j]))
            else:
                j = numbers.index(el[1]) # find by atom number
                charge_list.append(round(natural_charges[j], 5))
                num_element_list.append(str(numbers[j]+atoms[j]))

        #### Molecular Orbitals: HOMO/LUMO values
        homo_list, lumo_list = [], []
        if options.mo != False:
            homo, lumo = homo_lumo_values(file)
            homo_list.append(round(homo, 4)); lumo_list.append(round(lumo, 4))

        #### NMR tensor shielding values
        nmr_list = []
        if options.nmr != False:
            atoms, numbers, shielding_tensor = nmr_shielding(file)

            #### Scale Factor
            if options.scale != False:
                intercept, slope = float(options.scale.split(',')[0]), float(options.scale.split(',')[1])

            elements = options.nmr
            el = re.split('(\d+)', elements)
            if len(el) == 1:
                atom_i = [ i for i, a in enumerate(atoms) if a == elements ]
                for j in atom_i:
                    if options.scale != False:
                        chemical_shift = scale_chem_shift(slope, intercept, shielding_tensor[j])
                        nmr_list.append(round(chemical_shift, 4))
                        num_element_list.append(str(numbers[j]+atoms[j]))
                    else:
                        nmr_list.append(round(shielding_tensor[j], 4))
                        num_element_list.append(str(numbers[j]+atoms[j]))
            else:
                j = numbers.index(el[1])
                if options.scale != False: shielding_tensor = scale_chem_shift(slope, intercept, shielding_tensor[j])
                else: shielding_tensor = shielding_tensor[j]
                nmr_list.append(round(shielding_tensor, 4))
                num_element_list.append(str(numbers[j]+atoms[j]))


        #### Bond lengths
        bonded_atoms_list, distance_list = [], []
        if options.distance != False:
            elements = options.distance
            a, b = elements.split(',')

            if file.find('.xyz') > -1: atoms, coords = read_xyz_file(file)
            else: atoms, coords = par_xyz_data(file)
            bonded_atoms, bond_distances = gen_connectivity(coords, atoms, vdw_frac=0.5, cov_frac=1.1)

            for bond, distance in zip(bonded_atoms, bond_distances):
                i, j = bond
                if i.find(a) > -1 and j.find(b) > -1 or j.find(a) > -1 and i.find(b) > -1:
                    bonded_atoms_list.append(bond)
                    distance_list.append(distance)


        #### For object reference
        if options.charge:
            self.Chrg = charge_list
            self.Atom = num_element_list
            if options.verbose: print('  Atomic Charge: {} = {}'.format(self.Atom, self.Chrg))

        if options.mo:
            self.HOMO = homo_list
            self.LUMO = lumo_list
            if options.verbose: print('  HOMO and LUMO = {} and {}'.format(self.HOMO, self.LUMO))

        if options.nmr:
            self.NMR  = nmr_list
            self.Atom = num_element_list
            if options.verbose: print('  Atomic NMR Values: {} = {} '.format(self.Atom, self.NMR))

        if options.distance:
            self.Dist = distance_list
            self.Bonded_atoms = bonded_atoms_list
            if options.verbose: print('  Bond Distance: {} = {} '.format(self.Bonded_atoms, self.Dist))


class options_add:
	pass

def set_options(kwargs):
	#set default options and options provided
	options = options_add()
	#dictionary containing default values for options
	var_dict = {'verbose': ['verbose',False], 'v': ['verbose',False], 'charge': ['charge',False], 'nmr':['nmr',False],
                'scale': ['scale',False], 'mo':['mo',False], 'distance': ['distance',False], 'csv':['csv',False] }

	for key in var_dict:
		vars(options)[var_dict[key][0]] = var_dict[key][1]
	for key in kwargs:
		if key in var_dict:
			vars(options)[var_dict[key][0]] = kwargs[key]
		else:
			print("Warning! Option: [", key,":", kwargs[key],"] provided but no option exists, try -h to see available options.")

	return options


## Main code:
def main():
    files=[]
    parser = OptionParser(usage="Usage: %prog [options] <input1>.log <input2>.log ...")
    parser.add_option('--charge', dest="charge", action="store", default=False, help="Atomic charge at given element or atom indice.")
    parser.add_option('--nmr', dest="nmr", action="store", default=False, help="NMR shielding tensor at given element or atom indice.")
    parser.add_option('--scale', dest="scale", action="store", default=False, help="Scales isotripic values. Requires two values: first intercept then slope. (e.g., 31.8055,-1.0578)")
    parser.add_option('--mo', dest="mo", action="store_true", default=False, help="Returns the molecular orbitals for HOMO and LUMO.")
    parser.add_option('--distance', dest="distance", action="store", default=False, help="Gives bond length of given X atom in Angstrom.")
    parser.add_option('--csv', dest="csv", action="store", default=False, help="Prints out csv file of collected data.")
    parser.add_option('--verbose', dest="verbose", action="store_true", default=False, help="Request verbose print output")
    (options, args) = parser.parse_args()

    # Get input files from commandline
    if len(sys.argv) > 1:
        for elem in sys.argv[1:]:
            try:
                for file in glob(elem):
                    files.append(file)
            except IndexError: pass
    if len(files) == 0: sys.exit("    Please specify a valid input file and try again.")

    for file in files:
        if options.verbose != False: print('\no Running:', file)
        feats = disco(file, options=options)

        #### tabulate data
        csv_dict = {}
        if options.csv != False:
            csv_name = options.csv
            csv_dict['Name'] = [file]
            if options.charge != False:
                for a, b in zip(feats.Atom, feats.Chrg): csv_dict[a+' charge'] = [b]
            if options.mo != False:
                csv_dict['HOMO'] = [feats.HOMO[0]]
                csv_dict['LUMO'] = [feats.LUMO[0]]
            if options.nmr != False:
                for a, b in zip(feats.Atom, feats.NMR): csv_dict[a+' nmr'] = [b]
            if options.distance != False:
                for a, b in zip(feats.Bonded_atoms, feats.Dist): csv_dict[str(a[0]+'-'+a[1]+' length')] = [b]
        csv_df = pd.DataFrame(csv_dict)
        #### write to csv
        if options.csv != False:
            print('o Disco DFT features: \n', csv_df)
            csv_df.to_csv(f'{csv_name}_disco.csv', index=False)

if __name__ == "__main__":
    main()
