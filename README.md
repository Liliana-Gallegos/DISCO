![disco](DISCO_logo.png)

# DISCO
DISCO - distributing computed outputs.

Python script to parse through Gaussian NBO and GIAO outputs for atomistic and molecular properties.

Collects NBO atomic charges, NMR tensor values and/or chemical NMR shifts, HOMO, LUMO, and Bond distance values.

## Feature options
The following requires computed G16 outputs (*log)
* `--mo` - Highest occupied molecular orbital (HOMO) and lowest unoccupied molecular orbital (LUMO)
* `--charge [atom]` - Atomic natural charge
* `--nmr [atom]`    - NMR tensor values
* `--nmr [atom]` WITH `--scale [intercept, slope]` - Chemical NMR shift value in ppm

Bond distance utilizes molecular structures from: computed outputs (.log) or xyz coordinates (.xyz)
* `--distance [atom1, atom2]` - Bond length(s) of bonded atoms in Angstrom

## Other options
* `--verbose` - Prints outputs on the terminal.
* `--csv [filename]` - Saves outputs for single molecule in csv file. Requires a filename.

## Dependencies
Program requires pandas and numpy.


## Examples on Jupyter Notebook
See the Collecting_DFT-features.ipynb within the [Example_jupyter-notebook](https://github.com/Liliana-Gallegos/DISCO/tree/master/Example_jupyter-notebook) directory. 

To run in a Jupyter Notebook import: `import DISCO as cd` 
```
    import DISCO as cd
    
    #Create DISCO object
    mol = cd.disco(file, distance='atom1,atom2', charge=atom, nmr=atom, scale='intercept,slope', mo=True, verbose=True)
    
    #Grab DFT Parameters
    charge = mol.Chrg
    nmr    = mol.NMR
    distance = mol.Dist
    atom_index = mol.Atom
    homo = mol.HOMO
    lumo = mol.LUMO
   
```

## Examples on terminal command line
Examples for collecting DFT features using the terminal command line. Requires `--verbose` to display outputs on terminal.
1. Atomic charges for all Cl atoms:
```
>>>python DISCO.py mol1-Cl_NBO.log --charge Cl --verbose

   o Running: mol1-Cl_NBO.log
     Atomic Charge: ['11Cl', '13Cl', '14Cl'] = [0.03797, -0.81551, -0.81551]
```

2. Atomic charge for a specific indexed 3C atom:
```
>>>python DISCO.py mol1-Cl_NBO.log --charge 3 --verbose

   o Running: mol1-Cl_NBO.log
     Atomic Charge: ['3C'] = [0.0343]
```

3. NMR tensor value for N atom:
```
>>>python DISCO.py mol1-Cl_NMR.log --nmr N --verbose

   o Running: mol1-Cl_NMR.log
     Atomic NMR Values: ['6N'] = [-6.089]
```

4. NMR chemical shift values for all C atoms using the scaling factors obtained from [Cheshire NMR](http://cheshirenmr.info/ScalingFactors.htm) developed by the Tantillo group with mPW1PW91/6-311+G(2d,p) (giao, gas phase):
```
>>>python DISCO.py mol1-Cl_NMR.log --nmr C --scale '185.6667,-1.0303' --verbose

   o Running: mol1-Cl_NMR.log
     Atomic NMR Values: ['1C', '2C', '3C', '4C', '5C'] = [148.9686, 124.6071, 156.4963, 124.6073, 148.9704]
```

5. HOMO and LUMO molecular orbital values:
```
>>>python DISCO.py mol1-Cl_NBO.log --mo --verbose

   o Running: mol1-Cl_NBO.log
     HOMO and LUMO = [-0.3698] and [-0.017]
```

6. Bond distance for single C-Cl bond occurance:
```
>>>python DISCO.py mol1-Cl_NBO.log --distance 'C,Cl' --verbose

   o Running: mol1-Cl_NBO.log
     Bond Distance: [['3C', '11Cl']] = [1.7322] 
```

7. Bond distance for multiple N-H bond occurance: 
* Can also grab bond distances from XYZ files (*.xyz)
```
>>>python DISCO.py mol3-N_NMR.log mol4-2N_NMR.log --distance 'N,H' --verbose
   o Running: mol3-N_NMR.log
     Bond Distance: [['1N', '11H']] = [1.00521]
     
   o Running: mol4-2N_NMR.log
     Bond Distance: [['3N', '8H'], ['6N', '15H']] = [1.00876, 1.01189]


>>>python DISCO.py mol2-Br_NBO.xyz --distance 'Br,C' --verbose

   o Running: mol2-Br_NBO.xyz
     Bond Distance: [['1Br', '2C']] = [1.90708] 
```

8. Tabulate all data into dataframe (csv) by adding `--csv [filename]`.

* Example 1 - Specified by nitrogen atom `--charge N`: Tabulates by the first N atom.ex1_data.csv gets saved in the working folder.
```
>>>python DISCO.py mol*N_*NBO.log --charge N --verbose --csv ex1_data 

o Running: mol3-N_NBO.log
  Atomic Charge: ['1N', '16N'] = [-0.60228, -0.53661]

o Running: mol4-2N_NBO.log
  Atomic Charge: ['3N', '6N'] = [-0.64365, -0.68712]
o Disco DFT features: 
               Name N atom index  N charge
0   mol3-N_NBO.log           1N  -0.60228
1  mol4-2N_NBO.log           3N  -0.64365
```

* Example 2 - Specified by indexed atom `--charge 3`: Tabulates by the specified indexed atom. ex2_data.csv gets saved in the working folder.
```
>>>python DISCO.py mol1-Cl_NBO.log mol2-Br_NBO.log --charge 3 --verbose --csv ex2_data 

o Running: mol1-Cl_NBO.log
  Atomic Charge: ['3C'] = [0.0343]

o Running: mol2-Br_NBO.log
  Atomic Charge: ['3C'] = [-0.24698]
o Disco DFT features: 
               Name 3 atom index  3 charge
0  mol1-Cl_NBO.log           3C   0.03430
1  mol2-Br_NBO.log           3C  -0.24698
```

* Example 3 - Specified by indexed atom `--charge 2,3,4,5,Br`: Tabulates by multiple indexed specied by atom type or index. ex3_data.csv gets saved in the working folder.

<img src="https://github.com/Liliana-Gallegos/DISCO/blob/master/Example_jupyter-notebook/Ex3_mols.png" style="margin:auto" width="500"/>

```
>>>python DISCO.py mol6-Br_NBO.log mol2-Br_NBO.log --charge 2,3,4,5,Br --csv ex3_data

o Disco DFT features: 
               Name 2 atom index  2 charge 3 atom index  3 charge 4 atom index  4 charge 5 atom index  5 charge Br atom index  Br charge
0  mol6-Br_NBO.log           2C  -0.08636           3C  -0.23470           4C  -0.16334           5C  -0.06443           1Br    0.04563
1  mol2-Br_NBO.log           2C  -0.07746           3C  -0.24698           4C  -0.18143           5C  -0.21489           1Br    0.03987
```
