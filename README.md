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

Bond distance utilizes molecular structures from: computed outputs (*.log) or xyz coordinates (*.xyz)
* `--distance [atom1, atom2]` - Bond length(s) of bonded atoms in Angstrom

## Other options
* `--verbose` - Prints outputs on the terminal.
* `--csv` - Saves outputs for single molecule in csv file.


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

2. Atomic charge for the ipso-Carbon atom bonded to Cl with atom number, 3:
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

4. NMR tensor value for all C atoms:
```
>>>python DISCO.py mol1-Cl_NMR.log --nmr C --verbose

   o Running: mol1-Cl_NMR.log
     Atomic NMR Values: ['1C', '2C', '3C', '4C', '5C'] = [32.1843, 57.284, 24.4286, 57.2838, 32.1825] 
```

5. NMR chemical shift values for all C atoms using the scaling factors obtained from [Cheshire NMR](http://cheshirenmr.info/ScalingFactors.htm) developed by the Tantillo group with mPW1PW91/6-311+G(2d,p) (giao, gas phase):
```
>>>python DISCO.py mol1-Cl_NMR.log --nmr C --scale '185.6667,-1.0303' --verbose

   o Running: mol1-Cl_NMR.log
     Atomic NMR Values: ['1C', '2C', '3C', '4C', '5C'] = [152.5072, 126.647, 160.4979, 126.6472, 152.5091]
```

4. HOMO and LUMO molecular orbital values:
```
>>>python DISCO.py mol1-Cl_NBO.log --mo --verbose

   o Running: mol1-Cl_NBO.log
     HOMO and LUMO = [-0.3698] and [-0.017]
```

5. Bond distance for single C-Cl bond occurance:
```
>>>python DISCO.py mol1-Cl_NBO.log --distance 'C,Cl' --verbose

   o Running: mol1-Cl_NBO.log
     Bond Distance: [['3C', '11Cl']] = [1.7322] 
```

6. Bond distance for multiple N-H bond occurance:
```
>>>python DISCO.py mol3-N_NMR.log mol4-2N_NMR.log --distance 'N,H' --verbose
   o Running: mol3-N_NMR.log
     Bond Distance: [['1N', '11H']] = [1.00521]
     
   o Running: mol4-2N_NMR.log
     Bond Distance: [['3N', '8H'], ['6N', '15H']] = [1.00876, 1.01189]
```


7. Bond distance on .xyz file:
```
>>>python DISCO.py mol2-Br_NBO.xyz --distance 'Br,C' --verbose

   o Running: mol2-Br_NBO.xyz
     Bond Distance: [['1Br', '2C']] = [1.90708] 
```

