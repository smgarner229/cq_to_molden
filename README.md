# CQ To Molden

A simple python script which takes a CQ **output** file and produces molden
output files (via PySCF) which can be viewed with 
[JMol](https://jmol.sourceforge.net)

## Usage
### CQ
This script _should_ be useable with any (including NEO) single slater
calcualtions run with the 
[Chronus Quantum](https://urania.chem.washington.edu/chronusq/chronusq_public/-/wikis/home)
package.

Importantly, the `PRINTMOS` option must be set in the CQ input file
```
[SCF]
PRINTMOS=9
```

### Python
This script is used as
```
python3 cq_to_molden.py CQ_OUTPUT_FILE.out
```
This auto-detects the geometry, basis, and molecular orbitals and dumps them to
a file named `CQ_OUTPUT_FILE_electronic_MOs.molden`.
Additionally, if a NEO calculation is run, the protonic orbitals are dumped to
`CQ_OUTPUT_FILE_protonic_MOs.molden`.


## Dependencies
- numpy
- PySCF
- h5py

### PySCF Installation Instructions
With Python3 installed and utilizing the `pip` package manager, PySCF can be 
simply installed with 
```
pip install pyscf
```
This also will install numpy and scipy as dependencies

### PySCF Virtual Environments
To isolate an individual PySCF install, you can install it and it's dependencies
in a python virtual environment with the following commands
```
python3 -m venv PYSCF_ENV
source PYSCF_ENV/bin/activate
pip install --upgrade pip
pip install pyscf
```
Note that this creates a folder called `PYSCF_ENV` in your current working 
directory, which is where PySCF, Numpy, and h5py are installed to.

To exit the virtual environment, simply use the command
```
deactivate
```

## NOTE:
There's a bug in CQ when writing the MO coefficients, in that for protonic
coefficients the first coeff is attributed (incorrectly) to a non-quantum atom,
althought this is not always the case (and I have yet to figure out why)

Practically, the two lines on 289-290 
```
bas_per_atom[1]-=1
bas_per_atom.pop(0)
```
might need commented in / out of the scrip to account for this
(the lines should be kept in if the bug appears, and commented out if it
doesn't, as is the case in the quantum h2o example)
