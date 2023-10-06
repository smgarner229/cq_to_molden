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
