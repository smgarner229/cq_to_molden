from pyscf import gto
from pyscf.tools import molden

def write_electronic_MOs(mol,basename,extension="_electronic_MOs.molden"):
  filename = basename + extension
  # Make a molecule with the classical positions to plot electronic orbitals
  pyscfmol = gto.Mole(atom=mol.geometry,charge=mol.charge,spin=mol.mult)
  if mol.basis_name is not None:
    pyscfmol.basis = mol.basis_name
  else:
    pyscfmol.basis = mol.ebasis.basis
  pyscfmol.build()
  with open(filename,'w') as f1:
    molden.header(pyscfmol,f1)
    molden.orbital_coeff(pyscfmol,f1,mol.mo_coeff,ene=mol.mo_energy)

def write_protonic_MOs(mol,basename,extension="_protonic_MOs.molden"):
  # Make a dummy molecule to plot the protonic orbitals
  # probably need to pad the protonic molecular orbitals with 0's and ghost
  # basis centers on the rest of the atoms
  filename = basename+extension
  pyscfpmol = gto.Mole(atom=mol.p_geometry(),charge=mol.charge,spin=mol.mult)
  pyscfpmol.basis = mol.pbasis.basis
  pyscfpmol.build()
  with open(filename,'w') as f1:
    molden.header(pyscfpmol,f1)
    molden.orbital_coeff(pyscfpmol,f1,mol.pmo_coeff,ene=mol.pmo_energy)

