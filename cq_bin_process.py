import sys
import copy
import numpy as np
import h5py
from geom import Atom,Molecule

# Eventually abstract some of this away into a base class, deal with 
# redundacies for now
from cq_output_process import process_molecule

def read_bin_mos(h5file,mol):
  mo1 = np.asarray(h5file['SCF/MO1'])
  try:
    mo2 = np.asarray(h5file['SCF/MO2'])
  except:
    # Probably print some status messages for finding or not finding MO2
    pass
  # Transpose since CQ and PySCF have different standards
  mol.set_mo_coeff(mo1.T)
  # MO Energies (Hacky again for now)
  mol.set_mo_energy(range(mo1.shape[0]))
  return

def geom_helper_temp(cqinputfile):
  """
  Hacky for now, but need to feed the basis and the geometry into the PySCF
  object, and since these are not saved in the bin file, we'll temporarily
  use the input file to parse for them
  """
  with open(cqinputfile,'r') as f1:
    for line in f1:
      if "[Molecule]" in line:
        mol = process_molecule(f1)
      elif "[BASIS]" in line:
        basis_line = next(f1)
        basis=basis_line.split()[-1]
        mol.set_basis_name(basis)
  return mol

def skim_cq_bin(binfile,inputfilename=None):
  """
  TODO: Much of this can just be overwritten Gerardo's more intelligent work
  to be able to handle 2C/4C outputs
  """
  if inputfilename is None:
    inputfilename = binfile.split('.')[0]+".inp"
  mol = geom_helper_temp(inputfilename)
  # binfile
  bf = h5py.File(binfile,'r')
  read_bin_mos(bf,mol) 
  return mol

if __name__=="__main__":
  from molden import write_electronic_MOs,write_protonic_MOs
  basename = sys.argv[1].split('.')[0]
  mol = skim_cq_bin(sys.argv[1])
  write_electronic_MOs(mol,basename)

