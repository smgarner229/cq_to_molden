import sys
import copy
import numpy as np
from pyscf import gto
from geom import Atom,Molecule

ANGCOUNT = {0:1,
            1:3,
            2:5,
            3:7,
            4:9,
            5:11,
            6:13}

def factor(l,alpha):
  return pow(2*alpha/np.pi,3/4)*pow(pow(8*alpha,l),1/2)

class basis_info:
  def __init__(self,nbasis_,nprim_,nshell_,max_prim_,max_l_):
    self.__nbasis   = nbasis_
    self.__nprim    = nprim_
    self.__nshell   = nshell_
    self.__max_prim = max_prim_
    self.__max_l    = max_l_
  
  @property
  def nbasis(self):
    return self.__nbasis

class Basis:
  def __init__(self,bas_,bas_info_):
    self.__bas      = bas_
    self.__bas_info = bas_info_
    self.__func_count = None

  @property
  def basis(self):
    return self.__bas

  @property
  def nbasis(self):
    return self.__bas_info.nbasis

  def create(self,atoms,count_):
    """
    Parses the stored basis and assign blocks of basis functions to their
    respective atoms
    """
    basis_dict={}
    self.__func_count = count_
    last = 0
    for atom,count in zip(atoms,count_):
      n = 0
      for i,func in enumerate(self.__bas[last:]):
        n+=ANGCOUNT[func[0]]
        if n == count:
          basis_dict[atom.symbol] = self.__bas[last:last+i+1]
          last+=i+1
          break
    self.__bas = basis_dict        


def process_molecule(iterator):
  m = next(iterator)
  # Get the charge and multiplicity
  while "geom" not in m:
    if "charge" in m:
      charge = int(m.split()[-1])
    if "mult" in m:
      mult = int(m.split()[-1])
    m = next(iterator)
  mol = Molecule(charge,mult)
  m = next(iterator)
  # Read the geometry in
  while len(m.split()) > 2:
    q = False
    m_s = m.split()
    if 'Q' in m_s[-1]:
      q = True
      m_s = m_s[:-1]
    mol.add_atom(Atom(*m_s,q))
    m = next(iterator)
  print("")
  print("Detected geometry: ")
  print(mol)
  print("")
  return mol
  return gto.Mole(atom=g,charge=charge,mult=mult)

def fix_normalization(basis):
  new_basis = []
  for i,func in enumerate(basis):
    this_func = []
    ang = func[0]
    this_func.append(ang)
    for exp,cont in func[1:]:
      this_func.append((exp,cont/factor(ang,exp)))
    new_basis.append(this_func)
  return new_basis

def process_basis_detail(iterator,nbasis):
  basis = []
  cur = []
  m = next(iterator)
  while "=" not in m:
    # New basis function
    if len(m.split())==4 and "#" not in m:
      # Add angular momentum indicator
      cur.append(int(m.split()[1]))
      cur.append((float(m.split()[2]),float(m.split()[3])))
    # Current basis function
    elif len(m.split())==2 and "Shell" not in m:
      cur.append((float(m.split()[0]),float(m.split()[1])))
    # Add to the list of basis functions
    elif len(m.split())==0 and len(basis)<nbasis and len(cur)>0:
      basis.append(cur)
      cur=[]
    m = next(iterator) 
  basis=fix_normalization(basis)
  return basis

def process_basis(iterator):
  [next(iterator) for _ in range(2)]
  bas_info = basis_info(*[int(next(iterator).split()[-1]) for _ in range(5)])
  bas=process_basis_detail(iterator,bas_info.nbasis)
  return Basis(bas,bas_info)

def read_mos(iterator,mol,basis,NEO=False):
  nbasis = basis.nbasis
  num_blocks = nbasis // 4
  [next(iterator) for _ in range(2)]
  m = next(iterator)
  OrbEs=[]
  Orbs=[]
  bas_per_atom = []
  bas_assigned = False
  MOs=None
  # Big positive number
  ncol = 100000000
  while "--------" not in m:
    if "EigV" in m:
      tmp=[float(i) for i in m.split()[2:]]
      if len(Orbs)==ncol*nbasis:
        if MOs is None:
          MOs=np.reshape(np.asarray(Orbs),(nbasis,ncol))
          # Because there is a bug
          if NEO:
            pass
            #bas_per_atom[1]-=1
            #bas_per_atom.pop(0)
          bas_assigned = True
        else:
          MOs = np.hstack((MOs,np.reshape(np.asarray(Orbs),(nbasis,ncol))))
        Orbs=[]
      OrbEs+=tmp
      ncol = len(tmp)
    elif len(m.split())>ncol+1:
      Orbs+=[float(i) for i in m.split()[-ncol:]]
    # Detect the basis assignments
    if len(m.split())>ncol+3 and not bas_assigned:
      if(m.split()[2].count('-')):
        bas_per_atom.append(len(Orbs)//ncol)
    m=next(iterator)
  bas_per_atom.append(len(Orbs)//ncol+1)
  MOs = np.hstack((MOs,np.reshape(np.asarray(Orbs),(nbasis,ncol))))
  bas_func = [j-bas_per_atom[i-1] for i,j in enumerate(bas_per_atom)]
  bas_func.pop(0)
  if NEO:
    mol.p_basis_assign(bas_func)
    mol.set_pmo_coeff(MOs)
    mol.set_pmo_energy(OrbEs)
  else:
    mol.e_basis_assign(bas_func)
    mol.set_mo_coeff(MOs)
    mol.set_mo_energy(OrbEs)
  return 

def skim_cq_out(filename):
  """
  Skims a CQ output file for [TODO: Fill in this docstring]
  """

  basename = filename.split('.')[0]
  
  NEO  = False

  with open(filename,'r') as f1:
    for line in f1:
      if "[Molecule]" in line:
        mol = process_molecule(f1)
      elif "Basis Set Information" in line:
        if mol.ebasis is None:
          mol.add_ebasis(process_basis(f1))
        else:
          mol.add_pbasis(process_basis(f1))
          mol.set_NEO()
      elif "Canonical Molecular Orbital Coefficients (Alpha)" in line:
        # First read the proton MOs if NEO Calculation
        if mol.is_NEO and mol.pmo_coeff is None:
          read_mos(f1,mol,mol.pbasis,True)
        else:
          read_mos(f1,mol,mol.ebasis)
  return mol
 
if __name__ == "__main__":
  from molden import write_electronic_MOs,write_protonic_MOs
  basename = sys.argv[1].split('.')[0]
  mol = skim_cq_out(sys.argv[1])
  write_electronic_MOs(mol,basename)
  if mol.is_NEO:
    write_protonic_MOs(mol,basename)
