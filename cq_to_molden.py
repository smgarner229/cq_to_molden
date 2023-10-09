import sys
import copy
import numpy as np
from pyscf import gto,scf
from pyscf.tools import molden

ANGCOUNT = {0:1,
            1:3,
            2:5,
            3:7,
            4:9,
            5:11,
            6:13}

def factor(l,alpha):
  return pow(2*alpha/np.pi,3/4)*pow(pow(8*alpha,l),1/2)

class Atom:
  def __init__(self,sym_,x_,y_,z_,quantum_=False):
    if "GH" in sym_:
      sym_=sym_.replace("GH","X")
    self.__sym=sym_
    self.__x=float(x_)
    self.__y=float(y_)
    self.__z=float(z_)
    self.__quantum=quantum_

  def __str__(self):
    return "{0} {1:8f} {2:8f} {3:8f}\n".format(self.__sym,self.__x,self.__y,self.__z)

  @property
  def symbol(self):
    return self.__sym

  def set_symbol(self,symbol):
    self.__sym = symbol

  def add_index(self,num):
    self.__sym+=str(num)

  @property
  def is_quantum(self):
    return self.__quantum

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

class Molecule:
  def __init__(self,charge_,mult_):
    self.__atoms = []
    self.__qatoms = []
    self.__charge = charge_
    self.__mult = int((mult_-1)/2)
    self.__ebasis = None
    self.__pbasis = None
    self.__mo_coeff = None
    self.__pmo_coeff = None
    self.__atom_count = 0
    self.__is_NEO = False

  def add_atom(self,atom,quantum=False):
    atom.add_index(self.__atom_count)
    self.__atoms.append(atom)
    self.__qatoms.append(quantum)
    self.__atom_count += 1

  def __str__(self):
    g=""
    for i in self.__atoms:
      g += str(i)
    return g

  def set_NEO(self):
    self.__is_NEO = True

  @property
  def is_NEO(self):
    return self.__is_NEO

  @property
  def geometry(self):
    return str(self)

  @property
  def charge(self):
    return self.__charge

  @property
  def mult(self):
    return self.__mult

  @property
  def ebasis(self):
    return self.__ebasis

  def add_ebasis(self,ebas_):
    self.__ebasis = ebas_

  def e_basis_assign(self,count_):
    self.__ebasis.create(self.__atoms,count_)

  @property
  def mo_coeff(self):
    return self.__mo_coeff

  @property
  def mo_energy(self):
    return self.__mo_energy

  def set_mo_coeff(self,MOs):
    self.__mo_coeff = MOs

  def set_mo_energy(self,e):
    self.__mo_energy = e

  @property
  def pbasis(self):
    return self.__pbasis

  def add_pbasis(self,pbas_):
    self.__pbasis = pbas_

  def p_basis_assign(self,count_):
    q_atoms = [i for i in self.__atoms if i.is_quantum]
#    for atom in q_atoms:
#      atom.set_symbol(atom.symbol.replace('H','X'))
    self.__pbasis.create(q_atoms,count_)
#    for atom in q_atoms:
#      atom.set_symbol(atom.symbol.replace('X','H'))

  @property
  def pmo_coeff(self):
    return self.__pmo_coeff

  @property
  def pmo_energy(self):
    return self.__pmo_energy

  def set_pmo_coeff(self,MOs):
    self.__pmo_coeff = MOs

  def set_pmo_energy(self,e):
    self.__pmo_energy = e

  @property
  def atoms(self):
    return [i.symbol.upper() for i in self.__atoms]

  def get_unique(self):
    return set(self.atoms)

  def p_geometry(self):
#    for atom in self.__atoms:
#      if atom.is_quantum:
#        atom.set_symbol(atom.symbol.replace('H','X'))
    geom = copy.deepcopy(self.geometry)
#    for atom in self.__atoms:
#      if atom.is_quantum:
#        atom.set_symbol(atom.symbol.replace('X','H'))
    return geom

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
  
def write_electronic_MOs(mol,basename,extension="_electronic_MOs.molden"):
  filename = basename + extension
  # Make a molecule with the classical positions to plot electronic orbitals
  pyscfmol = gto.Mole(atom=mol.geometry,charge=mol.charge,spin=mol.mult)
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
 

if __name__ == "__main__":
  basename = sys.argv[1].split('.')[0]
  mol = skim_cq_out(sys.argv[1])
  write_electronic_MOs(mol,basename)
  if mol.is_NEO:
    write_protonic_MOs(mol,basename)
