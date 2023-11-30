
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

