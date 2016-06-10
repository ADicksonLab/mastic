from copy import copy

import networkx as nx

from mast.datastructures import *

print("Making TrackedMember")
tm = TrackedMember(idx=0, ids={'name': 'test'})
print(tm)
tm_c = copy(tm)
tms = [TrackedMember(idx=i) for i in range(5)]

print("Making TrackedList")
tl = TrackedList(tms)
print(tl)
tl_c = copy(tl)

print("Making Selections")
intsel = Selection(container=tl, sel=0, idx=0)
print(intsel)
intsel_c = copy(intsel)

selslice = slice(0,3)
slicesel = Selection(container=tl, sel=selslice)
print(slicesel)
sellist = [0,1]
listsel = Selection(container=tl, sel=sellist)
print(listsel)

print("Making SelectionList")
selections = SelectionList([intsel, slicesel, listsel])
print(selections)

print("Making AssociationType")
association_type = AssociationType()
print(association_type)

print("Making Association")
association = Association(members=selections, association=AssociationType())
print(association)

from mast.molecule import *

print("Making Atom")
atom = Atom(idx=0, coordinate=(1.0,1.0,1.0), velocity=(0.0,0.0,0.0))
print(atom)
atom_c = copy(atom)

print("Making AtomList")
atoms = AtomList([Atom(idx=i, coordinate=(float(i),float(i),float(i))) for i in range(2)])
print(atoms)
atoms_c = copy(atoms)

print("Making Bond from atoms")
bond = Bond(atoms, 0, 1, idx=0)
print(bond)
bond_c = copy(bond)

print("Making BondList from bond")
bonds = BondList([bond])
print(bonds)
bonds_c = copy(bonds)

print("Making Angle from None")
none_angle = Angle()
print(none_angle)

print("Making AngleList from none_angle")
none_angles = AngleList(members=[none_angle])

print("Making MoleculeTopology")
top = MoleculeTopology(bonds)
print(top)

print("Making bad MoleculeTopology")
atoms2 = AtomList([Atom(idx=i) for i in range(2)])
bond2 = Bond(atoms2, 0, 1)
bonds2 = BondList([bond, bond2])
try:
    top2 = MoleculeTopology(bonds2)
except ValueError:
    print("Bad MoleculeTopology ValueError caught")


print("Making Molecule")
molatoms = atoms + atoms2
print("atoms")
print(molatoms)
molbonds = BondList([Bond(molatoms, 0, 1), Bond(molatoms, 1, 2), Bond(molatoms, 2, 3)])
print("bonds")
print(bonds)
print("test MoleculeTopology")
moltop = MoleculeTopology(bonds)
print(moltop)

print("Molecule constructor")
mol = Molecule(atoms=molatoms, bonds=molbonds, angles=none_angles, idx=0)
print(mol)

print("copying Molecule, got type:", type(copy(mol)))

print("testing molecule overlap")
print("overlapping molecules:", mol.overlaps(copy(mol)))
mol2 = Molecule(atoms=AtomList([Atom(coordinate=(3.0,3.0,3.0))]))
print("non-overlapping molecules:", mol.overlaps(mol2))

print("Making MoleculeList")
mollist = MoleculeList([mol])
mols = MoleculeList([mol, mol2])

print("Importing the system module")
from mast.system import *

print("Single molecule System constructor")
sys1 = System(mol)

print("overlapping molecules")
try:
    overlap_sys = System([mol, copy(mol)])
except ValueError:
    print("Overlapping Molecules caught in constructor")

print("multi-molecule system, non-overlapping")
sys2 = System([mol, mol2])
print("get molecules")
print(sys2[:])
print(sys2.molecules)
      
print("making an empty SystemAssociation")
sys_assoc = SystemAssociation(members=None, association=None, system=None)
print("making empty SystemAssociation of sys1")
sys1_empty_assoc = SystemAssociation(members=None, association=None, system=sys1)
print("making SystemAssociation of sys2 molecules, with type AssociationType")
sys2_mol_assoc = SystemAssociation(members=sys2.molecules, association=AssociationType(), system=sys2)


print("importing the interaction.py module")
from mast.interaction import Interaction

print("Interaction constructor")
int1 = Interaction()
