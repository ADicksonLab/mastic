from copy import copy

import networkx as nx

from mast.TrackedStructures import *

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
intsel = Selection(container=tl, sel=0)
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


from mast.Molecule import *

print("Making Atom")
atom = Atom(idx=0)
print(atom)
atom_c = copy(atom)

print("Making AtomList")
atoms = AtomList([Atom(idx=i) for i in range(2)])
print(atoms)
atoms_c = copy(atoms)

print("Making Bond from atoms")
bond = Bond(atoms, 0, 1)
print(bond)
bond_c = copy(bond)

print("Making BondList from bond")
bonds = BondList([bond])
print(bonds)
bonds_c = copy(bonds)

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
mol = Molecule(atoms=molatoms, bonds=molbonds)
