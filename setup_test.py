from copy import copy

from mast.TrackedStructures import *

print("Making TrackedMember")
tm = TrackedMember(idx=0, ids={'name': 'test'})
tm
tm_c = copy(tm)
tms = [TrackedMember(idx=i) for i in range(5)]

print("Making TrackedList")
tl = TrackedList(tms)
tm
tl_c = copy(tl)

print("Making Selections")
intsel = Selection(container=tl, sel=0)
intsel
intsel_c = copy(intsel)

selslice = slice(0,3)
slicesel = Selection(container=tl, sel=selslice)

sellist = [0,1]
listsel = Selection(container=tl, sel=sellist)


from mast.Molecule import *

print("Making Atom")
atom = Atom(idx=0)
atom
print(atom)
atom_c = copy(atom)

print("Making AtomList")
atoms = AtomList([Atom(idx=i) for i in range(5)])
atoms
print(atoms)
atoms_c = copy(atoms)

print("Making Bond from atoms")
bond = Bond(atoms, 0, 1)
bond
print(bond)
bond_c = copy(bond)
