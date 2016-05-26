from copy import copy

from mast.TrackedStructures import *

tm = TrackedMember(idx=0, ids={'name': 'test'})

tms = [TrackedMember(idx=i) for i in range(5)]

tl = TrackedList(tms)

intsel = Selection(container=tl, sel=0)

selslice = slice(0,3)
slicesel = Selection(container=tl, sel=selslice)

sellist = [0,1]
listsel = Selection(container=tl, sel=sellist)



from mast.Molecule import *

atom = Atom(idx=0)

atoms = AtomList([Atom(idx=i) for i in range(5)])

atomsel = Selection(container=atoms, sel=0)

selbond = Bond(atomsel)
