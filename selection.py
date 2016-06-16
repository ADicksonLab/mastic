import numpy as np
import collections as col
import typing as ty
import collections.abc as colabc
from copy import copy


selection_registry = {}
sel_reg_counter = 0

class SelectionMember(object):
    def __init__(self, member):
        super().__init__()
        self.member = member
        # {selection_registry_id : id_in_selection}
        self.registry = {}

    def __repr__(self):
        return self.member.__repr__()

    def get_selections(self):
        global selection_registry
        return [selection_registry[sel_id] for sel_id in self.registry.keys()]

class CoordArray(SelectionMember, np.ndarray):
    def __init__(self, array):
        super().__init__()

Selected = ty.TypeVar('Selected')
SelectionIDs = ty.TypeVar('SelectionIDs')
class GenericSelection(SelectionMember, col.UserDict, ty.Mapping[SelectionIDs, Selected]):
    def __init__(self, container: ty.Container):
        super().__init__(self)

        assert '__getitem__' in dir(container), \
            "container must implement `__getitem__`, {} does not".format(
                container)
        assert container, "container must have at least one SelectionMember element"
        assert issubclass(type(container[0]), SelectionMember), \
            "container members must be a subclass of SelectionMember, not {}".format(
                type(container[0]))

        self.container = container
        self.sel_ids = SelectionIDs

    def __repr__(self):
        return "{0}[{1}]".format(self.container, self.sel_ids)

class IndexedSelection(GenericSelection[int, Selected]):
    def __init__(self, container: ty.Sequence, sel: ty.Sequence[int]):
        super().__init__(container)
        # register this selection
        global sel_reg_counter
        sel_reg_counter += 1
        self.sel_reg_id = sel_reg_counter
        global selection_registry
        selection_registry[self.sel_reg_id] = self

        # make the selections from container
        self.sel_ids = sel
        for sel_idx in sel:
            self[sel_idx] = container[sel_idx]
            # set this selection in the SelectionMember registry
            self[sel_idx].registry[self.sel_reg_id] = sel_idx

    def __repr__(self):
        return str(dict(self))

class CoordArraySelection(GenericSelection[int, np.ndarray]):
    def __init__(self, array: CoordArray, sel: ty.Sequence[int]):
        super().__init__(array)
        for sel_idx in sel:
            # slices like this are views into the array
            self[sel_idx] = self.container[sel_idx:(sel_idx+1)]

class Point(CoordArraySelection):
    def __init__(self, coords, coord_array=None):
        pass


class Atom(object):
    def __init__(self, element=None, coords=None):
        self.element = element
        self.coords = coords

    def __repr__(self):
        return "Atom<{0}>{1}".format(str(self.element), self.coords)

AtomList = ty.List[Atom]

class Bond(object):
    def __init__(self, atom1, atom2):
        self.atom1 = atom1
        self.atom2 = atom2

BondList = ty.List[Bond]

class Angle(object):
    def __init__(self, atom1, atom2, atom3):
        self.atom1 = atom1
        self.atom2 = atom2
        self.atom3 = atom3

class AtomSelection(IndexedSelection[int, Atom]):
    def __init__(self, atoms: ty.Sequence[Atom], sel: ty.Sequence[int]):
        super().__init__(atoms, sel)

class BondSelection(IndexedSelection[int, Bond]):
    def __init__(self, bonds: ty.Sequence[Bond], sel: ty.Sequence[int]):
        super().__init__(bonds, sel)


class Molecule(col.UserDict):
    def __init__(self, atoms: ty.Sequence[Atom],
                       bonds: ty.Sequence[Bond],
                       angles: ty.Sequence[Angle],
                       idx: int=None) -> None:
        super().__init__()
        self['atoms'] = atoms
        self['bonds'] = bonds
        self['angles'] = angles
        self.atoms = self['atoms']


if __name__ == "__main__":
    # test GenericSelection
    gensel = GenericSelection([SelectionMember(None)])
    print(gensel)

    # test SelectionMember
    string_selmember = SelectionMember('a')
    print(string_selmember)

    # test IndexedSelection
    selmembers = [SelectionMember(i) for i in [0,1,2]]
    print("selmembers[0] is a part of these selections", selmembers[0].get_selections())
    idxsel = IndexedSelection(selmembers, [0,2])
    print("idxsel", idxsel)
    print("selmembers[0] is a part of these selections", selmembers[0].get_selections())
    idxsel2 = IndexedSelection(selmembers, [0,1])
    print("selmembers[0] is a part of these selections", selmembers[0].get_selections())

    # test a selection from a list of SelectionMembers as a new type
    strings = [SelectionMember('a'), SelectionMember('b'), SelectionMember('c')]
    class StrSelection(IndexedSelection[int, str]):
        def __init__(self, strings: ty.Sequence[str], sel: ty.Sequence[int]):
            super().__init__(strings, sel)

    stringsel = StrSelection(strings, [0])
    print(stringsel)

    # make some atoms
    atoms = [Atom('C', (0,1,0)), Atom('B'), Atom('C')]
    print(atoms[0])
    # make a selection of those atoms
    atomsel = AtomSelection(atoms, [0])
    print(atomsel)
    bonds = [Bond(atoms[0], atoms[1]), Bond(atoms[1], atoms[2])]
    print(bonds)
    angles = [Angle(atoms[0], atoms[1], atoms[2])]
    print(angles)
    mol = Molecule(atoms, bonds, angles, idx=0)
    print(mol)
