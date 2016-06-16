import numpy as np

from mast.selection import CoordArray, Point, IndexedSelection

DIM_NUM_3D = 3

class AtomArray(CoordArray):
    def __init__(self, coord_array=None):

        super().__init__(coord_array)

class Atom(Point):
    def __init__(self, coords, atom_array=None, element=None):
        assert coords.shape[-1] == DIM_NUM_3D, \
            "coords must have 3-dimensions, not {}".format(
                coords.shape[-1])

        if atom_array:
            assert atom_array.shape[-1] == DIM_NUM_3D, \
                "atom_array must have 3-dimensions, not {}".format(
                    atom_array.shape[-1])

        super().__init__(coords, coord_array=atom_array)
        self.element = element

    def __repr__(self):
        return "Atom<{0}>{1}".format(str(self.element), self.coords)

# class Molecule(col.UserDict):
#     def __init__(self, atoms: ty.Sequence[Atom],
#                        bonds: ty.Sequence[Bond],
#                        angles: ty.Sequence[Angle],
#                        idx: int=None) -> None:
#         super().__init__()
#         self['atoms'] = atoms
#         self['bonds'] = bonds
#         self['angles'] = angles
#         self.atoms = self['atoms']


if __name__ == "__main__":

    print("making an CoordArray for atoms")
    array = np.array([[0,0,0], [0,0,1], [1,0,0]])
    atom_array = CoordArray(array)
    print(atom_array)

    print("making atoms")
    atom1 = Atom(np.array([5,5,5]))
    print(atom1)
    print(atom1.coords)
    atom2 = Atom(np.array([6,6,6]), element='C', atom_array=atom_array)
    print(atom2)
    print(atom2.coords)
    atoms = [atom1, atom2, Atom(np.array([0,1,0]))]
    # # make a selection of those atoms
    atomsel = IndexedSelection(atoms, [0,1])
    print(atomsel)

    # mol = Molecule(atoms, bonds, angles, idx=0)
    # print(mol)
