import numpy as np
import collections as col
from functools import reduce

from mast.selection import CoordArray, CoordArraySelection, \
    Point, IndexedSelection, SelectionDict, SelectionList

DIM_NUM_3D = 3

class Atom(Point):
    def __init__(self, coords=None, atom_array=None, array_idx=None, element=None):

        if not coords is None:
            assert coords.shape[-1] == DIM_NUM_3D, \
                "coords must have 3-dimensions, not {}".format(
                    coords.shape[-1])

        if atom_array:
            assert atom_array.shape[-1] == DIM_NUM_3D, \
                "atom_array must have 3-dimensions, not {}".format(
                    atom_array.shape[-1])


        super().__init__(coords=coords, coord_array=atom_array, array_idx=array_idx)
        self.element = element

        # stubs
        self.force_field = None


    def __repr__(self):
        return "Atom<{0}>{1}".format(str(self.element), self.coords)


class MoleculeType(object):
    def __init__(self):
        self._molecule = None
        self.features = None

class RDKitMoleculeType(MoleculeType):
    def __init__(self, rdkit_molecule):
        super().__init__()
        self.molecule = rdkit_molecule


    @property
    def atoms(self):
        return [atom for atom in self.molecule.GetAtoms()]

    @property
    def bonds(self):
        return [bond for bond in self.molecule.GetBonds()]

    def to_molecule(self, conformer_idx):
        """Construct a Molecule using a coordinates from a conformer of this
   rdchem.Mol.

        """
        assert self.molecule.GetNumConformers() > 0, \
            "{0} has no conformers".format(self)

        conformer = self.GetConformer(conformer_idx)
        atom_idxs = range(self.molecule.GetNumAtoms())
        # make the CoordArray
        positions = []
        for atom_idx in atom_idxs:
            position = conformer.GetAtomPosition(atom_idx)
            position = np.array([position.x, position.y, position.z])
        positions = reduce(lambda x,y: np.concatenate((x,y), axis=0), list(self.values()))
        coord_array = CoordArray(positions)
        # Make atoms out of the coord array
        atoms = []
        for atom_idx in atom_idxs:
            # TODO additional encapsualtion of atom data needed
            atom = Atom(atom_array=coord_array, array_idx=atom_idx)
            atoms.append(atom)

        return Molecule(atoms)


# class Molecule(SelectionDict):
#     def __init__(self, mol_input, *args):

#         if issubclass(type(mol_input), MoleculeType):
#             Molecule.type_constructor(mol_input, *args, **kwargs)
#         elif issubclass(type(mol_input), col.Sequence):
#             Molecule.atoms_constructor(mol_input, *args, **kwargs)
#         else:
#             raise TypeError("mol_input must be either a MoleculeType or a sequence of Atoms")


#     @classmethod
#     def type_constructor(self, mol_type):
#         pass


#     @classmethod
#     def atoms_constructor(self, atoms):
#         assert atoms, "atoms must exist, {}".format(atoms)
#         assert issubclass(type(atoms), col.Sequence), \
#             "atoms must be a subclass of collections.Sequence, not {}".format(
#                 type(atoms))
#         assert all([(lambda x: False if issubclass(type(x), Atom) else True)(atom)
#                     for atom in atoms]), \
#             "all elements in atoms must be a subclass of type Atom"

#         super().__init__(selection_dict=
#                          {'atoms' : atoms,
#                           'bonds' : bonds,
#                           'angles': angles})

#     @property
#     def molecule_type(self):
#         return self._molecule_type

#     @molecule_type.setter
#     def molecule_type(self, mol_type):
#         assert issubclass(type(mol_type), MoleculeType), \
#             "mol_type must be a subclass of MoleculeType, not {}".format(
#                 type(mol_type))
#         self._molecule_type = mol_type

#     @property
#     def atoms(self):
#         return self.data['atoms']

#     @property
#     def bonds(self):
#         return self.data['bonds']

#     @property
#     def angles(self):
#         return self.data['angles']



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
    atom3 = Atom(atom_array=atom_array, array_idx=0)
    print(atom3)
    print(atom3.coords)
    atoms = [atom1, atom2, Atom(np.array([0,1,0]))]
    # # make a selection of those atoms
    atomsel = IndexedSelection(atoms, [0,1])
    print(atomsel)

    # make a selection of atoms for bonds, and angle
    print("making a molecule")
    bonds = [IndexedSelection(atoms, [0,1]), IndexedSelection(atoms, [1,2])]
    angles = [IndexedSelection(atoms, [0,1,2])]
    mol = Molecule(atoms, bonds, angles)
    print(mol)

    from rdkit import Chem
    import os.path as osp
    tspo_dir = osp.expanduser("~/Dropbox/lab/tspo")
    PKA_pdb_path = osp.join(tspo_dir, "PKA.pdb")
    pka = Chem.MolFromPDBFile(PKA_pdb_path, removeHs=False)
    print(pka)

    pka_type = RDKitMoleculeType(pka)
    print(pka_type)
    pka_mol = Molecule(mol_type=pka_type)
    print(pka_mol)
    print(pka_mol.molecule_type)
