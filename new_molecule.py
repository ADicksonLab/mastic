import numpy as np
import collections as col
from functools import reduce
import os.path as osp

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

    # @classmethod
    # def create_molecule_type(cls, mol, mol_type=None):
    #     cls_name = str(super().__self_class__).strip("'>").split('.')[-1]
    #     return type(mol_type, (cls_name,), cls(mol).__dict__)


class RDKitMoleculeType(MoleculeType):
    def __init__(self, rdkit_molecule, mol_type=None):
        super().__init__()
        self.molecule = rdkit_molecule
        # TODO convert this attribute to a class factory
        assert isinstance(mol_type, str)
        self.molecule_type = mol_type
        self._features = None

    @property
    def atoms(self):
        return [atom for atom in self.molecule.GetAtoms()]

    @property
    def bonds(self):
        return [bond for bond in self.molecule.GetBonds()]

    @property
    def features(self):
        return self._features

    def find_features(self, fdef="BaseFeatures.fdef"):
        """Uses a feature definition (fdef) database to to find features in
the molecule.

        """
        from rdkit import RDConfig
        from rdkit.Chem import ChemicalFeatures

        assert isinstance(fdef, str)
        fdef_path = osp.join(RDConfig.RDDataDir, fdef)
        feature_factory = ChemicalFeatures.BuildFeatureFactory(fdef_path)
        factory_features = feature_factory.GetFeaturesForMol(self.molecule)
        features = {}
        for feature in factory_features:
            feature_info = {'family' : feature.GetFamily(),
                            'type' : feature.GetType(),
                            'atom_ids' : feature.GetAtomIds(),
                            'position' : feature.GetPos()}
            features[feature.GetId()] = feature_info

        return features

    def to_molecule(self, conformer_idx):
        """Construct a Molecule using a coordinates from a conformer of this
   rdchem.Mol.

        """
        assert self.molecule.GetNumConformers() > 0, \
            "{0} has no conformers".format(self)

        conformer = self.molecule.GetConformer(conformer_idx)
        atom_idxs = range(self.molecule.GetNumAtoms())
        # make the CoordArray
        positions = []
        for atom_idx in atom_idxs:
            position = conformer.GetAtomPosition(atom_idx)
            position = np.array([position.x, position.y, position.z])
            positions.append(position)
        positions = np.array(positions)
        coord_array = CoordArray(positions)
        # Make atoms out of the coord array
        atoms = []
        for atom_idx in atom_idxs:
            # TODO additional encapsualtion of atom data needed
            atom = Atom(atom_array=coord_array, array_idx=atom_idx)
            atoms.append(atom)

        # TODO handle bonds
        bonds = list(self.molecule.GetBonds())

        # TODO handle and create angles
        angles = None

        return Molecule(atoms, bonds, angles, mol_type=self.molecule_type, external_mol_rep=(RDKitMoleculeType, self))


class Molecule(SelectionDict):
    def __init__(self, mol_input, *args, **kwargs):

        if 'mol_type' not in kwargs.keys():
            mol_type = None
        else:
            mol_type = kwargs.pop('mol_type')

        if 'external_mol_rep' not in kwargs.keys():
            external_mol_rep = None
        else:
            assert isinstance(kwargs['external_mol_rep'], tuple), \
                "An external_mol_rep must be a tuple (external_type, external_mol), not {}".format(
                    kwargs['external_mol_rep'])
            external_mol_rep = kwargs.pop('external_mol_rep')

        if issubclass(type(mol_input), MoleculeType):
            molecule_dict = Molecule.type_constructor(mol_input, *args, **kwargs)
        elif issubclass(type(mol_input), col.Sequence):
            print("atoms constructor")
            molecule_dict = Molecule.atoms_constructor(mol_input, *args, **kwargs)
        else:
            raise TypeError("mol_input must be either a MoleculeType or a sequence of Atoms")

        super().__init__(selection_dict=molecule_dict)
        self._features = None
        self._molecule_type = mol_type
        # a dictionary of molecular representations from other libraries
        if external_mol_rep:
            self._external_mol_reps = {external_mol_rep[0] : external_mol_rep[1]}

    @classmethod
    def type_constructor(cls, mol_type):
        raise NotImplementedError


    @classmethod
    def atoms_constructor(cls, atoms, bonds, angles):
        assert atoms, "atoms must exist, {}".format(atoms)
        assert issubclass(type(atoms), col.Sequence), \
            "atoms must be a subclass of collections.Sequence, not {}".format(
                type(atoms))
        assert all([(lambda x: True if issubclass(type(x), Atom) else False)(atom)
                    for atom in atoms]), \
            "all elements in atoms must be a subclass of type Atom"
        molecule_dict= {'atoms' : atoms, 'bonds' : bonds, 'angles': angles}
        return molecule_dict

    @property
    def molecule_type(self):
        return self._molecule_type

    @molecule_type.setter
    def molecule_type(self, mol_type):
        assert issubclass(type(mol_type), MoleculeType), \
            "mol_type must be a subclass of MoleculeType, not {}".format(
                type(mol_type))
        self._molecule_type = mol_type

    @property
    def atoms(self):
        return self.data['atoms']

    @property
    def bonds(self):
        return self.data['bonds']

    @property
    def angles(self):
        return self.data['angles']

    @property
    def external_mol_reps(self):
        return self._external_mol_reps

    @property
    def features(self):
        return self._features

    # TODO should be in the Type class for the molecule not the
    # molecule itself
    def find_features(self):
        features = self._external_mol_reps[RDKitMoleculeType].find_features()
        self._features = features


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


    from rdkit import Chem
    import os.path as osp
    tspo_dir = osp.expanduser("~/Dropbox/lab/tspo")
    PKA_pdb_path = osp.join(tspo_dir, "PKA.pdb")
    pka = Chem.MolFromPDBFile(PKA_pdb_path, removeHs=False)
    print(pka)


    # PKA = RDKitMoleculeType.create_molecule_type(pka, mol_type='PKA')

    # PKA = RDKitMoleculeType.create_molecule_type(
    pka_type = RDKitMoleculeType(pka, mol_type="PKA")
    print(pka_type)

    # make a selection of atoms for bonds, and angle
    print("making a molecule")
    bonds = [IndexedSelection(atoms, [0,1]), IndexedSelection(atoms, [1,2])]
    angles = [IndexedSelection(atoms, [0,1,2])]
    mol = Molecule(atoms, bonds, angles)
    print(mol)

    print("Making a mast.Molecule from the RDKitMoleculeType")
    pka_mol = pka_type.to_molecule(0)
    # pka_mol = Molecule(mol_type=pka_type, coords=pka_coords)
    print(pka_mol)
    print(pka_mol.molecule_type)

    print("Getting features for the molecule from it's RDKit representation")
    pka_mol_rdmol = pka_mol.external_mol_reps[RDKitMoleculeType]
    pka_rdfeatures = pka_mol_rdmol.find_features()
    print(pka_rdfeatures)

    print("finding features using amst.molecule method")
    pka_mol.find_features()
    print(pka_mol.features)
