import numpy as np
import collections as col
from functools import reduce
import os.path as osp
from itertools import product
import math

import mast
import mast.interactions
from mast.selection import CoordArray, CoordArraySelection, \
    Point, IndexedSelection, SelectionDict, SelectionList, \
    SELECTION_REGISTRY, sel_reg_counter, \
    SelectionType, SelectionTypeLibrary

from mast.interactions import InteractionType

__all__ = ['Atom', 'AtomTypeLibrary', 'AtomType', 'Bond', 'MoleculeType', 'MoleculeTypeLibrary',
           'RDKitMoleculeType', 'Molecule', ]

DIM_NUM_3D = 3

class Atom(Point):
    def __init__(self, coords=None, atom_array=None, array_idx=None, atom_type=None):

        if coords is None:
            coords = np.array([np.nan, np.nan, np.nan])
        else:
            assert coords.shape[-1] == DIM_NUM_3D, \
                "coords must have 3-dimensions, not {}".format(
                    coords.shape[-1])

        if atom_array:
            assert atom_array.shape[-1] == DIM_NUM_3D, \
                "atom_array must have 3-dimensions, not {}".format(
                    atom_array.shape[-1])


        super().__init__(coords=coords, coord_array=atom_array, array_idx=array_idx)

        if atom_type is None:
            atom_type = AtomType({})
        self._atom_type = atom_type
        self._in_molecule = False
        self._in_bond = False


    def __repr__(self):
        return "Atom{0}{1}".format(self.atom_type, self.coords)

    @property
    def atom_type(self):
        return self._atom_type

    @property
    def molecule(self):
        if self._in_molecule is False:
            return None
        else:
            molecule = next((sel for sel in self.get_selections()
                             if isinstance(sel, Molecule)),
                            None)
            assert molecule
            return molecule

    @property
    def bonds(self):
        if self._in_bond is False:
            return None
        else:
            # TODO really weird bug where isinstance(Molecule(), Bond)
            # evaluates to True!!!! Need to find out why that is. For
            # now just test to make sure it is not a Molecule
            bonds = [sel for sel in self.get_selections() \
                     if (isinstance(sel, Bond) and not isinstance(sel, Molecule))]

            assert bonds
            return bonds


    @property
    def adjacent_atoms(self):

        # collect adjacent atoms
        adjacent_atoms = []
        for bond in self.bonds:
            # this atom's index in the bond
            curr_atom_idx = self.registry[bond.sel_reg_id]
            curr_atom = bond[curr_atom_idx]
            other_atom_idx = next((idx for idx in bond.keys() if idx != curr_atom_idx), None)
            other_atom = bond[other_atom_idx]
            assert curr_atom is not other_atom, "Cannot be bound to itself"
            adjacent_atoms.append(other_atom)
        return adjacent_atoms

class AtomTypeLibrary(SelectionTypeLibrary):
    def __init__(self):
        super().__init__()

    def add_type(self, atom_type, atom_name, rename=False):
        # type check that input types are AtomTypes
        assert isinstance(atom_type, AtomType), \
            "atom_type must be a subclass of AtomType, not {}".format(
                type(atom_type))

        return super().add_type(atom_type, atom_name, rename=rename)

class AtomType(SelectionType):
    def __init__(self, atom_attrs=None):
        super().__init__(attr_dict=atom_attrs)

class Bond(IndexedSelection):
    def __init__(self, atom_container=None, atom_ids=None):
        if atom_ids is not None:
            assert isinstance(atom_ids, tuple), \
                "atom_ids must be a length 2 tuple, not type{}".format(
                    type(atom_ids))
            assert len(atom_ids) == 2, \
                "atom_ids must be a length 2 tuple, not len={}".format(
                    len(atom_ids))
            assert all([(lambda x: isinstance(x, int))(i) for i in atom_ids]), \
                "atom_ids must be a length 2 tuple of ints"

        if atom_container is not None:
            assert issubclass(type(atom_container), col.Sequence), \
                "atom_container must be a subclass of collections.Sequence, not {}".format(
                    type(atom_container))

        super().__init__(atom_container, atom_ids)
        for atom in self.values():
            atom._in_bond = True

    @property
    def atoms(self):
        return tuple(self.values())

class MoleculeType(SelectionType):
    def __init__(self, mol_attrs=None):
        super().__init__(attr_dict=mol_attrs)

class MoleculeTypeLibrary(SelectionTypeLibrary):
    def __init__(self):
        super().__init__()

    def add_type(self, mol_type, mol_name, rename=False):
        # type check that input types are MoleculeTypes
        assert issubclass(type(mol_type)), \
            "mol_type must be a subclass of MoleculeType, not {}".format(
                type(mol_type))
        return super().add_type(mol_type, mol_name, rename=rename)

class RDKitMoleculeType(MoleculeType):
    def __init__(self, rdkit_molecule, mol_name=None):
        super().__init__()
        self.molecule = rdkit_molecule
        self.rdkit_mol = self.molecule
        self.name = mol_name
        self._features = None
        self._atom_type_library = None

    @property
    def atoms(self):
        return [atom for atom in self.molecule.GetAtoms()]

    @property
    def bonds(self):
        return [bond for bond in self.molecule.GetBonds()]

    @property
    def features(self):
        return self._features

    def atom_data(self, atom_idx):

        """Extracts useful information about an atom and returns it as a
dictionary.

        """
        atom = self.atoms[atom_idx]
        atom_dict = {}
        atom_dict['atomic_num'] = atom.GetAtomicNum()
        atom_dict['bond_degree_no_Hs'] = atom.GetDegree()
        # same but want a convenience attribute
        atom_dict['bond_degree'] = atom.GetDegree()
        atom_dict['bond_degree_with_Hs'] = atom.GetTotalDegree()
        # same but want a convenience attribute
        atom_dict['total_bond_degree'] = atom.GetTotalDegree()
        atom_dict['explicit_valence'] = atom.GetExplicitValence()
        atom_dict['implicit_valence'] = atom.GetImplicitValence()
        atom_dict['total_valence'] = atom.GetTotalValence()
        atom_dict['formal_charge'] = atom.GetFormalCharge()
        atom_dict['hybridization'] = atom.GetHybridization()

        atom_dict['is_aromatic'] = atom.GetIsAromatic()
        atom_dict['in_ring'] = atom.IsInRing()
        atom_dict['isotope'] = atom.GetIsotope()
        atom_dict['mass'] = atom.GetMass()
        atom_dict['num_radical_electrons'] = atom.GetNumRadicalElectrons()
        atom_dict['element'] = atom.GetSymbol()
        atom_dict['num_Hs'] = atom.GetTotalNumHs()
        monomer_info = atom.GetMonomerInfo()
        if monomer_info:
            atom_dict['pdb_name'] = monomer_info.GetName().strip()
            atom_dict['pdb_occupancy'] = monomer_info.GetOccupancy()
            atom_dict['pdb_residue_name'] = monomer_info.GetResidueName()
            atom_dict['pdb_temp_factor'] = monomer_info.GetTempFactor()

        atom_dict['rdkit_mol_idx'] = atom.GetIdx()
        atom_dict['name'] = atom_dict['rdkit_mol_idx']
        return atom_dict

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

        self._features = features


    def make_atom_type_library(self, type_name='name'):

        # initialize the library and names record
        atom_lib = AtomTypeLibrary()
        for atom_idx in range(self.molecule.GetNumAtoms()):
            # make a new atom type
            atom_type = AtomType(self.atom_data(atom_idx))
            atom_name = atom_type.__dict__[type_name]
            atom_lib.add_type(atom_type, atom_name=atom_name, rename=True)

        return atom_lib

    @property
    def atom_type_library(self):
        if self._atom_type_library is None:
            self._atom_type_library = self.make_atom_type_library()
        return self._atom_type_library

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
        self.make_atom_type_library()
        atoms = []
        for atom_idx in atom_idxs:
            atom_type = self.atom_type_library[atom_idx]
            atom = Atom(atom_array=coord_array, array_idx=atom_idx, atom_type=atom_type)
            atoms.append(atom)

        # handle bonds
        bonds = []
        for bond in self.molecule.GetBonds():
            bonds.append(Bond(atoms, (bond.GetBeginAtomIdx(), bond.GetEndAtomIdx())))

        # TODO handle and create angles
        angles = None

        return Molecule(atoms, bonds, angles, mol_type=self, external_mol_rep=(RDKitMoleculeType, self))


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
        self._in_system = False
        self._system = None
        self._features = None
        self._feature_families = None
        self._feature_family_selections = None
        self._feature_types = None
        self._feature_type_selections = None
        if mol_type is None:
            mol_type = MoleculeType({})
        self._molecule_type = mol_type
        self._atom_types = AtomTypeLibrary()
        # a dictionary of molecular representations from other libraries
        if external_mol_rep:
            self._external_mol_reps = {external_mol_rep[0] : external_mol_rep[1]}

        # set that each atom is in a molecule now
        for atom in self.atoms:
            atom._in_molecule = True
            self._atom_types.add_type(atom.atom_type, atom.atom_type.name, rename=True)


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
        assert not all([atom._in_molecule for atom in atoms]), \
            "all atoms must not be part of another molecule"



        molecule_dict= {'atoms' : atoms, 'bonds' : bonds, 'angles': angles}
        return molecule_dict

    @property
    def atom_types(self):
        return self._atom_types

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
    def system(self):
        return self._system

    @property
    def external_mol_reps(self):
        return self._external_mol_reps

    @property
    def features(self):
        return self.external_mol_reps[RDKitMoleculeType].features

    # TODO allow for tolerance
    def overlaps(self, other):
        """Check whether this molecule overlaps with another.
        Checks whether any two atoms in each molecule have the same coordinates.

        bool : returns True if any overlaps detected

        """
        assert isinstance(other, Molecule), \
            "Other must be type Molecule, not {}".format(type(other))

        pairs = product(self.atoms, other.atoms)
        try:
            pair = next(pairs)
        # if it is empty no overlaps
        except StopIteration:
            return False
        flag = True
        while flag:
            overlaps = np.isclose(pair[0].coords, pair[1].coords)
            if np.all(overlaps):
                return (pair[0], pair[1])
            else:
                try:
                    pair = next(pairs)
                except StopIteration:
                    flag = False
        return False

    def find_features(self):
        # find the features
        self._external_mol_reps[RDKitMoleculeType].find_features()
        self._feature_families = col.defaultdict(list)
        self._feature_types = col.defaultdict(list)
        for idx, info in self.features.items():
            self._feature_families[info['family']].append(idx)
            self._feature_types[info['type']].append(idx)

        # make selections out of these to each feature
        # for families
        self._feature_family_selections = SelectionDict()
        for family, idxs in self._feature_families.items():
            atom_idxs = self._feature_families[family]
            self._feature_family_selections[family] = IndexedSelection(self.atoms, atom_idxs)

        # for types
        self._feature_type_selections = SelectionDict()
        for ftype, idxs in self._feature_types.items():
            atom_idxs = self._feature_types[ftype]
            self._feature_type_selections[ftype] = IndexedSelection(self.atoms, atom_idxs)

    @property
    def family_selections(self):
        return self._feature_family_selections

    @property
    def type_selections(self):
        return self._feature_type_selections

    @property
    def feature_families(self):
        return set(self._feature_families.keys())

    @property
    def feature_types(self):
        return set(self._feature_types.keys())

    @property
    def feature_dataframe(self):
        import pandas as pd
        return pd.DataFrame(self.features, orient='index')

    def profile_interactions(self, interaction_types):
        assert all([issubclass(itype, InteractionType) for itype in interaction_types]), \
                   "All interaction_types must be a subclass of InteractionType"
        # go through each interaction_type and check for hits
        interactions = {}
        for interaction_type in interaction_types:
            # collect the specific features for each family
            family_features = {}
            for family in interaction_type.feature_families:
                try:
                    family_features[family] = self.family_selections[family].values()
                except KeyError:
                    print("No {0} features in {1} for profiling {2}".format(
                        family, self, interaction_type))
                    return None

            # pass these to the check class method of the InteractionType
            interactions[str(interaction_type)] = interaction_type.find_hits(**family_features)

        return interactions

if __name__ == "__main__":

    print("making an CoordArray for atoms")
    array = np.array([[0,0,0], [0,0,1], [1,0,0]])
    atom_array = CoordArray(array)
    print(atom_array)

    print("making atoms")
    atom1 = Atom(np.array([5,5,5]))
    print(atom1)
    print(atom1.coords)
    print("making AtomType")
    atom2_type = AtomType({'element':'C'})

    atom2 = Atom(np.array([6,6,6]), atom_array=atom_array, atom_type=atom2_type)
    print(atom2)
    print(atom2.coords)
    print("testing overlap of two atoms")
    print(atom2.overlaps(atom1))
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
    pka_type = RDKitMoleculeType(pka, mol_name="PKA")
    print(pka_type)

    print("getting positions objects from rdchem.Mol for atoms")
    conf = list(pka_type.molecule.GetConformers()).pop()
    positions = {idx : conf.GetAtomPosition(idx) for idx in range(conf.GetNumAtoms())}
    print(positions[0])
    print("Making an AtomTypeLibrary and Atom list for pka")
    atom_types = []
    atom_names = {}
    atoms = []
    pka_atom_type_library = AtomTypeLibrary()
    for atom_idx in range(len(pka_type.atoms)):
        atom_type = AtomType(pka_type.atom_data(atom_idx))
        atom_types.append(atom_type)

        atom_name = atom_type.pdb_name
        if atom_name in atom_names.keys() and \
           not pka_atom_type_library.attributes_match(atom_type):
            atom_names[atom_name] += 1
        elif atom_name not in atom_names.keys():
            atom_names[atom_name] = 0

        if atom_names[atom_name] > 0:
            pka_atom_type_library.add_type(atom_type, atom_name + str(atom_names[atom_name]) )
        else:
            pka_atom_type_library.add_type(atom_type, atom_name)

        coord = np.array([positions[atom_idx].x, positions[atom_idx].y, positions[atom_idx].z])
        atoms.append(Atom(coords=coord, atom_type=atom_type))

    print(pka_atom_type_library)
    print(atoms)

    # make a selection of atoms for bonds, and angle

    # making a bond between each atom and the next one in the index,
    # just to test
    print("making up bonds")
    idx_a = range(len(atoms))[:-1]
    idx_b = [a+1 for a in idx_a]
    idx = zip(idx_a, idx_b)
    bonds = [Bond(atoms, bond_idx) for bond_idx in idx]

    print("accessing bonds from an atom")
    print("Is an atom in a bond?")
    print(atoms[0]._in_bond)
    print("how many bonds is it in")
    print("first atom", len(atoms[0].bonds))
    print("second atom", len(atoms[1].bonds))
    print("get the bonds themselves")
    print(atoms[0].bonds)
    print("get the other atom in the bond")
    bond = atoms[0].bonds[0]
    curr_atom_idx = atoms[0].registry[bond.sel_reg_id]
    curr_atom = bond[curr_atom_idx]
    other_atom_idx = next((idx for idx in bond.keys() if idx != curr_atom_idx), None)
    other_atom = bond[other_atom_idx]
    print("idx:", other_atom_idx, "atom:", other_atom)
    print("same atom", curr_atom is atoms[curr_atom_idx])
    print("same thing using the method")
    print(atoms[0].adjacent_atoms)
    # angles still stubbed out for now
    angles = [IndexedSelection(atoms, [0,1,2])]
    print("making a MoleculeType")
    moltype = MoleculeType({'name':'test_molecule'})
    print("making a molecule")
    mol = Molecule(atoms, bonds, angles, mol_type=moltype)
    print(mol)
    print("atom_types in mol")
    print(mol.atom_types)

    print("Making a mast.Molecule from the RDKitMoleculeType")
    pka_mol = pka_type.to_molecule(0)
    # pka_mol = Molecule(mol_type=pka_type, coords=pka_coords)
    print(pka_mol)
    print(pka_mol.molecule_type)

    print("testing overlap of two molecules")
    print(pka_mol.overlaps(mol))

    print("finding features using mast.molecule method")
    pka_mol.find_features()
    print(pka_mol.features)
