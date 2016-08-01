import numpy as np
import collections as col
import os.path as osp
from itertools import product

from mast.selection import CoordArray, CoordArraySelection, \
    Point, IndexedSelection, SelectionDict, SelectionList, \
    SelectionType, SelectionTypeLibrary, Container
from mast.interactions import InteractionType

__all__ = ['Atom', 'AtomTypeLibrary', 'AtomType', 'Bond', 'MoleculeType', 'MoleculeTypeLibrary',
           'RDKitMoleculeType', 'Molecule', ]

DIM_NUM_3D = 3

ATOM_ATTRIBUTES = ['name',
                  'atomic_num',
                  'bond_degree_no_Hs',
                  'bond_degree_with_Hs',
                  'total_bond_degree',
                  'explicit_valence',
                  'implicit_valence',
                  'total_valence',
                  'formal_charge',
                  'hybridization',
                  'is_aromatic',
                  'in_ring',
                  'isotope',
                  'mass',
                  'num_radical_electrons',
                  'element',
                  'num_Hs',
                  'monomer_type',
                  'pdb_name',
                  'pdb_insertion_code',
                  'pdb_occupancy'
                  'pdb_residue_name',
                  'pdb_residue_number',
                  'pdb_serial_number',
                  'pdb_temp_factor',]

BOND_ATTRIBUTES = ['name',
                   'atom_types',
                   'bond_type',
                   'bond_type_number',
                   'is_aromatic',
                   'in_ring',
                   'stereo',
                   'is_conjugated',]

MOLECULE_ATTRIBUTES = ['name',
                       'pdb_name',
                       'atom_types',
                       'bond_types',
                       'angle_types',
                       'features']

def _atom_type_factory(atom_type_name, **atom_attrs):

    # simply keep track of which attributes the input did not provide
    for attr in AtomType.attributes:
        try:
            assert attr in atom_attrs.keys()
        except AssertionError:
            # logging
            print("Attribute {0} not found in atom input.".format(attr))

    # add the attributes that it has into the class
    attributes = {attr : None for attr in AtomType.attributes}
    for attr, value in atom_attrs.items():
        # make sure AtomType has an attribute that matches
        try:
            assert attr in AtomType.attributes
        # if it doesn't then report it
        except AssertionError:
            # logging
            print("Input attribute {0} not in AtomType attributes, ignoring.".format(attr))
        # if no error then add it
        else:
            attributes[attr] = value

    return type(atom_type_name, (AtomType,), attributes)

AtomType = type('AtomType', (object,), {attr : None for attr in ATOM_ATTRIBUTES})
AtomType.attributes = ATOM_ATTRIBUTES
AtomType.factory = _atom_type_factory


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
            atom_type = AtomType()
        self._atom_type = atom_type
        self._in_molecule = False
        self._in_bond = False


    def __repr__(self):
        return str(self.__class__)

    @property
    def atom_type(self):
        return self._atom_type

    @property
    def molecule(self):
        if self._in_molecule is False:
            return None
        else:
            molecule = next((sel for key, sel in self.registry
                             if isinstance(sel, Molecule)),
                            None)
            assert molecule
            return molecule

    @property
    def system(self):
        from mast.system import System

        if self._in_system is False:
            return None
        elif self._in_molecule is False:
            system = next((sel for key, sel in self.registry
                           if isinstance(sel, System)),
                          None)
            assert system
            return system
        else:
            system = self.molecule.system
            assert system
            return system

    @property
    def bonds(self):
        if self._in_bond is False:
            return None
        else:
            bonds = []
            for key, sel in self.registry:
                if isinstance(sel, Bond):
                    if not isinstance(sel, Molecule):
                        bonds.append(sel)

            assert bonds
            return bonds


    @property
    def adjacent_atoms(self):

        # collect adjacent atoms
        adjacent_atoms = []
        for bond in self.bonds:
            other_atom = next((a for a in bond.atoms if a is not self), None)
            assert self is not other_atom, "{} cannot be bound to itself".format(self)
            adjacent_atoms.append(other_atom)
        return adjacent_atoms

class BondType(object):
    attributes = BOND_ATTRIBUTES

    def __init__(self):
        pass

    @classmethod
    def factory(cls, bond_type_name, **bond_attrs):

        # must pass in atom types
        assert 'atom_types' in bond_attrs, \
            "'atom_types' must be in the bond_attrs"
        # should only pass in two atom types
        assert len(bond_attrs['atom_types']) == 2, \
            "atom_ids must be a length 2 tuple, not len={}".format(
                len(bond_attrs['atom_types']))
        assert all([(lambda x: issubclass(type(x), AtomType))(i)
                    for i in bond_attrs['atom_types']]), \
            "atom_ids must be a length 2 tuple of AtomType subclasses"

        # simply keep track of which attributes the input did not provide
        for attr in BondType.attributes:
            try:
                assert attr in bond_attrs.keys()
            except AssertionError:
                # logging
                print("Attribute {0} not found in bond input.".format(attr))

        # add the attributes that it has into the class
        attributes = {attr : None for attr in BondType.attributes}
        for attr, value in bond_attrs.items():
            # make sure BondType has an attribute that matches
            try:
                assert attr in BondType.attributes
            # if it doesn't then report it
            except AssertionError:
                # logging
                print("Input attribute {0} not in BondType attributes, ignoring.".format(attr))
            # if no error then add it
            else:
                attributes[attr] = value

        return type(bond_type_name, (BondType,), attributes)




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



class MoleculeType(object):

    attributes = MOLECULE_ATTRIBUTES

    def __init__(self, name=None, atom_types=None, ):
        if atom_types:
            self._atom_type_library = set(atom_types)
        else:
            self._atom_type_library = set()


        self.name = name
        self._atom_type_sequence = atom_types
        self._bonds = None

        # these must be manually set due to heavy computation time
        self._features = None
        self._feature_families = None
        self._feature_types = None


    @property
    def atom_type_library(self):
        return list(self._atom_type_library)

    @property
    def features(self):
        return self._features

    @features.setter
    def features(self, features):
        self._features = features

    @property
    def feature_families_map(self):
        """A dictionary mapping the feature families to the indices of the
        feature."""

        families = col.defaultdict(list)
        for idx, info in self.features.items():
            families[info['family']].append(idx)
        return families

    @property
    def feature_types_map(self):
        """A dictionary mapping the feature types to the indices of the
        feature."""

        types = col.defaultdict(list)
        for idx, info in self.features.items():
            types[info['type']].append(idx)

        return types

    @property
    def feature_families(self):
        return set(self.feature_families_map.keys())

    @property
    def feature_types(self):
        return set(self.feature_types_map.keys())

    @property
    def atom_types(self):
        return self._atom_type_sequence

    @property
    def bonds(self):
        return self._bonds

    def to_molecule(self, coords):
        """ Construct a Molecule using input coordinates with mapped indices"""

        coord_array = CoordArray(coords)

        # Make atoms out of the coord array
        atom_idxs = range(len(self.atom_types))
        atoms = []
        for atom_idx, atom_type in enumerate(self.atom_types):
            atom = Atom(atom_array=coord_array, array_idx=atom_idx, atom_type=atom_type)
            atoms.append(atom)

        # TODO handle bonds
        bonds = None

        # TODO handle and create angles
        angles = None

        return Molecule(atoms, bonds, angles, mol_type=self)

    @classmethod
    def factory(cls, mol_type_name, **molecule_attrs):

        # simply keep track of which attributes the input did not provide
        for attr in MoleculeType.attributes:
            try:
                assert attr in molecule_attrs.keys()
            except AssertionError:
                # logging
                print("Attribute {0} not found in molecule input.".format(attr))

        # add the attributes that it has into the class
        attributes = {attr : None for attr in MoleculeType.attributes}
        for attr, value in molecule_attrs.items():
            # make sure MoleculeType has an attribute that matches
            try:
                assert attr in MoleculeType.attributes
            # if it doesn't then report it
            except AssertionError:
                # logging
                print("Input attribute {0} not in MoleculeType attributes, ignoring.".format(attr))
            # if no error then add it
            else:
                attributes[attr] = value

        return type(mol_type_name, (MoleculeType,), attributes)

class Molecule(Container):
    def __init__(self, mol_input, *args, **kwargs):

        if 'mol_type' not in kwargs.keys():
            mol_type = None
        else:
            mol_type = kwargs.pop('mol_type')

        # if 'external_mol_rep' not in kwargs.keys():
        #     external_mol_rep = None
        # else:
        #     assert isinstance(kwargs['external_mol_rep'], tuple), \
        #         "An external_mol_rep must be a tuple (external_type, external_mol), not {}".format(
        #             kwargs['external_mol_rep'])
        #     external_mol_rep = kwargs.pop('external_mol_rep')

        # check to see which constructor to use
        if issubclass(type(mol_input), MoleculeType):
            molecule_dict = Molecule.type_constructor(mol_input, *args, **kwargs)
        elif issubclass(type(mol_input), col.Sequence):
            molecule_dict = Molecule.atoms_constructor(mol_input, *args, **kwargs)
        else:
            raise TypeError("mol_input must be either a MoleculeType or a sequence of Atoms")

        # call to parent class
        super().__init__()

        # set the atoms, bonds, and angles into this object
        self.atoms = molecule_dict['atoms']
        self.bonds = molecule_dict['bonds']
        self.angles = molecule_dict['angles']

        # set the molecule_type
        if mol_type is None:
            mol_type = MoleculeType()
        self._molecule_type = mol_type

        # initialize flags
        if 'system' in kwargs.keys():
            self._in_system = True
        else:
            self._in_system = False

        # set that each atom is in a molecule now
        for atom in self.atoms:
            atom._in_molecule = True
            if self._in_system is True:
                atom._in_system = True

        # attributes must explicitly be called due to computation time
        self._feature_family_selections = {}
        self._feature_type_selections = None
        self._internal_interactions = None

        self._external_mol_reps = []

        # an optional dictionary of molecular representations from
        # other libraries
        if "external_mol_rep" in kwargs:
            self._external_mol_reps.append(kwargs["external_mol_reps"])

    # the alternate constructors
    @classmethod
    def type_constructor(cls, mol_type, coords=None):
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

        molecule_dict = {'atoms' : atoms, 'bonds' : bonds, 'angles': angles}
        return molecule_dict

    # properties
    @property
    def atom_types(self):
        pass

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
    def system(self):
        from mast.system import System
        if self._in_system is False:
            return None
        else:
            system = next((sel for key, sel in self.registry
                           if isinstance(sel, System)),
                          None)
            assert system
            return system

    @property
    def atom_coords(self):
        coords = np.array([atom.coords for atom in self.atoms])
        return coords

    @property
    def external_mol_reps(self):
        return self._external_mol_reps

    @property
    def features(self):
        return self.molecule_type.features

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

    def make_feature_selections(self):
        family_selections = col.defaultdict(list)
        type_selections = col.defaultdict(list)
        for idx, feature in self.features.items():
            atom_idxs = list(feature['atom_ids'])
            # make the selection
            feature_selection = IndexedSelection(self.atoms, atom_idxs)
            # add it to it's families selections
            family_selections[feature['family']].append(feature_selection)
            # add it to it's type's selections
            type_selections[feature['type']].append(feature_selection)

        self._feature_family_selections = SelectionDict(family_selections)
        self._feature_type_selections = SelectionDict(type_selections)

    @property
    def family_selections(self):
        return self._feature_family_selections

    @property
    def type_selections(self):
        return self._feature_type_selections

    @property
    def feature_dataframe(self):
        import pandas as pd
        return pd.DataFrame(self.features, orient='index')

    @property
    def internal_interactions(self):
        return self._internal_interactions

    def profile_interactions(self, interaction_types):
        assert all([issubclass(itype, InteractionType) for itype in interaction_types]), \
                   "All interaction_types must be a subclass of InteractionType"
        # go through each interaction_type and check for hits
        interactions = {}
        for interaction_type in interaction_types:
            # collect the specific feature selections for each family
            family_feature_sels = {}
            for family in interaction_type.feature_families:
                # get the features from the molecule that are of the family
                try:
                    family_feature_sels[family] = self.family_selections[family]
                except KeyError:
                    # if there is none of a feature family then the
                    # interaction will not exist
                    print("No {0} features in {1} for profiling {2}".format(
                        family, self, interaction_type))
                    return None

            # pass these to the find_hits method of the InteractionType
            interactions[interaction_type] = interaction_type.find_hits(**family_feature_sels)

        self._internal_interactions = interactions

class RDKitMoleculeWrapper(object):
    def __init__(self, rdkit_molecule, mol_name=None):
        self.rdkit_molecule = rdkit_molecule
        self.mol_name = mol_name

    @property
    def atoms(self):
        return [atom for atom in self.rdkit_molecule.GetAtoms()]

    @property
    def bonds(self):
        return [bond for bond in self.rdkit_molecule.GetBonds()]

    @property
    def num_atoms(self):
        return self.rdkit_molecule.GetNumAtoms()

    @property
    def num_bonds(self):
        return self.rdkit_molecule.GetNumBonds()

    @property
    def num_conformers(self):
        return self.rdkit_molecule.GetNumConformers()

    def atoms_data(self):
        atoms_data = []
        for atom in self.rdkit_molecule.GetAtoms():
            atom_data = self.atom_data(atom.GetIdx())
            atoms_data.append(atom_data)

        return atoms_data

    def atom_data(self, atom_idx):
        """Extracts useful rdkit information about an atom and returns it as a
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
            atom_dict['monomer_type'] = monomer_info.GetMonomerType()
            atom_dict['pdb_name'] = monomer_info.GetName().strip()
            # atom_dict['pdb_chain_id'] = monomer_info.GetChainID()
            atom_dict['pdb_insertion_code'] = monomer_info.GetInsertionCode()
            # atom_dict['pdb_heteroatom'] = monomer_info.IsHeteroAtom()
            atom_dict['pdb_occupancy'] = monomer_info.GetOccupancy()
            atom_dict['pdb_residue_name'] = monomer_info.GetResidueName()
            atom_dict['pdb_residue_number'] = monomer_info.GetResidueNumber()
            atom_dict['pdb_serial_number'] = monomer_info.GetSerialNumber()
            # atom_dict['pdb_segment_number'] = monomer_info.GetSegmentNumber()
            atom_dict['pdb_temp_factor'] = monomer_info.GetTempFactor()

        atom_dict['rdkit_mol_idx'] = atom.GetIdx()
        atom_dict['name'] = atom_dict['rdkit_mol_idx']

        return atom_dict

    def bonds_data(self):
        bonds_data = []
        for bond_idx, bond in enumerate(self.bonds):
            bond_data = self.bond_data(bond_idx)
            bonds_data.append(bond_data)

        return bonds_data

    def bond_data(self, bond_idx):
        """Extracts useful rdkit information about an atom and returns it as a
        dictionary.

        """
        bond = self.bonds[bond_idx]
        bond_dict = {}
        bond_dict['bond_type'] = str(bond.GetBondType())
        bond_dict['bond_type_number'] = str(bond.GetBondTypeAsDouble())
        bond_dict['is_aromatic'] = bond.GetIsAromatic()
        bond_dict['in_ring'] = bond.IsInRing()
        bond_dict['stereo'] = str(bond.GetStereo())
        bond_dict['is_conjugated'] = bond.GetIsConjugated()
        atom1_idx = bond.GetBeginAtomIdx()
        atom2_idx = bond.GetEndAtomIdx()
        bond_dict['rdkit_atom_idxs'] = (atom1_idx, atom2_idx)
        bond_dict['rdkit_mol_idx'] = bond.GetIdx()
        bond_dict['name'] = bond_dict['rdkit_mol_idx']

        return bond_dict

    def molecule_data(self):
        """Extracts useful rdkit information about an atom and returns it as a
        dictionary.

        """
        molecule_dict = {}
        molecule_dict['name'] = self.mol_name
        ring_info = self.rdkit_molecule.GetRingInfo()
        molecule_dict['num_rings'] = ring_info.NumRings()
        molecule_dict['num_atoms'] = self.rdkit_molecule.GetNumAtoms()
        molecule_dict['num_bonds'] = self.rdkit_molecule.GetNumBonds()
        molecule_dict['num_heavy_atoms'] = self.rdkit_molecule.GetNumHeavyAtoms()

        return molecule_dict

    def find_features(self, fdef="BaseFeatures.fdef"):
        """Uses a feature definition (fdef) database to to find features in
        the molecule.

        Returns a tuple (dict : features, list : families, list : types)

        """
        from rdkit import RDConfig
        from rdkit.Chem import ChemicalFeatures

        assert isinstance(fdef, str)
        fdef_path = osp.join(RDConfig.RDDataDir, fdef)
        feature_factory = ChemicalFeatures.BuildFeatureFactory(fdef_path)
        factory_features = feature_factory.GetFeaturesForMol(self.rdkit_molecule)
        features = {}
        for feature in factory_features:
            feature_info = {'family' : feature.GetFamily(),
                            'type' : feature.GetType(),
                            'atom_ids' : feature.GetAtomIds(),
                            'position' : feature.GetPos()}
            features[feature.GetId()] = feature_info


        families = col.defaultdict(list)
        types = col.defaultdict(list)
        for idx, info in features.items():
            families[info['family']].append(idx)
            types[info['type']].append(idx)

        return features

    def get_conformer_coords(self, conf_idx):
        assert self.rdkit_molecule.GetNumConformers() > 0, \
            "{0} has no conformers".format(self)

        conformer = self.rdkit_molecule.GetConformer(conf_idx)
        atom_idxs = range(self.rdkit_molecule.GetNumAtoms())
        # make the CoordArray
        coords = []
        for atom_idx in atom_idxs:
            coord = conformer.GetAtomPosition(atom_idx)
            coord = np.array([coord.x, coord.y, coord.z])
            coords.append(coord)
        coords = np.array(coords)
        return coords

if __name__ == "__main__":
    pass
