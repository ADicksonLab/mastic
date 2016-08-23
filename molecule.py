import numpy as np
import collections as col
import os.path as osp
from itertools import product


from mast.selection import CoordArray, CoordArraySelection, \
    Point, IndexedSelection, SelectionsDict, SelectionsList, \
    Selection


import mast.config.molecule as mastmolconfig

__all__ = ['AtomType', 'BondType', 'MoleculeType',]

class AtomType(object):
    """Class for generating specific atom type classes with the factory
    method.

    Examples
    --------
    >>> CarbonAtomType = AtomType.factory("CarbonAtomType", **{"element" : 'C'})
    >>> CarbonAtomType.element
    'C'

    """
    attributes = mastmolconfig.ATOM_ATTRIBUTES

    def __init__(self, atom_type_name, **atom_attrs):
        # keep track of which attributes the input did not provide
        for attr in AtomType.attributes:
            try:
                assert attr in atom_attrs.keys()
            except AssertionError:
                # LOGGING
                pass
                # print("Attribute {0} not found in atom input.".format(attr))

        # add the attributes into the class
        attributes = {attr : None for attr in AtomType.attributes}
        for attr, value in atom_attrs.items():
            # Log the compliance of the attributes
            # if the attribute isn't declared in the attributes log it
            try:
                assert attr in AtomType.attributes
            except AssertionError:
                # LOGGING
                pass
                # print("Input attribute {0} not in AtomType attributes.".format(attr))

            # add it to the attributes
            attributes[attr] = value

        self.name = atom_type_name
        self.attributes_data = attributes
        self.__dict__.update(attributes)

    def to_atom(self, coords, flags=None):
        """Substantiate this AtomType with coordinates"""
        assert len(coords) == 3, \
            "coords must be length 3, not {}".format(len(coords))
        assert all([(lambda x: isinstance(x, float))(i)
                    for i in coords]), \
            "coords must be 3 floats, not {}".format(coords)
        coords = np.array(coords)
        atom = Atom(coords, atom_type=self, flags=flags)

        return atom

    def __str__(self):
        return self.name

    @property
    def record(self):
        # define the Record named tuple
        record_fields = ['AtomType'] + list(self.attributes_data.keys())
        AtomTypeRecord = col.namedtuple('AtomTypeRecord', record_fields)
        # build the values for it for this Type
        record_attr = {'AtomType' : self.name}
        record_attr.update(self.attributes_data)
        # make and return
        return AtomTypeRecord(**record_attr)

class BondType(object):
    """Class for generating specific bond type classes with the factory
    method.

    Examples
    --------
    First make the components of a bond:
    >>> carbon_attributes = {'element':'C', 'bond_degree':3}
    >>> oxygen_attributes = {'element':'O', 'bond_degree':3}
    >>> COCarbonAtomType = AtomType.factory("COCarbonAtomType", **carbon_attributes)
    >>> COOxygenAtomType = AtomType.factory("COOxygenAtomType", **oxygen_attributes)

    >>> CO_atoms = (COCarbonAtomType, COOxygenAtomType)
    >>> CO_attributes = {"bond_order":3}

    Then put them together:
    >>> COBondType = BondType.factory("COBondType", atom_types=CO_atoms, **CO_attributes)
    """
    attributes = mastmolconfig.BOND_ATTRIBUTES
    """Domain specific properties of the BondType"""

    def __init__(self, bond_type_name, atom_types=None, **bond_attrs):
        """Static method for generating bond types dynamically given a type
        name (which will be the class name) and a domain specific dictionary
        of bond attributes.

        See mast.config.molecule for standard BondType attributes.
        See class docstring for examples.
        """

        # must pass in atom types
        assert atom_types, \
            "'atom_types' must be provided"
        # should only pass in two atom types
        assert len(atom_types) == 2, \
            "atom_types must be a length 2 tuple, not len={}".format(
                len(atom_types))
        assert all([(lambda x: isinstance(x, AtomType))(i)
                    for i in atom_types]), \
            "atom_types must be a length 2 tuple of AtomType subclasses"

        # log keep track of which attributes the input did not provide
        for attr in BondType.attributes:
            try:
                assert attr in bond_attrs.keys()
            except AssertionError:
                # LOGGING
                pass
                # print("Attribute {0} not found in bond input.".format(attr))

        # add the attributes into the class
        attributes = {attr : None for attr in BondType.attributes}
        for attr, value in bond_attrs.items():
            try:
                assert attr in BondType.attributes
            # if it doesn't then log it
            except AssertionError:
                # LOGGING
                pass
                # print("Input attribute {0} not in BondType attributes.".format(attr))
            # add it regardless
            attributes[attr] = value

        self.name = bond_type_name
        self.attributes_data = attributes
        self.__dict__.update(attributes)
        self.atom_types = atom_types

    def to_bond(self, atom1_coords, atom2_coords):
        """Substantiate this Bond with given coordinates for each atom type."""
        # test the inputs
        for coords in [atom1_coords, atom2_coords]:
            assert len(coords) == 3, \
                "coords must be length 3, not {}".format(len(coords))
            assert all([(lambda x: isinstance(x, int) or isinstance(x, float))(i)
                        for i in coords]), \
                "coords must be 3 numbers, not {}".format(coords)

        atoms = []
        atoms.append(self.atom_types[0].to_atom(atom1_coords, flags='bond'))
        atoms.append(self.atom_types[1].to_atom(atom2_coords, flags='bond'))

        bond = Bond(atoms, (0,1), bond_type=self)
        return bond


    @property
    def record(self):
        # define the Record named tuple
        record_fields = ['BondType', 'AtomAType', 'AtomBType'] + \
                        list(self.attributes_data.keys())
        BondTypeRecord = col.namedtuple('BondTypeRecord', record_fields)
        # build the values for it for this Type
        record_attr = {'BondType' : self.name}
        record_attr['AtomAType'] = self.atom_types[0]
        record_attr['AtomBType'] = self.atom_types[1]
        record_attr.update(self.attributes_data)
        # make and return
        return BondTypeRecord(**record_attr)

class MoleculeType(object):
    """Class for generating specific bond type classes with the factory
    method.

    Examples
    --------

    We will make a carbon-monoxide molecule from scratch:

    >>> carbon_attributes = {'element':'C', 'bond_degree':3}
    >>> oxygen_attributes = {'element':'O', 'bond_degree':3}
    >>> COCarbonAtomType = AtomType.factory("COCarbonAtomType", **carbon_attributes)
    >>> COOxygenAtomType = AtomType.factory("COOxygenAtomType", **oxygen_attributes)

    >>> CO_atoms = (COCarbonAtomType, COOxygenAtomType)
    >>> CO_attributes = {"bond_order":3}
    >>> COBondType = BondType.factory("COBondType", atom_types=CO_atoms, **CO_attributes)

    Map the bonds to the atomtypes in a dict:

    >>> atom_types = [COCarbonAtomType, COOxygenAtomType]
    >>> bond_types = [COBondType]
    >>> bond_map = {0 : (0, 1)}

    Define domain specific information about the MoleculeType:

    >>> CO_attributes = {"name" : "carbon-monoxide", "toxic" : True}

    Put it together:
    >>> COType = MoleculeType.factory("COType", atom_types=atom_types, bond_types=bond_types, bond_map=bond_map, **CO_attributes)
    """
    attributes = mastmolconfig.MOLECULE_ATTRIBUTES
    """Domain specific properties of the MoleculeType"""

    def __init__(self, mol_type_name, atom_types=None,
                bond_types=None, bond_map=None,
                angle_types=None, angle_map=None, # stubbed out
                **molecule_attrs):
        """Static method for generating molecule types dynamically given a type
        name (which will be the class name) and a domain specific dictionary
        of molecule attributes.

        See mast.config.molecule for standard MoleculeType attributes.
        See class docstring for examples.
        """


        # check required inputs for validity
        assert atom_types, "atom_types must be provided"
        assert bond_types, "bond_types must be provided"
        assert bond_map, "bond_map must be provided"
        # assert angle_types, "angle_types must be provided"
        assert all([(lambda x: isinstance(x, AtomType))(i)
                    for i in atom_types]), \
            "atom_types must contain only AtomType subclasses"

        assert all([(lambda x: isinstance(x, BondType))(i)
                    for i in bond_types]), \
            "bond_types must contain only BondType subclasses"

        # check that the bond_map is valid
        for bond_idx, atom_idxs in bond_map.items():
            assert bond_idx < len(bond_types), \
                "bond_idx must be in range {0}, not {1}".format(
                    len(bond_types), bond_idx)
            assert len(atom_idxs) == 2, "bond {0} must map to a 2-tuple of ints"
            assert all([(lambda x: isinstance(x, int))(i)
                        for i in atom_idxs]), \
                "bond_map values must be 2-tuples of ints, not {}".format(atom_idxs)
            assert atom_idxs[0] < len(atom_types), \
                "bond ({2}) atom_type 0 index must be in range of atom_types {1}, not {0}".format(
                    atom_idxs[0], len(atom_types), bond_idx)
            assert atom_idxs[1] < len(atom_types), \
                "bond ({2}) atom_type 1 index must be in range of atom_types {1}, not {0}".format(
                    atom_idxs[1], len(atom_types), bond_idx)

        # check angle input
        # assert all([(lambda x: isinstance(x, AngleType))(i)
        #             for i in angle_types]), \
        #     "angle_types must contain only AngleType subclasses"

        # keep track of which attributes the input did not provide
        # compared to the config file
        for attr in MoleculeType.attributes:
            try:
                assert attr in molecule_attrs.keys()
            except AssertionError:
                pass
                # LOGGING
                # print("Attribute {0} not found in molecule input.".format(attr))

        # add the attributes into the class
        attributes = {attr : None for attr in MoleculeType.attributes}
        for attr, value in molecule_attrs.items():
            try:
                assert attr in MoleculeType.attributes
            # if it doesn't then log
            except AssertionError:
                # LOGGING
                pass
                # print("Input attribute {0} not in MoleculeType attributes.".format(attr))
            # add it regardless
            attributes[attr] = value

        # create the class with the domain specific attributes
        self.name = mol_type_name
        self.attributes_data = attributes
        self.__dict__.update(attributes)
        # add core attributes
        self.atom_types = atom_types
        self.bond_types = bond_types
        self.bond_map = bond_map
        # self.angle_types = angle_types
        # self.angle_map = angle_map
        self._atom_type_library = set(atom_types)
        self._bond_type_library = set(bond_types)
        # self.angle_type_library = set(angle_types)
        self._feature_types = {}

    @property
    def atom_type_library(self):
        """The unique AtomTypes in this MoleculeType."""
        return list(self._atom_type_library)

    @property
    def bond_type_library(self):
        """The unique BondTypes in this MoleculeType."""
        return list(self._bond_type_library)

    @property
    def feature_type_library(self):
        return set(self.feature_types.values())

    @property
    def feature_types(self):
        """The chemical features of this MoleculeType."""
        return self._feature_types

    def add_feature_type(self, key, feature_type):
        from mast.features import FeatureType
        assert isinstance(feature_type, FeatureType), \
            "feature must be a subclass of a mast.features.FeatureType class, not {}".format(
                feature_type)
        assert feature_type.molecule_type is self, \
            "The molecule of the feature must be this molecule, not {}".format(
                feature_type.molecule_type)
        self._feature_types[key] = feature_type

    def to_molecule(self, coords):
        """Substantiate this MoleculeType with the given coordinates for each
        atom.
        """

        # make one CoordArray to put everything in
        coord_array = CoordArray(coords)

        # Make atoms out of the coord array
        atoms = []
        for atom_idx, atom_type in enumerate(self.atom_types):
            atom = Atom(atom_array=coord_array, array_idx=atom_idx, atom_type=atom_type)
            atoms.append(atom)

        # make the bonds from the atoms
        bond_map = self.bond_map
        bonds = []
        for bond_idx, bond_type in enumerate(self.bond_types):
            bond = Bond(atom_container=atoms, atom_ids=bond_map[bond_idx],
                        bond_type=bond_type)
            bonds.append(bond)

        # TODO handle and create angles
        angles = None

        return Molecule(atoms=atoms, bonds=bonds, angles=angles, mol_type=self)


    @property
    def atom_type_records(self):
        return [atom_type.record for atom_type in self.atom_types]

    @property
    def atom_type_df(self):
        import pandas as pd
        return pd.DataFrame(self.atom_type_records)

    @property
    def bond_type_records(self):
        return [bond_type.record for bond_type in self.bond_types]

    @property
    def bond_type_df(self):
        import pandas as pd
        return pd.DataFrame(self.bond_type_records)

    @property
    def feature_type_records(self):
        return [feature_type.record for feat_key, feature_type in
                self.feature_types.items()]

    @property
    def feature_type_df(self):
        import pandas as pd
        return pd.DataFrame(self.feature_type_records)

    @property
    def record(self):
        # define the Record named tuple
        record_fields = ['MoleculeType'] + list(self.attributes_data.keys())
        MoleculeTypeRecord = col.namedtuple('MoleculeTypeRecord', record_fields)
        # build the values for it for this Type
        record_attr = {'MoleculeType' : self.name}
        record_attr.update(self.attributes_data)
        # make and return
        return MoleculeTypeRecord(**record_attr)


    def atom_types_by_attr(self, attr_field, attr_value):
        return [atom_type for atom_type in self.atom_types if
                atom_type.attributes_data[attr_field] == attr_value]

class Atom(Point):
    """The coordinate substantiation of an AtomType.

    Examples
    --------

    >>> import numpy as np

    Make an AtomType:

    >>> CarbonAtomType = AtomType.factory("CarbonAtomType", **{"element" : 'C'})

    Coordinates for it:

    >>> coords = np.array((0.0, 0.0, 0.0))

    Construct an Atom from it:

    >>> Atom(coords=coords, atom_type=CarbonAtomType)
    <class 'mast.molecule.Atom'>

    or

    >>> CarbonAtomType.to_atom(coords)
    <class 'mast.molecule.Atom'>
    """

    def __init__(self, coords=None, atom_array=None, array_idx=None, atom_type=None, flags=None):

        if coords is None:
            coords = np.array([np.nan, np.nan, np.nan])
        else:
            assert coords.shape[-1] == mastmolconfig.DIM_NUM_3D, \
                "coords must have 3-dimensions, not {}".format(
                    coords.shape[-1])

        assert atom_type is not None, \
            "An atom_type must be given."
        assert isinstance(atom_type, AtomType), \
            "atom_type must be a subclass of AtomType, not {}".format(atom_type)

        if atom_array:
            assert atom_array.shape[-1] == mastmolconfig.DIM_NUM_3D, \
                "atom_array must have 3-dimensions, not {}".format(
                    atom_array.shape[-1])


        super().__init__(coords=coords, coord_array=atom_array,
                         array_idx=array_idx, flags=flags)
        self._atom_type = atom_type


    @property
    def atom_type(self):
        """The AtomType subclass used to substantiate this Atom."""
        return self._atom_type

    @property
    def isin_molecule(self):
        """Boolean if this Atom is in a Molecule."""
        return 'molecule' in self.flags

    @property
    def molecule(self):
        """The Molecule this Atom is in else None."""
        if not self.isin_molecule:
            return None
        else:
            # to get the selection in the registry that contains this
            # SelectionMember search through them all and take the
            # first one that is a Molecule type
            molecule = next((sel for key, sel in self.registry
                             if isinstance(sel, Molecule)),
                            None)
            assert molecule
            return molecule

    @property
    def isin_system(self):
        """Boolean if this Atom is in a System."""
        return 'system' in self.flags

    @property
    def system(self):
        """The System this Atom is in else None."""
        from mast.system import System
        if not self.isin_system:
            return None
        else:
            # the atom may be alone in the system or in a molecule
            # if it is alone the system should select it directly
            # to get the selection in the registry that contains this
            # SelectionMember search through them all and take the
            # first one that is a System type
            if self.isin_molecule:
                system = next((sel for key, sel in self.molecule.registry
                               if isinstance(sel, System)),
                              None)
                assert system
                return system
            else:
                system = next((sel for key, sel in self.registry
                               if isinstance(sel, System)),
                              None)
                assert system
                return system

    @property
    def isin_bond(self):
        """Boolean if this is in a Bond."""
        return 'bond' in self.flags

    @property
    def bonds(self):
        """The Bonds this atom is in else None."""
        if not self.isin_bond:
            return None
        else:
            # go through and collect each bond it is a part of
            bonds = []
            for key, sel in self.registry:
                if isinstance(sel, Bond):
                    if not isinstance(sel, Molecule):
                        bonds.append(sel)

            assert bonds
            return bonds


    @property
    def adjacent_atoms(self):
        """The Atoms that are bonded to this Atom."""
        if self.isin_bond:
            # collect adjacent atoms
            adjacent_atoms = []
            for bond in self.bonds:
                other_atom = next((a for a in bond.atoms if a is not self), None)
                assert self is not other_atom, "{} cannot be bound to itself".format(self)
                adjacent_atoms.append(other_atom)
            return adjacent_atoms
        else:
            return None

    @property
    def record(self):
        # define the Record named tuple
        record_fields = ['AtomType', 'x', 'y', 'z']
        AtomRecord = col.namedtuple('AtomRecord', record_fields)
        # build the values for it for this Type
        record_attr = {'AtomType' : self.atom_type.name}
        record_attr['x'] = self.coords[0]
        record_attr['y'] = self.coords[1]
        record_attr['z'] = self.coords[2]
        # make and return
        return AtomRecord(**record_attr)

    def distance(self, other):
        if issubclass(type(other), Point):
            assert self.system is other.system, \
                "the two entities must be in the same system"

        return super().distance(other)

    def atoms_within_distance(self, distance, metric='euclidean'):
        """Return atoms in this atom's system that are within a distance."""
        from scipy.spatial.distance import cdist

        dists = cdist([self.coords], [atom.coords for atom in
                                      self.system.all_atoms], metric=metric)[0]
        close_bool_arr = (dists < distance) & (dists > 0.0)
        close_idxs = close_bool_arr.nonzero()[0]
        atoms = [atom for i, atom in enumerate(self.system.all_atoms)
                 if i in close_idxs]
        close_dists = dists[close_bool_arr]

        return (atoms, close_dists)

    def molecule_within_distance(self, distance, metric='euclidean'):
        pass

class Bond(IndexedSelection):
    """The coordinate substantiation of a BondType.

    Examples
    --------

    >>> carbon_attributes = {'element':'C', 'bond_degree':3}
    >>> oxygen_attributes = {'element':'O', 'bond_degree':3}
    >>> COCarbonAtomType = AtomType.factory("COCarbonAtomType", **carbon_attributes)
    >>> COOxygenAtomType = AtomType.factory("COOxygenAtomType", **oxygen_attributes)
    >>> CO_atoms = (COCarbonAtomType, COOxygenAtomType)
    >>> CO_attributes = {"bond_order":3}
    >>> COBondType = BondType.factory("COBondType", atom_types=CO_atoms, **CO_attributes)
    """
    def __init__(self, atom_container=None, atom_ids=None,
                 bond_array=None, array_idxs=None, bond_type=None):

        # check the bond_type
        assert bond_type is not None, \
            "A BondType subclass must be given"
        assert isinstance(bond_type, BondType), \
            "bond_type must be a subclass of BondType, not {}".format(bond_type)

        # if the atoms are passed in on their own, i.e. with coordinates already
        if atom_container is not None:
            assert issubclass(type(atom_container), col.Sequence), \
                "atom_container must be a subclass of collections.Sequence, not {}".format(
                    type(atom_container))
            assert len(atom_container) >= 2, \
                "atom_container must have at least 2 atoms, not {}".format(len(atom_container))

            if len(atom_container) == 2:
                atom_ids = (0, 1)
            elif atom_ids is not None:
                assert isinstance(atom_ids, tuple), \
                    "atom_ids must be a length 2 tuple, not type{}".format(
                        type(atom_ids))
                assert len(atom_ids) == 2, \
                    "atom_ids must be a length 2 tuple, not len={}".format(
                        len(atom_ids))
                assert all([(lambda x: isinstance(x, int))(i) for i in atom_ids]), \
                    "atom_ids must be a length 2 tuple of ints"

            # if no atoms are passed but an array and indices
            elif bond_array is not None and array_idxs is not None:
                raise NotImplementedError


        super().__init__(atom_container, atom_ids, flags=['bond'])
        self._bond_type = bond_type

    @property
    def atoms(self):
        """The Atoms in this Bond."""
        return tuple(self.values())

    @property
    def atom_types(self):
        """The AtomType subclasses of this Bond's Atoms."""
        return tuple([atom.atom_type for atom in self.atoms])

    @property
    def coords(self):
        """The coordinates of the two Atoms in this Bond."""
        atom1_coords = self.atoms[0].coords
        atom2_coords = self.atoms[1].coords
        return (atom1_coords, atom2_coords)

    @property
    def bond_type(self):
        """The BondType subclass of this Bond."""
        return self._bond_type

    @property
    def isin_molecule(self):
        """Boolean if this Bond is in a Molecule."""
        return 'molecule' in self.flags

    @property
    def molecule(self):
        """The Molecule this Bond is in else None."""
        if self.isin_molecule is False:
            return None
        else:
            # to get the selection in the registry that contains this
            # SelectionMember search through them all and take the
            # first one that is a Molecule type
            molecule = next((sel for key, sel in self.registry
                             if isinstance(sel, Molecule)),
                            None)
            assert molecule
            return molecule

    @property
    def isin_system(self):
        """Boolean if this Bond is in a System."""
        return 'system' in self.flags

    @property
    def system(self):
        """The System this Bond is in else None."""
        from mast.system import System
        if self.isin_system is False:
            return None
        else:
            # to get the selection in the registry that contains this
            # SelectionMember search through them all and take the
            # first one that is a System type
            system = next((sel for key, sel in self.registry
                           if isinstance(sel, System)),
                          None)
            assert system
            return system

    @property
    def record(self):
        # define the Record named tuple
        record_fields = ['BondType',
                         'AtomAType', 'atom_A_idx',
                         'x_A', 'y_A', 'z_A',
                         'AtomBType', 'atom_B_idx',
                         'x_B', 'y_B', 'z_B']

        BondRecord = col.namedtuple('BondRecord', record_fields)

        # build the values for it for this Type
        atom_idxs = list(self.data.keys())
        record_attr = {'BondType' : self.bond_type.name}
        # first atom
        record_attr['AtomAType'] = self.atom_types[0].name
        record_attr['atom_A_idx'] = atom_idxs[0]
        record_attr['x_A'] = self.atoms[0].coords[0]
        record_attr['y_A'] = self.atoms[0].coords[1]
        record_attr['z_A'] = self.atoms[0].coords[2]
        # second atom
        record_attr['AtomBType'] = self.atom_types[1].name
        record_attr['atom_B_idx'] = atom_idxs[1]
        record_attr['x_B'] = self.atoms[1].coords[0]
        record_attr['y_B'] = self.atoms[1].coords[1]
        record_attr['z_B'] = self.atoms[1].coords[2]

        return BondRecord(**record_attr)

class Molecule(SelectionsDict):
    """The coordinate substantiation of a MoleculeType.
    A Molecule has a MoleculeType subtype, coordinates for each atom,
    and may be part of coordinate systems.

    The easiest way to obtain a Molecule object is to use the
    mast.molecule.MoleculeType.to_molecule(coords) function.

    Examples
    --------

    Go through the steps to make a MoleculeType
    >>> carbon_attributes = {'element':'C', 'bond_degree':3}
    >>> oxygen_attributes = {'element':'O', 'bond_degree':3}
    >>> COCarbonAtomType = AtomType("COCarbonAtomType", **carbon_attributes)
    >>> COOxygenAtomType = AtomType("COOxygenAtomType", **oxygen_attributes)
    >>> CO_atoms = (COCarbonAtomType, COOxygenAtomType)
    >>> CO_attributes = {"bond_order":3}
    >>> COBondType = BondType("COBondType", atom_types=CO_atoms, **CO_attributes)
    >>> atom_types = [COCarbonAtomType, COOxygenAtomType]
    >>> bond_types = [COBondType]
    >>> bond_map = {0 : (0, 1)}
    >>> CO_attributes = {"name" : "carbon-monoxide", "toxic" : True}
    >>> COType = MoleculeType.factory("COType", atom_types=atom_types, bond_types=bond_types, bond_map=bond_map, **CO_attributes)

    And then just make some coordinates and use them:
    >>> C_coords = np.array((0.0, 0.0, 0.0))
    >>> O_coords = np.array((0.0, 0.0, 1.0))
    >>> coords = np.array([C_coords, O_coords])

    >>> COType.to_molecule(coords)
    <class 'mast.molecule.Molecule'>


    """

    def __init__(self, atoms=None, bonds=None, angles=None, mol_type=None):

        assert mol_type is not None, \
            "A MoleculeType subclass must be given"

        # we need either bonds or atoms or both but not either to make
        # a molecule
        assert not (atoms is None and bonds is None)

        # check atoms input for correctness
        if atoms is not None:
            # check that the atoms are correct inputs
            assert atoms, "atoms must exist, {}".format(atoms)
            assert issubclass(type(atoms), col.Sequence), \
                "atoms must be a subclass of collections.Sequence, not {}".format(
                    type(atoms))
            assert all([(lambda x: True if issubclass(type(x), Atom) else False)(atom)
                        for atom in atoms]), \
                "all elements in atoms must be a subclass of type Atom"
            assert not all([atom.isin_molecule for atom in atoms]), \
                "all atoms must not be part of another molecule"

        # check bonds input for correctness
        if bonds is not None:
            # if bonds were given check that they are proper
            assert bonds, "atoms must exist, {}".format(bonds)
            assert issubclass(type(bonds), col.Sequence), \
                "bonds must be a subclass of collections.Sequence, not {}".format(
                    type(bonds))
            assert all([(lambda x: True if isinstance(x, Bond) else False)(bond)
                        for bond in bonds]), \
                "all elements in bonds must be a subclass of type Bond, not {}".format(
                    [type(bond) for bond in bonds])
            assert all(['molecule' not in bond.flags for bond in bonds]), \
                "all bonds must not be part of another molecule"

            # check to make sure all the bonds are connected
            # TODO

        # if both were given make sure that the atoms and bonds are
        # the same
        if atoms is not None and bonds is not None:
            # a reduce here might help?
            bonds_atoms = set()
            for bond in bonds:
                a = set([atom for atom in bond.atoms])
                bonds_atoms.update(a)

            assert (set(atoms).issuperset(bonds_atoms) and set(atoms).issubset(bonds_atoms)), \
                "all atoms in bonds must also be a part of the atoms collection"

        elif atoms is None and bonds is not None:
            # if there are only bonds collect the atoms
            raise NotImplementedError
        elif atoms is None and bonds is not None:
            # if there are only atoms infer the bonds naively
            raise NotImplementedError

        # call to parent class
        selections_dict = {'atoms' : atoms, 'bonds' : bonds}
        super().__init__(selection_dict=selections_dict, flags=['molecule'])

        # set the atoms, bonds, and angles into this object
        self.atoms = atoms
        self.bonds = bonds
        # self.angles = molecule_dict['angles']

        self._molecule_type = mol_type

        self._features = self.make_features()

        # attributes must explicitly be called due to computation time
        self._internal_interactions = None

    @property
    def atom_types(self):
        """The AtomType subclasses for all the atoms in this Molecule."""
        return [atom.atom_type for atom in self.atoms]

    @property
    def bond_types(self):
        """The BondType subclasses for all the bonds in this Molecule."""
        return [bond.bond_type for bond in self.bonds]

    @property
    def molecule_type(self):
        """The MoleculeType subclass of this Molecule."""
        return self._molecule_type

    @property
    def isin_system(self):
        """Boolean if this Molecule is in a System."""
        return 'system' in self.flags

    @property
    def system(self):
        """The System this Molecule is in else None."""
        from mast.system import System
        if self.isin_system is False:
            return None
        else:
            system = next((sel for key, sel in self.registry
                           if isinstance(sel, System)),
                          None)
            assert system
            return system

    @property
    def atom_coords(self):
        """The coordinates from each Atom in this Molecule."""
        coords = np.array([atom.coords for atom in self.atoms])
        return coords

    @property
    def atom_records(self):
        return [atom.record for atom in self.atoms]

    @property
    def atom_df(self):
        import pandas as pd
        return pd.DataFrame(self.atom_records)

    @property
    def bond_records(self):
        return [bond.record for bond in self.bonds]

    @property
    def bond_df(self):
        import pandas as pd
        return pd.DataFrame(self.bond_records)

    @property
    def feature_types(self):
        """The chemical features of this Molecule's MoleculeType subclass."""
        return self.molecule_type.feature_types

    def atoms_by_attr(self, attr_field, attr_value):
        return [atom for atom in self.atoms if
                atom.atom_type.attributes_data[attr_field] == attr_value]

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

    def distance_to(self, other):
        from scipy.spatial.distance import cdist

        assert self.system is other.system, \
            "entities must be in the same system"

        if isinstance(other, Atom):
            return cdist(self.atom_coords, [other])[0][0]
        elif isinstance(other, Molecule):
            return cdist(self.atom_coords, other.atom_coords)

    def atoms_within_distance(self, distance, atom_idxs=None, metric='euclidean'):
        from scipy.spatial.distance import cdist
        # if atom_idxs were given we only want to use these from the
        # molecule
        if atom_idxs:
            ref_atom_coords = [atom.coords for i, atom in
                               enumerate(self.atoms) if i in atom_idxs]
        # otherwise use the whole molecule
        else:
            ref_atom_coords = self.atom_coords

        all_other_close_atoms = []

        # free atoms
        if self.system.atoms:
            # get the atoms that are not this one
            other_atoms = [atom for atom in self.system.atoms if
                             atom is not self]
            # calculate distances between them
            dists = cdist(ref_atom_coords, [atom.coords for atom in
                                             other_atoms], metric=metric)
            # 1) get a boolean array of the distances that meet the
            # criteria (dists < distance)
            # 2) get the indices of those (nonzero)
            # 3) we only need the indices from axis 1, which was the
            # other atoms list, which we only want the uniqe ones of.
            other_atom_idxs = set((dists < distance).nonzero()[1])
            close_other_atoms = [atom for i, atom in
                                 enumerate(other_atoms) if i in other_atom_idxs]
            all_other_close_atoms.extend(close_other_atoms)

        # molecules
        if self.system.molecules:
            other_molecules = [molecule for molecule in self.system.molecules if
                             molecule is not self]

            for molecule in other_molecules:
                dists = cdist(ref_atom_coords, molecule.atom_coords,
                              metric=metric)
                # 1) get a boolean array of the distances that meet the
                # criteria (dists < distance)
                # 2) get the indices of those (nonzero)
                # 3) we only need the indices from axis 1, which was the
                # other molecule atoms list, which we only want the uniqe ones of.
                mol_atoms_close_idxs = set((dists < distance).nonzero()[1])
                mol_atoms_close = [atom for i, atom in
                                   enumerate(molecule.atoms) if i in
                                   mol_atoms_close_idxs]

                all_other_close_atoms.extend(mol_atoms_close)

        return all_other_close_atoms

    def molecules_within_distance(self, distance, distances=False, selection='any'):
        pass

    def make_features(self):
        """Using the features attribute make IndexedSelections of those
        Atoms. Sets these to the feature_family_selections and
        feature_type_selections attributes of this object.

        Examples
        --------

        > molecule.make_feature_selections()

        Each feature is held here
        > molecule.feature_family_selections[0]
        IndexedSelection
        The selection is the feature's atoms
        > molecule.features[0][0]
        Atom

        """

        features = {}
        for key, feature_type in self.molecule_type.feature_types.items():
            features[key] = feature_type.to_feature(self)

        return features


    @property
    def features(self):
        return self._features

    @property
    def feature_records(self):
        return [feature.record for feature in self.features.values()]

    @property
    def feature_df(self):
        """Export a pandas.DataFrame of the features."""
        import pandas as pd
        return pd.DataFrame(self.feature_records)

    @property
    def internal_interactions(self):
        """The intramolecular interactions of this Molecule."""
        return self._internal_interactions

    def profile_interactions(self, interaction_types):
        """Given different InteractionTypes profile this Molecule for those
        interactions which are set to the internal_interactions
        attribute.

        """
        from mast.interactions import InteractionType
        assert all([isinstance(itype, InteractionType) for itype in interaction_types]), \
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

    @property
    def record(self):
        pass

if __name__ == "__main__":
    pass
