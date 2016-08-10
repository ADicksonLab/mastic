"""Features of molecules are defined by substructures of molecules
that match a signature in a database and are known to have some
physical property.

This class defines a FeatureType class for deriving new types of
features through either manual subclassing or through the factory
method for dynamic generation.

"""
import collections as col

import mast.config.features as mastfeatconfig

import mast.molecule as mastmol
import mast.selection as mastsel

class FeatureType(object):
    """Class for generating specific system type classes with the factory
    method.

    Examples
    --------


    """

    attributes = mastfeatconfig.FEATURE_ATTRIBUTES

    def __init__(self):
        pass

    @classmethod
    def to_feature(cls, molecule):
        """Substantiate a Feature by specifying a substantiated Molecule from
        which to make selections from.

        """
        return Feature(molecule=molecule, feature_type=cls)

    @classmethod
    def factory(cls, feature_type_name,
                molecule_type=None,
                atom_idxs=None, bond_idxs=None,
                **feature_attrs):
        """Static method for generating feature types dynamically given a type
        name (which will be the class name) and a domain specific dictionary
        of feature attributes.

        See mast.config.features for standard FeatureType attributes.
        See class docstring for examples.
        """

        # validate input
        assert molecule_type, "molecule_type must be given"
        assert issubclass(molecule_type, mastmol.MoleculeType), \
            "molecule_type must be a subclass of MoleculeType, not {}".format(
                molecule_type)

        assert atom_idxs or bond_idxs, \
            "At least one atom or bond must be specified to specify a feature"
        # TODO check to make sure the selections are valid
        if atom_idxs:
            pass
        if bond_idxs:
            pass

        # keep track of which attributes the input did not provide
        for attr in FeatureType.attributes:
            try:
                assert attr in feature_attrs.keys()
            except AssertionError:
                # LOGGING
                pass
                # print("Attribute {0} not found in feature input.".format(attr))

        # add the attributes into the class
        attributes = {attr : None for attr in FeatureType.attributes}
        for attr, value in feature_attrs.items():
            # Log the compliance of the attributes
            # if the attribute isn't declared in the attributes log it
            try:
                assert attr in FeatureType.attributes
            except AssertionError:
                # LOGGING
                pass
                # print("Input attribute {0} not in FeatureType attributes.".format(attr))

            # add it to the attributes
            attributes[attr] = value

        feature_type = type(feature_type_name, (FeatureType,), attributes)
        feature_type.name = feature_type_name
        feature_type.attributes_data = attributes
        # add core attributes
        feature_type.molecule_type = molecule_type
        if atom_idxs:
            feature_type.atom_idxs = atom_idxs
        else:
            feature_type.atom_idxs = []
        if bond_idxs:
            feature_type.bond_idxs = bond_idxs
        else:
            feature_type.bond_idxs = []


        # feature_type.angle_idxs = angle_idxs

        return feature_type

    @classmethod
    def atom_types(cls):
        """The AtomTypes in the FeatureType."""
        return [atom_type for i, atom_type in
                enumerate(cls.molecule_type.atom_types) if i in cls.atom_idxs]

    @classmethod
    def bond_types(cls):
        """The BondTypes in the FeatureType."""
        return [bond_type for i, bond_type in
                enumerate(cls.molecule_type.bond_types) if i in cls.bond_idxs]

    @classmethod
    def member_types(cls):
        return atom_types + bond_types

    @classmethod
    def record(cls):
        # define the Record namedtuple
        record_fields = ['FeatureType', 'MoleculeType',
                         'AtomTypes', 'atom_idxs',
                         'BondTypes', 'bond_idxs'] + \
                         list(cls.attributes_data.keys())

        FeatureTypeRecord = col.namedtuple('FeatureTypeRecord', record_fields)
        # build the values for it for this Type
        record_attr = {'FeatureType' : cls.name}
        record_attr['MoleculeType'] = cls.molecule_type().name
        record_attr['AtomTypes'] = cls.atom_types()
        record_attr['atom_idxs'] = cls.atom_idxs
        record_attr['BondTypes'] = cls.bond_types()
        record_attr['bond_idxs'] = cls.bond_idxs
        record_attr.update(cls.attributes_data)

        return FeatureTypeRecord(**record_attr)

class Feature(mastsel.SelectionsDict):
    """Feature, which is a collection of Atoms and Bonds selected from a
    single Molecule and are described by a domain specific set of
    attributes.

    Examples
    --------

    """
    def __init__(self, molecule=None, feature_type=None):

        assert isinstance(molecule, mastmol.Molecule), \
            "molecule must be a mast.molecule.Molecule instance, not {}".format(
                type(molecule))

        assert issubclass(feature_type, FeatureType), \
            "feature_type must be a subclass of FeatureType, not {}".format(
                feature_type)

        # make selections on the Molecule
        selections = {'atoms' : [molecule.atoms[i] for i in feature_type.atom_idxs],
                      'bonds' : [molecule.bonds[i] for i in feature_type.bond_idxs]}

        super().__init__(selection_dict=selections, flags=['feature'])
        self._feature_type = feature_type
        self._molecule = molecule

    @property
    def atoms(self):
        """Atoms in the feature"""
        return self['atoms']

    @property
    def bonds(self):
        """Bonds in the feature"""
        return self['bonds']

    @property
    def feature_type(self):
        """FeatureType this Feature substantiated."""
        return self._feature_type

    @property
    def molecule(self):
        """Molecule this Feature makes selections of."""
        return self._molecule

    @property
    def system(self):
        """The system this Feature's molecule is in."""
        return self.molecule.system

    @property
    def record(self):
        # would be too much to have all the coordinates
        return self.feature_type.record()
