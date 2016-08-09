"""Features of molecules are defined by substructures of molecules
that match a signature in a database and are known to have some
physical property.

"""

import mast.config.features as mastfeatconfig

import mast.molecule as mastmol
import mast.selection as mastsel

class FeatureType(object):

    attributes = mastfeatconfig.FEATURE_ATTRIBUTES

    def __init__(self):
        pass

    @classmethod
    def to_feature(cls, molecule):
        return Feature(molecule=molecule, feature_type=cls)

    @classmethod
    def factory(cls, feature_type_name,
                molecule_type=None,
                atom_idxs=None, bond_idxs=None,
                **feature_attrs):

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
        # add the attributes as a dict
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
        return [atom_type for i, atom_type in
                enumerate(cls.molecule_type.atom_types) if i in cls.atom_idxs]

    @classmethod
    def bond_types(cls):
        return [bond_type for i, bond_type in
                enumerate(cls.molecule_type.bond_types) if i in cls.bond_idxs]


class Feature(mastsel.SelectionsDict):

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
        return self['atoms']

    @property
    def bonds(self):
        return self['bonds']

    @property
    def feature_type(self):
        return self._feature_type

    @property
    def molecule(self):
        return self._molecule
