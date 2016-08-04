""" The system module. """
import collections as col

from mast.selection import SelectionsList, IndexedSelection
from mast.molecule import Atom, Bond, Molecule, AtomType, BondType, MoleculeType
import mast.config.system as mastsysconfig

__all__ = ['overlaps', 'SystemType', 'System', 'SystemAssociation']


class SystemType(object):
    attributes = mastsysconfig.SYSTEM_ATTRIBUTES
    def __init__(self):
        pass

    @classmethod
    def to_system(cls, members_coords):
        members = []
        for member_idx, member_coords in enumerate(members_coords):
            # create each member using the coordinates
            member_type = cls.member_types[member_idx]
            if issubclass(member_type, AtomType):
                members.append(member_type.to_atom(member_coords))
            elif issubclass(member_type, MoleculeType):
                members.append(member_type.to_molecule(member_coords))

        system = System(members, system_type=cls)
        return system

    @classmethod
    def molecule_types(cls):
        return [member_type for member_type in cls.member_types if
                issubclass(member_type, MoleculeType)]

    @classmethod
    def atom_types(cls):
        return [member_type for member_type in cls.member_types if

                issubclass(member_type, AtomType)]

    @staticmethod
    def factory(system_type_name, member_types=None,
                **system_attrs):
        assert member_types, "molecule_types must be provided"
        for member_type in member_types:
            assert (issubclass(member_type, MoleculeType) or
                    issubclass(member_type, AtomType)), \
                    "molecule_types must contain only MoleculeType or"\
                    " AtomType subclasses, not {}".format(
                        type(member_type))

        # keep track of which attributes the input did not provide
        # compared to the config file
        for attr in SystemType.attributes:
            try:
                assert attr in system_attrs.keys()
            except AssertionError:
                pass
                # LOGGING
                # print("Attribute {0} not found in system input.".format(attr))

        # add the attributes into the class
        attributes = {attr : None for attr in SystemType.attributes}
        for attr, value in system_attrs.items():
            try:
                assert attr in SystemType.attributes
            # if it doesn't then log
            except AssertionError:
                # LOGGING
                pass
                # print("Input attribute {0} not in SystemType attributes.".format(attr))
            # add it regardless
            attributes[attr] = value

        system_type = type(system_type_name, (SystemType,), attributes)
        system_type.member_types = member_types
        system_type.member_type_library = set(member_types)
        system_type.association_types = None

        return system_type


class System(SelectionsList):
    """System that contains non-overlapping molecules, assumed to be in
the same coordinate system.

    molecules : the molecules in the system

    """

    def __init__(self, members=None, system_type=None):

        assert issubclass(type(members), col.Sequence), \
            "members must be a subclass of collections.Sequence, not {}".format(
                type(members))

        for member in members:
            assert (issubclass(type(member), Atom) or \
                    issubclass(type(member), Molecule)), \
                    "all elements must be a subclass of type Atom or Molecule, not {}".format(
                        type(member))

        # check to make sure none of the atoms are overlapping
        assert not overlaps(members), \
            "molecule system members cannot be overlapping"

        if system_type:
            assert issubclass(system_type, SystemType), \
                "system_type must be a subclass of SystemType, not {}".format(
                    type(system_type))

        super().__init__(selection_list=members, flags=['system'])
        self._system_type = system_type
        self._associations = None

    def __repr__(self):
        return str(self.__class__)

    @property
    def system_type(self):
        return self._system_type

    @property
    def atom_types(self):
        return self.system_type.atom_types()

    @property
    def molecule_types(self):
        return self.system_type.molecule_types()

    @property
    def atoms(self):
        atoms = [member for  member in self if issubclass(type(member), Atom)]
        return atoms

    @property
    def molecules(self):
        molecules = [member for  member in self if issubclass(type(member), Molecule)]
        return molecules

    # is this even relevant?
    def molecules_sel(self):
        mol_indices = [i for i, member in enumerate(self) if issubclass(type(member), Molecule)]
        return IndexedSelection(self, mol_indices)

    @property
    def associations(self):
        return self._associations

    @property
    def associations_types(self):
        return self.system_type.association_types

    def make_feature_selections(self):
        """Make feature selections for all current features in the system's molecules"""
        for mol in self.molecules:
            mol.make_feature_selections()

    def overlaps(members):
        """Checks whether the members given overlap anything in this system."""
        pass

if __name__ == "__main__":
    pass
