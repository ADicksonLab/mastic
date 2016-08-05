""" The system module. """
import collections as col

from mast.selection import SelectionsList, IndexedSelection
from mast.molecule import Atom, Bond, Molecule, AtomType, BondType, MoleculeType


import mast.config.system as mastsysconfig

__all__ = ['overlaps', 'SystemType', 'System']

def overlaps(members):
    """Check to see if any iterable of substantiated members' coordinates
    overlap.

    """

    from itertools import combinations
    pairs = combinations(members, 2)
    try:
        pair = next(pairs)
    # if it is empty no overlaps
    except StopIteration:
        return False
    flag = True
    while flag:
        overlaps = pair[0].overlaps(pair[1])
        if overlaps:
            return overlaps
        else:
            try:
                pair = next(pairs)
            except StopIteration:
                flag = False
    return False


class SystemType(object):
    """Class for generating specific system type classes with the factory
    method.

    Examples
    --------

    Build some type to put into the system:

    >>> carbon_attributes = {'element':'C', 'bond_degree':3}
    >>> oxygen_attributes = {'element':'O', 'bond_degree':3}
    >>> COCarbonAtomType = AtomType.factory("COCarbonAtomType", **carbon_attributes)
    >>> COOxygenAtomType = AtomType.factory("COOxygenAtomType", **oxygen_attributes)
    >>> CO_atoms = (COCarbonAtomType, COOxygenAtomType)
    >>> CO_attributes = {"bond_order":3}
    >>> COBondType = BondType.factory("COBondType", atom_types=CO_atoms, **CO_attributes)
    >>> atom_types = [COCarbonAtomType, COOxygenAtomType]
    >>> bond_types = [COBondType]
    >>> bond_map = {0 : (0, 1)}
    >>> CO_attributes = {"name" : "carbon-monoxide", "toxic" : True}
    >>> COMoleculeType = MoleculeType.factory("COType", atom_types=atom_types, bond_types=bond_types, bond_map=bond_map, **CO_attributes)

    Make a SystemType that contains one COMoleculeType:

    >>> system_attrs = {'name' : 'carbon-monoxide-system'}
    >>> COSystemType = SystemType.factory("COSystemType", member_types=[COMoleculeType], **system_attrs)

    """
    attributes = mastsysconfig.SYSTEM_ATTRIBUTES
    def __init__(self):
        pass

    @classmethod
    def to_system(cls, members_coords):
        """Substantiate a System using input coordinates in the order of the
        members in the system.

        """
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
        """The MoleculeTypes of all Molecule system members."""
        return [member_type for member_type in cls.member_types if
                issubclass(member_type, MoleculeType)]

    @classmethod
    def atom_types(cls):
        """The AtomTypes of all Atom system members."""
        return [member_type for member_type in cls.member_types if

                issubclass(member_type, AtomType)]

    @classmethod
    def make_member_association_type(self, member_idxs, association_type=None):
        """Match an AssociationType to members of the SystemType"""
        raise NotImplementedError

    @classmethod
    def association_types(cls):
        return cls._association_types

    @classmethod
    def add_association_type(cls, association_type):
        # check to make sure that it's selection types are in this
        # SystemType
        from mast.interactions import AssociationType
        assert issubclass(association_type, AssociationType), \
            "association_type must be a subclass of mast.interactions.Association,"\
            " not {}".format(association_type)

        # there must be one member_type for each member_type in the association_type
        sys_member_count = {member_type :  cls.member_types.count(member_type)
                            for member_type in cls.member_types}
        assoc_member_count = {member_type : association_type.member_types.count(member_type)
                              for member_type in association_type.member_types}
        for member_type, assoc_count in assoc_member_count.items():
            assert assoc_count <= sys_member_count[member_type], \
                "There must be greater than or equal the number of {0}"\
                " in the system {1} as are in the association {2}".format(
                    member_type, cls, association_type)

        cls._association_types.append(association_type)

    @staticmethod
    def factory(system_type_name, member_types=None,
                **system_attrs):
        """Static method for generating system types dynamically given a type
        name (which will be the class name) and a domain specific dictionary
        of system attributes.

        See mast.config.molecule for standard SystemType attributes.
        See class docstring for examples.
        """

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
        system_type._association_types = []

        return system_type


class System(SelectionsList):
    """System that contains non-overlapping molecules, assumed to be in
the same coordinate system.

    molecules : the molecules in the system

    >>> import numpy as np
    >>> carbon_attributes = {'element':'C', 'bond_degree':3}
    >>> oxygen_attributes = {'element':'O', 'bond_degree':3}
    >>> COCarbonAtomType = AtomType.factory("COCarbonAtomType", **carbon_attributes)
    >>> COOxygenAtomType = AtomType.factory("COOxygenAtomType", **oxygen_attributes)
    >>> CO_atoms = (COCarbonAtomType, COOxygenAtomType)
    >>> CO_attributes = {"bond_order":3}
    >>> COBondType = BondType.factory("COBondType", atom_types=CO_atoms, **CO_attributes)
    >>> atom_types = [COCarbonAtomType, COOxygenAtomType]
    >>> bond_types = [COBondType]
    >>> bond_map = {0 : (0, 1)}
    >>> CO_attributes = {"name" : "carbon-monoxide", "toxic" : True}
    >>> COMoleculeType = MoleculeType.factory("COType", atom_types=atom_types, bond_types=bond_types, bond_map=bond_map, **CO_attributes)
    >>> system_attrs = {'name' : 'carbon-monoxide-system'}
    >>> COSystemType = SystemType.factory("COSystemType", member_types=[COMoleculeType], **system_attrs)

    Get coordinates for the things in the system
    >>> CO_coords = np.array([[0.0, 0.0, 0.0], [0.0, 0.0, 1.0]])
    >>> member_coords = [CO_coords]

    Substantiate a system from those:
    >>> COSystemType.to_system(member_coords)
    <class 'mast.system.System'>

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
        # substantiate the Associations in this System
        self._associations = []
        for association_type in self._system_type.association_types():
            self._associations.append(association_type.to_association(self))

    def __repr__(self):
        return str(self.__class__)

    @property
    def system_type(self):
        """The SystemType that substantiated this System."""
        return self._system_type

    @property
    def members(self):
        """The system members, of every type."""
        return self.data

    @property
    def atom_types(self):
        """The AtomTypes of every Atom system member."""
        return self.system_type.atom_types()

    @property
    def molecule_types(self):
        """The MoleculeTypes of every Molecule system member"""
        return self.system_type.molecule_types()

    @property
    def atoms(self):
        """The Atom system members."""
        atoms = [member for  member in self if issubclass(type(member), Atom)]
        return atoms

    @property
    def molecules(self):
        """The Molecule system members."""
        molecules = [member for  member in self if issubclass(type(member), Molecule)]
        return molecules

    # YAGNI?
    def molecules_sel(self):
        """Returns a selection on the system of just the Molecule system
        members.

        """
        mol_indices = [i for i, member in enumerate(self) if issubclass(type(member), Molecule)]
        return IndexedSelection(self, mol_indices)

    # YAGNI?
    def atoms_sel(self):
        """Returns a selection on the system of just the Atom system
        members.

        """
        atom_indices = [i for i, member in enumerate(self) if issubclass(type(member), Atom)]
        return IndexedSelection(self, atom_indices)

    @property
    def associations(self):
        return self._associations

    @property
    def associations_types(self):
        return self.system_type.association_types

    def make_feature_selections(self):
        """Make feature selections for all current features in the system's
        molecules.

        """
        for mol in self.molecules:
            mol.make_feature_selections()

    def overlaps(members):
        """Checks whether the members given overlap anything in this system."""
        for member in members:
            for sys_member in self.members:
                if overlaps([member, sys_member]):
                    return True
        # if none overlap
        return False


if __name__ == "__main__":
    pass
