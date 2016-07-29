""" The system module. """
import collections as col

from mast.selection import SelectionList, IndexedSelection, \
    SelectionType, SelectionTypeLibrary
from mast.molecule import Atom, Molecule, MoleculeTypeLibrary, MoleculeType

__all__ = ['overlaps', 'SystemType', 'System', 'SystemAssociation']

def overlaps(members):
    """Check to see if members overlap.

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

class SystemType(SelectionType):
    """Base type for systems, subclasses should implement a to_system
method.

    """
    def __init__(self, system_attrs=None):
        super().__init__(attr_dict=system_attrs)


        # TODO move functionality from System to here

class System(SelectionList):
    """System that contains non-overlapping molecules, assumed to be in
the same coordinate system.

    molecules : the molecules in the system

    """

    def __init__(self, members=None, system_type=None):

        assert issubclass(type(members), col.Sequence), \
            "members must be a subclass of collections.Sequence, not {}".format(
                type(members))

        type_test_func = lambda x: True if (issubclass(type(x), Atom) or
                                            issubclass(type(x), Molecule)) \
                                            else False
        assert all([(type_test_func)(member) for member in members]), \
            "all elements in atoms must be a subclass of type Atom"

        # check to make sure none of the atoms are overlapping
        assert not overlaps(members), \
            "molecule system members cannot be overlapping"

        if system_type:
            assert issubclass(type(system_type), SystemType), \
                "system_type must be a subclass of SystemType, not {}".format(
                    type(system_type))

        super().__init__(selection_list=members)
        for member in members:
            member._in_system = True
            member._system = self
        self._system_type = system_type
        self._molecule_types = MoleculeTypeLibrary()
        for member in members:
            self._molecule_types.add_type(member.molecule_type, mol_name=member.molecule_type.name)
        self._system_associations = None

    def __repr__(self):
        return str(self.__class__)

    @property
    def system_type(self):
        return self._system_type

    # TODO move to SystemType
    @property
    def molecule_types(self):
        return self._molecule_types

    # TODO move to SystemType
    def add_molecule_type(self, mol_type, mol_name=None):
        if not mol_name:
            mol_name = mol_type.name
        self._molecule_types.add_type(mol_type, mol_name)

    @property
    def molecules(self):
        molecules = [member for  member in self if issubclass(type(member), Molecule)]
        return molecules

    # TODO should this be a property or function?
    @property
    def molecules_sel(self):
        mol_indices = [i for i, member in enumerate(self) if issubclass(type(member), Molecule)]
        return IndexedSelection(self, mol_indices)

    # TODO move to SystemType
    @property
    def associations(self):
        return self._system_associations
    # TODO move to SystemType
    def find_features(self):
        """Find features in all members of the system. Currently only molecules."""
        for mol_type in self.molecule_types.values():
            mol_type.find_features()

    def make_feature_selections(self):
        """Make feature selections for all current features in the system's molecules"""
        for mol in self.molecules:
            mol.make_feature_selections()

if __name__ == "__main__":
    pass
