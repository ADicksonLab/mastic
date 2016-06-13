""" The system module. """
from itertools import product, combinations

from mast.datastructures import TrackedMember, TrackedList, Selection, SelectionList, Association
from mast.molecule import Molecule, MoleculeList

__all__ = ['System', 'SystemAssociation']

def overlaps(members):
    """Check to see if members overlap.

    STUB: just check to make sure that no two molecules have any atoms
    with the same coordinates.

    """

    pairs = combinations(members, 2)
    try:
        pair = next(pairs)
    # if it is empty no overlaps
    except StopIteration:
        return False
    flag = True
    while flag:
        if pair[0].overlaps(pair[1]):
            return True
        else:
            try:
                pair = next(pairs)
            except StopIteration:
                flag = False
    return False

class System(MoleculeList):
    """System that contains non-overlapping molecules, assumed to be in
the same coordinate system.

    molecules : the molecules in the system

    """

    def __init__(self, members=None):

        try:
            iter(members)
        except TypeError:
            members = [members]

        # check to make sure none of the atoms are overlapping
        if overlaps(members):
            raise ValueError("molecule system members cannot be overlapping")

        super().__init__(members=members)

class SystemAssociation(Association):
    """ Class for associations which only occur within a System. """

    def __init__(self, members=None, association=None, system=None, idx=None, ids=None):

        if not members or set(members) <= set(system.members):
            super().__init__(members=members, association=association)
        else:
            raise ValueError("Members of a SystemAssociation must all be in the same system")

        self._system = system

    @property
    def system(self):
        return self._system
