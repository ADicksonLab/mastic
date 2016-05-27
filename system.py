""" The system module. """

from mast.datastructures import TrackedMember, TrackedList, Selection, SelectionList
from mast.molecule import Molecule, MoleculeList

def overlap(members):
    """Check to see if members overlap.

    STUB: just check to make sure that no two molecules have any atoms
    with the same coordinates.

    """

    it1 = iter(members)
    it2 = iter(members)
    # if one is empty return False as they couldn't possibly overlap
    # TODO maybe should raise an error, SystemError
    try:
        member1 = next(it)
    except StopIteration:
        return False
    try:
        member2 = next(it)
    except StopIteration:
        return False

    flag = True
    while flag:
        while flag:
            
            if overlap([member1, member2]):
                return True
            else:
                try:
                    next(it2)
                except StopIteration:
                    flag = False
        try:
            next(it)
        except StopIteration:
            flag = False

        return True
    return False

class System(MoleculeList):
    """System that contains non-overlapping molecules, assumed to be in
the same coordinate system.

    molecules : the molecules in the system

    """

    def __init__(self, members=members):

        # check to make sure none of the atoms are overlapping
        if overlap(members):
            raise ValueError("molecule system members cannot be overlapping")

        super().__init__(members=members)


