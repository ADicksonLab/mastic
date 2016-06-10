""" The interactions module. """

from mast.datastructures import TrackedMember, TrackedList, Selection, SelectionList
from mast.molecule import Molecule, MoleculeList


class Interaction(SelectionList):
    """Base class for associating Selections from a SelectionList with
information about an about the interaction.

    """

    def __init__(self, members=None, interaction_type=None):
        super().__init__(members=members)
        self._interaction_type = interaction_type

    

