""" The interactions module. """

from mast.datastructures import AssociationType
from mast.system import SystemAssociation

class Interaction(SystemAssociation):
    """Base class for associating Selections from a SelectionList with
information about an about the interaction.

    """

    def __init__(self, members=None, interaction_type=None):
        super().__init__(members=members)
        self._interaction_type = interaction_type

    @property
    def interaction_type(self):
        return self._interaction_type


