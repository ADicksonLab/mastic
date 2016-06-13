""" The interactions module. """

from mast.datastructures import AssociationType
from mast.system import SystemAssociation

__all__ = ['Interaction', 'HydrogenBondInx',
           'InteractionType', 'HydrogenBondType']

class InteractionType(AssociationType):
    """ Prototype class for all intermolecular interactions."""
    def __init__(self, description=None):
        super().__init__(description=description)

class HydrogenBondType(InteractionType):
    """ Class for checking validity of a HydrogenBondInx."""

    @classmethod
    def donor(self, donor):
        ok = True

        return ok

    @classmethod
    def acceptor(self, donor):
        ok = True

        return ok

    @classmethod
    def distance(self, donor):
        ok = True

        return ok

    @classmethod
    def angle(self, donor):
        ok = True

        return ok

class Interaction(SystemAssociation):
    """Base class for associating Selections from a SelectionList with
information about an about the interaction.

    """

    def __init__(self, members=None, interaction_type=None):
        super().__init__(members=members)
        if not interaction_type:
            interaction_type = InteractionType()

        if not issubclass(type(interaction_type), InteractionType):
            raise TypeError(
                "interaction_type must be a subclass of InteractionType, not {}".format(
                type(interaction_type)))
        self._interaction_type = interaction_type

    @property
    def interaction_type(self):
        return self._interaction_type

class HydrogenBondInx(Interaction):
    """Interaction subclass that has HydrogenBondType type.

    """

    def __init__(self, donor=None, acceptor=None):
        description = "Hydrogen-Bond Interaction"

        # TODO get the distances
        distance=None
        angle=None

        super().__init__()
        # type checking of these using HydrogenBondType
        self._interaction_type = HydrogenBondType

        try:
            if HydrogenBondType.donor(donor):
                self._donor = donor
            else:
                raise ValueError("Donor, {}, is not valid".format(donor))
            if HydrogenBondType.acceptor(acceptor):
                self._acceptor = acceptor
            else:
                raise ValueError("Acceptor, {}, is not valid".format(acceptor))
            if HydrogenBondType.distance(distance):
                self._distance = distance
            else:
                raise ValueError("distance, {}, is not valid".format(distance))
            if HydrogenBondType.angle(angle):
                self._angle = angle
            else:
                raise ValueError("Angle, {}, is not valid".format(angle))
        except ValueError:
            print("not a HydrogenBondType Interaction")

    @property
    def donor(self):
        return self._donor

    @property
    def acceptor(self):
        return self._acceptor

    @property
    def distance(self):
        return self._distance

    @property
    def angle(self):
        return self._angle
