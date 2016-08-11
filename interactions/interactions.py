""" The interactions module. """
import itertools as it
from collections import defaultdict, Counter

import numpy as np
import numpy.linalg as la

from mast.selection import SelectionsList
from mast.system import System
import mast.selection as mastsel
import mast.molecule as mastmol
import mast.system as mastsys

import mast.config.interactions as mastinxconfig

__all__ = ['AssociationType', 'Association',
           'Interaction', 'HydrogenBondInx', 'NoHHydrogenBondInx'
           'InteractionType', 'HydrogenBondType', 'NoHHydrogenBondType']

class InteractionError(Exception):
    pass

class InteractionType(object):
    """Class for generating specific interaction type classes with the factory
    method.

    Examples
    --------

    """

    def __init__(self):
        pass

    @classmethod
    def factory(cls, interaction_type_name,
                feature_types=None,
                **interaction_attrs):

        assert feature_types, "feature_types must be given."
        for feature_type in feature_types:
            assert issubclass(feature_type, FeatureType), \
                "All feature_type members must be a subclass of FeatureType, " \
                "not, {}".format(feature_type)

        # keep track of which attributes the input did not provide
        # compared to the config file
        for attr in cls.attributes:
            try:
                assert attr in interaction_attrs.keys()
            except AssertionError:
                pass
                # LOGGING
                # print("Attribute {0} not found in interaction input.".format(attr))

        # add the attributes into the class
        attributes = {attr : None for attr in cls.attributes}
        for attr, value in interaction_attrs.items():
            try:
                assert attr in cls.attributes
            # if it doesn't then log
            except AssertionError:
                # LOGGING
                pass
                # print("Input attribute {0} not in InteractionType attributes.".format(attr))
            # add it regardless
            attributes[attr] = value

        interaction_type = type(interaction_type_name, (cls,), attributes)
        interaction_type.name = attributes
        interaction_type.attributes_data = attributes
        interaction_type.interaction_type = cls
        interaction_type.feature_types = feature_types

        return interaction_type

    @classmethod
    def check(cls, *args, **kwargs):
        """The principle class method for testing for the existence of an
        interaction from specified geometric constraint parameters.
        All subclasses of InteractionType should implement this for
        domain specificity.

        """
        pass

    @classmethod
    def find_hits(cls, member_a, member_b):
        """Returns all the 'hits' for interactions of features between two
        members (molecules, selections, etc.), where a hit is a
        feature set that satisfies the geometric constraints defined
        in the InteractionType.check function.

        """
        pass


class Interaction(SelectionsList):
    """Substantiates the InteractionType class by containing Feature
    objects.

    SelectionsList container for multiple features which together
    constitute an intermolecular interaction, e.g. donor and acceptor
    features in a hydrogen bond. Specific implementations of
    Interactions should inherit from this class, e.g. HydrogenBondInx
    substantiates the HydrogenBondType class, where the Inx suffix is
    short for interaction.

    """
    def __init__(self, features=None, system=None, interaction_type=None):

        assert interaction_type, "interaction_type must be given"
        assert issubclass(interaction_type, InteractionType), \
            "interaction_type must be a subclass of mast.interactions.InteractionType"

        for feature in features:
            assert feature.system is system, \
                "feature's system must be all the same"

        super().__init__(selection_list=features)
        self._interaction_type = interaction_type
        self._interaction_class = None

    @property
    def features(self):
        return self.data

    @property
    def feature_types(self):
        return [feature.feature_type for feature in self.features]

    @property
    def interaction_type(self):
        """The InteractionType this Interaction substantiates."""
        return self._interaction_type

    @property
    def interaction_class(self):
        return self._interaction_class

    @interaction_class.setter
    def interaction_class(self, interaction_class):
        assert issubclass(interaction_class, self.interaction_type), \
            "interaction_classes must be subclasses of this Interaction's" \
            "interaction_type, not {}".format(interaction_class)
        assert Counter(self.feature_types) == Counter(interaction_class.feature_types), \
            "the interaction_class must have the same number and types of features" \
            "as this Interaction, not {}".format(Counter(interaction_class.feature_types))

        self._interaction_class = interaction_class

if __name__ == "__main__":
    pass
