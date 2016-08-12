""" The interactions module. """
import itertools as it
from collections import defaultdict

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

    @staticmethod
    def factory(interaction_type_name,
                feature_types=None,
                **interaction_attrs):

        assert feature_types, "feature_types must be given."
        for feature_type in feature_types:
            assert issubclass(feature_type, FeatureType), \
                "All feature_type members must be a subclass of FeatureType, " \
                "not, {}".format(feature_type)

        
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

    @property
    def interaction_type(self):
        """The InteractionType this Interaction substantiates."""
        return self._interaction_type

if __name__ == "__main__":
    pass
