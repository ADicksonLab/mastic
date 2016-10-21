""" The interactions module. """
import itertools as it
from collections import defaultdict, Counter, namedtuple

import numpy as np
import numpy.linalg as la

from mast.selection import SelectionsList
from mast.system import System
import mast.selection as mastsel
import mast.molecule as mastmol
import mast.features as mastfeat
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

    interaction_name = "None"
    feature_keywords = None
    feature_order = None
    grouping_attribute = None
    # order is the number of features that participate in an interaction
    degree = 0

    # a string to use for formatting interaction class names
    inx_class_name_template = "{assoc_type}_{inx_name}_{idx}_InxClass"
    # This is used in interaction_classes for dynamically assigning
    # names to interaction classes in associations and defining the
    # interaction space of that association. Avoid making interaction
    # class names with this form unless you don't plan on using
    # interaction spaces

    def __init__(self, interaction_type_name,
                feature_types=None,
                association_type=None,
                assoc_member_pair_idxs=None,
                **interaction_attrs):

        assert feature_types, "feature_types must be given."
        for feature_type in feature_types:
            assert isinstance(feature_type, mastfeat.FeatureType), \
                "All feature_type members must be a subclass of FeatureType, " \
                "not, {}".format(feature_type)
        # keep track of which attributes the input did not provide
        # compared to the config file
        for attr in self.attributes:
            try:
                assert attr in interaction_attrs.keys()
            except AssertionError:
                pass
                # LOGGING
                # print("Attribute {0} not found in interaction input.".format(attr))

        # add the attributes into the class
        attributes = {attr : None for attr in self.attributes}
        for attr, value in interaction_attrs.items():
            try:
                assert attr in self.attributes
            # if it doesn't then log
            except AssertionError:
                # LOGGING
                pass
                # print("Input attribute {0} not in InteractionType attributes.".format(attr))
            # add it regardless
            attributes[attr] = value

        self.name = interaction_type_name
        self.attributes_data = attributes
        self.__dict__.update(attributes)
        self._feature_types = feature_types
        self.association_type = association_type
        self.assoc_member_pair_idxs = assoc_member_pair_idxs

    def __eq__(self, other):
        if not isinstance(other, type(self)):
            return False
        elif self.name != other.name:
            return False
        else:
            return True

    def __hash__(self):
        return self.name.__hash__()


    @property
    def feature_types(self):
        return self._feature_types

    @classmethod
    def interaction_classes(cls, association_type,
                            inx_class_name_template=None):
        """Receives an association and creates all of the possible interaction
        classes in the association. Interaction classes are simply
        instantiations of this InteractionType corresponding to
        topological representations of an interaction.

        Notice: This nomenclature may change, as it might make more
        sense to have this class (InteractionType) be the interaction
        class and the instantions be the interaction types to keep it
        more similar to what molecule types etc. are.

        The default way of naming interaction classes is given in the
        inx_class_name_template variable. Otherwise you can pass in an
        iterable

        """
        inx_class_name_template = cls.inx_class_name_template
        # for each member collect the features relevant to this
        # interaction type, in the feature_order, so initialize an
        # empty list for each
        members_features = [[] for feature_key in feature_order]
        # and go through each member_type
        for member_idx, member_type in enumerate(association_type.member_types):
            for feature_type in member_type.feature_types.values():
               # if the feature has one of the features in this interaction
                attr = cls.feature_inx_attribute(feature_type)
                if attr == cls.feature_order[member_idx]:
                    # and add it to the appropriate list
                    members_features[member_idx].append(feature_type)

        # for each of these combinations take the product of
        # compatible features
        feature_pairs = it.product(*members_features)

        # now that we have the pairs we will make interaction classes
        # for each of them
        inx_classes = []
        for inx_class_idx, feature_pair in enumerate(feature_pairs):
            # the name of the inx class
            inx_class_name = inx_class_name_template.format(
                assoc_type=association_type.name,
                inx_name=cls.interaction_name,
                idx=inx_class_idx)

            # stub
            inx_class_attributes = {}

            # create the interaction class
            inx_class = cls(inx_class_name,
                            feature_types=feature_pair,
                            association_type=association_type,
                            assoc_member_pair_idxs=association_type.member_idxs,
                            **inx_class_attributes)

            inx_classes.append(inx_class)

        return inx_classes

    @classmethod
    def feature_inx_attribute(cls, feature):
        """Check to see if this feature is part of the classes grouping
        attribute. Returns the feature attribute or None if it does
        not match.

        """
        feature_attribute = None
        for feature_key in cls.feature_keywords[cls.grouping_attribute]:
            if feature.attributes_data[cls.grouping_attribute] == feature_key:
                feature_attribute = feature_key

        return feature_attribute

    def check(self, *args, **kwargs):
        """The principle class method for testing for the existence of an
        interaction from specified geometric constraint parameters.
        All subclasses of InteractionType should implement this for
        domain specificity.

        """
        raise NotImplementedError

    def find_hits(self, member_a, member_b):
        """Returns all the 'hits' for interactions of features between two
        members (molecules, selections, etc.), where a hit is a
        feature set that satisfies the geometric constraints defined
        in the InteractionType.check function.

        """
        raise NotImplementedError


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
            "interaction_type must be a subclass of " \
            "mast.interactions.InteractionType, not {}".format(
                interaction_type)

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
        assert isinstance(interaction_class, self.interaction_type), \
            "interaction_classes must be subclasses of this Interaction's" \
            "interaction_type, not {}".format(interaction_class)
        assert Counter(self.feature_types) == Counter(interaction_class.feature_types), \
            "the interaction_class must have the same number and types of features" \
            "as this Interaction, not {}".format(Counter(interaction_class.feature_types))
        self._interaction_class = interaction_class


if __name__ == "__main__":
    pass
