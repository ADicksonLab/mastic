""" The interactions module. """
import itertools as it
import collections as col

from mast.selection import SelectionsList
import mast.features as mastfeat

# from mast.system import System
# import mast.selection as mastsel
# import mast.molecule as mastmol
# import mast.system as mastsys
# import mast.config.interactions as mastinxconfig
# import mast.config.features as mastfeatconfig

__all__ = ['Interaction', 'InteractionType', "InteractionError"]

class InteractionError(Exception):
    pass

class InteractionType(object):
    """Class for generating specific interaction type classes with the factory
    method.

    Examples
    --------

    """

    attributes = {}
    interaction_name = "None"

    feature_keys = None
    feature_classifiers = None

    # order is the number of features that participate in an interaction
    degree = 0
    commutative = True
    interaction_param_keys = []

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
        members_features = [[] for feature_key in cls.feature_keys]
        # and go through each member_type
        for member_idx, member_type in enumerate(association_type.member_types):
            for feature_type in member_type.feature_types.values():
                # get the classifiers for this feature (based on the
                # feature identification algorithms applied)
                # feature_classifiers = cls.feature_inx_attributes(feature_type)
                feature_classifiers = feature_type.feature_classifiers

                # if the feature has one of the classifiers for this member of the interaction
                inx_member_classifiers = cls.feature_classifiers[cls.feature_keys[member_idx]]
                if not set(feature_classifiers).isdisjoint(inx_member_classifiers):
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
    def feature_inx_attributes(cls, feature):
        """Check to see if this feature can be a participant in this
        interaction. Returns the classifiers that qualify it for
        participation in this interaction.

        """
        feature_inx_classifiers = []
        inx_classifiers = it.chain(*[classifiers for classifiers in
                                     cls.feature_classifiers.values()])
        if not set(feature.feature_classifiers).isdisjoint(inx_classifiers):
            feature_inx_classifiers.extend(feature.feature_classifiers)

        return feature_inx_classifiers

    @classmethod
    def check(cls, features, feature_tests):
        """The principle class method for testing for the existence of an
        interaction from specified geometric constraint parameters.
        All subclasses of InteractionType should implement this for
        domain specificity.

        """

        # initialize the list for storing every param value that will
        # be tested
        param_values = [None for test in feature_tests]
        # run the tests on the features
        for i, feature_test in enumerate(feature_tests):
            # run the test
            okay, param_value = feature_test(*features)
            # add the value good or bad to the param values
            param_values[i] = param_value
            # if it's bad return a bad okay and give the params
            if not okay:
                return False, tuple(param_values)

        # if they all passed, success!
        return True, tuple(param_values)


    @classmethod
    def new_find_hits(cls, members):
        pass
    @classmethod
    def find_hits(cls, members,
                  interaction_classes=interaction_classes,
                  return_feature_keys=False,
                  return_failed_hits=False):
        """Returns all the 'hits' for interactions of features between two
        members (molecules, selections, etc.), where a hit is a
        feature set that satisfies the geometric constraints defined
        in the InteractionType.check function.

        return_feature_keys will return the keys for the features
        instead of just the objects.

        return_failed_hits will return the parameters for each test in
        the check method for a potential interaction if at least one
        test was passed. This allows for debugging of the test
        functions in check.

        """

        # for each member collect the grouped features
        # initialize list of members (A, B, ...)
        members_features = [[] for i in members]
        for memb_idx, member in enumerate(members):
            # collect the features for this member/classifier
            for feature_key, feature in member.features.items():
                # get the classifiers this feature has
                feature_classifiers = feature.feature_type.feature_classifiers
                # if any match the classifiers for this interaction member keep it
                inx_member_classifiers = cls.feature_classifiers[cls.feature_keys[memb_idx]]
                if not set(feature_classifiers).isdisjoint(inx_member_classifiers):
                    members_features[memb_idx].append((feature_key, feature))

        # combine the features
        feature_tuples = it.product(members_features[0], members_features[1])

        if return_feature_keys:
            # initializing for the keys
            hit_pair_keys = []
        if return_failed_hits:
            failed_hits = []

        # initializing list for the actual Interaction objects
        inxs = []

        # for all the (feature_key, (features)) check if they are a
        # hit and make the Interaction object if they are
        for feature_key_pair in feature_tuples:
            if return_feature_keys:
                feature_keys = tuple([feature_key_tup[0] for feature_key_tup in feature_key_pair])
            features = tuple([feature_key_tup[1] for feature_key_tup in feature_key_pair])

            # call check for the InteractionType which checks to see
            # if the two features are in an interaction.
            # param_values should be in the same order as their labels in cls.interaction_param_keys
            okay, param_values = cls.check(*features)

            # if the feature pair did not pass do not construct an Interaction
            if not okay and return_failed_hits:
                if any([True for param in param_values if param is not None]):
                    failed_hits.append(param_values)
                continue
            elif not okay:
                continue

            # associate the parameter values with the names for them
            param_values = {param_name : param_val for param_name,
                            param_val in zip(cls.interaction_param_keys, param_values)}

            # otherwise make the Interaction from the features and the
            # values from check
            inx = cls.interaction_constructor(*features,
                                              interaction_class=None,
                                              check=False,
                                              **param_values)

            # classify the interaction if given an interaction space
            # of interaction classes
            if interaction_classes:
                interaction_class = match_inxclass(inx, interaction_classes)
                inx.interaction_class = interaction_class

            # add it to the list of Interactions
            inxs.append(inx)

            # and the feature keys to the feature key pairs
            if return_feature_keys:
                hit_pair_keys.append(feature_keys)

        if return_feature_keys:

            if return_failed_hits:
                return hit_pair_keys, failed_hits, inxs
            else:
                return hit_pair_keys, inxs
        elif return_failed_hits:
            return failed_hits, inxs
        else:
            return inxs

    @property
    def record(self):

        record_attr = {'InteractionType' : self.name}
        record_attr['AssociationType'] = self.association_type.name
        record_attr.update(self.attributes_data)

        return InteractionTypeRecord(**record_attr)

# InteractionTypeRecord
_interaction_type_record_fields = ['InteractionType', 'AssociationType']
InteractionTypeRecord = col.namedtuple('InteractionTypeRecord', _interaction_type_record_fields)

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
    def __init__(self, features=None, system=None,
                 interaction_type=None, interaction_class=None,
                 **param_values):

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
        self._interaction_class = interaction_class
        self._interaction_params = param_values
        # add the param values as individual attributes
        for param_name, param_value in param_values.items():
            self.__dict__[param_name] = param_value

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
        assert col.Counter(self.feature_types) == col.Counter(interaction_class.feature_types), \
            "the interaction_class must have the same number and types of features" \
            "as this Interaction, not {}".format(col.Counter(interaction_class.feature_types))
        self._interaction_class = interaction_class

    @property
    def interaction_params(self):
        return self._interaction_params

    @property
    def record(self):

        record_attr = {'InteractionType' : self.interaction_type.name}
        record_attr['InteractionClass'] = self.interaction_class

        return InteractionRecord(**record_attr)


# InteractionRecord
_interaction_record_fields = ['InteractionType', 'InteractionClass']
InteractionRecord = col.namedtuple('InteractionRecord', _interaction_record_fields)

def match_inxclass(inx, interaction_classes):
    """Given an Interaction object and a list of InteractionType objects
    (interaction classes) determines if the Interaction's features
    match a combination of features in one of the interaction classes.

    Returns the matching interaction class.

    """
    interaction_classes_it = iter(interaction_classes)
    found = False
    match = None
    # get the matching interaction class, throws error if no match
    while not found:

        # end the loop if no more interaction classes left
        try:
            inx_class = next(interaction_classes_it)
        except StopIteration:
            print("No matching interaction class given")
            break

        feature_pair = tuple([feature_type for feature_type
                        in inx_class.feature_types])
        # get the feature types to compare to the
        # feature pair in the inx class
        feature_types_tup = tuple([feature.feature_type for feature in
                                   inx.features])
        # if the interaction is not commutative the
        # order must be the same as the interaction class
        if not inx_class.commutative:
            if feature_pair == feature_types_tup:
                match = inx_class
                found = True
        # if the interaction is not commutative the
        # order might not be the same as the order in
        # the interaction class so we permute the
        # current inx feature types to check if any match
        else:
            for feature_pair_perm in it.permutations(feature_pair):
                if feature_pair_perm == feature_types_tup:
                    match = inx_class
                    found = True

    return match

if __name__ == "__main__":
    pass
