""" Module for storing interaction profiles of Systems and SystemTypes. """

import collections as col
import itertools as it

class ProfileError(Exception):
    pass

def get_inx_class_features(inx_class, system):

    features = []
    for i, feature_type in enumerate(inx_class.feature_types):
        # find the index of the feature in the member, order of features determines
        member_idx = inx_class.association_type.member_idxs[i]
        feat_idx = list(inx_class.association_type.member_types[i].feature_types.values())\
                   .index(feature_type)
        feature = list(system.members[member_idx].features.values())[feat_idx]


        feature = list(system.members[member_idx].features.values())[feat_idx]
        features.append(feature)

    return features

def profile_inx_class(inx_class, system):

    # get the features for this inx class
    features = get_inx_class_features(inx_class, system)

    # calculate inter-feature parameters and check against
    # constraints
    okay, param_values = inx_class.check(*features)

    if not okay:
        return None

    # if it passes we want to actually construct the interaction object
    # associate the parameter values with the names for them
    param_values = {param_name : param_val for param_name,
                    param_val in zip(inx_class.interaction_param_keys, param_values)}

    # construct the interaction
    inx = inx_class.interaction_constructor(*features,
                                            interaction_class=inx_class,
                                            check=False,
                                            **param_values)
    return inx

def profile_inx_classes(inx_classes, system):
    inxs = []
    # gather the substantiated features and check it
    for inx_class in inx_classes:
        inx = profile_inx_class(inx_class, system)
        # will add a None where the test fails to preserve position
        inxs.append(inx)

    return inxs

class InxSpaceProfiler(object):

    def __init__(self, interaction_space):
        self._inx_space = interaction_space

    def profile(self, system):
        return InxSpaceProfile(self._inx_space, system)

class InxSpaceProfile(object):

    def __init__(self, inx_space, system):
        self._system = system
        self._inx_space = inx_space

        self._association_idx = 0
        self._interaction_type_idx = 1
        self._subspace_map = {}
        self._inxs = []

        # do profiling for each subspaces interaction classes
        for subspace_tup, inx_class_idxs in self._inx_space.subspace_map.items():
            inx_classes = [self._inx_space[idx] for idx in inx_class_idxs]
            inx_hits = profile_inx_classes(inx_classes, self._system)

            # add these to the list of interactions keeping track of
            # their indices
            inx_idx = len(self._inxs)
            new_idxs = []
            for i, inx in enumerate(inx_hits):
                self._inxs.append(inx)
                new_idxs.append(inx_idx + i)

            self._subspace_map[subspace_tup] = new_idxs

    @property
    def n_inx_classes(self):
        return len(self._inxs)

    @property
    def inxs(self):
        return self._inxs

    @property
    def hit_idxs(self):
        return [i for i, inx in enumerate(self._inxs) if inx is not None]

    @property
    def hit_inxs(self):
        return [inx for inx in self._inxs if inx is not None]

    @property
    def vector(self):
        return [0 if inx is None else 1 for inx in self._inxs]

    @property
    def subspace_vector(self, association_type, interaction_type):
        key = (association_type, interaction_type)
        idxs = self._subspace_map[key]
        sel_inxs = [inx for i, inx in enumerate(self._inxs) if i in idxs]
        return [0 if inx is None else 1 for inx in sel_inxs]

    def hit_inx_records(self):
        return [inx.record for inx in self.hit_inxs]

    def hit_inx_df(self):
        import pandas as pd
        return pd.DataFrame(self.hit_inx_records())

    # TODO have a problem where a hit idx is assigned two inxs if they
    # are the opposite of each other in the association

    # def inx_type_hit_records(self, interaction_type):
    #     hit_records = []
    #     hit_idx_key = 'hit_idx'

    #     # we will want to make a new record for hits, so we get an
    #     # example record from the interaction
    #     inx_idxs = self.hits_by_inx_type(interaction_type)
    #     inx_record = self.inxs[inx_idxs[0]].record

    #     # add new hit_idx field
    #     record_fields = list(inx_record._fields)
    #     record_fields.append(hit_idx_key)
    #     # modify the name
    #     hit_record_name = "Hit" + type(inx_record).__name__
    #     # make the new namedtuple
    #     hit_record_type = col.namedtuple(hit_record_name, record_fields)

    #     # get the hits for this interaction type
    #     for hit_idx in inx_idxs:
    #         inx = self.inxs[hit_idx]
    #         # convert to a dictionary
    #         inx_dict_rec = inx.record._asdict()
    #         # add the hit index
    #         inx_dict_rec[hit_idx_key] = hit_idx
    #         # make the new record
    #         hit_record = hit_record_type(**inx_dict_rec)
    #         hit_records.append(hit_record)

    #     return hit_records

    # def association_hit_records(self, association_type):
    #     """Returns a dictionary of the hit records for AssociationType."""
    #     hit_records = []
    #     hit_idx_key = 'hit_idx'

    #     # we will want to make a new record for hits, so we get an
    #     # example record from the interaction
    #     inx_idxs = self.hits_by_association(association_type)
    #     inx_record = self.inxs[inx_idxs[0]].record

    #     # add new hit_idx field
    #     record_fields = list(inx_record._fields)
    #     record_fields.append(hit_idx_key)
    #     # modify the name
    #     hit_record_name = "Hit" + type(inx_record).__name__
    #     # make the new namedtuple
    #     hit_record_type = col.namedtuple(hit_record_name, record_fields)

    #     # get the hits for this interaction type
    #     for hit_idx in inx_idxs:
    #         inx = self.inxs[hit_idx]
    #         # convert to a dictionary
    #         inx_dict_rec = inx.record._asdict()
    #         # add the hit index
    #         inx_dict_rec[hit_idx_key] = hit_idx
    #         # make the new record
    #         hit_record = hit_record_type(**inx_dict_rec)
    #         hit_records.append(hit_record)

    #     return hit_records

    # def inx_type_hits_df(self, interaction_type):
    #     import pandas as pd
    #     return pd.DataFrame(self.inx_type_hit_records(interaction_type))

    @property
    def system(self):
        return self._system

    @property
    def interaction_space(self):
        return self._inx_space
    inx_space = interaction_space

    @property
    def interactions(self):
        return self._inxs
    inxs = interactions

    @property
    def subspace_map(self):
        return self._subspace_map

    def hits_by_inx_type(self, interaction_type):
        """Returns the indices of interactions matching the interaction_type"""

        return_idxs = []
        # for each subspace
        for assoc_inxtype_tup, idxs in self._subspace_map.items():
            # if the subspace involves the interaction type
            if interaction_type == assoc_inxtype_tup[self._interaction_type_idx]:
                # then we get the hit_idxs that match the ones in this subspace
                subspace_hit_idxs = [idx for idx in idxs if idx in self.hit_idxs]
                return_idxs.extend(subspace_hit_idxs)

        return return_idxs

    def hits_by_association(self, association_type):
        """Returns the indices of interactions matching the association_type"""
        return_idxs = []
        for assoc_inxtype_tup, idxs in self._subspace_map.items():
            if association_type == assoc_inxtype_tup[self._association_idx]:
                hit_idxs = [idx for idx in idxs if idx in self.hit_idxs]
                return_idxs.extend(hit_idxs)

        return return_idxs
