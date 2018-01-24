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
    hit_idxs = []
    # gather the substantiated features and check it
    for i, inx_class in enumerate(inx_classes):
        inx = profile_inx_class(inx_class, system)


        if inx is not None:
            hit_idxs.append(i)

            # DEBUG
            assert inx.interaction_class.name == inx_class.name, \
                "The interaction_class from the inx_class {0} does not match"\
                " the created Interaction {1} in hit {2}".format(inx_class.name,
                                                                 inx.interaction_class.name,
                                                                 i)
        # will add a None where the test fails to preserve position
        inxs.append(inx)


    return hit_idxs, inxs


def profiles_df(profiles, profile_ids=None):
    import pandas as pd
    hits_dfs = []
    for prof_idx, profile in profiles:
        hits_df = profile.hit_inx_df()
        if profile_ids is not None:
            prof_ids_col = [profile_ids[prof_idx]] * hits_df.shape[0]
        else:
            prof_ids_col = [prof_idx] * hits_df.shape[0]

        hits_df['profile_id'] = prof_ids_col

        hits_dfs.append(hits_df)

    master_df = pd.concat(hits_dfs)

    return master_df



class InxSpaceProfiler(object):

    def __init__(self, interaction_space):
        self._inx_space = interaction_space

    def profile(self, system):
        return InxSpaceProfile(self._inx_space, system)

class InxSpaceProfile(object):

    def __init__(self, inx_space, system):
        self._system = system
        self._inx_space = inx_space

        # profile by interaction class order
        self._hit_idxs, self._inxs = profile_inx_classes(self._inx_space, self._system)

    @property
    def n_inx_classes(self):
        return len(self._inx_space)

    @property
    def inxs(self):
        return self._inxs

    @property
    def hit_idxs(self):
        return self._hit_idxs
        #return [i for i, inx in enumerate(self._inxs) if inx is not None]

    @property
    def hit_inxs(self):
        return [inx for inx in self._inxs if inx is not None]

    @property
    def vector(self):
        return [0 if inx is None else 1 for inx in self._inxs]

    def hit_inx_records(self):
        return [inx.record for inx in self.hit_inxs]

    def hit_inx_dict(self):
        profile_dict = col.defaultdict(list)
        for inx in self.hit_inxs:
            for field, value in inx.record_dict.items():
                profile_dict[field].append(value)

        return profile_dict

    def hit_inx_df(self):
        import pandas as pd
        hit_df = pd.DataFrame(self.hit_inx_dict())
        hit_df['hit_idx'] = self.hit_idxs
        return hit_df

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
        return self._inx_space._subspace_map


    # @property
    # def subspace_vector(self, association_type, interaction_type):
    #     key = (association_type, interaction_type)
    #     idxs = self._subspace_map[key]
    #     sel_inxs = [inx for i, inx in enumerate(self._inxs) if i in idxs]
    #     return [0 if inx is None else 1 for inx in sel_inxs]

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


    # def hits_by_inx_type(self, interaction_type):
    #     """Returns the indices of interactions matching the interaction_type"""

    #     return_idxs = []
    #     # for each subspace
    #     for assoc_inxtype_tup, idxs in self._subspace_map.items():
    #         # if the subspace involves the interaction type
    #         if interaction_type == assoc_inxtype_tup[self._interaction_type_idx]:
    #             # then we get the hit_idxs that match the ones in this subspace
    #             subspace_hit_idxs = [idx for idx in idxs if idx in self.hit_idxs]
    #             return_idxs.extend(subspace_hit_idxs)

    #     return return_idxs

    # def hits_by_association(self, association_type):
    #     """Returns the indices of interactions matching the association_type"""
    #     return_idxs = []
    #     for assoc_inxtype_tup, idxs in self._subspace_map.items():
    #         if association_type == assoc_inxtype_tup[self._association_idx]:
    #             hit_idxs = [idx for idx in idxs if idx in self.hit_idxs]
    #             return_idxs.extend(hit_idxs)

    #     return return_idxs
