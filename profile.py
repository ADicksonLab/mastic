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
        try:
            feature = list(system.members[member_idx].features.values())[feat_idx]
        except IndexError:
            import ipdb; ipdb.set_trace()

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
        self._subspace_map = {}
        self._inxs = []
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


    # TODO
    @property
    def hits(self):
        return list(it.chain(*[prof.hits for prof in self.association_profiles.values()]))

    @property
    def records(self):
        return list(it.chain(*[prof.records for prof in self.association_profiles.values()]))

    @property
    def hits_df(self):
        import pandas as pd
        df = pd.DataFrame(self.records)
        df['hit_idx'] = self.hits
        return df

class AssociationProfile(object):

    def __init__(self, association):

        self._association = association
        self.association_type = association.association_type
        self._hits = []
        self._hit_records = []

    def profile_interactions(self, returns='hits'):

        inxs = []
        hits = []

        for inx_class_idx, inx_class in enumerate(self.association_type.interaction_subspace):
            # get the features for this interaction class from this association
            features = []
            for i, feature_type in enumerate(inx_class.feature_types):
                # find the index of the feature in the member, order of features determines
                feat_idx = list(self.association_type.member_types[i].feature_types.values())\
                           .index(feature_type)
                feature = list(self.members[i].features.values())[feat_idx]
                features.append(feature)

            okay, param_values = inx_class.check(*features)

            if not okay:
                # move onto next inx_class
                continue

            # associate the parameter values with the names for them
            param_values = {param_name : param_val for param_name,
                            param_val in zip(inx_class.interaction_param_keys, param_values)}

            inx = inx_class.interaction_constructor(*features,
                                                    interaction_class=inx_class,
                                                    check=False,
                                                    **param_values)
            inxs.append(inx)
            hits.append(inx_class_idx)

        if returns == 'inxs' or returns == 'interactions':
            return inxs
        elif returns == 'records':
            records = [inx.record for inx in inxs]
            return records
        elif returns == 'hits':
            records = [inx.record for inx in inxs]
            return hits, records

        hits, records = self._association.profile_interactions()
        self._hits.extend(hits)
        self._hit_records.extend(records)

    @property
    def hits(self):
        return self._hits

    @property
    def records(self):
        return self._hit_records

    @property
    def hits_df(self):
        import pandas as pd
        df = pd.DataFrame(self.records)
        df['hit_idx'] = self.hits
        return df
