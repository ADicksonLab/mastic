""" Module for storing interaction profiles of Systems and SystemTypes. """

import collections as col

class ProfileError(Exception):
    pass

class SystemTypeProfile(col.defaultdict):

    def __init__(self, system_type, system_profiles):
        # initialize the default dict as a list always
        super().__init__(list)

        self._system_type = system_type
        # store the values for the profile according to their InteractionType
        for profile in system_profiles:
            self[profile.interaction_type].append(profile)

    @property
    def system_type(self):
        return self._system_type

    @property
    def interaction_types(self):
        return list(self.keys())

    @property
    def profiles(self):
        return list(self.values())

class SystemProfileDataset(object):
    def __init__(self, system_type, system_profiles):
        self._system_type = system_type
        self._profiles = system_profiles

class SystemProfile(col.UserList):

    def __init__(self, interaction_space, inxs):

        # none of the inxs can be for outside the space, this could be
        # silently ignored but probably indicates another problem if so
        assert all([True if inx.interaction_class in interaction_space else False
                    for inx in inxs]), \
                    "All inxs must be of interaction classes in the interaction space."

        hits = []
        inx_records = []
        for inx in inxs:
            idx = interaction_space.index(inx.interaction_class)

            # there can only be one hit at most for each interaction class
            if idx in hits:
                raise ProfileError("can't have two interactions of the same interaction"
                                   "class in a single System")
            hits.append(idx)
            inx_records.append(inx.record)

        self._hits = hits
        self._hit_records = inx_records

    @property
    def hits(self):
        return self._hits

    @property
    def hit_records(self):
        return self._hit_records

    @property
    def hits_df(self):
        import pandas as pd
        return pd.DataFrame(self.hit_records)

class AssociationProfile(object):
    def __init__(self, interaction):
        pass
