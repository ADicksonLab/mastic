""" Module for storing interaction profiles of Systems and SystemTypes. """

import collections as col
import itertools as it

class ProfileError(Exception):
    pass


class SystemProfile(object):

    def __init__(self, system):
        self._system = system
        self._assoc_profiles = col.OrderedDict()

    def profile_associations(self):
        for assoc_term in self.system.system_type.assoc_member_idxs:
            self.profile_association(assoc_term)

    def profile_association(self, assoc_term):
        assoc_idx = self.system.system_type.unit_association_type_idxs[assoc_term]
        association = self.system.associations[assoc_idx]
        assoc_profile = AssociationProfile(association)
        assoc_profile.profile_interactions()
        self._assoc_profiles[assoc_term] = assoc_profile

    @property
    def system(self):
        return self._system

    @property
    def association_profiles(self):
        return self._assoc_profiles

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
        self._hits = []
        self._hit_records = []

    def profile_interactions(self):
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
