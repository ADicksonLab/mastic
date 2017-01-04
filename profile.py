""" Module for storing interaction profiles of Systems and SystemTypes. """

import collections as col

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

class SystemProfile(object):

    def __init__(self, interaction_type, inx_class_idxs, **param_values):
        self._interaction_type = interaction_type
        self._inx_class_idxs = inx_class_idxs
        self._param_values = param_values
        # add them to the namespace for exploratory purposes
        self.__dict__.update(param_values)

    @property
    def interaction_type(self):
        return self._interaction_type

    @property
    def inx_class_idxs(self):
        return self._inx_class_idxs

    @property
    def param_values(self):
        return self._param_values
