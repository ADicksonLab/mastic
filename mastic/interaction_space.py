import collections as col

class InteractionSpace(col.UserList):
    """An indexed list of InteractionType objects (interaction classes),
    which are accessible by their InteractionType."""

    def __init__(self, system_type):

        super().__init__()
        self._system_type = system_type

        # the order of the association+interaction_type mapping tuple key
        self._association_idx = 0
        self._interaction_type_idx = 1
        # the map
        self._subspace_map = {}

    @property
    def system_type(self):
        return self._system_type

    @property
    def subspace_map(self):
        return self._subspace_map

    def inxspace_dict(self):
        inxclasses_dict = col.defaultdict(list)
        for inxclass in self:
            for field, value in inxclass.record_dict.items():
                inxclasses_dict[field].append(value)

        return inxclasses_dict

    def to_dataframe(self):
        import pandas as pd
        return pd.DataFrame(self.inxspace_dict())
        #return pd.DataFrame([inx_class.record for inx_class in self])

    def add_association_subspace(self, association_type, interaction_type, return_idxs=False):
        inx_classes = interaction_type.interaction_classes(association_type)
        new_idxs = self.add_inx_classes(inx_classes)

        # TODO do we really need this??
        self._subspace_map[(association_type, interaction_type)] = new_idxs

        if return_idxs:
            return new_idxs

    def add_inx_classes(self, inx_classes):
        # get the index of the next would-be addition to the list
        inx_idx = len(self)

        new_idxs = []
        # add them and get the indices they would be in this space
        for i, inx_class in enumerate(inx_classes):
            self.append(inx_class)
            new_idxs.append(inx_idx + i)

        return new_idxs

    def by_inx_type(self, interaction_type):
        """Returns the indices of interaction classes matching the interaction_type"""

        return_idxs = []
        for assoc_inxtype_tup, idxs in self._subspace_map.items():
            if interaction_type == assoc_inxtype_tup[self._interaction_type_idx]:
                return_idxs.extend(idxs)

        return return_idxs

    def by_association(self, association_type):
        """Returns the indices of interaction classes matching the association_type"""
        return_idxs = []
        for assoc_inxtype_tup, idxs in self._subspace_map.items():
            if association_type == assoc_inxtype_tup[self._association_idx]:
                return_idxs.extend(idxs)

        return return_idxs

    def by_idxs(self, idxs):
        return [self[idx] for idx in idxs]


if __name__ == "__main__":
    pass
