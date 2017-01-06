import collections as col

class InteractionSpace(col.UserList):
    """An indexed list of InteractionType objects (interaction classes),
    which are accessible by their InteractionType."""

    def __init__(self, system_type, interaction_classes):
        # check that all of the interaction classes come from the same system
        assert all([inx_class.association_type in
                    system_type.association_types for inx_class in
                    interaction_classes]), \
                    "All interaction classes must come from the system_type"

        # add the interaction classes into this object
        super().__init__(interaction_classes)
        self._system_type = system_type

        # things that will be found and memoized when first called
        self._inx_type_idxs = None
        self._assoc_idxs = None

    @property
    def system_type(self):
        return self._system_type

    def add_inx_classes(self, inx_classes, return_idxs=False):
        # get the index of the next would-be addition to the list
        inx_idx = len(self)

        new_idxs = []
        # add them and get the indices thwy would be in this space
        for i, inx_class in enumerate(inx_classes):
            self.append(inx_class)
            new_idxs.append(inx_idx + i)

        if return_idxs:
            return new_idxs

    def by_inx_type(self, interaction_type):
        """Returns the indices of interaction classes matching the interaction_type"""

        return [i for i, inx_class in enumerate(self) if
                inx_class.interaction_type == interaction_type]

    def by_association(self, association_type):
        """Returns the indices of interaction classes matching the association_type"""

        return [i for i, inx_class in enumerate(self) if
                inx_class.association_type == association_type]

class InteractionSubSpace(col.UserList):
    """An indexed list of InteractionType objects (interaction classes),
    which are accessible by their InteractionType."""

    def __init__(self, association_type, interaction_classes):
        # check that all of the interaction classes come from the same system
        assert all([inx_class.association_type in association_type for
                    inx_class in interaction_classes]), \
                    "All interaction classes must come from the association_type"

        # add the interaction classes into this object
        super().__init__(interaction_classes)
        self._association_type = association_type

    @property
    def association_type(self):
        return self._association_type

    def by_inx_type(self, interaction_type):
        """Returns the indices of interaction classes matching the interaction_type"""
        return [i for i, inx_class in enumerate(self) if
                inx_class.interaction_type == interaction_type]


if __name__ == "__main__":
    pass
