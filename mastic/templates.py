class InxType(InteractionType):

    def __init__(self, inx_attr=None):
        super().__init__(attr_dict=inx_attrs)
    _feature_families = INX_FEATURE_FAMILIES
    feature_families = _feature_families
    _feature1_key = 'Key1'
    _feature2_key = 'Key2'
    _feature_types = INX_FEATURE_TYPES
    feature_type = _feature_types

    def __repr__(self):
        return str(self.__class__)

    @classmethod
    def find_hits(cls, **kwargs):
        """Takes in key-word arguments for the features. """

        from itertools import product
        # check that the keys ar okay in parent class
        super().find_hits(**kwargs)

        # inx specific stuff
        feature1 = kwargs[cls._feature1_key]
        feature2 = kwargs[cls._feature2_key]
        hits = []
        # make pairs of them to compare
        pairs = product(feature1, feature2)
        for pair in pairs:
            features = feature_dict(pair)
            # try to make a Inx object, which calls check
            try:
                # if it succeeds add it to the list of H-Bonds
                inx = Inx(**features)

            # else continue to the next pairing
            except InteractionError:
                continue

            hits.append(inx)

        return hits

    @classmethod
    def check(cls, feature1=None, feature2=None):
        """Check if the input features qualify for this type of
        interaction.

        returns (okay, param1, param2,...)"""
        param1 = calc_param1(feature1, feature2)
        if cls.check_param1(param1) is False:
            return (False, param1, None)

        param2 = calc_param2(feature1, feature2)
        if cls.check_param2(param2) is False:
            return (False, param1, param2)

        return (True, param1, param2)

    @classmethod
    def check_param1(cls, param1):
        if True:
            return True
        else:
            return False

    @classmethod
    def check_param2(cls, param2):
        if True:
            return True
        else:
            return False


class Inx(Interaction):
    """Interaction that has type InxType.

    """

    def __init__(self, feature1=None, feature2=None):

        okay, param1, param2 = InxType.check(feature1, feature2)
        if not okay:
            if param2 is None:
                raise InteractionError(
                    """feature1: {0}
feature2: {1}
                    param1 = {2} FAILED
                    param2 = not calculated""".format(feature1, feature2, param1))

            else:
                raise InteractionError(
                    """feature1: {0}
feature2: {1}
param1: {2}
param2 = {3} FAILED """.format(feature1, feature2, param1, param2))

        # success, finish creating interaction

        # TODO generalize
        system = feature1.system
        super().__init__(members=[feature1, feature2],
                         interaction_type=InxType,
                         system=system)
        self._interaction_type = InxType
        self._feature1 = feature1
        self._feature2 = feature2
        self._param1 = param1
        self._param2 = param2

    @property
    def feature1(self):
        return self._feature1

    @property
    def feature2(self):
        return self._feature2

    @property
    def param1(self):
        return self._param1

    @property
    def param2(self):
        return self._param2
