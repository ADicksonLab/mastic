class NoHHydrogenBondType(InteractionType):
    """Defines an InteractionType class for hydrogen bonds between
    Features without explicit hydrogens.

    """

    def __init__(self):
        pass

    _feature_families = mastinxconfig.HBOND_FEATURE_FAMILIES
    feature_families = _feature_families
    _donor_key = 'Donor'
    _acceptor_key = 'Acceptor'
    _grouping_attribute = 'rdkit_family'
    _feature_types = mastinxconfig.HBOND_FEATURE_TYPES
    feature_type = _feature_types

    def __repr__(self):
        return str(self.__class__)

    @classmethod
    def find_hits(cls, member_a, member_b):

        # check that the keys ar okay in parent class
        # super().find_hits(members_features)

        # for each member collect the grouped features
        # initialize list of members
        members_features = [{'donors':[], 'acceptors':[]} for member in [member_a, member_b]]
        for i, member in enumerate([member_a, member_b]):
            for feature_key, feature in member.features.items():
                # get groupby attribute to use as a key
                group_attribute = feature.feature_type.attributes_data[cls._grouping_attribute]

                if group_attribute == cls._acceptor_key:
                    acceptor_tup = (feature_key, feature)
                    members_features[i]['acceptors'].append(acceptor_tup)

                elif group_attribute == cls._donor_key:
                    donor_tup = (feature_key, feature)
                    members_features[i]['donors'].extend(donor_H_pairs_tup)

        donor_acceptor_pairs = []
        # pair the donors from the first with acceptors of the second
        donor_acceptor_pairs.extend(it.product(members_features[0]['donors'],
                                               members_features[1]['acceptors']))
        # pair the acceptors from the first with the donors of the second
        donor_acceptor_pairs.extend(it.product(members_features[1]['donors'],
                                               members_features[0]['acceptors']))

        # scan the pairs for hits
        hit_pair_keys = []
        hbonds = []
        for donor_tup, acceptor_tup in donor_acceptor_pairs:
            donor_feature_key = donor_tup[0]
            donor_feature = donor_tup[1]
            acceptor_feature_key = acceptor_tup[0]
            acceptor_feature = acceptor_tup[1]
            # try to make a HydrogenBondInx object, which calls check,
            #
            # OPTIMIZATION: otherwise we have to call check first then
            # the HydrogenBondInx constructor will re-call check to
            # get the angle and distance. If we allow passing and not
            # checking the angle and distance in the constructor then
            # it would be faster, however I am not going to allow that
            # in this 'safe' InteractionType, an unsafe optimized
            # version can be made separately if desired.
            try:
                hbond = NoHHydrogenBondInx(donor=donor_feature,
                                           acceptor=acceptor_feature)
            # else continue to the next pairing
            except InteractionError:
                continue
            # if it succeeds add it to the list of H-Bonds
            hbonds.append(hbond)
            # and the feature keys to the feature key pairs
            hit_pair_keys.append((donor_feature_key, acceptor_feature_key))

        return hit_pair_keys, hbonds

    @classmethod
    def check(cls, donor_atom, acceptor_atom):
        """Check if the 2 atoms qualify for a hydrogen bond."""

        from scipy.spatial.distance import cdist
        distance = cdist(np.array([donor_atom.coords]), np.array([acceptor_atom.coords]))[0,0]
        if cls.check_distance(distance) is False:
            return (False, distance, None)

        num_Hs = len([atom for atom in donor_atom.adjacent_atoms if atom.element == 'H'])
        total_bond_degree = len(donor_atom.adjacent_atoms)
        num_non_Hs = total_bond_degree - num_Hs

        # depending on what kind of donor it is assign the virtual H atom accordingly
        H_coords = None

        # primary amine
        if donor_atom.element == 'N' and num_non_Hs == 1 and num_Hs == 2:
            pass
        # secondary amine
        elif donor_atom.element == 'N' and num_non_Hs == 2 and num_Hs == 1:
            pass
        # charged primary amine
        elif donor_atom.element == 'N' and num_non_Hs == 1 and num_Hs == 3:
            pass
        # hydroxyl group
        elif donor_atom.element == 'O' and num_non_Hs == 1 and num_Hs == 1:
            pass
        # water
        elif donor_atom.element == 'O' and num_non_Hs == 0 and num_Hs == 2:
            pass

class NoHHydrogenBondInx(Interaction):
    """Substantiates NoHHydrogenBondType by selecting donor and acceptor
    features.

    """

    def __init__(self, donor=None, acceptor=None):

        donor_atom = donor.atoms[0]
        acceptor_atom = acceptor.atoms[0]
        okay, distance, angle = NoHHydrogenBondType.check(donor_atom, acceptor_atom)
        if not okay:
            if angle is None:
                raise InteractionError(
                    """donor: {0}
                    H: {1}
                    acceptor: {2}
                    distance = {3} FAILED
                    angle = not calculated""".format(donor_atom, H, acceptor_atom, distance))

            else:
                raise InteractionError(
                    """donor: {0}
                    H: {1}
                    acceptor: {2}
                    distance = {3}
                    angle = {4} FAILED""".format(donor_atom, H, acceptor_atom, distance, angle))

        # success, finish creating interaction
        atom_system = donor.system
        super().__init__(features=[donor, acceptor],
                         interaction_type=HydrogenBondType,
                         system=atom_system)
        self._donor = donor
        self._H = H
        self._acceptor = acceptor
        self._distance = distance
        self._angle = angle
