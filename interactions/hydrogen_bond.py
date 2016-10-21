"""The HydrogenBond module defines the HydrogenBondType and
HydrogenBondInx for explicit hydrogens.

"""
import itertools as it
from collections import namedtuple, defaultdict

import numpy as np
import numpy.linalg as la

import mast.config.interactions as mastinxconfig
from mast.interactions.interactions import InteractionType, Interaction, InteractionError


class HydrogenBondType(InteractionType):
    """Defines an InteractionType class for hydrogen bonds between members
    with explicit hydrogens.

    """

    interaction_name = "HydrogenBond"
    feature_keywords = mastinxconfig.HBOND_FEATURES
    donor_key = 'Donor'
    acceptor_key = 'Acceptor'
    feature_order = [donor_key, acceptor_key]
    grouping_attribute = 'rdkit_family'
    # order is the number of features that participate in an interaction
    degree = 2
    commutative = False

    def __init__(self, hydrogen_bond_type_name,
                 feature_types=None,
                 association_type=None,
                 assoc_member_pair_idxs=None,
                 **hydrogen_bond_attrs):

        super().__init__(hydrogen_bond_type_name,
                         feature_types=feature_types,
                         association_type=association_type,
                         assoc_member_pair_idxs=assoc_member_pair_idxs,
                         **hydrogen_bond_attrs)

        self.donor = feature_types[0]
        self.acceptor = feature_types[1]

    @classmethod
    def test_find_hits(cls, members,
                       interaction_classes=None,
                       distance_cutoff=mastinxconfig.HBOND_DIST_MAX,
                       angle_cutoff=mastinxconfig.HBOND_DON_ANGLE_MIN):

        # TODO checks

        # for each member collect the grouped features
        # initialize list of members (A, B, ...)
        members_features = tuple([[] for i in members])
        for i, member in enumerate(members):
            for feature_key, feature in member.features.items():
                # get groupby attribute to use as a key
                group_attribute = feature.feature_type.attributes_data[cls.grouping_attribute]


        ##### InteractionType specifc logic
                if group_attribute in cls.aromatic_keys:
                    aromatic_tup = (feature_key, feature)
                    members_features[i].append(aromatic_tup)



        # combine the features. this part may be improved by using the
        # degree and the symmetry of the interaction.
        # feature_tuples = it.product(members_features, repeat=cls.degree)

        # for now we rely on knowledge of the interaction
        feature_tuples = it.product(members_features[0],
                                    members_features[1])

        #####

        # scan the pairs for hits and assign interaction classes if given
        return super().find_hits(feature_tuples,
                                 # if there was an interaction space given
                                 interaction_classes=interaction_classes,
                                 # the parameters for the interaction existence
                                 distance_cutoff=mastinxconfig.HBOND_DIST_MAX,
                                 angle_cutoff=mastinxconfig.HBOND_DON_ANGLE_MIN)

    @classmethod
    def find_hits(cls, member_a, member_b,
                  interaction_classes=None,
                  distance_cutoff=mastinxconfig.HBOND_DIST_MAX,
                  angle_cutoff=mastinxconfig.HBOND_DON_ANGLE_MIN):

        # check that the keys ar okay in parent class
        # super().find_hits(members_features)

        # for each member collect the grouped features
        # initialize list of members
        members_features = [{'donors':[], 'acceptors':[]} for member in [member_a, member_b]]
        for i, member in enumerate([member_a, member_b]):
            for feature_key, feature in member.features.items():
                # get groupby attribute to use as a key
                group_attribute = feature.feature_type.attributes_data[cls.grouping_attribute]

                if group_attribute == cls.acceptor_key:
                    acceptor_tup = (feature_key, feature)
                    members_features[i]['acceptors'].append(acceptor_tup)

                elif group_attribute == cls.donor_key:
                    # get the donor-H pairs of atoms for this donor
                    donor_atom = feature.atoms[0]
                    donor_H_pairs = [(feature, atom) for atom in
                                     donor_atom.adjacent_atoms if
                                     atom.atom_type.element == 'H']
                    donor_H_pairs_tup = [(feature_key, donor_H_pair) for
                                         donor_H_pair in donor_H_pairs]
                    members_features[i]['donors'].extend(donor_H_pairs_tup)

        donor_acceptor_pairs = []
        # pair the donors from the first with acceptors of the second,
        # keeping track of which member it is in
        pairs = list(it.product(members_features[0]['donors'],
                           members_features[1]['acceptors']))
        # make tuples of (member order, donor-acceptors) ((0,1), (donor, acceptor))
        donor_acceptor_pairs.extend(zip([(0,1) for i in range(len(pairs))], pairs))

        # pair the acceptors from the first with the donors of the
        # second, keeping track of which member it is in
        pairs = list(it.product(members_features[1]['donors'],
                           members_features[0]['acceptors']))
        # make tuples of (member order, donor-acceptors) ((1,0), donor, acceptor)
        donor_acceptor_pairs.extend(zip([(1,0) for i in range(len(pairs))], pairs))

        # scan the pairs for hits
        hit_pair_keys = []
        hbonds = []
        for member_order, donor_acceptor_tup in donor_acceptor_pairs:
            donor_tup = donor_acceptor_tup[0]
            acceptor_tup = donor_acceptor_tup[1]
            donor_feature_key = donor_tup[0]
            donor_feature = donor_tup[1][0]
            h_atom = donor_tup[1][1]
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
                hbond = HydrogenBondInx(donor=donor_feature, H=h_atom,
                                        acceptor=acceptor_feature,
                                        distance_cutoff=distance_cutoff,
                                        angle_cutoff=angle_cutoff,
                                        interaction_class=None)
            # else continue to the next pairing
            except InteractionError:
                continue

            # classify the hbond if given classes
            interaction_class = None
            if interaction_classes:
                feature_pairs = [(inx_class.donor, inx_class.acceptor) for
                                 inx_class in interaction_classes]
                # get the matching interaction class, throws error if no match

                try:
                    interaction_classes_it = iter(interaction_classes)
                    found = False
                    while not found:
                        inx_class = next(interaction_classes_it)
                        feature_pair = (inx_class.donor, inx_class.acceptor)
                        if feature_pair == \
                           (donor_feature.feature_type, acceptor_feature.feature_type):
                            hbond.interaction_class = inx_class
                            found = True

                except StopIteration:
                    print("No matching interaction class given")

            # if it succeeds add it to the list of H-Bonds
            hbonds.append(hbond)
            # and the feature keys to the feature key pairs
            hit_pair_keys.append((member_order, (donor_feature_key, acceptor_feature_key,),))

        return hit_pair_keys, hbonds

    @classmethod
    def check(cls, donor_atom, h_atom, acceptor_atom,
              distance_cutoff=mastinxconfig.HBOND_DIST_MAX,
              angle_cutoff=mastinxconfig.HBOND_DON_ANGLE_MIN):
        """Checks if the 3 atoms qualify as a hydrogen bond. Returns a tuple
        (bool, float, float) where the first element is whether or not it
        qualified, the second and third are the distance and angle
        respectively. Angle may be None if distance failed to qualify.

        """
        from scipy.spatial.distance import cdist
        distance = cdist(np.array([donor_atom.coords]), np.array([acceptor_atom.coords]))[0,0]
        if cls.check_distance(distance, cutoff=distance_cutoff) is False:
            return (False, distance, None)

        v1 = donor_atom.coords - h_atom.coords
        v2 = acceptor_atom.coords - h_atom.coords
        try:
            angle = np.degrees(np.arccos(np.dot(v1, v2)/(la.norm(v1) * la.norm(v2))))
        except RuntimeWarning:
            print("v1: {0} \n"
                  "v2: {1}".format(v1, v2))
        if cls.check_angle(angle, cutoff=angle_cutoff) is False:
            return (False, distance, angle)

        return (True, distance, angle)

    @property
    def record(self):
        record_fields = ['interaction_class', 'interaction_type',
                         'association_type', 'assoc_member_pair_idxs',
                         'donor_feature_type', 'acceptor_feature_type'] + \
                         list(self.attributes_data.keys())
        HydrogenBondTypeRecord = namedtuple('HydrogenBondTypeRecord', record_fields)
        record_attr = {'interaction_class' : self.name}
        record_attr['interaction_type'] = self.interaction_name
        record_attr['association_type'] = self.association_type.name
        record_attr['assoc_member_pair_idxs'] = self.assoc_member_pair_idxs
        record_attr['donor_feature_type'] = self.feature_types[0].name
        record_attr['acceptor_feature_type'] = self.feature_types[1].name

        return HydrogenBondTypeRecord(**record_attr)

    def pdb_serial_output(self, inxs, path, delim=","):
        """Output the pdb serial numbers (index in pdb) of the donor and
        acceptor in each HBond to a file like:

        donor_1, acceptor_1
        donor_2, acceptor_2
        ...

        Notice: This will probably be removed in the future.

        """

        with open(path, 'w') as wf:
            for inx in inxs:
                wf.write("{0}{1}{2}\n".format(inx.donor.atom_type.pdb_serial_number,
                                              delim,
                                              inx.acceptor.atom_type.pdb_serial_number))

    #### Hydrogen Bond Specific methods
    # i.e. not necessarily found in other interaction types
    @classmethod
    def check_distance(cls, distance, cutoff=mastinxconfig.HBOND_DIST_MAX):
        """For a float distance checks if it is less than the configuration
        file HBOND_DIST_MAX value.

        """
        if distance < cutoff:
            return True
        else:
            return False

    @classmethod
    def check_angle(cls, angle, cutoff=mastinxconfig.HBOND_DON_ANGLE_MIN):
        """For a float distance checks if it is less than the configuration
        file HBOND_DON_ANGLE_MIN value.

        """

        if angle > cutoff:
            return True
        else:
            return False

    @classmethod
    def is_donor(cls, feature):
        if feature.attributes_data[cls.grouping_attribute] == cls.donor_key:
            return True
        else:
            return False

    @classmethod
    def is_acceptor(cls, feature):
        if feature.attributes_data[cls.grouping_attribute] == cls.acceptor_key:
            return True
        else:
            return False


class HydrogenBondInx(Interaction):
    """Substantiates HydrogenBondType by selecting donor and acceptor
    features, as well as the involved Hydrogen atom.

    """

    def __init__(self, donor=None, H=None, acceptor=None,
                 distance_cutoff=mastinxconfig.HBOND_DIST_MAX,
                 angle_cutoff=mastinxconfig.HBOND_DON_ANGLE_MIN,
                 interaction_class=None):

        donor_atom = donor.atoms[0]
        acceptor_atom = acceptor.atoms[0]
        okay, distance, angle = HydrogenBondType.check(donor_atom, H, acceptor_atom,
                                                       distance_cutoff=distance_cutoff,
                                                       angle_cutoff=angle_cutoff)
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
        self._interaction_class = interaction_class

    @property
    def donor(self):
        """The donor Feature in the hydrogen bond."""
        return self._donor

    @property
    def H(self):
        """The donated hydrogen Atom in the hydrogen bond."""
        return self._H

    @property
    def acceptor(self):
        """The acceptor Feature in the hydrogen bond."""
        return self._acceptor

    @property
    def distance(self):
        """The distance between the donor atom and the acceptor atom."""
        return self._distance

    @property
    def angle(self):
        """The angle (in degrees) between the donor atom, hydrogen atom, and
        acceptor atom with the hydrogen atom as the vertex.

        """
        return self._angle

    @property
    def interaction_class(self):
        return self._interaction_class

    @interaction_class.setter
    def interaction_class(self, inx_class):
        # TODO type checking
        self._interaction_class = inx_class

    @property
    def record(self):
        record_fields = ['interaction_class',
                         'donor_coords', 'acceptor_coords',
                         'distance', 'angle',
                         'H_coords']

        HydrogenBondInxRecord = namedtuple('HydrogenBondInxRecord', record_fields)
        record_attr = {'interaction_class' : self.interaction_class.name}
        record_attr['donor_coords'] = self.donor.atoms[0].coords
        record_attr['acceptor_coords'] = self.acceptor.atoms[0].coords
        record_attr['distance'] = self.distance
        record_attr['angle'] = self.angle
        record_attr['H_coords'] = self.H.coords

        return HydrogenBondInxRecord(**record_attr)

    def pp(self):

        string = ("{self_str}\n"
                  "donor: {donor_el}\n"
                  "       coords = {donor_coords}\n"
                  "       pdb_serial = {donor_serial}\n"
                  "       pdb_residue_name = {donor_resname}\n"
                  "       pdb_residue_index = {donor_resnum}\n"
                  "acceptor: {acceptor_el}\n"
                  "          coords = {acceptor_coords}\n"
                  "          pdb_serial = {acceptor_serial}\n"
                  "          pdb_residue_name = {acceptor_resname}\n"
                  "          pdb_residue_index = {acceptor_resnum}\n"
                  "distance: {distance}\n"
                  "angle: {angle}\n").format(
                      self_str=str(self),
                      donor_el=self.donor.atom_type.element,
                      donor_coords=tuple(self.donor.coords),
                      donor_mol_name=self.donor.molecule.molecule_type.name,
                      donor_serial=self.donor.atom_type.pdb_serial_number,
                      donor_resname=self.donor.atom_type.pdb_residue_name,
                      donor_resnum=self.donor.atom_type.pdb_residue_number,
                      acceptor_el=self.acceptor.atom_type.element,
                      acceptor_coords=tuple(self.acceptor.coords),
                      acceptor_mol_name=self.acceptor.molecule.molecule_type.name,
                      acceptor_serial=self.acceptor.atom_type.pdb_serial_number,
                      acceptor_resname=self.acceptor.atom_type.pdb_residue_name,
                      acceptor_resnum=self.acceptor.atom_type.pdb_residue_number,
                      distance=self.distance,
                      angle=self.angle,)

        print(string)

    def pickle(self, path):
        import sys
        sys.setrecursionlimit(10000)
        import pickle
        with open(path, 'wb') as wf:
            pickle.dump(self, wf)

