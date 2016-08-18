"""The PiStacking module defines the PiStackingType and
PiStackingInx for explicit hydrogens.

"""


import itertools as it
from collections import namedtuple

import numpy as np
import numpy.linalg as la

import mast.config.interactions as mastinxconfig
from mast.interactions.interactions import InteractionType, Interaction, InteractionError

class PiStackingType(InteractionType):
    """Defines an InteractionType class for hydrogen bonds between members
    with explicit hydrogens.

    """

    attributes = {}
    interaction_name = "PiStacking"
    feature_keywords = mastinxconfig.PISTACKING_FEATURES
    aromatic_keys = mastinxconfig.PISTACKING_FEATURES['rdkit_family']
    grouping_attribute = 'rdkit_family'

    def __init__(self, pi_stacking_type_name,
                 feature_types=None,
                 association_type=None,
                 assoc_member_pair_idxs=None,
                 **pi_stacking_attrs):

        super().__init__(pi_stacking_type_name,
                         feature_types=feature_types,
                         association_type=association_type,
                         assoc_member_pair_idxs=assoc_member_pair_idxs,
                         **pi_stacking_attrs)

    @classmethod
    def find_hits(cls, member_a, member_b):

        # check that the keys ar okay in parent class
        # super().find_hits(members_features)

        # for each member collect the grouped features
        # initialize list of members (A, B)
        members_features = ([], [])
        for i, member in enumerate([member_a, member_b]):
            for feature_key, feature in member.features.items():
                # get groupby attribute to use as a key
                group_attribute = feature.feature_type.attributes_data[cls.grouping_attribute]

                if group_attribute in cls.aromatic_keys:
                    aromatic_tup = (feature_key, feature)
                    members_features[i].append(aromatic_tup)

        # pair the aromatic features
        aromatic_pairs = []
        pi_stack_pairs = it.product(members_features[0],
                                    members_features[1])

        # scan the pairs for hits
        hit_pair_keys = []
        pistacks = []
        for arom_a_tup, arom_b_tup in pi_stack_pairs:
            arom_a_feature_key = atom_a_tup[0]
            arom_a_feature = donor_tup[1]
            arom_b_feature_key = arom_b_tup[0]
            arom_b_feature = arom_b_tup[1]
            # try to make a PiStackingInx object, which calls check,
            #
            # OPTIMIZATION: otherwise we have to call check first then
            # the PiStackingInx constructor will re-call check to
            # get the angle and distance. If we allow passing and not
            # checking the angle and distance in the constructor then
            # it would be faster, however I am not going to allow that
            # in this 'safe' InteractionType, an unsafe optimized
            # version can be made separately if desired.
            try:
                pistack = PiStackingInx(feature1=arom_a_feature,
                                        feature2=arom_b_feature)
            # else continue to the next pairing
            except InteractionError:
                continue
            # if it succeeds add it to the list of H-Bonds
            pistacks.append(pistack)
            # and the feature keys to the feature key pairs
            hit_pair_keys.append((arom_a_feature_key, arom_b_feature_key))

        return hit_pair_keys, pistacks

    @classmethod
    def check(cls, arom_a_atoms, arom_b_atoms):

        from scipy.spatial.distance import cdist

        arom_a_coords = np.array([atom.coords for atom in arom_a_atoms])
        arom_b_coords = np.array([atom.coords for atom in arom_b_atoms])
        arom_coords = [arom_a_coords, arom_b_coords]

        # 1) calculate the distance between centroids
        centroid_a = atom_a_coords.mean(axis=1)
        centroid_b = atom_b_coords.mean(axis=1)
        centroids = [centroid_a, centroid_b]
        centroid_distance = cdist(centroid_a, centroid_b)[0,0]
        # if this passes then move on
        if cls.check_centroid_distance(distance) is False:
            return (False, centroid_distance, None)

        # 2) determine whether it is parallel or perpendicular stacking

        # 2.1) calculate the normal vectors of the rings by using
        # vectors from the centroid to 2 different points on the ring
        arom_vectors = []
        arom_norms = []
        for i, atom_coords in enumerate(arom_coords):
            # choose the atoms
            a0 = atom_coords[0]
            if len(atom_coords) in [6,5]:
                a1 = 3
            else:
                raise InteractionError("aromatic rings without 5 or 6 atoms not supported")

            a0c = a0 - centroid[i]
            arom_vectors.append(a0c)
            a1c = a1 - centroid
            norm = a0c.cross(a1c)
            arom_norms.append(norm)

        # 2.2) calculate the angle between the normal vectors
        try:
            ring_normal_angle = np.degrees(np.arccos(
                np.dot(arom_norms[0], arom_norms[1])/(la.norm(
                    arom_norms[0]) * la.norm(arom_norms[1]))))
        except RuntimeWarning:
            print("v1: {0} \n"
                  "v2: {1}".format(arom_norms[0], arom_norms[1]))

        # 2.3) check the ring normal angles, we expect a string if it
        # passed or False if it failed
        ring_normal_result = cls.check_ring_normal_angle(ring_normal_angle)
        if ring_normal_result is False:
            return (False, distance, ring_normal_angle, )

        # 3) Depending on the type of pi-stacking make further checks
        # see Bahrach 2014 Computational Organic Chemistry for details

        # 3.p) The stacking is parallel
        elif result == 'parallel':
            # There are 3 different parallel possibilities for 6 member rings:
            #  - parallel stacked, yaw parallel
            #     77s from Bahrach
            #  - parallel stacked, yaw perpendicular
            #     77s' from Bahrach
            #  - parallel displaced, yaw parallel
            #     77pd from Bahrach

            # 3.p.1) calculate the angle betwee yaw vectors for each 6
            # member ring

            # calculate the projected vector from arom_b to arom_a plane
            arom_b_proj_vector = arom_vectors[0] * (arom_vectors[1].dot(arom_vectors[0]))
            # calculate the yaw angle
            yaw_angle = np.degrees(np.arccos(
                np.dot(arom_vectors[0], arom_b_proj_vector)/(la.norm(
                    arom_vectors[0]) * la.norm(arom_b_proj_vector))))

            return (True, distance, ring_normal_angle, yaw_angle, )
            # yaw_result = cls.check_yaw(yaw_angle)
            # if yaw_result is False:
            #     return (False, distance, ring_normal_angle, yaw_angle, )
            # else:
            #     return (True, distance, ring_normal_angle, yaw_angle, )
            # 3.p.2) for either parallel or perpendicular yaw
            # elif yaw_result == 'parallel-stacked':
            #     pass
            # elif yaw_result == 'parallel-displaced':
            #     pass
            # elif yaw_result == 'perpendicular':
            #     pass

        # 3.t) The stacking is perpendicular (T-stacking)
        elif result == 'perpendicular':
            # 3.t.1) find which ring bisects the other by getting the
            # distances from each atom of the rings to the centroid
            # the other ring
            arom_a_ds = cdist(arom_a_coords, [centroid_b])
            arom_b_ds = cdist(arom_b_coords, [centroid_a])
            # take the minimum from each comparison
            arom_a_min_d = min(arom_a_ds)
            arom_b_min_d = min(arom_b_ds)
            # the one that is closer is the bisecting aromatic ring
            mins = np.array([arom_a_min_d, arom_b_min_d])
            T_distance = mins.min()
            # if this meets the distance requirement we save the index
            # of (0 or 1) based on input, or return a False tuple
            T_distance_result = cls.check_T_distance(T_distance)
            if T_distance_result is False:
                return (False, distance, ring_normal_angle, T_distance, )
            if T_distance_result is True:
                return (True, distance, ring_normal_angle, T_distance, )
            # TODO further refinement
            else:
                bisecting_idx = mins.argmin()
                # There are 2 different perpendicular possibilities
                #  - rafter (the plane of the bisecting ring bisects across atoms)
                #     77t from Bahrach
                #  - perlin (the plane of the bisecting ring bisects across bonds)
                #     77t' from Bahrach

                # to detect this we project the line from the closest
                # atom to an adjacent atom onto the bisected ring and
                # calculate the yaw angle
                # TODO

        else:
            raise InteractionError("unknown result from check_ring_normal_angle")

    @classmethod
    def check_centroid_distance(cls, distance):
        """For a float distance checks if it is less than the configuration
        file HBOND_DIST_MAX value.

        """
        if distance < mastinxconfig.PISTACK_CENTROID_DISTANCE_MAX:
            return True
        else:
            return False

    @classmethod
    def check_ring_normal_angle(cls, angle):
        """For a float distance checks if it is less than the configuration
        file HBOND_DON_ANGLE_MIN value.

        """
        dev = mastinxconfig.PISTACK_ANGLE_DEVIATION
        if (angle > 180.0 - dev and angle < 180.0 + dev) or \
           (angle > 360.0 - dev and angle < 0.0 + dev):
            return 'parallel'
        elif (angle > 90.0 - dev and angle < 90.0 + dev):
            return 'perpendicular'
        else:
            return False


    @classmethod
    def check_T_distance(cls, distance):
        """For a float distance checks if it is less than the configuration
        file HBOND_DON_ANGLE_MIN value.

        """

        if distance > mastinxconfig.PISTACK_T_DIST:
            return True
        else:
            return False

    @classmethod
    def check_yaw(cls, angle):
        """For a float distance checks if it is less than the configuration
        file HBOND_DON_ANGLE_MIN value.

        """

        if angle > mastinxconfig.HBOND_DON_ANGLE_MIN:
            return True
        else:
            return False

    @property
    def record(self):
        record_fields = ['interaction_class', 'interaction_type',
                         'association_type', 'assoc_member_pair_idxs',
                         'donor_feature_type', 'acceptor_feature_type'] + \
                         list(self.attributes_data.keys())
        PiStackingTypeRecord = namedtuple('PiStackingTypeRecord', record_fields)
        record_attr = {'interaction_class' : self.name}
        record_attr['interaction_type'] = self.interaction_name
        record_attr['association_type'] = self.association_type.name
        record_attr['assoc_member_pair_idxs'] = self.assoc_member_pair_idxs
        record_attr['donor_feature_type'] = self.feature_types[0]
        record_attr['acceptor_feature_type'] = self.feature_types[1]

        return PiStackingTypeRecord(**record_attr)

    def pdb_serial_output(self, inxs, path, delim=","):
        """Output the pdb serial numbers (index in pdb) of the donor and
        acceptor in each HBond to a file like:

        donor_1, acceptor_1
        donor_2, acceptor_2
        ...

        """

        with open(path, 'w') as wf:
            for inx in inxs:
                wf.write("{0}{1}{2}\n".format(inx.donor.atom_type.pdb_serial_number,
                                              delim,
                                              inx.acceptor.atom_type.pdb_serial_number))


class PiStackingInx(Interaction):
    """Substantiates PiStackingType by selecting donor and acceptor
    features, as well as the involved Hydrogen atom.

    """

    def __init__(self, feature1=None, feature2=None):

        okay, distance, angle = PiStackingType.check(feature1.atoms, feature2.atoms)
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
                         interaction_type=PiStackingType,
                         system=atom_system)
        self._donor = donor
        self._H = H
        self._acceptor = acceptor
        self._distance = distance
        self._angle = angle

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
    def record(self):
        record_fields = ['interaction_class',
                         'donor_coords', 'acceptor_coords',
                         'distance', 'angle',
                         'H_coords']

        PiStackingInxRecord = namedtuple('PiStackingInxRecord', record_fields)
        record_attr = {'interaction_class' : self.interaction_class.name}
        record_attr['donor_coords'] = self.donor.atoms[0].coords
        record_attr['acceptor_coords'] = self.acceptor.atoms[0].coords
        record_attr['distance'] = self.distance
        record_attr['angle'] = self.angle
        record_attr['H_coords'] = self.H.coords

        return PiStackingInxRecord(**record_attr)

    def pickle(self, path):
        import sys
        sys.setrecursionlimit(10000)
        import pickle
        with open(path, 'wb') as wf:
            pickle.dump(self, wf)

