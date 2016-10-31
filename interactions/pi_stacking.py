"""The PiStacking module defines the PiStackingType and
PiStackingInx for explicit hydrogens.

"""


import itertools as it
from collections import namedtuple

import numpy as np
import numpy.linalg as la
from scipy.spatial.distance import cdist

import mast.config.interactions as mastinxconfig
from mast.interactions.interactions import InteractionType, Interaction, InteractionError

class PiStackingType(InteractionType):
    """Defines an InteractionType class for hydrogen bonds between members
    with explicit hydrogens.

    """

    attributes = {}
    interaction_name = "PiStacking"
    feature_keys = mastinxconfig.PISTACKING_FEATURE_KEYS
    feature_classifiers = mastinxconfig.PISTACKING_FEATURES
    degree = 2
    commutative = True
    interaction_param_keys = mastinxconfig.PISTACKING_PARAM_KEYS

    # parameters set from the config file
    centroid_max_distance = mastinxconfig.PISTACKING_CENTROID_DIST_MAX
    ring_normal_angle_deviation = mastinxconfig.PISTACKING_ANGLE_DEVIATION
    centroid_offset_max = mastinxconfig.PISTACKING_OFFSET_MAX
    heavy_atoms = mastinxconfig.PISTACKING_HEAVY_ATOMS_ELEMENT_SYMBOLS

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

    @staticmethod
    def interaction_constructor(*params, **kwargs):
        return PiStackingInx(*params, **kwargs)

    @classmethod
    def find_hits(cls, members,
                  interaction_classes=None,
                  return_feature_keys=False,
                  return_failed_hits=False):

        # TODO value checks

        # scan the pairs for hits and assign interaction classes if given
        return super().find_hits(members,
                                 interaction_classes=interaction_classes,
                                 return_feature_keys=return_feature_keys,
                                 return_failed_hits=return_failed_hits)

    @classmethod
    def check(cls, arom_a, arom_b):

        features = [arom_a, arom_b]
        feature_tests = [cls.test_features_centroid_distance,
                         cls.test_features_ring_normal_angle,
                         cls.test_features_centroid_offset,
                         cls.test_features_stacking_type]

        return super().check(features, feature_tests)



    @classmethod
    def check_centroid_distance(cls, distance):
        """For a float distance checks if it is less than the configuration
        file HBOND_DIST_MAX value.

        """
        if distance <= cls.centroid_max_distance:
            return True
        else:
            return False

    @classmethod
    def check_ring_normal_angle(cls, angle):
        """For a float distance checks if it is less than the configuration
        file HBOND_DON_ANGLE_MIN value.

        """
        dev = cls.ring_normal_angle_deviation
        if 0 <= angle <= dev or 90 - dev <= angle <= 90 + dev:
            return True
        else:
            return False

    @classmethod
    def check_centroid_offset_distance(cls, distance):
        """For a float distance checks if it is less than the configuration
        file HBOND_DIST_MAX value.

        """
        if distance <= cls.centroid_offset_max:
            return True
        else:
            return False

    @classmethod
    def check_stacking_type(cls, angle):
        dev = cls.ring_normal_angle_deviation
        if 0.0 <= angle <= dev:
            return 'parallel'
        elif 90 - dev <= angle <= 90 + dev:
            return 'perpendicular'
        else:
            return None

    @classmethod
    def test_features_centroid_distance(cls, arom_a, arom_b):
        arom_a_heavy_atom_coords = np.array([atom.coords for atom in arom_a.atoms if
                                             atom.atom_type.element in cls.heavy_atoms])
        arom_b_heavy_atom_coords = np.array([atom.coords for atom in arom_b.atoms if
                                             atom.atom_type.element in cls.heavy_atoms])

        centroid_distance = calc_centroid_distance(arom_a_heavy_atom_coords,
                                                   arom_b_heavy_atom_coords)

        # if this passes then move on
        if cls.check_centroid_distance(centroid_distance) is False:
            return False, centroid_distance
        else:
            return True, centroid_distance

    @classmethod
    def test_features_ring_normal_angle(cls, arom_a, arom_b):
        arom_a_heavy_atom_coords = np.array([atom.coords for atom in arom_a.atoms if
                                             atom.atom_type.element in cls.heavy_atoms])
        arom_b_heavy_atom_coords = np.array([atom.coords for atom in arom_b.atoms if
                                             atom.atom_type.element in cls.heavy_atoms])
        # calculate the normal angle
        ring_normal_angle = calc_arom_normal_angle(arom_a_heavy_atom_coords,
                                                   arom_b_heavy_atom_coords)

        if cls.check_ring_normal_angle(ring_normal_angle) is False:
            return False, ring_normal_angle
        else:
            return True, ring_normal_angle

    @classmethod
    def test_features_centroid_offset(cls, arom_a, arom_b):
        arom_a_heavy_atom_coords = np.array([atom.coords for atom in arom_a.atoms if
                                             atom.atom_type.element in cls.heavy_atoms])
        arom_b_heavy_atom_coords = np.array([atom.coords for atom in arom_b.atoms if
                                             atom.atom_type.element in cls.heavy_atoms])
        # calculate the centroid offset
        centroid_offset = calc_centroid_offset(arom_a_heavy_atom_coords,
                                               arom_b_heavy_atom_coords)

        if cls.check_centroid_offset_distance(centroid_offset) is False:
            return False, centroid_offset
        else:
            return True, centroid_offset

    @classmethod
    def test_features_stacking_type(cls, arom_a, arom_b):
        arom_a_heavy_atom_coords = np.array([atom.coords for atom in arom_a.atoms if
                                             atom.atom_type.element in cls.heavy_atoms])
        arom_b_heavy_atom_coords = np.array([atom.coords for atom in arom_b.atoms if
                                             atom.atom_type.element in cls.heavy_atoms])
        # Determine whether the stacking is parallel of perpendicular, as
        # a string
        ring_normal_angle = calc_arom_normal_angle(arom_a_heavy_atom_coords,
                                                   arom_b_heavy_atom_coords)
        stacking_type = cls.check_stacking_type(ring_normal_angle)
        if stacking_type is None:
            return False, stacking_type
        else:
            return True, stacking_type

    @property
    def record(self):
        record_fields = ['interaction_class', 'interaction_type',
                         'association_type', 'assoc_member_pair_idxs',
                         'arom_a_feature_type', 'arom_b_feature_type'] + \
                         list(self.attributes_data.keys())
        PiStackingTypeRecord = namedtuple('PiStackingTypeRecord', record_fields)
        record_attr = {'interaction_class' : self.name}
        record_attr['interaction_type'] = self.interaction_name
        record_attr['association_type'] = self.association_type.name
        record_attr['assoc_member_pair_idxs'] = self.assoc_member_pair_idxs
        record_attr['arom_a_feature_type'] = self.feature_types[0]
        record_attr['arom_b_feature_type'] = self.feature_types[1]

        return PiStackingTypeRecord(**record_attr)


def calc_centroid_offset(arom_a_coords, arom_b_coords):
    """Project the centers of each ring over each other and get
    the offset"""

    # get the centroid coordinates
    centroid_a, centroid_b = calc_centroids(arom_a_coords, arom_b_coords)
    # get the norms that are facing each other
    face_norm_a, face_norm_b = calc_arom_facing_norms(arom_a_coords, arom_b_coords)

    # the vector going from centroid a to b
    centroid_vec_a = centroid_b - centroid_a
    # vector from b to a
    centroid_vec_b = centroid_a - centroid_b

    # calculate the rejection of the centroid vector on the normal
    # face vector, which is the projection on the plane defined by the
    # normal vector in the direction of the centroid vector
    norm_a_proj = calc_vector_rejection(centroid_vec_a, face_norm_a)
    norm_b_proj = calc_vector_rejection(centroid_vec_b, face_norm_b)

    # compare the magnitudes of the two
    centroid_offset = min(la.norm(norm_a_proj), la.norm(norm_b_proj))

    return centroid_offset

def calc_vector_rejection(vec_a, vec_b):
    """Reject vec_a onto vec_b"""

    projection_vec = calc_vector_projection(vec_a, vec_b)
    # a_2 = a - a_1
    return vec_a - projection_vec

def calc_vector_projection(vec_a, vec_b):
    """Project vec_a onto vec_b."""

    # a_1 = (a dot b_norm) dot b_norm
    vec_b_norm = vec_b / la.norm(vec_b)
    return np.dot(np.dot(vec_a, vec_b_norm), vec_b_norm)

def calc_arom_facing_norms(arom_a_coords, arom_b_coords):
    """Given two aromatic rings get the normal vectors that face the other ring"""

    centroids = calc_centroids(arom_a_coords, arom_b_coords)
    arom_norms = calc_arom_norms(arom_a_coords, arom_b_coords)

    face_norms = []
    for i, arom_norm in enumerate(arom_norms):
        # get the index of the other arom
        j = 1 if i ==0 else 0
        norm_up = arom_norm
        norm_down = -1 * arom_norm
        # get the norm so that it points to the other ring
        d_up = cdist([norm_up + centroids[i]], [centroids[j]])
        d_down = cdist([norm_down + centroids[i]], [centroids[j]])
        norm = norm_up if d_up < d_down else norm_down
        face_norms.append(norm)

    return face_norms

def calc_centroids(arom_a_coords, arom_b_coords):
    centroid_a = arom_a_coords.mean(axis=0)
    centroid_b = arom_b_coords.mean(axis=0)
    centroids = [centroid_a, centroid_b]

    return centroid_a, centroid_b

def calc_centroid_distance(arom_a_coords, arom_b_coords):
    centroid_a, centroid_b = calc_centroids(arom_a_coords, arom_b_coords)
    centroid_distance = cdist([centroid_a], [centroid_b])[0,0]
    return centroid_distance

def calc_arom_norms(arom_a_coords, arom_b_coords):
    # Calculate and check the angle between the two ring normal vectors

    centroids = calc_centroids(arom_a_coords, arom_b_coords)

    arom_norms = []
    for i, arom_coords in enumerate([arom_a_coords, arom_b_coords]):
        # get the coordinates of two atoms on the ring
        a0 = arom_coords[0]
        if len(arom_coords) in [6,5]:
            a1 = arom_coords[2]
        else:
            raise InteractionError("aromatic rings without 5 or 6 atoms not supported")

        # get two plane vectors
        a0c = a0 - centroids[i]
        a1c = a1 - centroids[i]
        norm = np.cross(a0c, a1c)
        arom_norms.append(norm)

    return tuple(arom_norms)

def calc_arom_normal_angle(arom_a_coords, arom_b_coords):

    # get the normal vectors
    arom_norm_a, arom_norm_b = calc_arom_norms(arom_a_coords, arom_b_coords)

    try:
        # flip one of them because it is opposite the other
        ring_normal_angle = np.degrees(np.arccos(
            np.dot(arom_norm_a, -1 * arom_norm_b)/(la.norm(
                arom_norm_a) * la.norm(arom_norm_b))))
    except RuntimeWarning:
        print("v1: {0} \n"
              "v2: {1}".format(arom_norms[0], arom_norms[1]))

    # if normals happen to be opposite directions correct and get
    # the angle that is non-negative and smallest
    alt_angle = 180 - ring_normal_angle
    ring_normal_angle = min(ring_normal_angle, alt_angle) if not\
                        alt_angle < 0 else ring_normal_angle

    return ring_normal_angle



    # @classmethod
    # def check_T_distance(cls, distance):
    #     """For a float distance checks if it is less than the configuration
    #     file HBOND_DON_ANGLE_MIN value.

    #     """

    #     if distance > mastinxconfig.PISTACK_T_DIST:
    #         return True
    #     else:
    #         return False

    # @classmethod
    # def my_check_ring_normal_angle(cls, angle,
    #                             dev=ring_normal_angle_deviation):
    #     """For a float distance checks if it is less than the configuration
    #     file HBOND_DON_ANGLE_MIN value.

    #     """
    #     if (angle > 180.0 - dev and angle < 180.0 + dev) or \
    #        (angle > 360.0 - dev and angle < 0.0 + dev):
    #         return 'parallel'
    #     elif (angle > 90.0 - dev and angle < 90.0 + dev):
    #         return 'perpendicular'
    #     else:
    #         return False

    # @classmethod
    # def my_check(cls, arom_a_atoms, arom_b_atoms):

    #     # parameter initialization for return
    #     centroid_distance = None
    #     ring_normal_angle = None
    #     T_distance = None
    #     proj_centroid_distance = None

    #     # coordinates for atoms of aromatic rings (heavy atoms only)
    #     arom_a_coords = np.array([atom.coords for atom in arom_a_atoms])
    #     arom_b_coords = np.array([atom.coords for atom in arom_b_atoms])
    #     arom_coords = [arom_a_coords, arom_b_coords]

    #     # 1) calculate the distance between centroids
    #     centroid_a = atom_a_coords.mean(axis=1)
    #     centroid_b = atom_b_coords.mean(axis=1)
    #     centroids = [centroid_a, centroid_b]
    #     centroid_distance = cdist(centroid_a, centroid_b)[0,0]
    #     # if this passes then move on
    #     if cls.check_centroid_distance(distance) is False:
    #         return (False, centroid_distance, ring_normal_angle,
    #                 T_distance, proj_centroid_distance,)

    #     # 2) determine whether it is parallel or perpendicular stacking

    #     # 2.1) calculate the normal vectors of the rings by using
    #     # vectors from the centroid to 2 different points on the ring
    #     arom_plane_vectors = []
    #     arom_norms = []
    #     for i, atom_coords in enumerate(arom_coords):
    #         # choose the atoms
    #         a0 = atom_coords[0]
    #         if len(atom_coords) in [6,5]:
    #             a1 = 3
    #         else:
    #             raise InteractionError("aromatic rings without 5 or 6 atoms not supported")

    #         a0c = a0 - centroid[i]
    #         arom_plane_vectors.append(a0c)
    #         a1c = a1 - centroid
    #         norm = a0c.cross(a1c)
    #         arom_norms.append(norm)

    #     # 2.2) calculate the angle between the normal vectors
    #     try:
    #         ring_normal_angle = np.degrees(np.arccos(
    #             np.dot(arom_norms[0], arom_norms[1])/(la.norm(
    #                 arom_norms[0]) * la.norm(arom_norms[1]))))
    #     except RuntimeWarning:
    #         print("v1: {0} \n"
    #               "v2: {1}".format(arom_norms[0], arom_norms[1]))

    #     # 2.3) check the ring normal angles, we expect a string if it
    #     # passed or False if it failed
    #     ring_normal_result = cls.check_ring_normal_angle(ring_normal_angle)
    #     if ring_normal_result is False:
    #         return (False, centroid_distance, ring_normal_angle,
    #                 T_distance, proj_centroid_distance,)

    #     # 3) Project the closest carbon onto the other ring to see if
    #     # it is oriented correctly. A different atom will be use for each

    #     # 3.1) Get the atom to be projected for each condition.
    #     # 3.1.p) The stacking is parallel
    #     ref_arom_idx = None
    #     proj_arom_idx = None
    #     proj_atom_idx = None
    #     # there is an extra parameter that needs checked for the T
    #     T_distance = None
    #     elif result == 'parallel':
    #         # our choice doesn't matter here so arbitrarily choose
    #         # reference ring
    #         ref_arom_idx = 0
    #         proj_arom_idx = 1
    #         # we need to find which atom is closest to the centroid of
    #         # the reference
    #         proj_arom_coords = arom_coords[proj_arom_idx]
    #         ref_centroid = centroids[ref_arom_idx]
    #         proj_arom_ds = cdist(proj_arom_coords, [ref_centroid])
    #         proj_atom_idx = arom_a_ds.argmin()

    #     # 3.1.t) The stacking is perpendicular (T-stacking)
    #     elif result == 'perpendicular':
    #         # 3.1.t.1) find which ring bisects the other by getting the
    #         # distances from each atom of the rings to the centroid
    #         # the other ring
    #         arom_a_ds = cdist(arom_a_coords, [centroid_b])
    #         arom_b_ds = cdist(arom_b_coords, [centroid_a])
    #         arom_ds = [arom_a_ds, arom_b_ds]
    #         # take the minimum from each comparison
    #         arom_a_min_d = min(arom_a_ds)
    #         arom_b_min_d = min(arom_b_ds)
    #         # the one that is closer is the bisecting aromatic ring
    #         mins = np.array([arom_a_min_d, arom_b_min_d])
    #         T_distance = mins.min()
    #         # if this meets the distance requirement we save the index
    #         # of (0 or 1) based on input, or return a False tuple
    #         T_distance_result = cls.check_T_distance(T_distance)
    #         if T_distance_result is False:
    #             return (False, centroid_distance ring_normal_angle,
    #                     T_distance, proj_centroid_distance,)
    #         elif T_distance_result is True:
    #             # set the appropriate reference etc.
    #             ref_arom_idx = mins.argmin()
    #             proj_arom_idx = mins.argmax()
    #             proj_atom_idx = arom_ds[mins.argmin()].argmin()
    #         else:
    #             raise InteractionError("unknown result from check_T_distance")

    #     else:
    #         raise InteractionError("unknown result from check_ring_normal_angle")


    #     # 3.2) project the point to the reference ring plane
    #     proj_point = arom_plane_vectors[ref_arom_idx] * \
    #                  (arom_coords[proj_arom_idx][proj_atom_idx].dot(
    #                      arom_plane_vectors[ref_arom_idx]))
    #     proj_centroid_distance = cdist([proj_point], [centroids[ref_arom_idx]])
    #     offset_result = cls.check_offset_distance(proj_centroid_distance)
    #     if offset_result is False:
    #         return (False, centroid_distance ring_normal_angle,
    #                 T_distance, proj_centroid_distance,)
    #     elif offset_result is True:
    #         return (True, centroid_distance ring_normal_angle,
    #                 T_distance, proj_centroid_distance,)
    #     else:
    #         raise InteractionError("unknown result from check_projection_centroid_distance")

class PiStackingInx(Interaction):
    """Substantiates PiStackingType by selecting donor and acceptor
    features, as well as the involved Hydrogen atom.

    """

    interaction_type = PiStackingType

    def __init__(self, arom_a, arom_b,
                 check=True,
                 interaction_class=None,
                 **param_values):

        if check:
            # use the default settings for the interaction_type only
            # for implicit checks, the idea is that we want the user
            # to mutate the InteractionType to change the
            # classification criteria
            okay, param_values = self.interaction_type.check(arom_a.atoms,
                                                             arom_b.atoms,)

            if not okay:
                raise InteractionError

        # success, finish creating interaction
        atom_system = arom_a.system
        super().__init__(features=[arom_a, arom_b],
                         interaction_type=self.interaction_type,
                         system=atom_system,
                         **param_values)
        self._arom_a = arom_a
        self._arom_b = arom_b

    @property
    def arom_a(self):
        return self._arom_a
    @property
    def arom_b(self):
        return self._arom_b

    @property
    def record(self):
        record_fields = ['interaction_class'] + \
                        self.interaction_type.interaction_param_keys

        PiStackingInxRecord = namedtuple('PiStackingInxRecord', record_fields)
        record_attr = {'interaction_class' : self.interaction_class.name}

        return PiStackingInxRecord(**record_attr, **self.interaction_params)


#### parallel yaw calculations
            # # There are 3 different parallel possibilities for 6 member rings:
            # #  - parallel stacked, yaw parallel
            # #     77s from Bahrach
            # #  - parallel stacked, yaw perpendicular
            # #     77s' from Bahrach
            # #  - parallel displaced, yaw parallel
            # #     77pd from Bahrach

            # # 3.p.1) calculate the angle betwee yaw vectors for each 6
            # # member ring

            # # calculate the projected vector from arom_b to arom_a plane
            # arom_b_proj_vector = arom_vectors[0] * (arom_vectors[1].dot(arom_vectors[0]))
            # # calculate the yaw angle
            # yaw_angle = np.degrees(np.arccos(
            #     np.dot(arom_vectors[0], arom_b_proj_vector)/(la.norm(
            #         arom_vectors[0]) * la.norm(arom_b_proj_vector))))

            # return (True, distance, ring_normal_angle, yaw_angle, )
            # # yaw_result = cls.check_yaw(yaw_angle)
            # # if yaw_result is False:
            # #     return (False, distance, ring_normal_angle, yaw_angle, )
            # # else:
            # #     return (True, distance, ring_normal_angle, yaw_angle, )
            # # 3.p.2) for either parallel or perpendicular yaw
            # # elif yaw_result == 'parallel-stacked':
            # #     pass
            # # elif yaw_result == 'parallel-displaced':
            # #     pass
            # # elif yaw_result == 'perpendicular':
            # #     pass


    # @classmethod
    # def check_yaw(cls, angle):
    #     """For a float distance checks if it is less than the configuration
    #     file HBOND_DON_ANGLE_MIN value.

    #     """

    #     if angle > mastinxconfig.HBOND_DON_ANGLE_MIN:
    #         return True
    #     else:
    #         return False
