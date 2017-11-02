"""The PiCation module defines the PiCationType and
PiCationInx for explicit hydrogens.

"""


import itertools as it
from collections import namedtuple

import numpy as np
import numpy.linalg as la
from scipy.spatial.distance import cdist, euclidean

import mastic.config.interactions as masticinxconfig
from mastic.interactions.interactions import InteractionType, Interaction, InteractionError

class PiCationType(InteractionType):
    """Defines an InteractionType class for hydrogen bonds between members
    with explicit hydrogens.

    """

    attributes = {}
    interaction_name = "PiCation"
    feature_keys = masticinxconfig.PICATION_FEATURE_KEYS
    feature_classifiers = masticinxconfig.PICATION_FEATURES
    degree = masticinxconfig.PICATION_COMMUTATIVITY
    commutative = masticinxconfig.PICATION_COMMUTATIVITY
    interaction_param_keys = masticinxconfig.PICATION_PARAM_KEYS

    # parameters set from the config file
    centroid_max_distance = masticinxconfig.PICATION_CENTROID_DIST_MAX
    centroid_offset_max = masticinxconfig.PICATION_OFFSET_MAX
    amine_normal_angle_min = masticinxconfig.PICATION_AMINE_NORMAL_ANGLE_MIN
    tertiary_amine_feature_classifiers = masticinxconfig.PICATION_TERTIARY_AMINE_FEATURE
    heavy_atoms = masticinxconfig.PICATION_HEAVY_ATOMS_ELEMENT_SYMBOLS

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
        return PiCationInx(*params, **kwargs)

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
    def check(cls, arom, cation):

        features = [arom, cation]
        feature_tests = [cls.test_features_centroid_distance,
                         cls.test_features_centroid_offset,
                         cls.test_features_tertiary_amine]

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
    def check_centroid_offset_distance(cls, distance):
        """For a float distance checks if it is less than the configuration
        file HBOND_DIST_MAX value.

        """
        if distance <= cls.centroid_offset_max:
            return True
        else:
            return False

    @classmethod
    def test_features_centroid_distance(cls, arom, cation):
        arom_heavy_atom_coords = np.array([atom.coords for atom in arom.atoms if
                                             atom.atom_type.element in cls.heavy_atoms])
        arom_centroid = calc_centroid(arom_heavy_atom_coords)
        cation_coords = cation.atoms[0].coords
        centroid_distance = euclidean(arom_centroid, cation_coords)

        if cls.check_centroid_distance(centroid_distance) is False:
            return False, centroid_distance
        else:
            return True, centroid_distance

    @classmethod
    def test_features_centroid_offset(cls, arom, cation):
        arom_heavy_atom_coords = np.array([atom.coords for atom in arom.atoms if
                                           atom.atom_type.element in cls.heavy_atoms])
        cation_coord = cation.atoms[0].coords
        # calculate the centroid offset
        centroid_offset = calc_centroid_offset(arom_a_heavy_atom_coords,
                                               cation_coord)

        if cls.check_centroid_offset_distance(centroid_offset) is False:
            return False, centroid_offset
        else:
            return True, centroid_offset

    @classmethod
    def test_features_tertiary_amine(cls, arom, cation):
        # check for if the cation is a tertiary amine
        if True:
            return True, None

        tertiary_amine_angle = calc_tertiary_amine_angle(cation)

        if cls.check_tertiary_amine_angle(tertiary_amine_angle):
            return False, tertiary_amine_angle
        else:
            return True, tertiary_amine_angle

    @property
    def record(self):
        record_fields = ['interaction_class', 'interaction_type',
                         'association_type', 'assoc_member_pair_idxs',
                         'arom_a_feature_type', 'arom_b_feature_type'] + \
                         list(self.attributes_data.keys())
        PiCationTypeRecord = namedtuple('PiCationTypeRecord', record_fields)
        record_attr = {'interaction_class' : self.name}
        record_attr['interaction_type'] = self.interaction_name
        record_attr['association_type'] = self.association_type.name
        record_attr['assoc_member_pair_idxs'] = self.assoc_member_pair_idxs
        record_attr['arom_a_feature_type'] = self.feature_types[0]
        record_attr['arom_b_feature_type'] = self.feature_types[1]

        return PiCationTypeRecord(**record_attr)

def calc_tertiary_amine_angle(cation):
    angle = None
    return angle

def calc_centroid_offset(arom_coords, other_centroid):
    arom_centroid = calc_centroid(arom_a_coords)
    norm = calc_arom_norm(arom_coords)
    # get the direction facing the other centroid

def calc_aroms_centroid_offset(arom_a_coords, arom_b_coords):
    """Project the centers of each ring over each other and get
    the offset"""

    # get the centroid coordinates
    centroid_a, centroid_b = [calc_centroid(arom_coords) for
                              arom_coords in (arom_a_coords, arom_b_coords)]

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

def calc_facing_vector(vec_up, point):
    vec_down = -1 * vec_up

    d_up = euclidean(vec_up, point)
    d_down = euclidean(vec_down, point)

    face_vec = vec_up if d_up < d_down else vec_down
    return face_vec

def calc_arom_facing_norms(arom_a_coords, arom_b_coords):
    """Given two aromatic rings get the normal vectors that face the other ring"""

    centroids = [calc_centroid(arom_coords) for arom_coords in [arom_a_coords, arom_b_coords]]
    arom_norms = calc_arom_norms(arom_a_coords, arom_b_coords)

    face_norms = []
    for i, arom_norm in enumerate(arom_norms):
        # get the index of the other arom
        j = 1 if i ==0 else 0
        norm = calc_facing_vector(arom_norm + centroids[i], centroids[j])
        # norm_up = arom_norm
        # norm_down = -1 * arom_norm
        # # get the norm so that it points to the other ring
        # d_up = euclidean(norm_up + centroids[i], centroids[j])
        # d_down = cdist(norm_down + centroids[i], centroids[j])
        # norm = norm_up if d_up < d_down else norm_down
        face_norms.append(norm)

    return face_norms

def calc_centroid(atom_coords):
    return atom_coords.mean(axis=0)

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

def calc_arom_norm(arom_coords):
    centroid = calc_centroid(arom_coords)

    a0 = arom_coords[0]
    if len(arom_coords) in [6,5]:
        a1 = arom_coords[2]
    else:
        raise InteractionError("aromatic rings without 5 or 6 atoms not supported")

    a0c = a0 - centroid
    a1c = a1 - centroid
    norm = np.cross(a0c, a1c)
    return norm

def calc_angle(v1, v2):
    try:
        # flip one of them because it is opposite the other
        angle = np.degrees(np.arccos(
            np.dot(v1, v2)/(la.norm(
                v1) * la.norm(v2))))
    except RuntimeWarning:
        print("v1: {0} \n"
              "v2: {1}".format(v1, v2))

    return angle

def calc_arom_normal_angle(arom_a_coords, arom_b_coords):

    # get the normal vectors
    arom_norm_a, arom_norm_b = calc_arom_norms(arom_a_coords, arom_b_coords)
    arom_norm_b = -1 * arom_norm_b
    ring_normal_angle = calc_angle(arom_norm_a, arom_norm_b)

    # if normals happen to be opposite directions correct and get
    # the angle that is non-negative and smallest
    alt_angle = 180 - ring_normal_angle
    ring_normal_angle = min(ring_normal_angle, alt_angle) if not\
                        alt_angle < 0 else ring_normal_angle

    return ring_normal_angle


class PiCationInx(Interaction):
    """Substantiates PiCationType by selecting donor and acceptor
    features, as well as the involved Hydrogen atom.

    """

    interaction_type = PiCationType

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

        PiCationInxRecord = namedtuple('PiCationInxRecord', record_fields)
        record_attr = {'interaction_class' : self.interaction_class.name}

        return PiCationInxRecord(**record_attr, **self.interaction_params)
