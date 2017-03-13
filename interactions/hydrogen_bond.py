"""The HydrogenBond module defines the HydrogenBondType and
HydrogenBondInx for explicit hydrogens.

"""
import itertools as it
from collections import namedtuple, defaultdict

import numpy as np
import numpy.linalg as la
from scipy.spatial.distance import cdist

import mast.config.interactions as mastinxconfig
import mast.config.features as mastfeatconfig
from mast.interactions.interactions import InteractionType, Interaction, InteractionError


class HydrogenBondType(InteractionType):
    """Defines an InteractionType class for hydrogen bonds between members
    with explicit hydrogens.

    """

    ## class attributes that need to exist
    attributes = {}
    interaction_name = "HydrogenBond"
    feature_keys = mastinxconfig.HYDROGEN_BOND_FEATURE_KEYS
    feature_classifiers = mastinxconfig.HYDROGEN_BOND_FEATURES
    # degree is the number of features that participate in an interaction
    degree = mastinxconfig.HYDROGEN_BOND_DEGREE
    commutative = mastinxconfig.HYDROGEN_BOND_COMMUTATIVITY
    interaction_param_keys = mastinxconfig.HYDROGEN_BOND_PARAM_KEYS

    ## specific to this class parameters but make defaults easier and
    ## for writing other similar InteractionTypes
    distance_cutoff = mastinxconfig.HYDROGEN_BOND_DIST_MAX
    angle_cutoff = mastinxconfig.HYDROGEN_BOND_DON_ANGLE_MIN

    ## convenience class attributes particular to this class
    donor_key = feature_keys[0]
    acceptor_key = feature_keys[1]
    donor_feature_classifiers = feature_classifiers[donor_key]
    acceptor_feature_classifiers = feature_classifiers[acceptor_key]


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

        self.donor_type = feature_types[0]
        self.acceptor_type = feature_types[1]

    @staticmethod
    def interaction_constructor(*params, **kwargs):
        return HydrogenBondInx(*params, **kwargs)

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
    def check(cls, donor, acceptor):
        """Checks if the 3 atoms qualify as a hydrogen bond. Returns a tuple
        (bool, float, float) where the first element is whether or not it
        qualified, the second and third are the distance and angle
        respectively. Angle may be None if distance failed to qualify.

        Compatible with RDKit Acceptor and Donor features

        """

        # assemble the features and their tests
        features = [donor, acceptor]
        feature_tests = [cls.test_features_distance, cls.test_features_angle]
        # pass to parent function, this returns a results tuple of the
        # form: (okay, (param_values))
        return super().check(features, feature_tests)

    @classmethod
    def test_features_distance(cls, donor, acceptor):
        donor_atom = donor.atoms[0]
        acceptor_atom = acceptor.atoms[0]

        # calculate the distance
        distance = cls.calc_distance(donor_atom, acceptor_atom)

        # check it
        # if it fails return a false okay and the distance
        if cls.check_distance(distance) is False:
            return False, distance
        # otherwise return that it was okay
        else:
            return True, distance

    @classmethod
    def test_features_angle(cls, donor, acceptor):
        donor_atom = donor.atoms[0]
        acceptor_atom = acceptor.atoms[0]

        # if the distance passes we want to check the angle, which we
        # will need the coordinates of the adjacent hydrogens to the donor
        h_atoms = [atom for atom in donor_atom.adjacent_atoms
                   if atom.atom_type.element == 'H']

        # check to see if even 1 hydrogen atom satisfies the angle
        # constraint
        okay_angle = None
        h_atoms_iter = iter(h_atoms)
        h_atoms_angles = []
        # if it doesn't the loop will end and a false okay and the bad
        # angles will be returned
        while okay_angle is None:
            try:
                h_atom = next(h_atoms_iter)
            # none are found to meet the constraint
            except StopIteration:
                return False, tuple(h_atoms_angles)

            # calculate the angle for this hydrogen
            angle = cls.calc_angle(donor_atom, acceptor_atom, h_atom)

            # save the angle
            h_atoms_angles.append(angle)

            # check if the angle meets the constraints
            if cls.check_angle(angle) is True:
                okay_angle = angle

        # if it succeeds in finding a good angle return the first one
        return True, okay_angle

    #### Hydrogen Bond Specific methods
    # i.e. not necessarily found in other interaction types

    @classmethod
    def check_distance(cls, distance):
        """For a float distance checks if it is less than the configuration
        file HYDROGEN_BOND_DIST_MAX value.

        """
        if distance < cls.distance_cutoff:
            return True
        else:
            return False

    @classmethod
    def check_angle(cls, angle):
        """For a float distance checks if it is less than the configuration
        file HYDROGEN_BOND_DON_ANGLE_MIN value.

        """

        if angle > cls.angle_cutoff:
            return True
        else:
            return False

    @classmethod
    def is_donor(cls, feature):
        if feature.attributes_data[cls.grouping_attribute] in cls.donor_keys:
            return True
        else:
            return False

    @classmethod
    def is_acceptor(cls, feature):
        if feature.attributes_data[cls.grouping_attribute] == cls.acceptor_key:
            return True
        else:
            return False

    @staticmethod
    def calc_distance(donor_atom, acceptor_atom):
        return cdist(np.array([donor_atom.coords]), np.array([acceptor_atom.coords]))[0,0]

    @staticmethod
    def calc_angle(donor_atom, acceptor_atom, h_atom):
        v1 = donor_atom.coords - h_atom.coords
        v2 = acceptor_atom.coords - h_atom.coords
        return np.degrees(np.arccos(np.dot(v1, v2)/(la.norm(v1) * la.norm(v2))))

    @property
    def record(self):
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

# HydrogenBondTypeRecord
_hydrogen_bond_type_record_fields = ['interaction_class', 'interaction_type',
                                     'association_type', 'assoc_member_pair_idxs',
                                     'donor_feature_type', 'acceptor_feature_type']
HydrogenBondTypeRecord = namedtuple('HydrogenBondTypeRecord', _hydrogen_bond_type_record_fields)

class HydrogenBondInx(Interaction):
    """Substantiates HydrogenBondType by selecting donor and acceptor
    features, as well as the involved Hydrogen atom.

    """

    interaction_type = HydrogenBondType

    def __init__(self, donor, acceptor,
                 check=True,
                 interaction_class=None,
                 **param_values):

        donor_atom = donor.atoms[0]
        acceptor_atom = acceptor.atoms[0]
        if check:
            okay, param_values = self.interaction_type.check(donor_atom, acceptor_atom)

            if not okay:
                raise InteractionError

        # success, finish creating interaction
        atom_system = donor.system
        super().__init__(features=[donor, acceptor],
                         interaction_type=self.interaction_type,
                         system=atom_system,
                         interaction_class=interaction_class,
                         **param_values)
        self._donor = donor
        self._acceptor = acceptor

    @property
    def donor(self):
        """The donor Feature in the hydrogen bond."""
        return self._donor

    # TODO implement a way to find all the H atoms that satisfy the interaction
    @property
    def H(self):
        """The donated hydrogen Atom in the hydrogen bond."""
        raise NotImplementedError
        # return self._H

    @property
    def acceptor(self):
        """The acceptor Feature in the hydrogen bond."""
        return self._acceptor

    @property
    def record(self):
        record_attr = {'interaction_class' : self.interaction_class.name}
        record_attr['donor_coords'] = self.donor.atoms[0].coords
        record_attr['acceptor_coords'] = self.acceptor.atoms[0].coords
        # TODO because the H might be ambiguous, see the H property
        # record_attr['H_coords'] = self.H.coords

        return HydrogenBondInxRecord(**record_attr, **self.interaction_params)

# HydrogenBondInxRecord
_hydrogen_bond_inx_record_fields = ['interaction_class',
                                    'donor_coords', 'acceptor_coords'] + \
                                    HydrogenBondType.interaction_param_keys
                                    #'H_coords']
HydrogenBondInxRecord = namedtuple('HydrogenBondInxRecord', _hydrogen_bond_inx_record_fields)



# def pp(self):

#     string = ("{self_str}\n"
#               "donor: {donor_el}\n"
#               "       coords = {donor_coords}\n"
#               "       pdb_serial = {donor_serial}\n"
#               "       pdb_residue_name = {donor_resname}\n"
#               "       pdb_residue_index = {donor_resnum}\n"
#               "acceptor: {acceptor_el}\n"
#               "          coords = {acceptor_coords}\n"
#               "          pdb_serial = {acceptor_serial}\n"
#               "          pdb_residue_name = {acceptor_resname}\n"
#               "          pdb_residue_index = {acceptor_resnum}\n"
#               "distance: {distance}\n"
#               "angle: {angle}\n").format(
#                   self_str=str(self),
#                   donor_el=self.donor.atom_type.element,
#                   donor_coords=tuple(self.donor.coords),
#                   donor_mol_name=self.donor.molecule.molecule_type.name,
#                   donor_serial=self.donor.atom_type.pdb_serial_number,
#                   donor_resname=self.donor.atom_type.pdb_residue_name,
#                   donor_resnum=self.donor.atom_type.pdb_residue_number,
#                   acceptor_el=self.acceptor.atom_type.element,
#                   acceptor_coords=tuple(self.acceptor.coords),
#                   acceptor_mol_name=self.acceptor.molecule.molecule_type.name,
#                   acceptor_serial=self.acceptor.atom_type.pdb_serial_number,
#                   acceptor_resname=self.acceptor.atom_type.pdb_residue_name,
#                   acceptor_resnum=self.acceptor.atom_type.pdb_residue_number,
#                   distance=self.distance,
#                   angle=self.angle,)

#     print(string)
