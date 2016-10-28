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
                  distance_cutoff=distance_cutoff,
                  angle_cutoff=angle_cutoff):

        # TODO value checks

        # scan the pairs for hits and assign interaction classes if given
        return super().find_hits(members,
                                 interaction_classes=interaction_classes,
                                 # the parameters for the interaction existence
                                 distance_cutoff=distance_cutoff,
                                 angle_cutoff=angle_cutoff)

    @classmethod
    def check(cls, donor, acceptor,
              distance_cutoff=distance_cutoff,
              angle_cutoff=angle_cutoff):
        """Checks if the 3 atoms qualify as a hydrogen bond. Returns a tuple
        (bool, float, float) where the first element is whether or not it
        qualified, the second and third are the distance and angle
        respectively. Angle may be None if distance failed to qualify.

        Compatible with RDKit Acceptor and Donor features

        """


        donor_atom = donor.atoms[0]
        acceptor_atom = acceptor.atoms[0]

        # calculate the distance
        distance = cls.calc_distance(donor_atom, acceptor_atom)

        # check it
        if cls.check_distance(distance, cutoff=distance_cutoff) is False:
            return (False, (None, None))

        # if the distance passes we want to check the angle, which we
        # will need the coordinates of the adjacent hydrogens to the donor
        h_atoms = [atom for atom in donor_atom.adjacent_atoms
                   if atom.atom_type.element == 'H']

        # check to see if even 1 hydrogen atom satisfies the angle
        # constraint
        okay_angle = None
        h_atoms_iter = iter(h_atoms)
        try:
            while okay_angle is None:
                h_atom = next(h_atoms_iter)

                angle = cls.calc_angle(donor_atom, acceptor_atom, h_atom)

                # check if the angle meets the constraints
                if cls.check_angle(angle, cutoff=angle_cutoff) is False:
                    okay_angle = angle

        # none are found to meet the constraint
        except StopIteration:
            return (False, (distance, None))

        # return in the order of cls.interaction_params
        return (True, (distance, okay_angle))

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
    @staticmethod
    def calc_distance(donor_atom, acceptor_atom):
        return cdist(np.array([donor_atom.coords]), np.array([acceptor_atom.coords]))[0,0]

    @staticmethod
    def calc_angle(donor_atom, acceptor_atom, h_atom):
        v1 = donor_atom.coords - h_atom.coords
        v2 = acceptor_atom.coords - h_atom.coords
        return np.degrees(np.arccos(np.dot(v1, v2)/(la.norm(v1) * la.norm(v2))))

    @staticmethod
    def check_distance(distance, cutoff=distance_cutoff):
        """For a float distance checks if it is less than the configuration
        file HYDROGEN_BOND_DIST_MAX value.

        """
        if distance < cutoff:
            return True
        else:
            return False

    @staticmethod
    def check_angle(angle, cutoff=angle_cutoff):
        """For a float distance checks if it is less than the configuration
        file HYDROGEN_BOND_DON_ANGLE_MIN value.

        """

        if angle > cutoff:
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


class HydrogenBondInx(Interaction):
    """Substantiates HydrogenBondType by selecting donor and acceptor
    features, as well as the involved Hydrogen atom.

    """

    interaction_type = HydrogenBondType
    interaction_params = {key : None for key in interaction_type.interaction_param_keys}

    def __init__(self, donor, acceptor,
                 check=True,
                 interaction_class=None,
                 distance=None,
                 angle=None,
                 distance_cutoff=interaction_type.distance_cutoff,
                 angle_cutoff=interaction_type.angle_cutoff):

        donor_atom = donor.atoms[0]
        acceptor_atom = acceptor.atoms[0]
        if check:
            okay, param_values = self.interaction_type.check(donor_atom, acceptor_atom,
                                                           distance_cutoff=distance_cutoff,
                                                           angle_cutoff=angle_cutoff)

            distance, angle = param_values[0], param_values[1]
            if not okay:
                if angle is None:
                    raise InteractionError(
                        """donor: {0}
                        H: {1}
                        acceptor: {2}
                        distance = {3} FAILED
                        angle = not calculated""".format(donor_atom,
                                                         H,
                                                         acceptor_atom,
                                                         distance))

                else:
                    raise InteractionError(
                        """donor: {0}
                        H: {1}
                        acceptor: {2}
                        distance = {3}
                        angle = {4} FAILED""".format(donor_atom,
                                                     H,
                                                     acceptor_atom,
                                                     distance, angle))
            elif (distance is None) or (angle is None) and (check is False):
                raise ValueError("Must provide distance and angle if check=False is passed")

        # success, finish creating interaction
        atom_system = donor.system
        super().__init__(features=[donor, acceptor],
                         interaction_type=self.interaction_type,
                         system=atom_system)
        self._donor = donor
        self._acceptor = acceptor
        self._distance = distance
        self._angle = angle
        self._interaction_class = interaction_class

    @property
    def donor(self):
        """The donor Feature in the hydrogen bond."""
        return self._donor

    # TODO implement a way to find the H atoms that satisfy the interaction
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
                         'distance', 'angle',]
                         #'H_coords']

        HydrogenBondInxRecord = namedtuple('HydrogenBondInxRecord', record_fields)
        record_attr = {'interaction_class' : self.interaction_class.name}
        record_attr['donor_coords'] = self.donor.atoms[0].coords
        record_attr['acceptor_coords'] = self.acceptor.atoms[0].coords
        record_attr['distance'] = self.distance
        record_attr['angle'] = self.angle
        # TODO
        # record_attr['H_coords'] = self.H.coords

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
