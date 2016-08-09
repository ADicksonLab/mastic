""" The interactions module. """
import itertools as it
from collections import defaultdict

import numpy as np
import numpy.linalg as la

from mast.selection import SelectionsList
from mast.system import System
import mast.selection as mastsel
import mast.molecule as mastmol
import mast.system as mastsys

import mast.config.interactions as mastinxconfig

__all__ = ['AssociationType', 'Association',
           'Interaction', 'HydrogenBondInx', 'NoHHydrogenBondInx'
           'InteractionType', 'HydrogenBondType', 'NoHHydrogenBondType']

class InteractionError(Exception):
    pass

# basically just a named grouping of multiple selections from a system
# with some methods
class AssociationType(object):
    """Class for defining a relationship between two selections in a
    SystemType allowing for convenient introspection and calculations
    on only the selections. Usually contained within the SystemType
    and substantiated automatically when the system is substantiated
    with the to_system method.

    Examples
    --------
    >>> carbon_attributes = {'element':'C', 'bond_degree':3}
    >>> oxygen_attributes = {'element':'O', 'bond_degree':3}
    >>> COCarbonAtomType = AtomType.factory("COCarbonAtomType", **carbon_attributes)
    >>> COOxygenAtomType = AtomType.factory("COOxygenAtomType", **oxygen_attributes)
    >>> CO_atoms = (COCarbonAtomType, COOxygenAtomType)
    >>> CO_attributes = {"bond_order":3}
    >>> COBondType = BondType.factory("COBondType", atom_types=CO_atoms, **CO_attributes)
    >>> atom_types = [COCarbonAtomType, COOxygenAtomType]
    >>> bond_types = [COBondType]
    >>> bond_map = {0 : (0, 1)}
    >>> CO_attributes = {"name" : "carbon-monoxide", "toxic" : True}
    >>> COMoleculeType = MoleculeType.factory("COType", atom_types=atom_types, bond_types=bond_types, bond_map=bond_map, **CO_attributes)
    >>> system_attrs = {'name' : 'carbon-monoxide-system'}
    >>> member_types = [COMoleculeType, COMoleculeType]
    >>> COSystemType = SystemType.factory("COSystemType", member_types=member_types, **system_attrs)

    We can associate the two CO molecules in the SystemType by
    creating a mapping for which system member has what selection:

    >>> selection_map = {0 : ..., 1 : ...}
    >>> selection_types = [mastsel.Selection, mastsel.Selection]

    So for system members 0 and 1 we will make a mast.Selection of the
    whole member (...) when the AssociationType is substantiated.

    We can also store metadata about an AssociationType if you would
    like but we will keep it simple:

    >>> association_attrs = {'name' : 'carbon-monoxide-carbon-monoxide-association'}

    >>> COCOAssociationType = mastinx.AssociationType.factory("COCOAssociationType", system_type=COSystemType, selection_map=selection_map, selection_types=selection_types, **association_attrs)

    """
    attributes = mastinxconfig.ASSOCIATION_ATTRIBUTES

    def __init__(self):
        pass

    @classmethod
    def to_association(cls, system):
        """Substantiate the association by providing the System to make
        selections on.

        """
        assert system.system_type is cls.system_type, \
            "The system_type of system must be {0}, not {1}".format(
                cls.system_type, system.system_type)

        return Association(system=system, association_type=cls)

    @staticmethod
    def factory(association_type_name,
                system_type=None,
                selection_map=None, selection_types=None,
                **association_attrs):
        """Static method for generating association types dynamically given a type
        name (which will be the class name) and a domain specific dictionary
        of association attributes.

        association_type_name :: name for the generated class
        system_type :: SystemType
        selection_map :: {system_member_idx : member_selection_ids}
        selection_types :: list of classes inheriting from mast.selection.GenericSelection
        association_attrs :: domain specific metadata dictionary

        See mast.config.interactions for standard AssociationType attributes.
        See class docstring for examples.
        """

        # check that system_type is given and correct
        assert system_type, "system_type must be given"
        assert issubclass(system_type, mastsys.SystemType), \
            "system_type must be a subclass of SystemType, not {}}".format(
                system_type)

        # check that there are the same number of selection_map
        # records and selection_types
        assert len(selection_map) == len(selection_types)

        # check that the selection_map is correct
        for i, selmap in enumerate(selection_map.items()):
            member_idx, sel_ids = (selmap[0], selmap[1])
            # validate it indexes a system member
            assert member_idx < len(system_type.member_types), \
                "member index {0} in selection_map out of"\
                " range of {2}, length {1}".format(
                    member_idx, len(system_type.member_types), system_type)
            # validate the sel_ids for the member
            member = system_type.member_types[member_idx]

        # validate the selection_types
        # for selection_type in selection_types:
        #     assert issubclass(selection_type, mastsel.GenericSelection), \
        #         "The selection_type must be a subclass of" \
        #         " mast.selection.GenericSelection, not {}".format(
        #             selection_type)

        # keep track of which attributes the input did not provide
        for attr in AssociationType.attributes:
            try:
                assert attr in association_attrs.keys()
            except AssertionError:
                # LOGGING
                pass
                # print("Attribute {0} not found in association input.".format(attr))

        # add the attributes into the class
        attributes = {attr : None for attr in AssociationType.attributes}
        for attr, value in association_attrs.items():
            # Log the compliance of the attributes
            # if the attribute isn't declared in the attributes log it
            try:
                assert attr in AssociationType.attributes
            except AssertionError:
                # LOGGING
                pass
                # print("Input attribute {0} not in AssociationType attributes.".format(attr))

            # add it to the attributes
            attributes[attr] = value

        association_type = type(association_type_name, (AssociationType,), attributes)
        # add core attributes
        association_type.system_type = system_type
        association_type.selection_map = selection_map
        association_type.selection_types = selection_types

        return association_type


class Association(SelectionsList):
    def __init__(self, system=None, association_type=None):
        # TODO check to make sure that all the atoms are in the same system
        # print(system, id(system))
        # print([bool(member in system) for member in members])
        # print([(member.molecule.system, id(member.molecule.system)) for member in members])
        # if all([(member in system) for member in members]):
        #     super().__init__(association_list=members, association_type=association_type)
        # else:
        #     raise ValueError("Members of a SystemAssociation must all be in the same system")


        # check validity of association_type
        assert issubclass(association_type, AssociationType), \
            "association_type must be a subclass of AssociationType, not {}".format(
                association_type)

        # check validity of the system
        assert isinstance(system, mastsys.System), \
            "system must be of type System, not {}".format(type(system))
        assert system.system_type is association_type.system_type, \
            "the system must be of the system_type in the association_type, not {}".format(
                system.system_type)

        # for selection in selections:
        #     # check member_types to make sure they are in a system
        #     assert 'system' in selection.flags, \
        #         "member_types must be in a system, {0} flags are {1}".format(
        #             selection, selection.flags)
        # # make sure they are in the same system
        # systems = [selection.find_selections(selection_type=[SystemType]) for
        #            selection in selections]
        # assert all([system is systems[0] for system in _systems]), \
        #     "All selections must be of the same system"

        # make selections on the system
        selections = []
        for i, selmap in enumerate(association_type.selection_map.items()):
            member_idx, sel_ids = (selmap[0], selmap[1])
            member = system[member_idx]
            # if the selection_type is None do not make a selection of
            # the member, instead just save the whole member
            if association_type.selection_types[i] is None:
                selection = member
            # otherwise we will make a selection with the type
            elif issubclass(association_type.selection_types[i], mastsel.GenericSelection):
                selection = association_type.selection_types[i](member, sel_ids)
            else:
                raise TypeError("No handler for this type in making an Association selection")

            selections.append(selection)

        # create the SelectionsList
        super().__init__(selection_list=selections)
        self._association_type = association_type
        self._system = system
        self._interactions = None

    @property
    def members(self):
        return self.data

    @property
    def system(self):
        return self._system

    @property
    def system_type(self):
        return self._system.system_type

    @property
    def association_type(self):
        return self._association_type

    @property
    def interactions(self):
        return self._interactions

    def profile_interactions(self, interaction_types,
                             intramember_interactions=False):
        assert all([issubclass(itype, InteractionType) for itype in interaction_types]), \
                   "All interaction_types must be a subclass of InteractionType"

        # if intramember interactions is True make pairs of each
        # member to itself
        if intramember_interactions:
            member_pairs = it.combinations_with_replacement(self.members, 2)
            # the key to each pairing is a tuple of the members indices
            member_idx_pairs = list(it.combinations_with_replacement(
                range(len(self.members)), 2))
        # if intramember_interactions is False only get interactions between members
        else:
            member_pairs = it.combinations(self.members, 2)
            # the key to each pairing is a tuple of the members indices
            member_idx_pairs = list(it.combinations(range(len(self.members)), 2))

        # go through each interaction_type and check for hits
        interactions = {}
        inx_feature_key_pairs = {}
        for interaction_type in interaction_types:

            inx_hits = {}
            member_feature_key_pairs = {}
            # for each pair find the hits in this interaction_type
            for idx, member_pair in enumerate(member_pairs):
                member_a = member_pair[0]
                member_b = member_pair[1]
                feature_key_pairs, pair_hits = interaction_type.find_hits(member_a,
                                                       member_b)
                inx_hits[member_idx_pairs[idx]] = pair_hits
                member_feature_key_pairs[member_idx_pairs[idx]] = feature_key_pairs

            # save the results for this interaction for all member pairs
            interactions[interaction_type] = inx_hits
            inx_feature_key_pairs[interaction_type] = member_feature_key_pairs

        # set the interactions for only the intermember interactions
        return inx_feature_key_pairs, interactions


class InteractionType(object):
    """ Prototype class for all intermolecular interactions."""
    def __init__(self):
        pass

    @classmethod
    def check(cls, *args, **kwargs):
        raise NotImplementedError

    @classmethod
    def find_hits(cls, members_features):
        pass
        # # make sure all the necessary key word argument families were passed
        # for family in cls._feature_families:
        #     assert family in kwargs.keys(), \
        #         "{0} feature family must be in keyword arguments".format(
        #             family)
        # # make sure there are no extra key word argument families
        # for key in kwargs:
        #     assert key in cls._feature_families, \
        #         "{0} is not a feature needed for finding hits for {1}".format(
        #             key, cls)

class HydrogenBondType(InteractionType):
    """ Class for checking validity of a HydrogenBondInx."""

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
                    # get the acceptor atom
                    acceptor_tup = (feature_key, feature)
                    members_features[i]['acceptors'].append(acceptor_tup)

                elif group_attribute == cls._donor_key:
                    # get the donor-H pairs of atoms for this donor
                    donor_atom = feature.atoms[0]
                    donor_H_pairs = [(feature, atom) for atom in
                                     donor_atom.adjacent_atoms if
                                     atom.atom_type.element == 'H']
                    donor_H_pairs_tup = [(feature_key, donor_H_pair) for
                                         donor_H_pair in donor_H_pairs]
                    members_features[i]['donors'].extend(donor_H_pairs_tup)

        # member_a_acceptors = []
        # member_a_donors = []
        # for feature_key, feature in member_a.features.items():
        #     # get groupby attribute to use as a key
        #     group_attribute = feature.feature_type.attributes_data[cls._grouping_attribute]

        #     if group_attribute == cls._acceptor_key:
        #         # get the acceptor atom
        #         acceptor_tup = (feature_key, acceptor)
        #         member_a_acceptors.append(acceptor_tup)

        #     elif group_attribute == cls._donor_key:
        #         # get the donor-H pairs of atoms for this donor
        #         donor_atom = feature.atoms[0]
        #         donor_H_pairs = [(feature, atom) for atom in
        #                          donor_atom.adjacent_atoms if
        #                          atom.atom_type.element == 'H']
        #         donor_H_pairs_tup = [(feature_key, donor_H_pair) for
        #                              donor_H_pair in donor_H_pairs]
        #         member_a_donors.extend(donor_H_pairs_tup)

        # member_b_acceptors = []
        # member_b_donors = []
        # for feature_key, feature in member_b.features.items():
        #     # get groupby attribute to use as a key
        #     group_attribute = feature.feature_type.attributes_data[cls._grouping_attribute]

        #     if group_attribute == cls._acceptor_key:
        #         # get the acceptor atom
        #         acceptor_tup = (feature_key, feature)
        #         member_b_acceptors.append(acceptor_tup)

        #     elif group_attribute == cls._donor_key:
        #         # get the donor-H pairs of atoms for this donor
        #         donor_atom = feature.atoms[0]
        #         donor_H_pairs = [(feature, atom) for atom in
        #                          donor_atom.adjacent_atoms if
        #                          atom.atom_type.element == 'H']
        #         donor_H_pairs_tup = [(feature_key, donor_H_pair) for
        #                              donor_H_pair in donor_H_pairs]
        #         member_b_donors.extend(donor_H_pairs_tup)


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
            donor_feature = donor_tup[1][0]
            h_atom = donor_tup[1][1]
            acceptor_feature_key = acceptor_tup[0]
            acceptor_feature = acceptor_tup[1]
            # try to make a HydrogenBondInx object, which calls check,
            # OPTIMIZATION
            #
            # otherwise we have to call check first then the
            # HydrogenBondInx constructor will re-call check to get
            # the angle and distance. If we allow passing and not
            # checking the angle and distance in the constructor then
            # it would be faster, however I am not going to allow that
            # in this 'safe' InteractionType, an unsafe optimized
            # version can be made separately if desired.
            try:
                hbond = HydrogenBondInx(donor=donor_feature, H=h_atom, acceptor=acceptor_feature)
            # else continue to the next pairing
            except InteractionError:
                continue
            # if it succeeds add it to the list of H-Bonds
            hbonds.append(hbond)
            # and the feature keys to the feature key pairs
            hit_pair_keys.append((donor_feature_key, acceptor_feature_key))

        return hit_pair_keys, hbonds

    @classmethod
    def check(cls, donor_atom, h_atom, acceptor_atom):
        """Checks if the 3 atoms qualify as a hydrogen bond. Returns a tuple
        (bool, float, float) where the first element is whether or not it
        qualified, the second and third are the distance and angle
        respectively. Angle may be None if distance failed to qualify.

        """
        from scipy.spatial.distance import cdist
        distance = cdist(np.array([donor_atom.coords]), np.array([acceptor_atom.coords]))[0,0]
        if cls.check_distance(distance) is False:
            return (False, distance, None)

        v1 = donor_atom.coords - h_atom.coords
        v2 = acceptor_atom.coords - h_atom.coords
        try:
            angle = np.degrees(np.arccos(np.dot(v1,v2)/(la.norm(v1) * la.norm(v2))))
        except RuntimeWarning:
            print("v1: {0} \n"
                  "v2: {1}".format(v1, v2))
        if cls.check_angle(angle) is False:
            return (False, distance, angle)

        return (True, distance, angle)

    @classmethod
    def check_distance(cls, distance):
        if distance < mastinxconfig.HBOND_DIST_MAX:
            return True
        else:
            return False

    @classmethod
    def check_angle(cls, angle):
        if angle > mastinxconfig.HBOND_DON_ANGLE_MIN:
            return True
        else:
            return False

    @classmethod
    def pdb_serial_output(self, inxs, path, delim=","):
        """Output the pdb serial numbers (index in pdb) of the donor and
        acceptor in each HBond to a file:

        donor_1, acceptor_1
        donor_2, acceptor_2"""

        with open(path, 'w') as wf:
            for inx in inxs:
                wf.write("{0}{1}{2}\n".format(inx.donor.atom_type.pdb_serial_number,
                                              delim,
                                              inx.acceptor.atom_type.pdb_serial_number))

class Interaction(SelectionsList):

    def __init__(self, features=None, system=None, interaction_type=None):

        assert interaction_type, "interaction_type must be given"
        assert issubclass(interaction_type, InteractionType), \
            "interaction_type must be a subclass of mast.interactions.InteractionType"

        for feature in features:
            assert feature.system is system, \
                "feature's system must be all the same"

        super().__init__(selection_list=features)
        self._interaction_type = interaction_type

    @property
    def interaction_type(self):
        return self._interaction_type

class HydrogenBondInx(Interaction):
    """Interaction subclass that has HydrogenBondType type.

    """

    def __init__(self, donor=None, H=None, acceptor=None):

        donor_atom = donor.atoms[0]
        acceptor_atom = acceptor.atoms[0]
        okay, distance, angle = HydrogenBondType.check(donor_atom, H, acceptor_atom)
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

    @property
    def donor(self):
        return self._donor

    @property
    def H(self):
        return self._H

    @property
    def acceptor(self):
        return self._acceptor

    @property
    def distance(self):
        return self._distance

    @property
    def angle(self):
        return self._angle

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
                      self_str = str(self),
                      donor_el = self.donor.atom_type.element,
                      donor_coords = tuple(self.donor.coords),
                      donor_mol_name = self.donor.molecule.molecule_type.name,
                      donor_serial = self.donor.atom_type.pdb_serial_number,
                      donor_resname = self.donor.atom_type.pdb_residue_name,
                      donor_resnum = self.donor.atom_type.pdb_residue_number,
                      acceptor_el = self.acceptor.atom_type.element,
                      acceptor_coords = tuple(self.acceptor.coords),
                      acceptor_mol_name = self.acceptor.molecule.molecule_type.name,
                      acceptor_serial = self.acceptor.atom_type.pdb_serial_number,
                      acceptor_resname = self.acceptor.atom_type.pdb_residue_name,
                      acceptor_resnum = self.acceptor.atom_type.pdb_residue_number,
                      distance = self.distance,
                      angle = self.angle,
                 )

        print(string)

    def pickle(self, path):
        import sys
        sys.setrecursionlimit(10000)
        import pickle
        with open(path, 'wb') as wf:
            pickle.dump(self, wf)

if __name__ == "__main__":
    pass
