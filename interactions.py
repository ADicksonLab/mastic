""" The interactions module. """
import numpy as np
import numpy.linalg as la

from mast.selection import SelectionsList

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

    def profile_interactions(self, interaction_types, between=None):
        """Profiles all interactions of all features for everything in the
        association. the between flag sets the inter-selection interactions to be profiled.

        E.g.: association.profile_interactions([HydrogenBondType], between=Molecule)

        """
        from itertools import combinations
        from collections import defaultdict

        assert not between is None, "between must be provided"
        assert all([issubclass(itype, InteractionType) for itype in interaction_types]), \
                   "All interaction_types must be a subclass of InteractionType"
        # go through each interaction_type and check for hits
        interactions = {}
        for interaction_type in interaction_types:
            # collect the specific features for each family
            family_features = defaultdict(list)
            for family in interaction_type.feature_families:
                for member in self:
                    try:
                        family_features[family].extend(member.family_selections[family])
                    except KeyError:
                        print("No {0} in {1}".format(family, member))
                        continue

            # pass these to the check class method of the InteractionType
            all_inxs = interaction_type.find_hits(**family_features)
            # separate into intra- and inter-member interactions
            intramember_inxs = []
            intermember_inxs = []
            for inx in all_inxs:
                # get the selections of type `between`
                selections = []
                for member in inx:
                    member_sels = [member for key, member in member.registry]
                    for sel in member_sels:
                        if isinstance(sel, between):
                            selections.append(sel)
                # if they are all in the same selection they are in
                # the same member
                sel_pairs = combinations(selections, 2)
                intra = True
                for pair in sel_pairs:
                    if pair[0] is pair[1]:
                        pass
                    else:
                        intra = False
                if intra:
                    intramember_inxs.append(inx)
                else:
                    intermember_inxs.append(inx)

            interactions[interaction_type] = intermember_inxs

        # TODO handle the intramember_interactions

        # set the interactions for only the intermember interactions
        self._interactions = interactions


class InteractionType(object):
    """ Prototype class for all intermolecular interactions."""
    def __init__(self):
        pass

    @classmethod
    def check(cls, *args, **kwargs):
        raise NotImplementedError

    @classmethod
    def find_hits(cls, **kwargs):
        # make sure all the necessary key word argument families were passed
        for family in cls._feature_families:
            assert family in kwargs.keys(), \
                "{0} feature family must be in keyword arguments".format(
                    family)
        # make sure there are no extra key word argument families
        for key in kwargs:
            assert key in cls._feature_families, \
                "{0} is not a feature needed for finding hits for {1}".format(
                    key, cls)



class HydrogenBondType(InteractionType):
    """ Class for checking validity of a HydrogenBondInx."""

    def __init__(self):
        pass

    _feature_families = mastinxconfig.HBOND_FEATURE_FAMILIES
    feature_families = _feature_families
    _donor_key = 'Donor'
    _acceptor_key = 'Acceptor'
    _feature_types = mastinxconfig.HBOND_FEATURE_TYPES
    feature_type = _feature_types

    def __repr__(self):
        return str(self.__class__)

    @classmethod
    def find_hits(cls, **kwargs):
        """Takes in key-word arguments for the donors and acceptor atom
        IndexedSelections. As an interface find_hits must take in more
        generic selections.

        """
        from itertools import product
        # check that the keys ar okay in parent class
        super().find_hits(**kwargs)

        # Hbond specific stuff
        donors = kwargs[cls._donor_key]
        acceptors = kwargs[cls._acceptor_key]
        donor_Hs = []
        # find Hs for all donors and make Donor-H pairs
        for donor in donors:
            # donors are given in as INdexedSelections by the
            # interface must get the atom
            donor = list(donor.values())[0]
            Hs = [atom for atom in donor.adjacent_atoms if
                        atom.atom_type.element == 'H']
            for H in Hs:
                donor_Hs.append((donor, H))
        # make pairs of Donor-H and acceptors
        hits = []
        # acceptors are given in as IndexedSelections
        # make pairs of them to compare
        acceptors = [list(acceptor.values())[0] for acceptor in acceptors]
        pairs = product(donor_Hs, acceptors)
        for pair in pairs:
            donor_atom = pair[0][0]
            h_atom = pair[0][1]
            acceptor_atom = pair[1]
            # try to make a HydrogenBondInx object, which calls check
            try:
                hbond = HydrogenBondInx(donor=donor_atom, H=h_atom, acceptor=acceptor_atom)
            # else continue to the next pairing
            except InteractionError:
                continue
            # if it succeeds add it to the list of H-Bonds
            hits.append(hbond)

        return hits

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


class NoHHydrogenBondType(InteractionType):
    """Class for checking validity of a NoHHydrogenBondInx. Different
from HydrogenBondType in that having the hydrogens present is not necessary."""

    def __init__(self, hbond_attrs=None):
        super().__init__(attr_dict=hbond_attrs)

    _feature_families = mastinxconfig.HBOND_FEATURE_FAMILIES
    feature_families = _feature_families
    _donor_key = 'Donor'
    _acceptor_key = 'Acceptor'
    _feature_types = mastinxconfig.HBOND_FEATURE_TYPES
    feature_type = _feature_types

    def __repr__(self):
        return str(self.__class__)

    @classmethod
    def find_hits(cls, **kwargs):
        """Takes in key-word arguments for the donors and acceptor atom
IndexedSelections. As an interface find_hits must take in more generic selections."""
        from itertools import product
        # check that the keys ar okay in parent class
        super().find_hits(**kwargs)

        # Hbond specific stuff
        donors = kwargs[cls._donor_key]
        acceptors = kwargs[cls._acceptor_key]
        donor_Hs = []
        # determine the virtual Hydrogens the donor would have
        for donor in donors:
            # donors are given in an IndexedSelections by the
            # interface, must get the atom from this to use
            donor = list(donor.values())[0]
            # get the inferred number of hydrogen bonds
            # TODO dependent on RDKit setting the num_Hs attribute
            num_virtual_Hs = donor.num_Hs
            adj_atoms = donor.adjacent_atoms
            # make a pseudo atom for each virtual hydrogen
            H_coords = virtual_H_coords(donor, adj_atoms, num_virtual_Hs)
            # make pairs of the pseudo-atoms and donors
            for H_coord in H_coords:
                # create a pseudo-atom
                H = PseudoAtom(coords=H_coords)
                donor_Hs.append((donor, H))
        # make pairs of Donor-H and acceptors
        hits = []
        # acceptors are given in as IndexedSelections
        # make pairs of them to compare
        acceptors = [list(acceptor.values())[0] for acceptor in acceptors]
        pairs = product(donor_Hs, acceptors)
        for pair in pairs:
            donor_atom = pair[0][0]
            h_atom = pair[0][1]
            acceptor_atom = pair[1]
            # try to make a HydrogenBondInx object, which calls check
            try:
                hbond = NoHHydrogenBondInx(donor=donor_atom, pseudo_H=h_atom, acceptor=acceptor_atom)
            # else continue to the next pairing
            except InteractionError:
                continue
            # if it succeeds add it to the list of H-Bonds
            hits.append(hbond)

        return hits

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



class Interaction(Association):
    """Base class for associating Selections from a SelectionsList with
information about an about the interaction.

    """

    def __init__(self, members=None, interaction_type=None, system=None):
        super().__init__(members=members, system=system)
        if not interaction_type:
            interaction_type = InteractionType()

        if not issubclass(interaction_type, InteractionType):
            raise TypeError(
                "interaction_type must be a subclass of InteractionType, not {}".format(
                type(interaction_type)))
        self._interaction_type = interaction_type

    @property
    def interaction_type(self):
        return self._interaction_type

class HydrogenBondInx(Interaction):
    """Interaction subclass that has HydrogenBondType type.

    """

    def __init__(self, donor=None, H=None, acceptor=None):

        okay, distance, angle = HydrogenBondType.check(donor, H, acceptor)
        if not okay:
            if angle is None:
                raise InteractionError(
                    """donor: {0}
H: {1}
acceptor: {2}
distance = {3} FAILED
                    angle = not calculated""".format(donor, H, acceptor, distance))

            else:
                raise InteractionError(
                    """donor: {0}
H: {1}
acceptor: {2}
distance = {3}
                    angle = {4} FAILED""".format(donor, H, acceptor, distance, angle))

        # success, finish creating interaction

        atom_system = donor.molecule.system
        super().__init__(members=[donor,H,acceptor],
                         interaction_type=HydrogenBondType,
                         system=atom_system)
        self._interaction_type = HydrogenBondType
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
