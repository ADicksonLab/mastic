""" The system module. """
import collections as col
import itertools as it
from functools import reduce
import operator as op

from mast.selection import SelectionsList, IndexedSelection
from mast.molecule import Atom, Bond, Molecule, AtomType, BondType, MoleculeType
import mast.selection as mastsel


import mast.config.system as mastsysconfig

__all__ = ['overlaps', 'SystemType', 'System']

def overlaps(members):
    """Check to see if any iterable of substantiated members' coordinates
    overlap.

    """

    from itertools import combinations
    pairs = combinations(members, 2)
    try:
        pair = next(pairs)
    # if it is empty no overlaps
    except StopIteration:
        return False
    flag = True
    while flag:
        overlaps = pair[0].overlaps(pair[1])
        if overlaps:
            return overlaps
        else:
            try:
                pair = next(pairs)
            except StopIteration:
                flag = False
    return False

class SystemType(object):
    """Class for generating specific system type classes with the factory
    method.

    Examples
    --------

    Build some type to put into the system:

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

    Make a SystemType that contains one COMoleculeType:

    >>> system_attrs = {'name' : 'carbon-monoxide-system'}
    >>> COSystemType = SystemType.factory("COSystemType", member_types=[COMoleculeType], **system_attrs)

    """
    attributes = mastsysconfig.SYSTEM_ATTRIBUTES

    def __init__(self, system_type_name, member_types=None,
                **system_attrs):
        """Static method for generating system types dynamically given a type
        name (which will be the class name) and a domain specific dictionary
        of system attributes.

        See mast.config.molecule for standard SystemType attributes.
        See class docstring for examples.
        """

        assert member_types, "molecule_types must be provided"
        for member_type in member_types:
            assert (isinstance(member_type, MoleculeType) or
                    isinstance(member_type, AtomType)), \
                    "molecule_types must contain only MoleculeType or"\
                    " AtomType instancees, not {}".format(
                        type(member_type))

        # keep track of which attributes the input did not provide
        # compared to the config file
        for attr in SystemType.attributes:
            try:
                assert attr in system_attrs.keys()
            except AssertionError:
                pass
                # LOGGING
                # print("Attribute {0} not found in system input.".format(attr))

        # add the attributes into the class
        attributes = {attr : None for attr in SystemType.attributes}
        for attr, value in system_attrs.items():
            try:
                assert attr in SystemType.attributes
            # if it doesn't then log
            except AssertionError:
                # LOGGING
                pass
                # print("Input attribute {0} not in SystemType attributes.".format(attr))
            # add it regardless
            attributes[attr] = value

        self.name = system_type_name
        self.attributes_data = attributes
        self.__dict__.update(attributes)
        self.member_types = member_types
        self.member_type_library = set(member_types)
        self._association_types = []

    def to_system(self, members_coords):
        """Substantiate a System using input coordinates in the order of the
        members in the system.

        """
        # give the members coordinates
        members = []
        for member_idx, member_coords in enumerate(members_coords):
            # create each member using the coordinates
            member_type = self.member_types[member_idx]
            if isinstance(member_type, AtomType):
                members.append(member_type.to_atom(member_coords))
            elif isinstance(member_type, MoleculeType):
                members.append(member_type.to_molecule(member_coords))

        system = System(members, system_type=self)

        return system

    @property
    def molecule_types(self):
        """The MoleculeTypes of all Molecule system members."""
        return [member_type for member_type in self.member_types if
                isinstance(member_type, MoleculeType)]

    @property
    def molecule_type_records(self):
        return [molecule_type.record for molecule_type in self.molecule_types]

    @property
    def molecule_type_df(self):
        import pandas as pd
        return pd.DataFrame(self.molecule_type_records)

    @property
    def atom_types(self):
        """The AtomTypes of all Atom system members."""
        return [member_type for member_type in self.member_types if
                isinstance(member_type, AtomType)]

    @property
    def atom_type_records(self):
        return [atom_type.record for atom_type in self.atom_types]

    @property
    def atom_type_df(self):
        import pandas as pd
        return pd.DataFrame(self.atom_type_records)

    
    def make_member_association_type(self, member_idxs, association_type=None):
        """Match an AssociationType to members of the SystemType"""
        raise NotImplementedError

    @property
    def association_types(self):
        return self._association_types

    
    def add_association_type(self, association_type):
        # check to make sure that it's selection types are in this
        # SystemType
        assert isinstance(association_type, AssociationType), \
            "association_type must be a instance of mast.interactions.Association,"\
            " not {}".format(association_type)

        # check that it is an AssociationType of this SystemType
        assert association_type.system_type is self, \
            "The SystemType of the association_type must be {0}, not {1}".format(
                self, association_type.system_type)

        self._association_types.append(association_type)

    @property
    def association_type_records(self):
        return [assoc_type.record for assoc_type in self.association_types]

    @property
    def association_type_df(self):
        import pandas as pd
        return pd.DataFrame(self.association_type_records)


    @property
    def record(self):
        # define the Record named tuple
        record_fields = ['SystemType'] + list(self.attributes_data.keys())
        SystemTypeRecord = col.namedtuple('SystemTypeRecord', record_fields)
        # build the values for it for this Type
        record_attr = {'SystemType' : self.name}
        record_attr.update(self.attributes_data)
        # make and return
        return SystemTypeRecord(**record_attr)


class System(SelectionsList):
    """System that contains non-overlapping molecules, assumed to be in
the same coordinate system.

    molecules : the molecules in the system

    >>> import numpy as np
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
    >>> COSystemType = SystemType.factory("COSystemType", member_types=[COMoleculeType], **system_attrs)

    Get coordinates for the things in the system
    >>> CO_coords = np.array([[0.0, 0.0, 0.0], [0.0, 0.0, 1.0]])
    >>> member_coords = [CO_coords]

    Substantiate a system from those:
    >>> COSystemType.to_system(member_coords)
    <class 'mast.system.System'>

    """

    def __init__(self, members=None, system_type=None):

        assert issubclass(type(members), col.Sequence), \
            "members must be a subclass of collections.Sequence, not {}".format(
                type(members))

        for member in members:
            assert (issubclass(type(member), Atom) or \
                    issubclass(type(member), Molecule)), \
                    "all elements must be a subclass of type Atom or Molecule, not {}".format(
                        type(member))

        # check to make sure none of the atoms are overlapping
        assert not overlaps(members), \
            "molecule system members cannot be overlapping"

        if system_type:
            assert isinstance(system_type, SystemType), \
                "system_type must be a instance of SystemType, not {}".format(
                    type(system_type))

        super().__init__(selection_list=members, flags=['system'])
        self._system_type = system_type
        # substantiate the Associations in this System
        self._associations = []
        for i, association_type in enumerate(self._system_type.association_types):
            association_class_name = "Association{}".format(i)
            self._associations.append(
                association_type.to_association(self,
                                                association_name=association_class_name))

    @property
    def system_type(self):
        """The SystemType that substantiated this System."""
        return self._system_type

    @property
    def members(self):
        """The system members, of every type."""
        return self.data

    @property
    def atom_types(self):
        """The AtomTypes of every Atom system member."""
        return self.system_type.atom_types

    @property
    def molecule_types(self):
        """The MoleculeTypes of every Molecule system member"""
        return self.system_type.molecule_types

    @property
    def atoms(self):
        """The Atom system members."""
        atoms = [member for  member in self if issubclass(type(member), Atom)]
        return atoms

    @property
    def molecules(self):
        """The Molecule system members."""
        molecules = [member for  member in self if issubclass(type(member), Molecule)]
        return molecules

    @property
    def all_atoms(self):
        """All atoms in the system including those in molecules."""
        all_atoms = self.atoms
        for molecule in self.molecules:
            all_atoms.extend(molecule.atoms)

        return all_atoms

    # YAGNI?
    def molecules_sel(self):
        """Returns a selection on the system of just the Molecule system
        members.

        """
        mol_indices = [i for i, member in enumerate(self) if issubclass(type(member), Molecule)]
        return IndexedSelection(self, mol_indices)

    # YAGNI?
    def atoms_sel(self):
        """Returns a selection on the system of just the Atom system
        members.

        """
        atom_indices = [i for i, member in enumerate(self) if issubclass(type(member), Atom)]
        return IndexedSelection(self, atom_indices)

    @property
    def associations(self):
        return self._associations

    @property
    def associations_types(self):
        return self.system_type.association_types

    def make_feature_selections(self):
        """Make feature selections for all current features in the system's
        molecules.

        """
        for mol in self.molecules:
            mol.make_feature_selections()

    def overlaps(members):
        """Checks whether the members given overlap anything in this system."""
        for member in members:
            for sys_member in self.members:
                if overlaps([member, sys_member]):
                    return True
        # if none overlap
        return False

    @property
    def record(self):
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

    >>> COCOAssociationType = mastsys.AssociationType.factory("COCOAssociationType", system_type=COSystemType, selection_map=selection_map, selection_types=selection_types, **association_attrs)

    """
    attributes = mastsysconfig.ASSOCIATION_ATTRIBUTES

    def __init__(self, association_type_name,
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
        assert isinstance(system_type, SystemType), \
            "system_type must be a instance of SystemType, not {}}".format(
                system_type)

        # check that there are the same number of selection_map
        # records and selection_types
        assert len(selection_map) == len(selection_types)

        # check that the selection_map is correct
        for i, selmap in enumerate(selection_map):
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

        self.name = association_type_name
        self.attributes_data = attributes
        self.__dict__.update(attributes)
        # add core attributes
        self.system_type = system_type
        self.selection_map = selection_map
        self.selection_types = selection_types


    def to_association(self, system, association_name=None):
        """Substantiate the association by providing the System to make
        selections on.

        """
        assert system.system_type is self.system_type, \
            "The system_type of system must be {0}, not {1}".format(
                self.system_type, system.system_type)

        return Association(system=system, association_type=self,
                           association_name=association_name)



    @property
    def record(self):
        # define the Record named tuple
        record_fields = ['AssociationType'] + list(self.attributes_data.keys())
        AssociationTypeRecord = col.namedtuple('AssociationTypeRecord', record_fields)
        # build the values for it for this Type
        record_attr = {'AssociationType' : self.name}
        record_attr.update(self.attributes_data)
        # make and return
        return AssociationTypeRecord(**record_attr)



class Association(SelectionsList):

    def __init__(self, system=None, association_type=None, association_name=None):
        # TODO check to make sure that all the atoms are in the same system
        # print(system, id(system))
        # print([bool(member in system) for member in members])
        # print([(member.molecule.system, id(member.molecule.system)) for member in members])
        # if all([(member in system) for member in members]):
        #     super().__init__(association_list=members, association_type=association_type)
        # else:
        #     raise ValueError("Members of a SystemAssociation must all be in the same system")


        # check validity of association_type
        assert isinstance(association_type, AssociationType), \
            "association_type must be a instance of AssociationType, not {}".format(
                association_type)

        # check validity of the system
        assert isinstance(system, System), \
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
        for i, selmap in enumerate(association_type.selection_map):
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
        self.name = association_name
        self._association_type = association_type
        self._system = system
        self._interactions = None

    @property
    def members(self):
        """The members in the Association."""
        return self.data

    @property
    def system(self):
        """The System the Association selects on."""
        return self._system

    @property
    def system_type(self):
        """The SystemType of the System the Association selects on."""
        return self._system.system_type

    @property
    def association_type(self):
        """The AssociationType this Association substantiates."""
        return self._association_type

    @property
    def interactions(self):
        """Interactions that this Association contains between members. Must
        be set manually due to computation time associated with
        profiling interactions.

        """
        return self._interactions

    @property
    def record(self):
        pass

    def profile_interactions(self, interaction_types,
                             intramember_interactions=False):
        """Accepts any number of InteractionType instancees and identifies
        Interactions between the members of the association using the
        InteractionType.find_hits function.

        Examples
        --------

        """
        from mast.interactions.interactions import InteractionType

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

            member_inx_hits = {}
            member_feature_key_pairs = {}
            # for each pair find the hits in this interaction_type
            for idx, member_pair in enumerate(member_pairs):
                member_a = member_pair[0]
                member_b = member_pair[1]
                feature_key_pairs, pair_hits = interaction_type.find_hits(member_a,
                                                       member_b)
                member_inx_hits[member_idx_pairs[idx]] = pair_hits
                member_feature_key_pairs[member_idx_pairs[idx]] = feature_key_pairs


            # make interaction classes with interaction_type constructor
            for member_pair_idx, item in enumerate(member_feature_key_pairs.items()):
                # member_pair_idx is the pair of members
                # get the member and feature idx pairs (member_a_idx, member_b_idx)
                member_pair_idxs = item[0]
                members = (self.members[member_pair_idxs[0]], self.members[member_pair_idxs[1]])

                # the (member_order, feature_pair) s
                member_feature_pairs = item[1]

                # get the unique feature pairs
                unique_feature_pairs = set(member_feature_pairs)

                # make an instance of interaction_type for all of these
                for inx_class_idx, member_features_tup in enumerate(unique_feature_pairs):
                    member_order = member_features_tup[0]
                    # get the correct member by using the member order
                    # specified from find_hits (either 0 or 1 for the
                    # two members focused on here)
                    member_a = members[member_order[0]]
                    member_b = members[member_order[1]]
                    # feature pair indices for a and b
                    inx_pair_idxs = member_features_tup[1]
                    # inx_class_idx is the index of the unique pair of
                    # features (or the interaction class)
                    # inx_pair_idxs is a tuple of the feature idxs
                    # that make up the interaction class, for each
                    # member respectively

                    # collect the actual objects
                    feat_a_type = member_a.features[inx_pair_idxs[0]].feature_type
                    feat_b_type = member_b.features[inx_pair_idxs[1]].feature_type
                    feat_types = [feat_a_type, feat_b_type]
                    # make a name for the interaction class
                    inx_class_name = "{0}_{1}_{2}InxType".format(
                        interaction_type.interaction_name,
                        self.name,
                        inx_class_idx)
                    # TODO empty for now but could add stuff in the future
                    inx_class_attrs = {}
                    # make the interaction class
                    inx_class = interaction_type(inx_class_name,
                                                 feature_types=feat_types,
                                                 association_type=self,
                                                 assoc_member_pair_idxs=member_pair_idxs,
                                                 **inx_class_attrs)

                    # associate the inx hits for this inx class

                    feature_pairs_idxs = [pair[1] for pair in member_feature_pairs]
                    # get the indices of the hits that are the same
                    # features of this pair
                    class_inxs_idxs = [idx for idx, pair in
                                       enumerate(feature_pairs_idxs)
                                       if pair == inx_pair_idxs]

                    # then get the actual interaction hits (not just
                    # indices) for this member_pair
                    member_inxs = member_inx_hits[member_pair_idxs]
                    # then filter these inxs for those that match the
                    # feature paired indices
                    class_inxs = [inx for i, inx in
                                  enumerate(member_inxs) if i in class_inxs_idxs]
                    # set the interaction_class attribute in each of
                    # these to inx_class
                    for class_inx in class_inxs:
                        class_inx.interaction_class = inx_class

            # flatten the member_inx_hits dict
            all_inx_hits = reduce(op.add, [inxs for inxs in member_inx_hits.values()])
            # add them to the interaction type entry
            interactions[interaction_type] = all_inx_hits
            inx_feature_key_pairs[interaction_type] = member_feature_key_pairs

        # set the interactions for only the intermember interactions
        return inx_feature_key_pairs, interactions

if __name__ == "__main__":
    pass
