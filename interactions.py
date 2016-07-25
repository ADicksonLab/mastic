""" The interactions module. """
import numpy as np
import numpy.linalg as la

from mast.selection import SelectionList, SelectionType

__all__ = ['Interaction', 'HydrogenBondInx', 'NoHHydrogenBondInx'
           'InteractionType', 'HydrogenBondType', 'NoHHydrogenBondType']



# families and types from the RDKit feature definition database
# (.fdef)
HBOND_FEATURE_FAMILIES = ['Donor', 'Acceptor']
HBOND_FEATURE_TYPES = []

PISTACK_FEATURE_FAMILIES = ['Aromatic']
PISTACK_FEATURE_TYPES = []

HYDROPH_FEATURE_FAMILIES = ['Hydrophobe', 'LumpedHydrophobe']
HYDROPH_FEATURE_TYPES = []

PICATION_FEATURE_FAMILIES = ['Aromatic', 'PosIonizable']
PICATION_FEATURE_TYPES = []

SALTBRIDGE_FEATURE_FAMILIES = ['PosIonizable', 'NegIonizable']
SALTBRIDGE_FEATURE_TYPES = []

HALOGEN_FEATURE_FAMILIES = []
HALOGEN_FEATURE_TYPES = ['HalogenAcceptor']

WATER_BRIDGE_FEATURE_FAMILIES = []
WATER_BRIDGE_FEATURE_TYPES = []

METAL_FEATURE_FAMILIES = []
METAL_FEATURE_TYPES = []

# Determines allowed deviation from planarity in aromatic rings
AROMATIC_PLANARITY = 5.0 # /AA

# Max. distance between hydrogen bond donor and acceptor (Hubbard & Haider, 2001) + 0.6 A
HBOND_DIST_MAX = 4.1 # /AA
# Min. angle at the hydrogen bond donor (Hubbard & Haider, 2001) + 10
HBOND_DON_ANGLE_MIN = 100 # /AA


# Max. distance for parallel or offset pistacking (McGaughey, 1998)
PISTACK_DIST_MAX = 7.5 # /AA
# Max. Deviation from parallel or perpendicular orientation (in degrees)
PISTACK_ANG_DEV = 30 # /degrees
# Maximum offset of the two rings (corresponds to the radius of benzene + 0.5 A)
PISTACK_OFFSET_MAX = 2.0 # /AA

# Some distance thresholds were extended (max. 1.0A) if too restrictive too account for low-quality structures
# Distance cutoff for detection of hydrophobic contacts
HYDROPH_DIST_MAX = 4.0 # /AA

# Max. distance between charged atom and aromatic ring center (Gallivan and Dougherty, 1999)
PICATION_DIST_MAX = 6.0 # /AA

# Max. distance between centers of charge for salt bridges (Barlow and Thornton, 1983) + 1.5
SALTBRIDGE_DIST_MAX = 5.5 # /AA

# Max. distance between oxy. and halogen (Halogen bonds in biological molecules., Auffinger)+0.5
HALOGEN_DIST_MAX = 4.0 # /AA
# Optimal acceptor angle (Halogen bonds in biological molecules., Auffinger)
HALOGEN_ACC_ANGLE = 120 # /degrees
# Optimal donor angle (Halogen bonds in biological molecules., Auffinger)
HALOGEN_DON_ANGLE = 165 # /degrees
# Max. deviation from optimal angle
HALOGEN_ANGLE_DEV = 30 # /degrees

# Min. distance between water oxygen and polar atom (Jiang et al., 2005) -0.1
WATER_BRIDGE_MINDIST = 2.5 # /AA
# Max. distance between water oxygen and polar atom (Jiang et al., 2005) +0.4
WATER_BRIDGE_MAXDIST = 4.0 # /AA
# Min. angle between acceptor, water oxygen and donor hydrogen (Jiang et al., 2005) - 5
WATER_BRIDGE_OMEGA_MIN = 75 # /degrees
# Max. angle between acceptor, water oxygen and donor hydrogen (Jiang et al., 2005)
WATER_BRIDGE_OMEGA_MAX = 140 # /degrees
# Min. angle between water oxygen, donor hydrogen and donor atom (Jiang et al., 2005)
WATER_BRIDGE_THETA_MIN = 100 # /degrees

# Max. distance between metal ion and interacting atom (Harding, 2001)
METAL_DIST_MAX = 3.0 # /AA


class InteractionError(Exception):
    pass

class Association(SelectionList):
    def __init__(self, association_list=None, association_type=None):
        super().__init__(selection_list=association_list)
        self._association_type = association_type

    @property
    def association_type(self):
        return self._association_type


class AssociationType(SelectionType):
    def __init__(self, assoc_dict=None):
        super().__init__(attr_dict=assoc_dict)

class SystemAssociation(Association):
    def __init__(self, members=None, association_type=None, system=None):
        # TODO check to make sure that all the atoms are in the same system
        # print(system, id(system))
        # print([bool(member in system) for member in members])
        # print([(member.molecule.system, id(member.molecule.system)) for member in members])
        # if all([(member in system) for member in members]):
        #     super().__init__(association_list=members, association_type=association_type)
        # else:
        #     raise ValueError("Members of a SystemAssociation must all be in the same system")

        super().__init__(association_list=members, association_type=association_type)
        self._system = system
        self._interactions = None

    @property
    def system(self):
        return self._system

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


class InteractionType(AssociationType):
    """ Prototype class for all intermolecular interactions."""
    def __init__(self, inx_attrs=None):
        super().__init__(attr_dict=inx_attrs)
        self._feature_families = None
        self._feature_types = None

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

    def __init__(self, hbond_attrs=None):
        super().__init__(attr_dict=hbond_attrs)

    _feature_families = HBOND_FEATURE_FAMILIES
    feature_families = _feature_families
    _donor_key = 'Donor'
    _acceptor_key = 'Acceptor'
    _feature_types = HBOND_FEATURE_TYPES
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
        if distance < HBOND_DIST_MAX:
            return True
        else:
            return False

    @classmethod
    def check_angle(cls, angle):
        if angle > HBOND_DON_ANGLE_MIN:
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

    _feature_families = HBOND_FEATURE_FAMILIES
    feature_families = _feature_families
    _donor_key = 'Donor'
    _acceptor_key = 'Acceptor'
    _feature_types = HBOND_FEATURE_TYPES
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
        if distance < HBOND_DIST_MAX:
            return True
        else:
            return False

    @classmethod
    def check_angle(cls, angle):
        if angle > HBOND_DON_ANGLE_MIN:
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

class PiStackingType(InteractionType):

    def __init__(self, pi_stacking_attr=None):
        super().__init__(attr_dict=pi_stacking_attrs)
    _feature_families = PISTACK_FEATURE_FAMILIES
    feature_families = _feature_families
    _feature1_key = 'Pi1'
    _feature2_key = 'Pi2'
    _feature_types = PISTACK_FEATURE_TYPES
    feature_type = _feature_types

    def __repr__(self):
        return str(self.__class__)

    @classmethod
    def find_hits(cls, **kwargs):
        """Takes in key-word arguments for the features. """

        from itertools import product
        # check that the keys ar okay in parent class
        super().find_hits(**kwargs)

        # pi_stacking specific stuff
        feature1 = kwargs[cls._feature1_key]
        feature2 = kwargs[cls._feature2_key]
        hits = []
        # make pairs of them to compare
        pairs = product(feature1, feature2)
        for pair in pairs:
            features = feature_dict(pair)
            # try to make a PiStacking object, which calls check
            try:
                # if it succeeds add it to the list of H-Bonds
                pi_stacking = PiStacking(**features)

            # else continue to the next pairing
            except InteractionError:
                continue

            hits.append(pi_stacking)

        return hits

    @classmethod
    def check(cls, feature1=None, feature2=None):
        """Check if the input features qualify for this type of
        interaction.

        returns (okay, param1, param2,...)"""
        param1 = calc_param1(feature1, feature2)
        if cls.check_param1(param1) is False:
            return (False, param1, None)

        param2 = calc_param2(feature1, feature2)
        if cls.check_param2(param2) is False:
            return (False, param1, param2)

        return (True, param1, param2)

    @classmethod
    def check_param1(cls, param1):
        if True:
            return True
        else:
            return False

    @classmethod
    def check_param2(cls, param2):
        if True:
            return True
        else:
            return False


class Interaction(SystemAssociation):
    """Base class for associating Selections from a SelectionList with
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
    from rdkit import Chem
    from rdkit.Chem import AllChem
    import os.path as osp
    from copy import copy

    from mast.interactions import HydrogenBondType, InteractionType

    trypsin_dir = osp.expanduser("~/Dropbox/lab/trypsin")
    trypsin_pdb_path = osp.join(trypsin_dir,  "trypsin_Hs.pdb")
    trypsin = Chem.MolFromPDBFile(trypsin_pdb_path, removeHs=False, sanitize=False)
    ben_pdb_path = osp.join(trypsin_dir, "BEN_Hs.pdb")
    ben = Chem.MolFromPDBFile(ben_pdb_path, removeHs=False, sanitize=False)

    from mast.molecule import RDKitMoleculeType

    print("loading RDKit molecules")
    trypsin_type = RDKitMoleculeType(trypsin, mol_name="trypsin")
    ben_type = RDKitMoleculeType(ben, mol_name="BEN")

    print("finding features")
    ben_type.find_features()
    trypsin_type.find_features()

    print("loading into mast.Molecules")
    ben_mol = ben_type.to_molecule(0)
    trypsin_mol = trypsin_type.to_molecule(0)

    print( "making a system")
    from mast.system import System
    trypsys = System([ben_mol, trypsin_mol])

    print("finding system features")
    trypsys.find_features()

    print("setting feature selections")
    trypsys.make_feature_selections()

    print("finding Hbonds in BEN")
    ben_mol.profile_interactions([HydrogenBondType])
    print(ben_mol.internal_interactions)

    print("finding Hbonds in trypsin")
    trypsin_mol.profile_interactions([HydrogenBondType])
    print(trypsin_mol.internal_interactions)

    print("making SystemAssociation for receptor-ligand complex")
    rec_lig_attrs = {}
    rec_lig_attrs['ligand_type'] = ben_type
    rec_lig_attrs['receptor_type'] = trypsin_type
    rec_lig_attrs['name'] = 'trypsin-benzamidine-complex'
    rec_lig_type = AssociationType(rec_lig_attrs)
    rec_lig_assoc = SystemAssociation(members=[trypsys[0],trypsys[1]],
                                                     system=trypsys,
                                      association_type=rec_lig_type)
    rec_lig_assoc = SystemAssociation(members=[trypsys[0],trypsys[1]],
                                                     system=trypsys)


    from mast.molecule import Molecule
    print("testing Hbond interaction between molecules in the receptor ligand association")
    rec_lig_assoc.profile_interactions([HydrogenBondType], between=Molecule)
