""" The interactions module. """

from itertools import product

from scipy.spatial.distance import cdist
import numpy as np
import numpy.linalg as la

from mast.selection import SelectionList, SelectionType

__all__ = ['Interaction', 'HydrogenBondInx',
           'InteractionType', 'HydrogenBondType']
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

    def profile_interactions(self, interaction_types):
        """Profiles all interactions of all features for everything in the
association.

        """
        assert all([issubclass(itype, InteractionType) for itype in interaction_types]), \
                   "All interaction_types must be a subclass of InteractionType"
        # go through each interaction_type and check for hits
        interactions = {}
        for interaction_type in interaction_types:
            # collect the specific features for each family
            family_features = {}
            for family in interaction_type.feature_families:
                for member in self.members:
                    try:
                        family_features[family] = member.family_selections[family].values()
                    except KeyError:
                        print("No {0} features in {1} for profiling {2}".format(
                            family, self, interaction_type))
                        return None

            # pass these to the check class method of the InteractionType
            all_inxs = interaction_type.find_hits(**family_features)

            # separate into intra- and inter-member interactions
            for inx in all_inxs:
                inx.members

            intermember_interactions[str(interaction_type)]


        self._interactions = intermember_interactions




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


# from the RDKit feature definition database (.fdef)
HBOND_FEATURE_FAMILIES = ['Donor', 'Acceptor']
HBOND_FEATURE_TYPES = []
# Max. distance between hydrogen bond donor and acceptor (Hubbard & Haider, 2001) + 0.6 A
HBOND_DIST_MAX = 4.1
# Min. angle at the hydrogen bond donor (Hubbard & Haider, 2001) + 10
HBOND_DON_ANGLE_MIN = 100
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
        return "HydrogenBond"

    @classmethod
    def find_hits(cls, **kwargs):
        """ Takes in key-word arguments for the donors and acceptor atoms"""
        # check that the keys ar okay in parent class
        super().find_hits(**kwargs)

        # Hbond specific stuff
        donors = kwargs[cls._donor_key]
        acceptors = kwargs[cls._acceptor_key]
        donor_Hs = []
        # find Hs for all donors and make Donor-H pairs
        for donor in donors:
            Hs = [atom for atom in donor.adjacent_atoms if
                        atom.atom_type.element == 'H']
            for H in Hs:
                donor_Hs.append((donor, H))
        # make pairs of Donor-H and acceptors
        hits = []
        pairs = product(donor_Hs, acceptors)
        for pair in pairs:
            donor_atom = pair[0][0]
            h_atom = pair[0][1]
            acceptor_atom = pair[1]
            # try to make a HydrogenBondInx object, which calls check
            try:
                # if it succeeds add it to the list of H-Bonds
                hbond = HydrogenBondInx(donor=donor_atom, H=h_atom, acceptor=acceptor_atom)

            # else continue to the next pairing
            except InteractionError:
                continue

            hits.append(hbond)

        return hits

    @classmethod
    def check(cls, donor_atom, h_atom, acceptor_atom):
        """Checks if the 3 atoms qualify as a hydrogen bond. Returns a tuple
(bool, float, float) where the first element is whether or not it
qualified, the second and third are the distance and angle
respectively. Angle may be None if distance failed to qualify.

        """
        distance = cdist(np.array([donor_atom.coords]), np.array([acceptor_atom.coords]))[0,0]
        if cls.check_distance(distance) is False:
            return (False, distance, None)

        v1 = donor_atom.coords + h_atom.coords
        v2 = h_atom.coords + acceptor_atom.coords
        angle = np.arccos(np.dot(v1,v2)/(la.norm(v1) * la.norm(v2)))
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
        if angle < HBOND_DON_ANGLE_MIN:
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
    print("loading into mast.Molecules")
    ben_mol = ben_type.to_molecule(0)
    trypsin_mol = trypsin_type.to_molecule(0)

    from mast.system import System
    print( "making a system")
    trypsys = System([ben_mol, trypsin_mol])
    print("finding all features")
    trypsys.find_features()

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


    print("testing Hbond interaction between molecules in the receptor ligand association")
    rec_lig_inxs = rec_lig_assoc.profile_interactions([HydrogenBondType])
