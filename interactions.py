""" The interactions module. """

from scipy.spatial.distance import cdist
import numpy as np
import numpy.linalg as la

from mast.selection import AssociationType
from mast.system import SystemAssociation, System

__all__ = ['Interaction', 'HydrogenBondInx',
           'InteractionType', 'HydrogenBondType']

class InteractionType(AssociationType):
    """ Prototype class for all intermolecular interactions."""
    def __init__(self, description=None):
        super().__init__(description=description)
        self._feature_families = None
        self._feature_types = None

# Max. distance between hydrogen bond donor and acceptor (Hubbard & Haider, 2001) + 0.6 A
HBOND_DIST_MAX = 4.1
# Min. angle at the hydrogen bond donor (Hubbard & Haider, 2001) + 10
HBOND_DON_ANGLE_MIN = 100
class HydrogenBondType(InteractionType):
    """ Class for checking validity of a HydrogenBondInx."""

    @classmethod
    def check(cls, donor_atom, h_atom, acceptor_atom):
        distance = cdist(np.array([donor_atom.coords]), np.array([acceptor_atom.coords]))[0,0]
        print("distance:", distance)
        if not cls.check_distance(distance):
            return False

        v1 = donor_atom.coords + h_atom.coords
        v2 = h_atom.coords + acceptor_atom.coords
        angle = np.arccos(np.dot(v1,v2)/(la.norm(v1) * la.norm(v2)))
        print("angle:", angle)
        if not cls.check_angle(angle):
            return False

        return True

    @classmethod
    def check_distance(cls, distance):
        if distance < HBOND_DIST_MAX:
            print("distance HIT")
            return True
        else:
            return False

    @classmethod
    def check_angle(cls, angle):
        if angle < HBOND_DON_ANGLE_MIN:
            print("angle HIT")
            return True
        else:
            return False

class Interaction(SystemAssociation):
    """Base class for associating Selections from a SelectionList with
information about an about the interaction.

    """

    def __init__(self, members=None, interaction_type=None):
        super().__init__(members=members)
        if not interaction_type:
            interaction_type = InteractionType()

        if not issubclass(type(interaction_type), InteractionType):
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

    def __init__(self, donor=None, acceptor=None):
        description = "Hydrogen-Bond Interaction"

        # TODO get the distances
        distance=None
        angle=None

        super().__init__()
        # type checking of these using HydrogenBondType
        self._interaction_type = HydrogenBondType

        try:
            if HydrogenBondType.donor(donor):
                self._donor = donor
            else:
                raise ValueError("Donor, {}, is not valid".format(donor))
            if HydrogenBondType.acceptor(acceptor):
                self._acceptor = acceptor
            else:
                raise ValueError("Acceptor, {}, is not valid".format(acceptor))
            if HydrogenBondType.distance(distance):
                self._distance = distance
            else:
                raise ValueError("distance, {}, is not valid".format(distance))
            if HydrogenBondType.angle(angle):
                self._angle = angle
            else:
                raise ValueError("Angle, {}, is not valid".format(angle))
        except ValueError:
            print("not a HydrogenBondType Interaction")

    @property
    def donor(self):
        return self._donor

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

    trypsin_dir = osp.expanduser("~/Dropbox/lab/trypsin")
    trypsin_pdb_path = osp.join(trypsin_dir, "trypsin.pdb")
    trypsin = Chem.MolFromPDBFile(trypsin_pdb_path, removeHs=False)
    # trypsin = Chem.AddHs(trypsin)
    ben_pdb_path = osp.join(trypsin_dir, "BEN.pdb")
    ben = Chem.MolFromPDBFile(ben_pdb_path, removeHs=False)
    ben = Chem.AddHs(ben)
    ben_Hs_embed = copy(ben)
    # AllChem.EmbedMolecule(ben)
    # AllChem.UFFOptimizeMolecule(ben_Hs)

    from mast.molecule import RDKitMoleculeType

    print("loading RDKit molecules")
    trypsin_type = RDKitMoleculeType(trypsin, mol_type="trypsin")
    ben_type = RDKitMoleculeType(ben, mol_type="BEN")
    print("loading into mast.Molecules")
    ben_mol = ben_type.to_molecule(0)
    trypsin_mol = trypsin_type.to_molecule(0)

    print("making a system")
    tryp_sys = System([ben_mol, trypsin_mol])

    print("making SystemAssociation")
    rec_lig_assoc = SystemAssociation(members=[tryp_sys[0],tryp_sys[1]],
                                                     system=tryp_sys)

    print("finding BEN features")
    ben_mol.find_features()
    print("finding trypsin features")
    trypsin_mol.find_features()

    print("testing an Hbond interaction")
    donors = ben_mol.family_selections['Donor']
    donor_keys = list(donors.keys())
    acceptors = trypsin_mol.family_selections['Acceptor']
    acceptor_keys = list(acceptors.keys())

    from itertools import product
    # test all pairs for if they are hydrogen bonded
    h_bonds = []
    for pair in product(donors.values(), acceptors.values()):
        donor = pair[0]
        acceptor = pair[1]
        donor_Hs = [atom for atom in donor.adjacent_atoms if atom.atom_type.element == 'H']
        distance = cdist(np.array([donor.coords]), np.array([acceptor.coords]))[0,0]
        print(distance)
        if HydrogenBondType.check_distance(distance):
            print("HIT")

