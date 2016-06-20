""" The system module. """
from itertools import product, combinations
import collections as col

from mast.selection import SelectionList, Association, IndexedSelection
from mast.molecule import Atom, Molecule

def overlaps(members):
    """Check to see if members overlap.

    """

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

class System(SelectionList):
    """System that contains non-overlapping molecules, assumed to be in
the same coordinate system.

    molecules : the molecules in the system

    """

    def __init__(self, members=None):

        assert issubclass(type(members), col.Sequence), \
            "members must be a subclass of collections.Sequence, not {}".format(
                type(members))

        type_test_func = lambda x: True if (issubclass(type(x), Atom) or
                                            issubclass(type(x), Molecule)) \
                                            else False
        assert all([(type_test_func)(member) for member in members]), \
            "all elements in atoms must be a subclass of type Atom"

        # check to make sure none of the atoms are overlapping
        assert not overlaps(members), \
            "molecule system members cannot be overlapping"

        super().__init__(selection_list=members)

class SystemAssociation(Association):
    def __init__(self, members=None, association_type=None, system=None):

        if all([(mol in system) for mol in members]):
            super().__init__(association_list=members, association_type=association_type)
        else:
            raise ValueError("Members of a SystemAssociation must all be in the same system")

        self._system = system

    @property
    def system(self):
        return self._system

    def profile_interactions(self, interaction_types):
        for interaction_type in interaction_types:
            pass


if __name__ == "__main__":

    from rdkit import Chem
    import os.path as osp
    trypsin_dir = osp.expanduser("~/Dropbox/lab/trypsin")
    trypsin_pdb_path = osp.join(trypsin_dir, "trypsin.pdb")
    trypsin = Chem.MolFromPDBFile(trypsin_pdb_path, removeHs=False)
    ben_pdb_path = osp.join(trypsin_dir, "BEN.pdb")
    ben = Chem.MolFromPDBFile(ben_pdb_path, removeHs=False)

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


    ligand = IndexedSelection(rec_lig_assoc, [0])
    receptor = IndexedSelection(rec_lig_assoc, [1])

    print(rec_lig_assoc[0].registry)
    print(rec_lig_assoc[1].registry)
    print(rec_lig_assoc[0].get_selections())

    print("finding BEN features")
    ben_mol.find_features()
    print(ben_mol.feature_families)
    print(ben_mol.feature_types)
    print(ben_mol.family_selections)
    print(ben_mol.type_selections)
    # print(ben_mol.feature_dataframe)
    print("finding trypsin features")
    trypsin_mol.find_features()
