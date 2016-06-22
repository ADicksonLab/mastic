""" The system module. """
from itertools import product, combinations
import collections as col

from mast.selection import SelectionList, Association, IndexedSelection, \
    SelectionType, SelectionTypeLibrary
from mast.molecule import Atom, Molecule, MoleculeTypeLibrary, MoleculeType

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

class SystemType(SelectionType):
    """Base type for systems, subclasses should implement a to_system
method.

    """
    def __init__(self, system_attrs=None):
        super().__init__(attr_dict=system_attrs)


class System(SelectionList):
    """System that contains non-overlapping molecules, assumed to be in
the same coordinate system.

    molecules : the molecules in the system

    """

    def __init__(self, members=None, system_type=None):

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

        if system_type:
            assert issubclass(system_type, SystemType), \
                "system_type must be a subclass of SystemType, not {}".format(
                    type(system_type))

        super().__init__(selection_list=members)
        for member in member:
            member._in_system = True
        self._system_type = system_type
        self._molecule_types = MoleculeTypeLibrary()
        self._system_associations = None

    @property
    def system_type(self):
        return self._system_type

    @property
    def molecule_types(self):
        return self._molecule_types

    def add_molecule_type(self, mol_type, mol_name=None):
        if not mol_name:
            mol_name = mol_type.name
        self._molecule_types.add_type(mol_type, mol_name)

    @property
    def molecules(self):
        molecules = [member for  member in self if issubclass(type(member), Molecule)]
        return molecules

    # TODO should this be a property or function?
    @property
    def molecules_sel(self):
        mol_indices = [i for i, member in enumerate(self) if issubclass(type(member), Molecule)]
        return IndexedSelection(self, mol_indices)

    @property
    def associations(self):
        return self._system_associations

    def find_features(self):
        """Find features in all members of the system. Currently only molecules."""

        for mol in molecules:
            mol.find_features()

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
            interaction_type()


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

    print("making a SystemType")
    systype = SystemType({'name': 'trypsin-benzamidine-complex'})
    print("making a system")
    tryp_sys = System([ben_mol, trypsin_mol], system_type=systype)

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

