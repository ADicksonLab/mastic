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


from mast.molecule import Molecule
print("testing Hbond interaction between molecules in the receptor ligand association")
rec_lig_assoc.profile_interactions([HydrogenBondType], between=Molecule)
