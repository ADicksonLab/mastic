import os.path as osp

trypsin_dir = osp.expanduser("~/Dropbox/lab/trypsin")
trypsin_pdb_path = osp.join(trypsin_dir, "trypsin_Hs.pdb")
ben_pdb_path = osp.join(trypsin_dir, "BEN_Hs.pdb")

from rdkit import Chem

BEN_rdkit = Chem.MolFromPDBFile(ben_pdb_path, removeHs=False, sanitize=False)
trypsin_rdkit = Chem.MolFromPDBFile(trypsin_pdb_path, removeHs=False, sanitize=False)

from mast.interfaces.rdkit import RDKitMoleculeWrapper
BEN_rdkit_wrapper = RDKitMoleculeWrapper(BEN_rdkit, mol_name="BEN")
trypsin_rdkit_wrapper = RDKitMoleculeWrapper(trypsin_rdkit, mol_name="Trypsin")

print("making molecule type for trypsin")
TrypsinType = trypsin_rdkit_wrapper.make_molecule_type()
print("making molecule type for benzamidine")
BENType = BEN_rdkit_wrapper.make_molecule_type()

print("Finding features for trypsin")
TrypsinType.features = trypsin_rdkit_wrapper.find_features()
print("finding features for benzamidine")
BENType.features = BEN_rdkit_wrapper.find_features()


BEN_coords = BEN_rdkit_wrapper.get_conformer_coords(0)
trypsin_coords = trypsin_rdkit_wrapper.get_conformer_coords(0)
member_coords = [BEN_coords, trypsin_coords]

BEN_mol = BENType.to_molecule(BEN_coords)
trypsin_mol = TrypsinType.to_molecule(trypsin_coords)

member_types = [BENType, TrypsinType]
system_attrs = {'name' : 'trypsin-benzamidine-complex'}
from mast.system import System, SystemType
import ipdb; ipdb.set_trace()
TrypsinBenzamidineSystemType = SystemType.factory("TrypsinBenzamidineSystemType",
                                                  member_types=member_types,
                                                  **system_attrs)


trypsys = TrypsinBenzamidineSystemType.to_system(member_coords)
print(trypsys.system_type)
print(trypsys.molecule_types)
print(BEN_mol.isin_system)
print(trypsin_mol.isin_system)
mols = trypsys.molecules
print(mols)
mol_sels = trypsys.molecules_sel()
print(mol_sels)

from mast.interactions import AssociationType, SystemAssociation
print("making AssociationType")
rec_lig_attrs = {}
rec_lig_attrs['ligand_type'] = ben_type
rec_lig_attrs['receptor_type'] = trypsin_type
rec_lig_attrs['name'] = 'trypsin-benzamidine-complex'
rec_lig_type = AssociationType(**rec_lig_attrs)
print(rec_lig_type)
print("making SystemAssociation")
rec_lig_assoc = SystemAssociation(members=[trypsys[0],trypsys[1]],
                                                 system=trypsys,
                                  association_type=rec_lig_type)

print(rec_lig_assoc[0].registry)
print(rec_lig_assoc[1].registry)
