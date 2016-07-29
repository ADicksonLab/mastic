from rdkit import Chem
import os.path as osp

trypsin_dir = osp.expanduser("~/Dropbox/lab/trypsin")
trypsin_pdb_path = osp.join(trypsin_dir, "trypsin_Hs.pdb")
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

from mast.system import System, SystemType
print("making a SystemType")
systype = SystemType({'name' : 'trypsin-benzamidine-complex',
                      'trypsin_type' : trypsin_type,
                      'benzamidine_type' : ben_type,})
print("making a system")
trpsys = System([ben_mol, trypsin_mol], system_type=systype)
print(trpsys.system_type)
print(trpsys.molecule_types)
print(ben_mol._in_system)
print(ben_mol.system)
print(trypsin_mol._in_system)
print(trypsin_mol.system)
mols = trpsys.molecules
print(mols)
mol_sels = trpsys.molecules_sel
print(mol_sels)

from mast.interactions import AssociationType, SystemAssociation
print("making AssociationType")
rec_lig_attrs = {}
rec_lig_attrs['ligand_type'] = ben_type
rec_lig_attrs['receptor_type'] = trypsin_type
rec_lig_attrs['name'] = 'trypsin-benzamidine-complex'
rec_lig_type = AssociationType(rec_lig_attrs)
print(rec_lig_type)
print("making SystemAssociation")
rec_lig_assoc = SystemAssociation(members=[trpsys[0],trpsys[1]],
                                                 system=trpsys,
                                  association_type=rec_lig_type)

print(rec_lig_assoc[0].registry)
print(rec_lig_assoc[1].registry)
