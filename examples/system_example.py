import os.path as osp

from mast.system import System, SystemType, AssociationType, Association

from rdkit import Chem
from mast.interfaces.rdkit import RDKitMoleculeWrapper

trypsin_dir = osp.expanduser("~/Dropbox/lab/trypsin")
ben_pdb_path = osp.join(trypsin_dir, "BEN_Hs.pdb")

BEN_rdkit = Chem.MolFromPDBFile(ben_pdb_path, removeHs=False, sanitize=False)

BEN_rdkit_wrapper = RDKitMoleculeWrapper(BEN_rdkit, mol_name="BEN")

print("making molecule type for benzamidine")
BENType = BEN_rdkit_wrapper.make_molecule_type(find_features=True)

BEN_coords = BEN_rdkit_wrapper.get_conformer_coords(0)
member_coords = [BEN_coords]

member_types = [BENType]
system_attrs = {'name' : 'benzamidine-system'}

BenzamidineSystemType = SystemType("BenzamidineSystemType",
                                                  member_types=member_types,
                                                  **system_attrs)

MockAssocType = AssociationType('MockAssocType', system_type=BenzamidineSystemType,
                                selection_map={0 : None},
                                selection_types=[None],
                                name='uuuhhhh')
BenzamidineSystemType.add_association_type(MockAssocType)

bensys = BenzamidineSystemType.to_system(member_coords)

molecule = bensys.molecules[0]

association = bensys.associations[0]
