from rdkit import Chem
from rdkit.Chem import AllChem
import os.path as osp

import mast.selection as mastsel
import mast.molecule as mastmol
import mast.system as mastsys
import mast.interactions as mastinx

from mast.interfaces.rdkit import RDKitMoleculeWrapper

import mast.config.interactions as mastinxconfig

from rdkit import Chem

trypsin_dir = osp.expanduser("~/Dropbox/lab/trypsin")
trypsin_pdb_path = osp.join(trypsin_dir, "trypsin_Hs.pdb")
ben_pdb_path = osp.join(trypsin_dir, "BEN_Hs.pdb")

BEN_rdkit = Chem.MolFromPDBFile(ben_pdb_path, removeHs=False, sanitize=False)
trypsin_rdkit = Chem.MolFromPDBFile(trypsin_pdb_path, removeHs=False, sanitize=False)

BEN_rdkit_wrapper = RDKitMoleculeWrapper(BEN_rdkit, mol_name="BEN")
trypsin_rdkit_wrapper = RDKitMoleculeWrapper(trypsin_rdkit, mol_name="Trypsin")

print("making molecule type for trypsin")
TrypsinType = trypsin_rdkit_wrapper.make_molecule_type(find_features=True)
print("making molecule type for benzamidine")
BENType = BEN_rdkit_wrapper.make_molecule_type(find_features=True)

BEN_coords = BEN_rdkit_wrapper.get_conformer_coords(0)
trypsin_coords = trypsin_rdkit_wrapper.get_conformer_coords(0)
member_coords = [BEN_coords, trypsin_coords]

BEN_mol = BENType.to_molecule(BEN_coords)
trypsin_mol = TrypsinType.to_molecule(trypsin_coords)

member_types = [BENType, TrypsinType]
system_attrs = {'name' : 'trypsin-benzamidine-complex'}

TrypsinBenzamidineSystemType = mastsys.SystemType.factory("TrypsinBenzamidineSystemType",
                                                          member_types=member_types,
                                                          **system_attrs)

print("making an AssociationType of the receptor and ligand in the TrypsinBenzamidineSystemType")
rec_lig_attrs = {'name' : 'trypsin-benzamidine-complex'}
# rec_lig_attrs['ligand_type'] = ben_type
# rec_lig_attrs['receptor_type'] = trypsin_type
selection_map = {0 : None, 1 : None}
selection_types = [None, None]
TrypsinBenzamidineAssociationType = \
            mastsys.AssociationType.factory("TrypsinBenzamidineAssociationType",
                                            system_type=TrypsinBenzamidineSystemType,
                                            selection_map=selection_map,
                                            selection_types=selection_types,
                                            **rec_lig_attrs)

# add it to the SystemType
TrypsinBenzamidineSystemType.add_association_type(TrypsinBenzamidineAssociationType)

# now when we make the system the selections are put into an
# Association that can be profiled
trypsys = TrypsinBenzamidineSystemType.to_system(member_coords)

# from mast.molecule import Molecule
print("testing Hbond interaction between molecules in the receptor ligand association")
tryp_ben_assoc = trypsys.associations[0]


intermember_key_pairs, intermember_interactions = \
tryp_ben_assoc.profile_interactions([mastinx.HydrogenBondType])

intramember_key_pairs, intramember_interactions = \
tryp_ben_assoc.profile_interactions([mastinx.HydrogenBondType],
                                    intramember_interactions=True)
