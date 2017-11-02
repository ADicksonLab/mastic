import os.path as osp
import pickle
import sys
# must set this to pickle substantiated systems sometimes
sys.setrecursionlimit(100000)

# load the system type pickle in
system_pkl_path = osp.join(".", "Trypsin_Benzamidine_SystemType.pkl")
with open(system_pkl_path, 'rb') as rf:
    Trypsin_Benzamidine_System = pickle.load(rf)

from rdkit import Chem
import mastic.tests.data as masticdata

BEN_MOL_path = osp.join(".", "benzamidine.mol")
BEN_MOL_rdkit = Chem.MolFromMolFile(BEN_MOL_path, sanitize=True)
BEN_PDB_path = osp.join(".", "BEN+Hs_3ptb.pdb")
BEN_PDB_rdkit = Chem.MolFromPDBFile(BEN_PDB_path, removeHs=False, sanitize=False)
trypsin_PDB_path = osp.join(".", "trypsin+Hs_3ptb.pdb")
trypsin_rdkit = Chem.MolFromPDBFile(trypsin_PDB_path, removeHs=False, sanitize=False)

from mastic.interfaces.rdkit import AssignBondOrdersFromTemplate

BEN_rdkit = AssignBondOrdersFromTemplate(BEN_MOL_rdkit, BEN_PDB_rdkit)

from mastic.interfaces.rdkit import RDKitMoleculeWrapper

BEN_rdkit_wrapper = RDKitMoleculeWrapper(BEN_rdkit, mol_name="BEN")
trypsin_rdkit_wrapper = RDKitMoleculeWrapper(trypsin_rdkit, mol_name="Trypsin")

BEN_coords = BEN_rdkit_wrapper.get_conformer_coords(0)
trypsin_coords = trypsin_rdkit_wrapper.get_conformer_coords(0)
member_coords = [BEN_coords, trypsin_coords]

# substantiate the system with the new association
cryst_system = Trypsin_Benzamidine_System.to_system(member_coords)

# pickle it
system_cryst_pkl_path = osp.join(".", "Trypsin_Benzamidine_System_cryst.pkl")
with open(system_cryst_pkl_path, 'wb') as wf:
    pickle.dump(cryst_system, wf)
