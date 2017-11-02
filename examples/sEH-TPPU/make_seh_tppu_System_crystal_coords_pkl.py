import os.path as osp
import pickle
import sys
# must set this to pickle substantiated systems sometimes
sys.setrecursionlimit(100000)

# load the system type pickle in
system_pkl_path = osp.join(".", "sEH_TPPU_SystemType.pkl")
with open(system_pkl_path, 'rb') as rf:
    seh_tppu_System = pickle.load(rf)

from rdkit import Chem
import mastic.tests.data as masticdata

TPPU_MOL_path = osp.join(".", "TPPU.mol")
TPPU_MOL_rdkit = Chem.MolFromMolFile(TPPU_MOL_path, sanitize=True)
TPPU_PDB_path = osp.join(".", "TPPU.pdb")
TPPU_PDB_rdkit = Chem.MolFromPDBFile(TPPU_PDB_path, removeHs=False, sanitize=False)
seh_PDB_path = osp.join(".", "sEH.pdb")
seh_rdkit = Chem.MolFromPDBFile(seh_PDB_path, removeHs=False, sanitize=False)

from mastic.interfaces.rdkit import AssignBondOrdersFromTemplate

TPPU_rdkit = AssignBondOrdersFromTemplate(TPPU_MOL_rdkit, TPPU_PDB_rdkit)

from mastic.interfaces.rdkit import RDKitMoleculeWrapper

TPPU_rdkit_wrapper = RDKitMoleculeWrapper(TPPU_rdkit, mol_name="TPPU")
seh_rdkit_wrapper = RDKitMoleculeWrapper(seh_rdkit, mol_name="sEH")

TPPU_coords = TPPU_rdkit_wrapper.get_conformer_coords(0)
seh_coords = seh_rdkit_wrapper.get_conformer_coords(0)
member_coords = [TPPU_coords, seh_coords]

# substantiate the system with the new association
cryst_system = seh_tppu_System.to_system(member_coords)

# pickle it
system_cryst_pkl_path = osp.join(".", "sEH_TPPU_System_cryst.pkl")
with open(system_cryst_pkl_path, 'wb') as wf:
    pickle.dump(cryst_system, wf)
