import pickle
import os.path as osp

import numpy as np

from rdkit import Chem

inputs_path = osp.realpath("../../examples/sEH-TPPU")
TPPU_MOL_path = osp.join(inputs_path, "TPPU.mol")
TPPU_MOL_rdkit = Chem.MolFromMolFile(TPPU_MOL_path, sanitize=True)
TPPU_PDB_path = osp.join(inputs_path, "TPPU.pdb")
TPPU_PDB_rdkit = Chem.MolFromPDBFile(TPPU_PDB_path, removeHs=False, sanitize=False)
seh_PDB_path = osp.join(inputs_path, "sEH.pdb")
seh_rdkit = Chem.MolFromPDBFile(seh_PDB_path, removeHs=False, sanitize=False)

from mastic.interfaces.rdkit import AssignBondOrdersFromTemplate

TPPU_rdkit = AssignBondOrdersFromTemplate(TPPU_MOL_rdkit, TPPU_PDB_rdkit)

from mastic.interfaces.rdkit import RDKitMoleculeWrapper

TPPU_rdkit_wrapper = RDKitMoleculeWrapper(TPPU_rdkit, mol_name="TPPU")
seh_rdkit_wrapper = RDKitMoleculeWrapper(seh_rdkit, mol_name="sEH")

TPPU_coords = TPPU_rdkit_wrapper.get_conformer_coords(0)
seh_coords = seh_rdkit_wrapper.get_conformer_coords(0)
member_coords = [TPPU_coords, seh_coords]

# pickle everything
TPPU_rdkit_pkl_path = osp.join(".", "TPPU_rdkit.pkl")
seh_rdkit_pkl_path = osp.join(".", "sEH_rdkit.pkl")
TPPU_coords_path = osp.join(".", "TPPU.npy")
seh_coords_path = osp.join(".", "sEH.npy")
with open(TPPU_rdkit_pkl_path, 'wb') as wf:
    pickle.dump(TPPU_rdkit_wrapper, wf)
with open(seh_rdkit_pkl_path, 'wb') as wf:
    pickle.dump(seh_rdkit_wrapper, wf)
np.save(TPPU_coords_path, TPPU_coords)
np.save(seh_coords_path, seh_coords)
