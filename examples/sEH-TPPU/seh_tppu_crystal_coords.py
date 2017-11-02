import os.path as osp
import numpy as np

from rdkit import Chem
import mastic.tests.data as masticdata
from mastic.interfaces.rdkit import AssignBondOrdersFromTemplate
from mastic.interfaces.rdkit import RDKitMoleculeWrapper

TPPU_MOL_path = osp.join(".", "TPPU.mol")
TPPU_MOL_rdkit = Chem.MolFromMolFile(TPPU_MOL_path, sanitize=True)
TPPU_PDB_path = osp.join(".", "TPPU.pdb")
TPPU_PDB_rdkit = Chem.MolFromPDBFile(TPPU_PDB_path, removeHs=False, sanitize=False)
seh_PDB_path = osp.join(".", "sEH.pdb")
seh_rdkit = Chem.MolFromPDBFile(seh_PDB_path, removeHs=False, sanitize=False)

TPPU_rdkit = AssignBondOrdersFromTemplate(TPPU_MOL_rdkit, TPPU_PDB_rdkit)

TPPU_rdkit_wrapper = RDKitMoleculeWrapper(TPPU_rdkit, mol_name="TPPU")
seh_rdkit_wrapper = RDKitMoleculeWrapper(seh_rdkit, mol_name="sEH")

TPPU_coords = TPPU_rdkit_wrapper.get_conformer_coords(0)
seh_coords = seh_rdkit_wrapper.get_conformer_coords(0)

# write the coordinates out to a binary file
np.save("TPPU_coords.npy", TPPU_coords)
np.save("sEH_coords.npy", seh_coords)
