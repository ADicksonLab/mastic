import pickle
import os.path as osp

import mast.molecule as mastmol
import mast.system as mastsys
import mast.interactions as mastinx
from mast.interactions.hydrogen_bond import HydrogenBondType

from mast.interfaces.rdkit import RDKitMoleculeWrapper

import mast.config.interactions as mastinxconfig

import mast.tests.data as mastdata

from rdkit import Chem
from rdkit.Chem import Draw
from rdkit.Chem import AllChem

mastdir = "/home/salotz/Dropbox/devel/mast"
masttmpdir = osp.join(mastdir, "tmp")

# BEN no Hs pdb
BEN_pdb_path = osp.join(mastdir, "tests/data/3ptb_BEN.pdb")
BEN_pdb_rdkit = Chem.MolFromPDBFile(BEN_pdb_path, removeHs=False, sanitize=True)
ben_pkl_path = mastdata.BEN_mastmol_path
BEN_rdkit_wrapper = RDKitMoleculeWrapper(BEN_pdb_rdkit, mol_name="BEN")
BEN_Molecule = BEN_rdkit_wrapper.make_molecule_type(find_features=True)
AllChem.Compute2DCoords(BEN_pdb_rdkit)

Chem.Draw.MolToFile(BEN_pdb_rdkit, osp.join(masttmpdir, "BEN_pdb.png"))

# BEN Hs pdb
BEN_Hs_pdb_path = osp.join(mastdir, "tests/data/BEN+Hs_3ptb.pdb")
BEN_Hs_pdb_rdkit = Chem.MolFromPDBFile(BEN_Hs_pdb_path, removeHs=False, sanitize=True)
BEN_Hs_pkl_path = mastdata.BEN_Hs_mastmol_path
BEN_Hs_rdkit_wrapper = RDKitMoleculeWrapper(BEN_Hs_pdb_rdkit, mol_name="BEN")
BEN_Hs_Molecule = BEN_Hs_rdkit_wrapper.make_molecule_type(find_features=True)
AllChem.Compute2DCoords(BEN_Hs_pdb_rdkit)

Chem.Draw.MolToFile(BEN_Hs_pdb_rdkit, osp.join(masttmpdir, "BEN_Hs_pdb.png"))

# BEN Hs mol
BEN_mol_path = "/home/salotz/Dropbox/devel/mast/tests/data/benzamidine.mol"
BEN_mol_rdkit = Chem.MolFromMolFile(BEN_mol_path, removeHs=False, sanitize=True)
BEN_mol_rdkit_wrapper = RDKitMoleculeWrapper(BEN_mol_rdkit, mol_name="BEN")
BEN_mol_Molecule = BEN_mol_rdkit_wrapper.make_molecule_type(find_features=True)
AllChem.Compute2DCoords(BEN_mol_rdkit)

Chem.Draw.MolToFile(BEN_mol_rdkit, osp.join(masttmpdir, "BEN_mol.png"))


print("making TPPU")
TPPU_pdb_path = "/home/salotz/Dropbox/devel/mast/tests/data/TPPU_Hs.pdb"
TPPU_Hs_pdb_rdkit = Chem.MolFromPDBFile(TPPU_pdb_path, removeHs=False, sanitize=True)
TPPU_Hs_pdb_rdkit_wrapper = RDKitMoleculeWrapper(TPPU_Hs_pdb_rdkit, mol_name="TPPU")
TPPU_Hs_pdb_Molecule = TPPU_Hs_pdb_rdkit_wrapper.make_molecule_type(find_features=True)

TPPU_mol_path = "/home/salotz/Dropbox/devel/mast/tests/data/TPPU.mol"
TPPU_mol_rdkit = Chem.MolFromMolFile(TPPU_mol_path, removeHs=False, sanitize=True)
TPPU_mol_rdkit_wrapper = RDKitMoleculeWrapper(TPPU_mol_rdkit, mol_name="TPPU")
TPPU_mol_Molecule = TPPU_mol_rdkit_wrapper.make_molecule_type(find_features=True)


AllChem.Compute2DCoords(TPPU_Hs_pdb_rdkit)
AllChem.Compute2DCoords(TPPU_mol_rdkit)

Chem.Draw.MolToFile(TPPU_Hs_pdb_rdkit, osp.join(masttmpdir, "TPPU_pdb.png"))
Chem.Draw.MolToFile(TPPU_mol_rdkit, osp.join(masttmpdir, "TPPU_mol.png"))

# TPPU_Hs_pkl_path = "TPPU+Hs+features.pkl"
# with open(TPPU_Hs_pkl_path, 'wb') as pkl_wf:
#     pickle.dump(TPPU_Hs_Molecule, pkl_wf)

# print("making sEH")
# sEH_path = "/home/salotz/Dropbox/devel/mast/tests/data/sEH_Hs.pdb"
# sEH_Hs_rdkit = Chem.MolFromPDBFile(sEH_path)
# sEH_Hs_rdkit_wrapper = RDKitMoleculeWrapper(sEH_Hs_rdkit, mol_name="sEH")
# sEH_Hs_Molecule = sEH_Hs_rdkit_wrapper.make_molecule_type(find_features=True)

# sEH_Hs_pkl_path = "sEH+Hs+features.pkl"
# with open(sEH_Hs_pkl_path, 'wb') as pkl_wf:
#     pickle.dump(sEH_Hs_Molecule, pkl_wf)
