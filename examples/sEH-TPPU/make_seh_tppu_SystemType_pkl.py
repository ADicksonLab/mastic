"""The starting point for doing anything with mastic. We first need to
load into mastic datastructures and use a feature identifier (here we
use rdkit and is currently the only implementation for feature
identification) to find the relevant features for doing profiling
on. Then we use these to create templates for molecules and systems
which will then be filled with coordinate data and profiled.

"""
import os.path as osp
import pickle

from rdkit import Chem
from mastic.interfaces.rdkit import AssignBondOrdersFromTemplate
from mastic.interfaces.rdkit import RDKitMoleculeWrapper

from mastic.system import SystemType


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

TPPU_Molecule = TPPU_rdkit_wrapper.make_molecule_type(find_features=True)

seh_Molecule = seh_rdkit_wrapper.make_molecule_type(find_features=True)

seh_pkl_path = osp.join(".", "sEHMoleculeType.pkl")
with open(seh_pkl_path, 'wb') as wf:
    pickle.dump(seh_Molecule, wf)

member_types = [TPPU_Molecule, seh_Molecule]
# just some extra metadata so we know where the features come from in
# the future
system_attrs = {'molecule_source' : 'rdkit'}
sEH_TPPU_SystemType = SystemType("sEH_TPPU_System",
                                     member_types=member_types,
                                     **system_attrs)

system_pkl_path = osp.join(".", "sEH_TPPU_SystemType.pkl")
with open(system_pkl_path, 'wb') as wf:
    pickle.dump(sEH_TPPU_SystemType, wf)
