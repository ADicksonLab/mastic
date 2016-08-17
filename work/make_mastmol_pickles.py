import pickle

from rdkit import Chem

import mast.molecule as mastmol
import mast.system as mastsys
import mast.interactions as mastinx
from mast.interactions.hydrogen_bond import HydrogenBondType

from mast.interfaces.rdkit import RDKitMoleculeWrapper

import mast.config.interactions as mastinxconfig

import mast.tests.data as mastdata

from rdkit import Chem

# trypsin no Hs
print("running for trypsin with no Hs")
try:
    trypsin_rdkit = Chem.MolFromPDBBlock(mastdata.trypsin_3ptb, removeHs=False, sanitize=True)
    trypsin_pkl_path = mastdata.trypsin_mastmol_path
    trypsin_rdkit_wrapper = RDKitMoleculeWrapper(trypsin_rdkit, mol_name="Trypsin")
    Trypsin_Molecule = trypsin_rdkit_wrapper.make_molecule_type(find_features=True)
    with open(trypsin_pkl_path, 'wb') as pkl_wf:
        pickle.dump(Trypsin_Molecule, pkl_wf)
except Exception:
    print("trypsin with no Hs was not successful")

# BEN no Hs from .mol file
print("running for BEN with no Hs")
try:
    BEN_rdkit = Chem.MolFromMolBlock(mastdata.benzamidine, removeHs=False, sanitize=True)
    ben_pkl_path = mastdata.BEN_mastmol_path
    BEN_rdkit_wrapper = RDKitMoleculeWrapper(BEN_rdkit, mol_name="BEN")
    BEN_Molecule = BEN_rdkit_wrapper.make_molecule_type(find_features=True)
    with open(ben_pkl_path, 'wb') as pkl_wf:
        pickle.dump(BEN_Molecule, pkl_wf)

except Exception:
    print("BEN with no Hs was not successful")

# trypsin Hs
print("running for trypsin with Hs")
try:
    trypsin_Hs_rdkit = Chem.MolFromPDBBlock(mastdata.trypsin_Hs_3ptb, removeHs=False, sanitize=True)
    trypsin_Hs_pkl_path = mastdata.trypsin_Hs_mastmol_path
    trypsin_Hs_rdkit_wrapper = RDKitMoleculeWrapper(trypsin_Hs_rdkit, mol_name="Trypsin")
    Trypsin_Hs_Molecule = trypsin_Hs_rdkit_wrapper.make_molecule_type(find_features=True)
    with open(trypsin_Hs_pkl_path, 'wb') as pkl_wf:
        pickle.dump(Trypsin_Hs_Molecule, pkl_wf)

except Exception:
    print("trypsin with Hs was not successful")

# BEN Hs from pdb
print("running for BEN with Hs")
try:
    BEN_Hs_rdkit = Chem.MolFromPDBBlock(mastdata.BEN_Hs_3ptb, removeHs=False, sanitize=True)
    BEN_Hs_pkl_path = mastdata.BEN_Hs_mastmol_path
    BEN_Hs_rdkit_wrapper = RDKitMoleculeWrapper(BEN_Hs_rdkit, mol_name="BEN")
    BEN_Hs_Molecule = BEN_Hs_rdkit_wrapper.make_molecule_type(find_features=True)
    with open(BEN_Hs_pkl_path, 'wb') as pkl_wf:
        pickle.dump(BEN_Hs_Molecule, pkl_wf)

except Exception:
    print("BEN with Hs was not successful")

# TPPU
# print("running for TPPU with Hs")
# try:
#     TPPU_Hs_rdkit = Chem.MolFromPDBBlock(mastdata.TPPU_Hs)
#     TPPU_Hs_pkl_path = mastdata.TPPU_Hs_mastmol_path
#     TPPU_Hs_rdkit_wrapper = RDKitMoleculeWrapper(TPPU_Hs_rdkit, mol_name="TPPU")
#     TPPU_Hs_Molecule = TPPU_Hs_rdkit_wrapper.make_molecule_type(find_features=True)
#     with open(TPPU_Hs_pkl_path, 'wb') as pkl_wf:
#         pickle.dump(TPPU_Hs_Molecule, pkl_wf)
# except Exception:
#     print("TPPU with Hs not succesful")

# # sEH
# print("running for sEH with Hs")
# try:
#     sEH_Hs_rdkit = Chem.MolFromPDBBlock(mastdata.sEH_Hs)
#     sEH_Hs_pkl_path = mastdata.sEH_Hs_mastmol_path
#     sEH_Hs_rdkit_wrapper = RDKitMoleculeWrapper(sEH_Hs_rdkit, mol_name="sEH")
#     sEH_Hs_Molecule = sEH_Hs_rdkit_wrapper.make_molecule_type(find_features=True)
#     with open(sEH_Hs_pkl_path, 'wb') as pkl_wf:
#         pickle.dump(sEH_Hs_Molecule, pkl_wf)
# except Exception:
#     print("sEH with Hs not succesful")

# sulfasalazine

# difluoromethylornithine
