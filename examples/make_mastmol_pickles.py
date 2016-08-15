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
except Except:
    print("trypsin with no Hs was not successful")

# BEN no Hs
print("running for BEN with no Hs")
try:
    BEN_rdkit = Chem.MolFromPDBBlock(mastdata.BEN_3ptb, removeHs=False, sanitize=True)
    ben_pkl_path = mastdata.BEN_mastmol_path
    BEN_rdkit_wrapper = RDKitMoleculeWrapper(BEN_rdkit, mol_name="BEN")
    BEN_Molecule = BEN_rdkit_wrapper.make_molecule_type(find_features=True)
    with open(ben_pkl_path, 'wb') as pkl_wf:
        pickle.dump(BEN_Molecule, pkl_wf)

except Except:
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

except Except:
    print("trypsin with Hs was not successful")

# BEN Hs
print("running for BEN with Hs")
try:
    BEN_Hs_rdkit = Chem.MolFromPDBBlock(mastdata.BEN_Hs_3ptb, removeHs=False, sanitize=True)
    ben_Hs_pkl_path = mastdata.BEN_Hs_mastmol_path
    BEN_Hs_rdkit_wrapper = RDKitMoleculeWrapper(BEN_Hs_rdkit, mol_name="BEN")
    BEN_Hs_Molecule = BEN_Hs_rdkit_wrapper.make_molecule_type(find_features=True)
    with open(ben_Hs_pkl_path, 'wb') as pkl_wf:
        pickle.dump(BEN_Hs_Molecule, pkl_wf)

except Except:
    print("BEN with Hs was not successful")
