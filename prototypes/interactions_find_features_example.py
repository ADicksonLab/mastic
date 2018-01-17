import pandas as pd

from rdkit import Chem
from rdkit.Chem import AllChem
import os.path as osp

import mastic.selection as masticsel
import mastic.molecule as masticmol
import mastic.system as masticsys
import mastic.interactions as masticinx
from mastic.interactions.hydrogen_bond import HydrogenBondType

from mastic.interfaces.rdkit import RDKitMoleculeWrapper

import mastic.config.interactions as masticinxconfig

import mastic.tests.data as masticdata

# without Hs straight from pdb
BEN_rdkit = Chem.MolFromPDBBlock(masticdata.BEN_3ptb, removeHs=False, sanitize=True)
trypsin_rdkit = Chem.MolFromPDBBlock(masticdata.trypsin_3ptb, removeHs=False, sanitize=True)

# with hydrogens added files
BEN_Hs_rdkit = Chem.MolFromPDBBlock(masticdata.BEN_Hs_3ptb, removeHs=False, sanitize=True)
trypsin_Hs_rdkit = Chem.MolFromPDBBlock(masticdata.trypsin_Hs_3ptb, removeHs=False, sanitize=True)

# the MoleculeType pickle paths
trypsin_pkl_path = masticdata.trypsin_masticmol_path
ben_pkl_path = masticdata.BEN_masticmol_path

trypsin_Hs_pkl_path = masticdata.trypsin_Hs_masticmol_path
ben_Hs_pkl_path = masticdata.BEN_Hs_masticmol_path



# make wrappers
BEN_rdkit_wrapper = RDKitMoleculeWrapper(BEN_rdkit, mol_name="BEN")
trypsin_rdkit_wrapper = RDKitMoleculeWrapper(trypsin_rdkit, mol_name="Trypsin")

BEN_Hs_rdkit_wrapper = RDKitMoleculeWrapper(BEN_Hs_rdkit, mol_name="BEN")
trypsin_Hs_rdkit_wrapper = RDKitMoleculeWrapper(trypsin_Hs_rdkit, mol_name="Trypsin")

print("making molecule type for trypsin")
Trypsin_Molecule = trypsin_rdkit_wrapper.make_molecule_type(find_features=True)

Trypsin_Hs_Molecule = trypsin_Hs_rdkit_wrapper.make_molecule_type(find_features=True)
print("making molecule type for benzamidine")
BEN_Molecule = BEN_rdkit_wrapper.make_molecule_type(find_features=True)
BEN_Hs_Molecule = BEN_Hs_rdkit_wrapper.make_molecule_type(find_features=True)

print("pickling the molecules for later use")
import pickle

with open(trypsin_pkl_path, 'wb') as pkl_wf:
    pickle.dump(Trypsin_Molecule, pkl_wf)
with open(ben_pkl_path, 'wb') as pkl_wf:
    pickle.dump(BEN_Molecule, pkl_wf)

with open(ben_Hs_pkl_path, 'wb') as pkl_wf:
    pickle.dump(BEN_Hs_Molecule, pkl_wf)
with open(trypsin_Hs_pkl_path, 'wb') as pkl_wf:
    pickle.dump(Trypsin_Hs_Molecule, pkl_wf)

# making substantiated molecules
BEN_coords = BEN_rdkit_wrapper.get_conformer_coords(0)
trypsin_coords = trypsin_rdkit_wrapper.get_conformer_coords(0)
member_coords = [BEN_coords, trypsin_coords]

BEN_mol = BEN_Molecule.to_molecule(BEN_coords)
trypsin_mol = Trypsin_Molecule.to_molecule(trypsin_coords)

member_types = [BEN_Molecule, Trypsin_Molecule]
system_attrs = {'name' : 'trypsin-benzamidine-complex'}

Trypsin_Benzamidine_System = masticsys.SystemType("Trypsin_Benzamidine_System",
                                                          member_types=member_types,
                                                          **system_attrs)

print("making an AssociationType of the receptor and ligand in the Trypsin_Benzamidine_System")
rec_lig_attrs = {'name' : 'trypsin-benzamidine-complex'}
# rec_lig_attrs['ligand_type'] = ben_type
# rec_lig_attrs['receptor_type'] = trypsin_type
selection_map = {0 : None, 1 : None}
selection_types = [None, None]
Trypsin_Benzamidine_Association = \
            masticsys.AssociationType("Trypsin_Benzamidine_Association",
                                            system_type=Trypsin_Benzamidine_System,
                                            selection_map=selection_map,
                                            selection_types=selection_types,
                                            **rec_lig_attrs)

# add it to the SystemType
Trypsin_Benzamidine_System.add_association_type(Trypsin_Benzamidine_Association)


# now when we make the system the selections are put into an
# Association that can be profiled
trypsys = Trypsin_Benzamidine_System.to_system(member_coords)

# from mastic.molecule import Molecules
print("testing Hbond interaction between molecules in the receptor ligand association")
tryp_ben_assoc = trypsys.associations[0]

intermember_key_pairs, intermember_interactions = \
tryp_ben_assoc.profile_interactions([HydrogenBondType])

intramember_key_pairs, intramember_interactions = \
tryp_ben_assoc.profile_interactions([HydrogenBondType],
                                    intramember_interactions=True)

intermember_inx_class_df = pd.DataFrame([inx.interaction_class.record() for
                                 inx in intermember_interactions[HydrogenBondType]])
intramember_inx_class_df = pd.DataFrame([inx.interaction_class.record() for
                                 inx in intramember_interactions[HydrogenBondType]])

intermember_inx_df = pd.DataFrame([inx.record for inx in
                           intermember_interactions[HydrogenBondType]])
intramember_inx_df = pd.DataFrame([inx.record for inx in
                           intramember_interactions[HydrogenBondType]])
