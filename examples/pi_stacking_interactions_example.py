import pandas as pd

from rdkit import Chem
from rdkit.Chem import AllChem, rdMolAlign
import os.path as osp

import mast.selection as mastsel
import mast.molecule as mastmol
import mast.system as mastsys
import mast.interactions as mastinx
from mast.interactions.pi_stacking import PiStackingType
from mast.interfaces.rdkit import RDKitMoleculeWrapper, AssignBondOrdersFromTemplate
import mast.config.interactions as mastinxconfig
import mast.tests.data as mastdata

# to get the aromaticity for the Molecule_Type we need to load the
# .mol file and then we will associate this with the coordinates of
# the PDB file

# load the necessary files with rdkit
BEN_MOL_rdkit = Chem.MolFromMolBlock(mastdata.benzamidine_MOL, sanitize=True)
BEN_PDB_rdkit = Chem.MolFromPDBBlock(mastdata.BEN_Hs_3ptb, removeHs=False, sanitize=True)
trypsin_rdkit = Chem.MolFromPDBBlock(mastdata.trypsin_3ptb, removeHs=False, sanitize=True)

# assign the bond orders from the .mol template for benzamidine
print("assigning bond orders")
BEN_rdkit = AssignBondOrdersFromTemplate(BEN_MOL_rdkit, BEN_PDB_rdkit)

# put it in the mast.interfaces.rdkit wrapper
print("making wrappers")
BEN_rdkit_wrapper = RDKitMoleculeWrapper(BEN_rdkit, mol_name="BEN")
trypsin_rdkit_wrapper = RDKitMoleculeWrapper(trypsin_rdkit, mol_name="Trypsin")

# get the coordinates
BEN_coords = BEN_rdkit_wrapper.get_conformer_coords(0)
trypsin_coords = trypsin_rdkit_wrapper.get_conformer_coords(0)
member_coords = [BEN_coords, trypsin_coords]

# make a MoleculeType
print("making MoleculeType")
BEN_Molecule = BEN_rdkit_wrapper.make_molecule_type(find_features=True)
# Trypsin_Molecule = trypsin_rdkit_wrapper.make_molecule_type(find_features=True)
Trypsin_Molecule = mastdata.Trypsin_Hs_Molecule
# substantiate a Molecule (don't need this)
# BEN_mol = BEN_Molecule.to_molecule(BEN_coords)
# trypsin_mol = Trypsin_Molecule.to_molecule(trypsin_coords)

# make the SystemType
print("making the SystemType")
member_types = [BEN_Molecule, Trypsin_Molecule]
system_attrs = {'name' : 'trypsin-benzamidine-complex'}

Trypsin_Benzamidine_System = mastsys.SystemType("Trypsin_Benzamidine_System",
                                                member_types=member_types,
                                                **system_attrs)

# make an AssociationType for the receptor ligand complex
print("making an AssociationType of the receptor and ligand in the Trypsin_Benzamidine_System")
rec_lig_attrs = {'name' : 'trypsin-benzamidine-complex'}
# rec_lig_attrs['ligand_type'] = ben_type
# rec_lig_attrs['receptor_type'] = trypsin_type
selection_map = {0 : None, 1 : None}
selection_types = [None, None]
Trypsin_Benzamidine_Association = \
            mastsys.AssociationType("Trypsin_Benzamidine_Association",
                                    system_type=Trypsin_Benzamidine_System,
                                    selection_map=selection_map,
                                    selection_types=selection_types,
                                    **rec_lig_attrs)

# add it to the SystemType
Trypsin_Benzamidine_System.add_association_type(Trypsin_Benzamidine_Association)

# now when we make the system the selections are put into an
# Association that can be profiled
trypsys = Trypsin_Benzamidine_System.to_system(member_coords)


print("testing Hbond interaction between molecules in the receptor ligand association")
tryp_ben_assoc = trypsys.associations[0]

# profile just the intermember(molecular) interactions
intermember_key_pairs, intermember_interactions = \
tryp_ben_assoc.profile_interactions([PiStackingType])

# profile the intermolecular interactions plus the intramolecular
# interactions
intramember_key_pairs, intramember_interactions = \
tryp_ben_assoc.profile_interactions([PiStackingType],
                                    intramember_interactions=True)

# present the data in dataframes
intermember_inx_class_df = pd.DataFrame([inx.interaction_class.record for
                                 inx in intermember_interactions[HydrogenBondType]])
intramember_inx_class_df = pd.DataFrame([inx.interaction_class.record for
                                 inx in intramember_interactions[HydrogenBondType]])

intermember_inx_df = pd.DataFrame([inx.record for inx in
                           intermember_interactions[HydrogenBondType]])
intramember_inx_df = pd.DataFrame([inx.record for inx in
                           intramember_interactions[HydrogenBondType]])
