import pandas as pd

from rdkit import Chem
from rdkit.Chem import AllChem
import os.path as osp

import mast.selection as mastsel
import mast.molecule as mastmol
import mast.system as mastsys
import mast.interactions as mastinx
from mast.interactions.hydrogen_bond import HydrogenBondType

from mast.interfaces.rdkit import RDKitMoleculeWrapper

import mast.config.interactions as mastinxconfig

import mast.tests.data as mastdata

from rdkit import Chem


# without Hs straight from pdb
BEN_rdkit = Chem.MolFromPDBBlock(mastdata.BEN_3ptb, removeHs=False, sanitize=True)
trypsin_rdkit = Chem.MolFromPDBBlock(mastdata.trypsin_3ptb, removeHs=False, sanitize=True)

BEN_rdkit_wrapper = RDKitMoleculeWrapper(BEN_rdkit, mol_name="BEN")
trypsin_rdkit_wrapper = RDKitMoleculeWrapper(trypsin_rdkit, mol_name="Trypsin")

BEN_coords = BEN_rdkit_wrapper.get_conformer_coords(0)
trypsin_coords = trypsin_rdkit_wrapper.get_conformer_coords(0)
member_coords = [BEN_coords, trypsin_coords]

BEN_Molecule = BEN_rdkit_wrapper.make_molecule_type(find_features=True)
# Trypsin_Molecule = 

BEN_mol = BEN_Molecule.to_molecule(BEN_coords)
trypsin_mol = Trypsin_Molecule.to_molecule(trypsin_coords)

member_types = [BEN_Molecule, Trypsin_Molecule]
system_attrs = {'name' : 'trypsin-benzamidine-complex'}

Trypsin_Benzamidine_System = mastsys.SystemType("Trypsin_Benzamidine_System",
                                                member_types=member_types,
                                                **system_attrs)

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

# from mast.molecule import Molecules
print("testing Hbond interaction between molecules in the receptor ligand association")
tryp_ben_assoc = trypsys.associations[0]

intermember_key_pairs, intermember_interactions = \
tryp_ben_assoc.profile_interactions([HydrogenBondType])

intramember_key_pairs, intramember_interactions = \
tryp_ben_assoc.profile_interactions([HydrogenBondType],
                                    intramember_interactions=True)

intermember_inx_class_df = pd.DataFrame([inx.interaction_class.record for
                                 inx in intermember_interactions[HydrogenBondType]])
intramember_inx_class_df = pd.DataFrame([inx.interaction_class.record for
                                 inx in intramember_interactions[HydrogenBondType]])

intermember_inx_df = pd.DataFrame([inx.record for inx in
                           intermember_interactions[HydrogenBondType]])
intramember_inx_df = pd.DataFrame([inx.record for inx in
                           intramember_interactions[HydrogenBondType]])
