import os.path as osp
import pickle

import mastic.interactions.pi_stacking as pinx
import mastic.interactions.hydrogen_bond as hinx

work_dir = "/home/salotz/Dropbox/devel/mastic/work/pi_stacking"

# load the SystemType
benzene_system_pkl_path = osp.join(work_dir, "Benzene_Benzene_SystemType.pkl")
with open(benzene_system_pkl_path, 'rb') as rf:
    Benzene_Benzene_SystemType = pickle.load(rf)

# load the coordinates for the reference benzene
ref_benzene_PDB_path = osp.join(work_dir, "ref_benzene.pdb")

from rdkit import Chem

ref_benzene_rdkit = Chem.MolFromPDBFile(ref_benzene_PDB_path, removeHs=False, sanitize=False)

from mastic.interfaces.rdkit import RDKitMoleculeWrapper

benzene_rdkit_wrapper = RDKitMoleculeWrapper(ref_benzene_rdkit, mol_name="benzene")

ref_benzene_coords = benzene_rdkit_wrapper.get_conformer_coords(0)

from mastic.interactions.pi_stacking import PiStackingType

# get the interaction space for pi-stacking
pistack_inx_classes = Benzene_Benzene_SystemType.interaction_space([(0,1)], PiStackingType)[(0,1)]

# profile the stacked one that should qualify
stacked_member_coords = [ref_benzene_coords, test_benzenes['stacked']]
stacked_system = Benzene_Benzene_SystemType.to_system(stacked_member_coords)

# profile the interactions between the two rings
stacked_inxs = stacked_system.associations[0].\
               profile_interactions([PiStackingType],
                            interaction_classes=pistack_inx_classes)\
                            [PiStackingType]

# substantiate the systems and profile each one
test_inxs = {}
test_failed_hits = {}
for test_name, test_benzene in test_benzenes.items():
    member_coords = [ref_benzene_coords, test_benzene]
    system = Benzene_Benzene_SystemType.to_system(member_coords)

    # profile the interactions between the two rings
    failed_hits, all_inxs = system.associations[0].\
           profile_interactions([PiStackingType],
                                interaction_classes=pistack_inx_classes,
                                return_failed_hits=True)
    inxs = all_inxs[PiStackingType]
    test_failed_hits[test_name] = failed_hits
    test_inxs[test_name] = inxs
