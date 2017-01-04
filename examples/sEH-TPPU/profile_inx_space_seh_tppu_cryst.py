import os.path as osp
import pickle
import numpy as np

from mast.interactions.hydrogen_bond import HydrogenBondType

# load the system type pickle in
system_pkl_path = osp.join(".", "sEH_TPPU_SystemType.pkl")
with open(system_pkl_path, 'rb') as rf:
    sEH_TPPU_SystemType = pickle.load(rf)

# we generate interaction space for the system for all combinations of
# members for a degree 2 interaction (e.g. hydrogen bonds)
assoc_terms = sEH_TPPU_SystemType.association_polynomial(
    degree=2,
    permute=True,
    replace=True,
    return_idxs=True)

# this gives them to you organized by which association they fall under
hbond_inx_class_idxs = sEH_TPPU_SystemType.generate_interaction_space(
    assoc_terms, HydrogenBondType)

# load the coordinates for the members
member_coords = [np.load('TPPU_coords.npy'), np.load("sEH_coords.npy")]

# substantiate the system
system = sEH_TPPU_SystemType.to_system(member_coords)

# profile the associations you are interested in
lig_rec_idx = assoc_terms.index((0,1))
rec_lig_idx = assoc_terms.index((1,0))

lig_rec_inxs = system.associations[lig_rec_idx].profile_interactions(
    [HydrogenBondType])[HydrogenBondType]

rec_lig_inxs = system.associations[rec_lig_idx].profile_interactions(
    [HydrogenBondType])[HydrogenBondType]
