import os.path as osp
import pickle

from mast.interactions.hydrogen_bond import HydrogenBondType

# load the system type pickle in
system_pkl_path = osp.join(".", "Trypsin_Benzamidine_SystemType.pkl")
with open(system_pkl_path, 'rb') as rf:
    Trypsin_Benzamidine_SystemType = pickle.load(rf)

# load the crystal structure substantiated system
system_cryst_pkl_path = osp.join(".", "Trypsin_Benzamidine_System_cryst.pkl")
with open(system_cryst_pkl_path, 'rb') as pkl_rf:
    Trypsin_Benzamidine_System_cryst = pickle.load(pkl_rf)

# use the association polynomial function of the system
assoc_terms = Trypsin_Benzamidine_SystemType.association_polynomial(
    # input the degree of the interaction
    interaction_degree=HydrogenBondType.degree,
    # return the indices of the system members instead of the members
    # themselves
    return_idxs=True,
    # whether or not the interaction is symmetric or not
    commutative=False)

# this gives them to you organized by which association they fall under
hbond_inx_classes_assocs = Trypsin_Benzamidine_SystemType.interaction_space(
                                     assoc_terms, HydrogenBondType)

# if you want the whole collection of interaction classes in one list
from itertools import chain
hbond_inx_classes = list(chain(*[inx_classes for inx_classes in
                                 hbond_inx_classes_assocs.values()]))

# however we are only interested in one association
rec_lig_association = Trypsin_Benzamidine_System_cryst.associations[1]
rec_lig_association_type = rec_lig_association.association_type
rec_lig_member_idxs = rec_lig_association_type.member_idxs
rec_lig_inx_classes = hbond_inx_classes_assocs[rec_lig_member_idxs]

# profile that association
import ipdb; ipdb.set_trace()
rec_lig_inxs = rec_lig_association.profile_interactions(
    [HydrogenBondType],
    interaction_classes=rec_lig_inx_classes)[HydrogenBondType]
