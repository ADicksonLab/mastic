import os.path as osp
import pickle
import itertools as it

from mastic.interactions.hydrogen_bond import HydrogenBondType
from mastic.interactions.pi_stacking import PiStackingType

import sys
sys.recursion_depth

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
hbond_inx_classes = list(it.chain(*[inx_classes for inx_classes in
                                 hbond_inx_classes_assocs.values()]))

# however we are only interested in one association
rec_lig_association = Trypsin_Benzamidine_System_cryst.associations[1]
rec_lig_association_type = rec_lig_association.association_type
rec_lig_member_idxs = rec_lig_association_type.member_idxs
rec_lig_inx_classes = hbond_inx_classes_assocs[rec_lig_member_idxs]

# profile that association
rec_lig_inxs = rec_lig_association.profile_interactions(
    [HydrogenBondType],
    interaction_classes=rec_lig_inx_classes)[HydrogenBondType]

import ipdb; ipdb.set_trace()
# now for pi-stacking, we only want intraprotein and ligand-protein
# interactions, and it is commutative so we only need one inter- term
pistack_inx_classes_assocs = Trypsin_Benzamidine_SystemType.interaction_space(
                                     [(1,1), (0,1)], PiStackingType)

intraprotein_pistack_inxs = Trypsin_Benzamidine_System_cryst.associations[2].profile_interactions(
    [PiStackingType],
    interaction_classes=pistack_inx_classes_assocs[(1,1)])

PL_pistack_inxs = Trypsin_Benzamidine_System_cryst.associations[2].profile_interactions(
    [PiStackingType],
    interaction_classes=pistack_inx_classes_assocs[(1,1)])
