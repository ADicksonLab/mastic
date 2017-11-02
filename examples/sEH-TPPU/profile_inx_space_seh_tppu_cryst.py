import os.path as osp
import pickle
import numpy as np

from mastic.interactions.hydrogen_bond import HydrogenBondType
import mastic.profile as masticprof
from mastic.interaction_space import InteractionSpace

# load the system type pickle in
system_pkl_path = osp.join(".", "sEH_TPPU_SystemType.pkl")
with open(system_pkl_path, 'rb') as rf:
    sEH_TPPU_SystemType = pickle.load(rf)

# generate the interaction space we will be profiling
inx_space = InteractionSpace(sEH_TPPU_SystemType)

# we want associations for all combinations of members for a degree 2
# interaction (e.g. hydrogen bonds)
# so we could use this method I've commented out or just define it ourselves
# assoc_terms = sEH_TPPU_SystemType.association_polynomial(
#     degree=2,
#     permute=True,
#     replace=True,
#     return_idxs=True)

assoc_terms = [(0,1), (1,0)]
# make the unit associations, interaction classes, and add to interaction space
for assoc_term in assoc_terms:
    # make the unit AssociationTypes
    assoc_idx = sEH_TPPU_SystemType.make_unit_association_type(assoc_term)
    association_type = sEH_TPPU_SystemType.association_types[assoc_idx]

    # make HydrogenBondType interaction classes for this association
    # in the inx_space
    inx_space.add_association_subspace(association_type, HydrogenBondType)

# make a Profiler for the inx space
profiler = masticprof.InxSpaceProfiler(inx_space)


# load the coordinates for the members
member_coords = [np.load('TPPU_coords.npy'), np.load("sEH_coords.npy")]
# substantiate the system
system = sEH_TPPU_SystemType.to_system(member_coords)

# profile the interaction space over the system
system_profile = profiler.profile(system)
