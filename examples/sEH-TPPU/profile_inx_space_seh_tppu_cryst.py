"""Script to do a simple profiling of intermolecular hydrogen bonds."""
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

# the (unit) association terms are used to simply build associations (which
# are directional, thus terms must be permuted) between whole
# molecules. All possible assocition terms for a two molecule system
# is: [(0,1), (1,0), (1,1), (0,0)]. Where (0,1) and (1,0) cover
# assymetric interactions (e.g. donor-acceptor, as opposed to
# symmetric interactions like pi-stacking) between two molecules and
# (1,1) and (0,0) cover intramolecular interactions

# we do not want intramolecular interactions (intraprotein interaction
# space is very large and takes a long time and is probably better
# done with some careful selections and associations)
assoc_terms = [(0,1), (1,0)]

# make the unit associations, interaction classes, and add to interaction space
for assoc_term in assoc_terms:
    # make the unit AssociationTypes from the association terms we defined
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

# the profile has the hit indexes for the interaction space (i.e. the
# indices in the interaction space where interactions were observed)
# and the actual interaction objects, which contains the actual
# interaction parameter values and the features, which contain the
# atoms and their actual coordinates

# Output results as dataframe
hits_df = system_profile.hit_inx_df()
