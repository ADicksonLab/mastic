"""Example of how to perform profiling of a subset of a full
system. There is multiple steps:

1. Determine a selection of atoms from each of the molecules that you
   want to profile interactions between.

2. Create an AssociationType between these selections of atoms (from
   the unsubstantiated SystemType).

3. Generate interaction classes (instances of the InteractionType
   e.g. HydrogenBondType) from the AssociationType.

4. Create an InteractionSpace from only these interaction classes.

5. Create a Profiler from this interaction space.

6. Substantiate the system(s) you want to profile with the SystemType
and coordinates.

7. Profile interactions of the system(s) with the profiler.

8. Save the interactions for analysis.

"""

import os.path as osp
import pickle

import numpy as np

from mastic.system import AssociationType
from mastic.molecule import MoleculeAtomSelection
from mastic.interaction_space import InteractionSpace
from mastic.interactions.hydrogen_bond import HydrogenBondType
from mastic.profile import InxSpaceProfiler

# load the SystemType we will add associations to
system_type_pkl_path = osp.join(".", "sEH_TPPU_SystemType.pkl")
with open(system_type_pkl_path, 'rb') as rf:
    sEH_TPPU_SystemType = pickle.load(rf)

# load a system we can find what the binding site indices are from
system_pkl_path = osp.join(".", "sEH_TPPU_System_cryst.pkl")
with open(system_pkl_path, 'rb') as rf:
    cryst_system = pickle.load(rf)

# find atoms in the binding site using a cutoff distance from the
# ligand
binding_site_cutoff_dist = 4 #in Angstroms \AA

# find the atoms within this distance
binding_site_atoms = cryst_system.molecules[0].atoms_within_distance(
    binding_site_cutoff_dist)

# get the indices of these atoms to define the AssociationType
binding_site_atom_idxs = [cryst_system.molecules[1].atoms.index(atom) for
                          atom in binding_site_atoms]

# you might also want to get the pdb serial numbers so you can
# visually check to see where these atoms are
binding_site_atom_serials = [atom.atom_type.pdb_serial_number for atom
                             in binding_site_atoms]

# the selection map tells the association the index of the member and
# the indices of the atoms to include as one component of the
# association. By selection None as the indices no selection will be
# made and the whole molecule will be a component
selection_map = [(1, binding_site_atom_idxs), (0, None)]

# The selection types correspond to the elements in the selection map
# and tell the AssociationType what kind of selection to make on the
# molecule. Setting one of them to None should mean the selection map
# also had no indices selected and it should use the whole system
# member. The MoleculeAtomSelection allows for selection of atoms in a
# Molecule or MoelculeType.
selection_types = [MoleculeAtomSelection, None]


# make the actual association between the selections
sehBS_tppu_assoc = AssociationType("sEHBS-TPPU",
                                   system_type=sEH_TPPU_SystemType,
                                   selection_map=selection_map,
                                   selection_types=selection_types
                                       )

# create only these interaction classes in this association
inx_classes = HydrogenBondType.interaction_classes(sehBS_tppu_assoc)

# create an InteractionSpace for the system
inxspace = InteractionSpace(sEH_TPPU_SystemType)

# add these inxclasses to the inxspace
inxspace.add_inx_classes(inx_classes)

# We probably want to save a record of this interaction space so that
# we can find out more information about the hits later once we have
# found them from profiling. This includes the interaction type,
# association type, the indices of the members for each association
# member, and the feature types.
inxspace_df = inxspace.to_dataframe()


# we now have just the interaction space that we want for profiling,
# so we make a profiler with this inxspace
profiler = InxSpaceProfiler(inxspace)


# substantiate the updated SystemType with the new association
# load the coordinates for the members
member_coords = [np.load('TPPU_coords.npy'), np.load("sEH_coords.npy")]
system = sEH_TPPU_SystemType.to_system(member_coords)

# profile the system with our profiler
profile = profiler.profile(system)


# Now we need to get these interactions into a format that we can save
# to disk and analyze later. We note that pickling is not really
# efficient here, due to the potentially large amount of data this
# could produce from large datasets, and (unfortunately) due to the
# current (flawed) implementation of mastic that makes pickling
# molecules difficult (but doable if you increase python's recursion
# limit)

# Assuming you have a record of all the interaction classes in your
# interaction space you would only need to save the hits and maybe the
# parameters by which they satisfied the constraints of the
# interaction type.

# For this you can use a convenient output to a pandas Dataframe
hits_df = profile.hit_inx_df()

# or you can get just the records themselves (which are used to make
# the dataframe anyhow)
hit_records = profile.hit_inx_records()

# there are some other methods and attributes though that may be
# useful including getting just a list of the hit indices for this
# profile
hit_idxs = profile.hit_idxs

# the interactions that were hits
hit_inxs = profile.hit_inxs

# or a feature vector which is the length of the interaction space
# with 1s for the hits and 0s elsewhere
inxspace_vec = profile.vector

# there is other information you may need to do a thorough
# analysis. For instance to get the index of an atom in an interaction
# from the SystemType or MoleculeType you would need information from
# the FeatureType itself. However, given the feature type from the
# inxspace data that was save above we can find it in the SystemType
# (MoleculeType) that we already have saved as a pickle this
# information can be found.
