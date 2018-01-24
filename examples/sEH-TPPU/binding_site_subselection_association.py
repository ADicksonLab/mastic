import os.path as osp
import pickle

import numpy as np

from mastic.system import AssociationType
from mastic.molecule import MoleculeTypeAtomSelection

# load the SystemType we will add associations to
system_type_pkl_path = osp.join(".", "sEH_TPPU_SystemType.pkl")
with open(system_type_pkl_path, 'rb') as rf:
    sEH_TPPU_SystemType = pickle.load(rf)

# substantiate a crystal structure system so we can figure out where
# the bidning site is
# load the coordinates for the members
member_coords = [np.load(osp.realpath(osp.join('.', 'TPPU_coords.npy'))),
                 np.load(osp.realpath(osp.join('.', 'sEH_coords.npy')))]

# substantiate the system
cryst_system = sEH_TPPU_SystemType.to_system(member_coords)

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
# member. The MoleculeTypeAtomSelection allows for selection of atoms in a
# Molecule or MoelculeType.
selection_types = [MoleculeTypeAtomSelection, None]


# make the actual association
sehBS_tppu_assoc = AssociationType("sEHBS-TPPU",
                                   system_type=sEH_TPPU_SystemType,
                                   selection_map=selection_map,
                                   selection_types=selection_types
                                       )

# now you can add it to the original SystemType if you want, or you
# can use it to generate interaction classes for a particular
# interaction type (see other examples)
sEH_TPPU_SystemType.add_association_type(sehBS_tppu_assoc)


# then if our interaction is assymetric (which is the case for
# HydrogenBondType) we need to do it the other way around.
selection_map = [(0, None), (1, binding_site_atom_idxs)]
selection_types = [None, MoleculeTypeAtomSelection]
tppu_sehBS_assoc = AssociationType("TPPU-sEHBS",
                                   system_type=sEH_TPPU_SystemType,
                                   selection_map=selection_map,
                                   selection_types=selection_types
                                       )
sEH_TPPU_SystemType.add_association_type(tppu_sehBS_assoc)

