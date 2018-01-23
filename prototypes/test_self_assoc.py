import os.path as osp
import pickle

import numpy as np

from mastic.system import AssociationType
from mastic.molecule import MoleculeTypeAtomSelection
from mastic.interactions.hydrogen_bond import HydrogenBondType

inputs_path = "../examples/sEH-TPPU"
# load the SystemType we will add associations to
system_type_pkl_path = osp.realpath(osp.join(inputs_path, "sEH_TPPU_SystemType.pkl"))
with open(system_type_pkl_path, 'rb') as rf:
    sEH_TPPU_SystemType = pickle.load(rf)


# substantiate a crystal structure system so we can figure out where
# the bidning site is
# load the coordinates for the members
member_coords = [np.load(osp.realpath(osp.join(inputs_path, 'TPPU_coords.npy'))),
                 np.load(osp.realpath(osp.join(inputs_path, 'sEH_coords.npy')))]

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
selection_map = [(1, binding_site_atom_idxs), (1, binding_site_atom_idxs)]

# The selection types correspond to the elements in the selection map
# and tell the AssociationType what kind of selection to make on the
# molecule. Setting one of them to None should mean the selection map
# also had no indices selected and it should use the whole system
# member. The MoleculeTypeAtomSelection allows for selection of atoms in a
# Molecule or MoelculeType.
selection_types = [MoleculeTypeAtomSelection, MoleculeTypeAtomSelection]


# make the actual association
sehBS_sehBS_assoc = AssociationType("sEHBS-sEHBS",
                                   system_type=sEH_TPPU_SystemType,
                                   selection_map=selection_map,
                                   selection_types=selection_types
                                       )

# make inxclasses from this
BS_only_inxclasses = HydrogenBondType.interaction_classes(sehBS_sehBS_assoc)


# for comparison make inxclasses for the whole protein-protein


selection_map = [(1, None), (1, None)]

# The selection types correspond to the elements in the selection map
# and tell the AssociationType what kind of selection to make on the
# molecule. Setting one of them to None should mean the selection map
# also had no indices selected and it should use the whole system
# member. The MoleculeTypeAtomSelection allows for selection of atoms in a
# Molecule or MoelculeType.
selection_types = [None, None]


# make the actual association
seh_seh_assoc = AssociationType("sEH-sEH",
                                   system_type=sEH_TPPU_SystemType,
                                   selection_map=selection_map,
                                   selection_types=selection_types
                                       )

# make inxclasses from this
self_inxclasses = HydrogenBondType.interaction_classes(seh_seh_assoc)
