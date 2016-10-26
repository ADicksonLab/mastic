import os.path as osp

from rdkit import Chem
import mast.tests.data as mastdata

TPPU_MOL_path = osp.join(".", "TPPU.mol")
TPPU_MOL_rdkit = Chem.MolFromMolFile(TPPU_MOL_path, sanitize=True)
TPPU_PDB_path = osp.join(".", "TPPU.pdb")
TPPU_PDB_rdkit = Chem.MolFromPDBFile(TPPU_PDB_path, removeHs=False, sanitize=False)
seh_PDB_path = osp.join(".", "sEH.pdb")
seh_rdkit = Chem.MolFromPDBFile(seh_PDB_path, removeHs=False, sanitize=False)

from mast.interfaces.rdkit import AssignBondOrdersFromTemplate

TPPU_rdkit = AssignBondOrdersFromTemplate(TPPU_MOL_rdkit, TPPU_PDB_rdkit)

from mast.interfaces.rdkit import RDKitMoleculeWrapper

TPPU_rdkit_wrapper = RDKitMoleculeWrapper(TPPU_rdkit, mol_name="TPPU")
seh_rdkit_wrapper = RDKitMoleculeWrapper(seh_rdkit, mol_name="sEH")

TPPU_coords = TPPU_rdkit_wrapper.get_conformer_coords(0)
seh_coords = seh_rdkit_wrapper.get_conformer_coords(0)

TPPU_Molecule = TPPU_rdkit_wrapper.make_molecule_type(find_features=True)

seh_Molecule = seh_rdkit_wrapper.make_molecule_type(find_features=True)

import os.path as osp
import pickle

seh_pkl_path = osp.join(".", "sEHMoleculeType.pkl")
with open(seh_pkl_path, 'wb') as wf:
    pickle.dump(seh_Molecule, wf)

import mast.system as mastsys

member_types = [TPPU_Molecule, seh_Molecule]
system_attrs = {'molecule_source' : 'rdkit'}
seh_tppu_System = mastsys.SystemType("sEH_TPPU_System",
                                                member_types=member_types,
                                                **system_attrs)

# when we make associations for assymmetric interactions we need to
# define an association of A -> B and B -> A so we define the receptor
# -> ligand interactions and ligand -> receptor interactions, this
# really only means the donors -> acceptors from the members.

# these variables are used to define selections of only part of system
# members, which we will ignore for now
selection_map_BA = [(1, None), (0, None)]
selection_types = [None, None]

rec_lig_attrs = {'info' : 'receptor-ligand'}
seh_tppu_Association = \
            mastsys.AssociationType("sEH_TPPU_Association",
                                    system_type=seh_tppu_System,
                                    selection_map=selection_map_BA,
                                    selection_types=selection_types,
                                    **rec_lig_attrs)

seh_tppu_System.add_association_type(seh_tppu_Association)

selection_map_AB = selection_map_BA[::-1]
lig_rec_attrs = {'info' : 'ligand-receptor'}
tppu_seh_Association = \
            mastsys.AssociationType("TPPU_sEH_Association",
                                    system_type=seh_tppu_System,
                                    selection_map=selection_map_AB,
                                    selection_types=selection_types,
                                    **lig_rec_attrs)

seh_tppu_System.add_association_type(tppu_seh_Association)


# put them together in the order they are as system members
member_coords = [TPPU_coords, seh_coords]

# substantiate the system
system = seh_tppu_System.to_system(member_coords)

binding_site_cutoff_dist = 4 #in Angstroms \AA

# find the atoms within this distance
binding_site_atoms = system.molecules[0].atoms_within_distance(
    binding_site_cutoff_dist)

# get the indices of these atoms to define the AssociationType
binding_site_atom_idxs = [system.molecules[1].atoms.index(atom) for
                          atom in binding_site_atoms]

# you might also want to get the pdb serial numbers so you can
# visually check to see where these atoms are
binding_site_atom_serials = [atom.atom_type.pdb_serial_number for atom
                             in binding_site_atoms]

import mast.molecule as mastmol

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
selection_types = [mastmol.MoleculeAtomSelection, None]

# instantiate the association
sehBS_tppu_assoc = mastsys.AssociationType("sEHBS-TPPU",
                                         system_type=seh_tppu_System,
                                         selection_map=selection_map,
                                         selection_types=selection_types)

# add it to the system
seh_tppu_System.add_association_type(sehBS_tppu_assoc)

import os.path as osp
import pickle

system_pkl_path = osp.join(".", "sEH_TPPU_SystemType.pkl")
with open(system_pkl_path, 'wb') as wf:
    pickle.dump(seh_tppu_System, wf)
