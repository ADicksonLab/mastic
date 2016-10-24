from rdkit import Chem
import mast.tests.data as mastdata

BEN_MOL_rdkit = Chem.MolFromMolBlock(mastdata.benzamidine_MOL, sanitize=True)
BEN_PDB_rdkit = Chem.MolFromPDBBlock(mastdata.BEN_Hs_3ptb, removeHs=False, sanitize=True)
trypsin_rdkit = Chem.MolFromPDBBlock(mastdata.trypsin_3ptb, removeHs=False, sanitize=True)

from mast.interfaces.rdkit import AssignBondOrdersFromTemplate

BEN_rdkit = AssignBondOrdersFromTemplate(BEN_MOL_rdkit, BEN_PDB_rdkit)

from mast.interfaces.rdkit import RDKitMoleculeWrapper

BEN_rdkit_wrapper = RDKitMoleculeWrapper(BEN_rdkit, mol_name="BEN")
trypsin_rdkit_wrapper = RDKitMoleculeWrapper(trypsin_rdkit, mol_name="Trypsin")

BEN_coords = BEN_rdkit_wrapper.get_conformer_coords(0)
trypsin_coords = trypsin_rdkit_wrapper.get_conformer_coords(0)

BEN_Molecule = BEN_rdkit_wrapper.make_molecule_type(find_features=True)

Trypsin_Molecule = trypsin_rdkit_wrapper.make_molecule_type(find_features=True)

import os.path as osp
import pickle

trypsin_pkl_path = osp.join(".", "TrypsinMoleculeType.pkl")
with open(trypsin_pkl_path, 'wb') as wf:
    pickle.dump(Trypsin_Molecule, wf)

import mast.system as mastsys

member_types = [BEN_Molecule, Trypsin_Molecule]
system_attrs = {'molecule_source' : 'rdkit'}
Trypsin_Benzamidine_System = mastsys.SystemType("Trypsin_Benzamidine_System",
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
Trypsin_Benzamidine_Association = \
            mastsys.AssociationType("Trypsin_Benzamidine_Association",
                                    system_type=Trypsin_Benzamidine_System,
                                    selection_map=selection_map_BA,
                                    selection_types=selection_types,
                                    **rec_lig_attrs)

Trypsin_Benzamidine_System.add_association_type(Trypsin_Benzamidine_Association)

selection_map_AB = selection_map_BA[::-1]
lig_rec_attrs = {'info' : 'ligand-receptor'}
Benzamidine_Trypsin_Association = \
            mastsys.AssociationType("Benzamidine_Trypsin_Association",
                                    system_type=Trypsin_Benzamidine_System,
                                    selection_map=selection_map_AB,
                                    selection_types=selection_types,
                                    **lig_rec_attrs)

Trypsin_Benzamidine_System.add_association_type(Benzamidine_Trypsin_Association)


# put them together in the order they are as system members
member_coords = [BEN_coords, trypsin_coords]

# substantiate the system
system = Trypsin_Benzamidine_System.to_system(member_coords)

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
TrypsinBS_Benzamidine_assoc = mastsys.AssociationType("TrypsinBS-Benzamidine",
                                         system_type=Trypsin_Benzamidine_System,
                                         selection_map=selection_map,
                                         selection_types=selection_types)

# add it to the system
Trypsin_Benzamidine_System.add_association_type(TrypsinBS_Benzamidine_assoc)

import os.path as osp
import pickle

system_pkl_path = osp.join(".", "Trypsin_Benzamidine_SystemType.pkl")
with open(system_pkl_path, 'wb') as wf:
    pickle.dump(Trypsin_Benzamidine_System, wf)
