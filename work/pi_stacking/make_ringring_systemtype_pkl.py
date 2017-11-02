import os.path as osp

work_dir = "/home/salotz/Dropbox/devel/mastic/work/pi_stacking"

ref_ring_path = osp.join(work_dir, "benzene_hex.pdb")

from rdkit import Chem

ref_ring_PDB_rdkit = Chem.MolFromPDBFile(ref_ring_path, removeHs=False, sanitize=False)

from mastic.interfaces.rdkit import RDKitMoleculeWrapper

ring_rdkit_wrapper = RDKitMoleculeWrapper(ref_ring_PDB_rdkit, mol_name="6ring")

ref_ring_coords = ring_rdkit_wrapper.get_conformer_coords(0)

ring_Molecule = ring_rdkit_wrapper.make_molecule_type(find_features=True)

import pickle

import mastic.system as masticsys

member_types = [ring_Molecule, ring_Molecule]
system_attrs = {'molecule_source' : 'rdkit'}
Ring_Ring_System = masticsys.SystemType("Ring_Ring_System",
                                                member_types=member_types,
                                                **system_attrs)

# when we make associations for assymmetric interactions we need to
# define an association of A -> B and B -> A so we define the receptor
# -> ligand interactions and ligand -> receptor interactions, this
# really only means the donors -> acceptors from the members.


selection_map_AB = [(0, None), (1, None)]
selection_types = [None, None]
assoc1_attrs = {'info' : 'ring1-ring2'}
Ring1_Ring2_Association = \
            masticsys.AssociationType("Ring1_Ring2_Association",
                                    system_type=Ring_Ring_System,
                                    selection_map=selection_map_AB,
                                    selection_types=selection_types,
                                    **assoc1_attrs)
Ring_Ring_System.add_association_type(Ring1_Ring2_Association)

selection_map_BA = selection_map_AB[::-1]
assoc2_attrs = {'info' : 'ring2-ring1'}
Ring2_Ring1_Association = \
            masticsys.AssociationType("Ring2_Ring1_Association",
                                    system_type=Ring_Ring_System,
                                    selection_map=selection_map_BA,
                                    selection_types=selection_types,
                                    **assoc2_attrs)
Ring_Ring_System.add_association_type(Ring2_Ring1_Association)

import pickle

system_pkl_path = osp.join(".", "Ring_Ring_SystemType.pkl")
with open(system_pkl_path, 'wb') as wf:
    pickle.dump(Ring_Ring_System, wf)
