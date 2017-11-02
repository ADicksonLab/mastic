import os.path as osp

work_dir = "/home/salotz/Dropbox/devel/mastic/work/pi_stacking"

ref_benzene_PDB_path = osp.join(work_dir, "ref_benzene.pdb")
ref_benzene_MOL_path = osp.join(work_dir, "benzene.mol")

from rdkit import Chem

ref_benzene_PDB_rdkit = Chem.MolFromPDBFile(ref_benzene_PDB_path, removeHs=False, sanitize=False)
ref_benzene_MOL_rdkit = Chem.MolFromMolFile(ref_benzene_MOL_path, sanitize=True)

from mastic.interfaces.rdkit import AssignBondOrdersFromTemplate

ref_benzene_rdkit = AssignBondOrdersFromTemplate(ref_benzene_MOL_rdkit, ref_benzene_PDB_rdkit)

from mastic.interfaces.rdkit import RDKitMoleculeWrapper

benzene_rdkit_wrapper = RDKitMoleculeWrapper(ref_benzene_rdkit, mol_name="benzene")

ref_benzene_coords = benzene_rdkit_wrapper.get_conformer_coords(0)

Benzene_Molecule = benzene_rdkit_wrapper.make_molecule_type(find_features=True)

import os.path as osp

import mastic.system as masticsys

member_types = [Benzene_Molecule, Benzene_Molecule]
system_attrs = {'molecule_source' : 'rdkit'}
Benzene_Benzene_System = masticsys.SystemType("Benzene_Benzene_System",
                                                member_types=member_types,
                                                **system_attrs)

# when we make associations for assymmetric interactions we need to
# define an association of A -> B and B -> A so we define the receptor
# -> ligand interactions and ligand -> receptor interactions, this
# really only means the donors -> acceptors from the members.


selection_map_AB = [(0, None), (1, None)]
selection_types = [None, None]
assoc1_attrs = {'info' : 'benzene1-benzene2'}
Benzene1_Benzene2_Association = \
            masticsys.AssociationType("Benzene1_Benzene2_Association",
                                    system_type=Benzene_Benzene_System,
                                    selection_map=selection_map_AB,
                                    selection_types=selection_types,
                                    **assoc1_attrs)
Benzene_Benzene_System.add_association_type(Benzene1_Benzene2_Association)

selection_map_BA = selection_map_AB[::-1]
assoc2_attrs = {'info' : 'benzene2-benzene1'}
Benzene2_Benzene1_Association = \
            masticsys.AssociationType("Benzene2_Benzene1_Association",
                                    system_type=Benzene_Benzene_System,
                                    selection_map=selection_map_BA,
                                    selection_types=selection_types,
                                    **assoc2_attrs)
Benzene_Benzene_System.add_association_type(Benzene2_Benzene1_Association)

import pickle

system_pkl_path = osp.join(".", "Benzene_Benzene_SystemType.pkl")
with open(system_pkl_path, 'wb') as wf:
    pickle.dump(Benzene_Benzene_System, wf)
