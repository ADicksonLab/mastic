import os.path as osp
import rdkit
from rdkit.Chem import MolFromPDBFile
trypsin_dir = osp.expanduser("~/Dropbox/lab/trypsin/tryp_net")
trypsin_pdb_path = osp.join(trypsin_dir,  "trypsin.pdb")
trypsin = MolFromPDBFile(trypsin_pdb_path, removeHs=False, sanitize=False)
ben_pdb_path = osp.join(trypsin_dir, "BEN.pdb")
ben = MolFromPDBFile(ben_pdb_path, removeHs=False, sanitize=False)

# set the path so we can import my module
import sys
sys.path.append("/home/salotz/Dropbox/devel")
from mastic.molecule import RDKitMoleculeType

print("loading RDKit molecules")
trypsin_type = RDKitMoleculeType(trypsin, mol_name="trypsin")
ben_type = RDKitMoleculeType(ben, mol_name="BEN")

print("finding features")
ben_type.find_features()
trypsin_type.find_features()

print("loading trajectories")
import mdtraj as mdj

traj_path = osp.join(trypsin_dir, "allframes.dcd")
traj_top = osp.join(trypsin_dir, "frame_0.pdb")
ben_indices = range(0,18)
trypsin_indices = range(18, 3247)

ben_traj = mdj.load_dcd(traj_path, top=traj_top, atom_indices=ben_indices)
trypsin_traj = mdj.load_dcd(traj_path, top=traj_top, atom_indices=trypsin_indices)


# slice only the frames we want
ben_coords = ben_traj.xyz[:]
trypsin_coords = trypsin_traj.xyz[:]

# the units are different in the trajectory
ben_coords = ben_coords * 10
trypsin_coords = trypsin_coords * 10

print("making the system")
from mastic.system import SystemType, System

sys_type = SystemType({'name' : 'trypsin-benzamidine-complex',
                          'trypsin_type' : trypsin_type,
                          'benzamidine_type' : ben_type})

from mastic.interactions import AssociationType, SystemAssociation, HydrogenBondType

rec_lig_attr = {'ligand_type' : ben_type,
                'receptor_type' : trypsin_type,
                'name' : 'trypsin_benzamidine-complex'}
rec_lig_type = AssociationType(rec_lig_attr)

import sys

def profile_coords(mol_types, coords, sys_type, assoc_sel_idxs, assoc_type, inx_type):

    # make the molecules
    molecules = []
    for i, mol_type in enumerate(mol_types):
        molecule = mol_type.to_molecule_from_coords(coords[i])
        molecule.make_feature_selections()
        molecules.append(molecule)

    # make the system
    system = System(molecules, system_type=sys_type)

    # make the associations
    members = [system[sel_idx] for sel_idx in assoc_sel_idxs]
    assoc = SystemAssociation(members=members,
                              system=system,
                              association_type=assoc_type)

    # profile intermolecular interactions
    assoc.profile_interactions([HydrogenBondType], between=Molecule)

    return assoc.interactions, system, assoc


def frame_profile(frame_idx, serial_data_path, pickle_path,
                  mol_types, coords, sys_type, assoc_sel_idxs, assoc_type, inx_type):

                
    inxs, system, assoc = profile_coords(mol_types, coords, sys_type, assoc_sel_idxs, assoc_type, inx_type)

    # data output
    inx_type.pdb_serial_output(inxs[inx_type], serial_data_path, delim=" ")

    # persistent storage
    with open(pickle_path, 'wb') as f:
        pickle.dump(inxs, f)

    print("--------------------------------------------------------------------------------")
    print("frame", frame_idx)
    print("----------------------------------------")
    print("size of inxs {}".format(sys.getsizeof(inxs)))
    print("size of system {}".format(sys.getsizeof(system)))
    print("size of assoc {}".format(sys.getsizeof(assoc)))
    if len(inxs[inx_type]) > 0:
        print(len(inxs[inx_type]), "intermolecular hydrogen bonds")
        for inx in inxs[inx_type]:
            inx.pp()
    else:
        print(0, "intermolecular hydrogen bonds")

from mastic.molecule import Molecule
import mastic.interactions as minx
from mastic.interactions import HydrogenBondType
import mastic.selection as masticsel
import pickle
import sys
import os
sys.setrecursionlimit(100000)

print("profile each frame")
print("parameters:")
print("Max distance:", minx.HBOND_DIST_MAX)
print("minimum angle:", minx.HBOND_DON_ANGLE_MIN)

output_dir = osp.join(trypsin_dir, "all_frames_parallel")
# attempt to make it and except if it already exists
try:
    os.mkdir(output_dir)
except OSError:
    pass
frame_idxs = range(ben_coords.shape[0])
inx_type = HydrogenBondType
mol_types = [ben_type, trypsin_type]
sys_type = sys_type
assoc_sel_idxs = [0,1]
assoc_type = rec_lig_type
pickle_paths = [osp.join(output_dir, "all_frames_lig_{}.pkl".format(frame_idx)) for
                frame_idx in frame_idxs]
serial_paths = [osp.join(output_dir, "all_frames_lig_{}.dat".format(frame_idx)) for
                frame_idx in frame_idxs]

from joblib import delayed
# make the generator for the interactions
print("MAKING GENERATOR")
gen = (delayed(frame_profile)(frame_idx, serial_paths[frame_idx], pickle_paths[frame_idx],
                              mol_types, [ben_coords[frame_idx], trypsin_coords[frame_idx]],
                               sys_type, assoc_sel_idxs, assoc_type, inx_type)
         for frame_idx in frame_idxs)

from joblib import Parallel
jobs = 3
print("Running {} jobs".format(jobs))
Parallel(n_jobs=jobs, verbose=5)(gen)
