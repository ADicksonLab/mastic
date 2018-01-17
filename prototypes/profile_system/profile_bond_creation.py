import os.path as osp
import pickle
import time

import numpy as np

from mastic.selection import CoordArray
from mastic.molecule import Atom, Bond

# must set this to pickle substantiated systems sometimes
#sys.setrecursionlimit(100000)

# load the system type pickle in
inputs_path = osp.realpath("../../examples/sEH-TPPU")
system_pkl_path = osp.join(inputs_path, "sEH_TPPU_SystemType.pkl")
with open(system_pkl_path, 'rb') as rf:
    seh_tppu_System = pickle.load(rf)

# load the rdkit wrapper pickles and the coordinates
TPPU_rdkit_pkl_path = osp.join(".", "TPPU_rdkit.pkl")
seh_rdkit_pkl_path = osp.join(".", "sEH_rdkit.pkl")
TPPU_coords_path = osp.join(".", "TPPU.npy")
seh_coords_path = osp.join(".", "sEH.npy")
with open(TPPU_rdkit_pkl_path, 'rb') as rf:
    TPPU_rdkit_wrapper = pickle.load(rf)
with open(seh_rdkit_pkl_path, 'rb') as rf:
    seh_rdkit_wrapper = pickle.load(rf)
TPPU_coords = np.load(TPPU_coords_path)
seh_coords = np.load(seh_coords_path)
member_coords = [TPPU_coords, seh_coords]

seh_type = seh_tppu_System.member_types[1]

# make atoms for the whole molecule
coord_array = CoordArray(seh_coords)
atoms = []
for atom_idx, atom_type in enumerate(seh_type.atom_types):
    atom = Atom(atom_array=coord_array, array_idx=atom_idx, atom_type=atom_type)
    atoms.append(atom)

bond_idx = 0
bond_type = seh_type.bond_types[bond_idx]
atom_ids = seh_type.bond_map[bond_idx]

print("Making the Bond")
start = time.time()
bond = Bond(atom_container=atoms, atom_ids=atom_ids,
                        bond_type=bond_type)
end = time.time()

print("Time for creating an sEH bond was {} seconds".format(end-start))
