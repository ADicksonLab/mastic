import os.path as osp
import pickle
import time

import numpy as np
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


start = time.time()
system = seh_tppu_System.to_system(member_coords)
end = time.time()

print("Time for SystemType.to_system was {} seconds".format(end-start))

