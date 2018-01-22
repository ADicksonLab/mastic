"""Simple example showing how to take a SystemType template object and
creating a System from it (called substantiation, because we are
giving it substance, the coordinates namely). Systems can be used to
calculate geometric properties of a conformation while still
containing all the coordinate independent information in the single
SystemType.

"""
import os.path as osp
import pickle
import sys
import numpy as np

import mastic.profile as masticprof

sys.setrecursionlimit(100000)

# load the system type pickle in
system_pkl_path = osp.join(".", "sEH_TPPU_SystemType.pkl")
with open(system_pkl_path, 'rb') as rf:
    sEH_TPPU_SystemType = pickle.load(rf)

# load the coordinates for the members
member_coords = [np.load('TPPU_coords.npy'), np.load("sEH_coords.npy")]

# substantiate the system
system = sEH_TPPU_SystemType.to_system(member_coords)

with open("sEH_TPPU_System_cryst.pkl", 'wb') as wf:
    pickle.dump(system, wf)
