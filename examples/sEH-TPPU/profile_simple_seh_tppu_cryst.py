import os.path as osp
import pickle

from mast.interactions.hydrogen_bond import HydrogenBondType

system_cryst_pkl_path = osp.join(".", "sEH_TPPU_System_cryst.pkl")
with open(system_cryst_pkl_path, 'rb') as pkl_rf:
    seh_tppu_System_cryst = pickle.load(pkl_rf)

tryp_tppu_prof_results_0 = seh_tppu_System_cryst.associations[0].test_profile_interactions([HydrogenBondType])
tryp_tppu_prof_results_1 = seh_tppu_System_cryst.associations[1].test_profile_interactions([HydrogenBondType])
