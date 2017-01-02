import os.path as osp
import pickle

from mast.interactions.hydrogen_bond import HydrogenBondType

system_cryst_pkl_path = osp.join(".", "Trypsin_Benzamidine_System_cryst.pkl")
with open(system_cryst_pkl_path, 'rb') as pkl_rf:
    Trypsin_Benzamidine_System_cryst = pickle.load(pkl_rf)

tryp_ben_prof_results_0 = Trypsin_Benzamidine_System_cryst.associations[0].profile_interactions([HydrogenBondType])
tryp_ben_prof_results_1 = Trypsin_Benzamidine_System_cryst.associations[1].profile_interactions([HydrogenBondType])
