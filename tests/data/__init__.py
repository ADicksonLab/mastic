from io import StringIO
import sys
import os.path as osp
import pickle

# get the path to this module directory
data_dir = osp.dirname(sys.modules[__name__].__file__)

BEN_file_name = "3ptb_BEN.pdb"
trypsin_file_name = "3ptb_trypsin.pdb"
waters_3ptb_file_name = "3ptb_water.pdb"

BEN_Hs_file_name = "BEN+Hs_3ptb.pdb"
trypsin_Hs_file_name = "trypsin+Hs_3ptb.pdb"

benzamidine_file_name = "benzamidine.mol"

# sEH_Hs_file_name = "sEH_Hs.pdb"
# TPPU_Hs_file_name = "TPPU_Hs.pdb"

Top7_file_name = "Top7_1qys.pdb"
chignolin_file_name = "chignolin_5awl.pdb"



ligand_structure_files = [BEN_file_name, BEN_Hs_file_name]
protein_structure_files = [trypsin_file_name, trypsin_Hs_file_name,
                           Top7_file_name, chignolin_file_name]
solvent_files = [waters_3ptb_file_name]

structure_files = ligand_structure_files + protein_structure_files + solvent_files

# load each data example as a string
BEN_path = osp.join(data_dir, BEN_file_name)
with open(BEN_path, 'r') as rf:
    BEN_3ptb = rf.read()

BEN_Hs_path = osp.join(data_dir, BEN_Hs_file_name)
with open(BEN_Hs_path, 'r') as rf:
    BEN_Hs_3ptb = rf.read()

benzamidine_MOL_path = osp.join(data_dir, benzamidine_file_name)
with open(benzamidine_MOL_path, 'r') as rf:
    benzamidine_MOL = rf.read()

trypsin_path = osp.join(data_dir, trypsin_file_name)
with open(trypsin_path, 'r') as rf:
    trypsin_3ptb = rf.read()

trypsin_Hs_path = osp.join(data_dir, trypsin_Hs_file_name)
with open(trypsin_Hs_path, 'r') as rf:
    trypsin_Hs_3ptb = rf.read()


waters_3ptb_path = osp.join(data_dir, waters_3ptb_file_name)
with open(waters_3ptb_path, 'r') as rf:
    waters_3ptb = rf.read()


# sEH_Hs_path = osp.join(data_dir, sEH_Hs_file_name)
# with open(sEH_Hs_path, 'r') as rf:
#     sEH_Hs_3ptb = rf.read()

# TPPU_Hs_path = osp.join(data_dir, TPPU_Hs_file_name)
# with open(TPPU_Hs_path, 'r') as rf:
#     TPPU_Hs_3ptb = rf.read()

Top7_path = osp.join(data_dir, Top7_file_name)
with open(Top7_path, 'r') as rf:
    Top7_1qys = rf.read()

chignolin_path = osp.join(data_dir, chignolin_file_name)
with open(chignolin_path, 'r') as rf:
    chignolin_5awl = rf.read()


# precomputed molecules with features already detected

# from rdkit feature detection
trypsin_mastmol_path = osp.join(data_dir, "trypsin+features_mastmol.pkl")
with open(trypsin_mastmol_path, 'rb') as pkl_rf:
    Trypsin_Molecule = pickle.load(pkl_rf)

# BEN_mastmol_path = osp.join(data_dir, "BEN+features_mastmol.pkl")
# with open(BEN_mastmol_path, 'rb') as pkl_rf:
#     BEN_Molecule = pickle.load(pkl_rf)

trypsin_Hs_mastmol_path = osp.join(data_dir, "trypsin+Hs+features_mastmol.pkl")
with open(trypsin_mastmol_path, 'rb') as pkl_rf:
    Trypsin_Hs_Molecule = pickle.load(pkl_rf)

# BEN_Hs_mastmol_path = osp.join(data_dir, "BEN+Hs+features_mastmol.pkl")
# with open(BEN_Hs_mastmol_path, 'rb') as pkl_rf:
#     BEN_Hs_Molecule = pickle.load(pkl_rf)

# sEH_Hs_mastmol_path = osp.join(data_dir, "sEH+Hs+features_mastmol.pkl")
# with open(sEH_Hs_mastmol_path, 'rb') as pkl_rf:
#     sEH_Hs_Molecule = pickle.load(pkl_rf)

# TPPU_Hs_mastmol_path = osp.join(data_dir, "TPPU+Hs+features_mastmol.pkl")
# with open(TPPU_Hs_mastmol_path, 'rb') as pkl_rf:
#     TPPU_Hs_Molecule = pickle.load(pkl_rf)
