from io import StringIO
import sys
import os.path as osp
import pickle

# get the path to this module directory
data_dir = osp.dirname(sys.modules[__name__].__file__)

BEN_file_name = "BEN_3ptb.pdb"
trypsin_file_name = "trypsin_3ptb.pdb"
Top7_file_name = "Top7_1qys.pdb"
chignolin_file_name = "chignolin_5awl.pdb"

ligand_structure_files = [BEN_file_name]
protein_structure_files = [trypsin_file_name, Top7_file_name, chignolin_file_name]
structure_files = ligand_structure_files + protein_structure_files

# load each data example as a string
BEN_path = osp.join(data_dir, BEN_file_name)
with open(BEN_path, 'r') as rf:
    BEN_3ptb = rf.read()

trypsin_path = osp.join(data_dir, trypsin_file_name)
with open(trypsin_path, 'r') as rf:
    trypsin_3ptb = rf.read()

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

BEN_mastmol_path = osp.join(data_dir, "BEN+features_mastmol.pkl")
with open(BEN_mastmol_path, 'rb') as pkl_rf:
    BEN_Molecule = pickle.load(pkl_rf)
