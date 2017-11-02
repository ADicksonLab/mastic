import os.path as osp

from mastic.tests.hexagons import hexagon3d, write_hexagon_pdb, \
    pdb_row, end_row, benzene_bond_length, benzene_H_bond_length, \
    benzene3d, write_benzene_pdb

work_dir = "/home/salotz/Dropbox/devel/mastic/work/pi_stacking"
write_benzene_pdb(benzene3d(), osp.join(work_dir, "ref_benzene.pdb"))
