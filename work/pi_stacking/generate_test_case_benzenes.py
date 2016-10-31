import os.path as osp
import numpy as np
import math

import mast.config.interactions as mastinxconfig

from mast.tests.hexagons import benzene3d, write_benzene_pdb, \
    pdb_row, end_row, benzene_bond_length

test_benzenes = {}

# test directly stacked benzenes
stacked_centroid = [0.0, 0.0, 3.0]
max_stacked_centroid = [0.0, 0.0, mastinxconfig.PISTACKING_CENTROID_DIST_MAX]
test_benzenes['stacked'] = benzene3d(centroid=stacked_centroid)
# with a 90 degree turn
test_benzenes['stacked_90'] = benzene3d(centroid=stacked_centroid, z_theta=0.5*math.pi)

# give a parallel example at the maximum centroid distance
test_benzenes['pll_max'] = benzene3d(centroid=max_stacked_centroid)

# give a parallel example at maximum angular deviation
test_benzenes['pll_ang'] = benzene3d(centroid=max_stacked_centroid,
                    x_theta=np.radians(mastinxconfig.PISTACKING_ANGLE_DEVIATION))

# # give a parallel example for different yaw categories
# # give an optimal parallel displaced example
# test_benzenes['pll_x_disp'] = benzene3d(centroid=[1.0, 0.0, mastinxconfig.PISTACKING_CENTROID_DIST_MAX])
# test_benzenes['pll_y_disp'] = benzene3d(centroid=[0.0, 1.0, mastinxconfig.PISTACKING_CENTROID_DIST_MAX])
# test_benzenes['pll_xy_disp'] = benzene3d(centroid=[1.0, 1.0, mastinxconfig.PISTACKING_CENTROID_DIST_MAX])

# give a perpendicular example with centroid at maximum distance
test_benzenes['perp_max'] = benzene3d(centroid=[0.0, 0.0, mastinxconfig.PISTACKING_CENTROID_DIST_MAX],
                     y_theta=math.pi*0.5)

# give a perpendicular example with closest bisecting atom near just
# under the maximum distance
test_benzenes['perp_close_max'] = benzene3d(centroid=[0.0, 0.0,
                                     mastinxconfig.PISTACKING_CENTROID_DIST_MAX -
                                     benzene_bond_length],
                           y_theta=math.pi*0.5)

# give a perpendicular example at maximum angular deviation
test_benzenes['perp_close_max_dev'] = benzene3d(centroid=[0.0, 0.0,
                                         mastinxconfig.PISTACKING_CENTROID_DIST_MAX
                                         - benzene_bond_length],
                               y_theta=math.pi*0.5 +
                               np.radians(mastinxconfig.PISTACKING_ANGLE_DEVIATION))

# give a perpendicular example for different twist categories
test_benzenes['perp_close_max_twist'] = benzene3d(centroid=[0.0, 0.0,
                                           mastinxconfig.PISTACKING_CENTROID_DIST_MAX -
                                           benzene_bond_length],
                                                  y_theta=math.pi*0.5,
                                                  z_theta=math.pi*0.5)

# # example of close but not projected onto the other benzene, parallel
# test_benzenes['pll_displaced'] = benzene3d(centroid=[benzene_bond_length,
#                                                      0.0, 5.0])

# test_benzenes['pll_off_center'] = benzene3d(centroid=[benzene_bond_length +
#                                                       mastinxconfig.PISTACKING_OFFSET_MAX,
#                                                       0.0, 5.0])

for test_name, test_benzene in test_benzenes.items():
    file_name = "{}.pdb".format(test_name)
    write_benzene_pdb(test_benzene, osp.join(work_dir, file_name))
