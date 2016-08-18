import numpy as np
import math

from mast.tests.hexagons import hexagon3d, write_hexagon_pdb
import mast.config.interactions as mastinxconfig

# give a parallel example at the maximum centroid distance
pll_max = hexagon3d(centroid=[0.0, 0.0, mastinxconfig.PISTACK_CENTROID_DIST_MAX])
write_hexagon_pdb(pll_max, "tmp/pll_max.pdb")

# give a parallel example but translate along the plane to end of
# centroid distance

# give a parallel example at maximum angular deviation
pll_ang = hexagon3d(centroid=[0.0, 0.0, mastinxconfig.PISTACK_CENTROID_DIST_MAX],
                    x_theta=np.radians(mastinxconfig.PISTACK_ANGLE_DEVIATION))
write_hexagon_pdb(pll_ang, "tmp/pll_ang.pdb")
# give a parallel example for different yaw categories
# give an optimal parallel displaced example
pll_x_disp = hexagon3d(centroid=[1.0, 0.0, mastinxconfig.PISTACK_CENTROID_DIST_MAX])
write_hexagon_pdb(pll_x_disp, "tmp/pll_x_disp.pdb")
pll_y_disp = hexagon3d(centroid=[0.0, 1.0, mastinxconfig.PISTACK_CENTROID_DIST_MAX])
write_hexagon_pdb(pll_y_disp, "tmp/pll_y_disp.pdb")
pll_xy_disp = hexagon3d(centroid=[1.0, 1.0, mastinxconfig.PISTACK_CENTROID_DIST_MAX])
write_hexagon_pdb(pll_xy_disp, "tmp/pll_xy_disp.pdb")

# give a perpendicular example with centroid at maximum distance

# give a perpendicular example with closest bisecting atom at maximum distance

# give a perpendicular example at maximum angular deviation

# give a perpendicular example for different twist categories
