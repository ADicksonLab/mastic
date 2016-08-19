import numpy as np
import math

from mast.tests.hexagons import hexagon3d, write_hexagon_pdb, \
    pdb_row, end_row, benzene_bond_length
import mast.config.interactions as mastinxconfig

# unit stuff
with open("tmp/origin.pdb", 'w') as wf:
    wf.write(pdb_row(0, 0.0, 0.0, 0.0, el='S'))
    wf.write(end_row)

write_hexagon_pdb(hexagon3d(), "tmp/benzene_hex.pdb")
stacked_centroid = [0,0,1]
write_hexagon_pdb(hexagon3d(centroid=stacked_centroid),
                  "tmp/stacked_hex.pdb",
                  centroid=stacked_centroid)
write_hexagon_pdb(hexagon3d(centroid=stacked_centroid, z_theta=0.5*math.pi),
                  "tmp/stacked_90.pdb",
                  centroid=stacked_centroid)
write_hexagon_pdb(hexagon3d(x_theta=0.5*math.pi),
                  "tmp/unit_xrot90.pdb")
write_hexagon_pdb(hexagon3d(y_theta=0.5*math.pi),
                  "tmp/unit_yrot90.pdb")
write_hexagon_pdb(hexagon3d(z_theta=0.5*math.pi),
                  "tmp/unit_zrot90.pdb")
write_hexagon_pdb(hexagon3d(x_theta=-0.5*math.pi),
                  "tmp/unit_xrot-90.pdb")
write_hexagon_pdb(hexagon3d(y_theta=-0.5*math.pi),
                  "tmp/unit_yrot-90.pdb")
write_hexagon_pdb(hexagon3d(z_theta=-0.5*math.pi),
                  "tmp/unit_zrot-90.pdb")

# give a parallel example at the maximum centroid distance
pll_max = hexagon3d(centroid=[0.0, 0.0, mastinxconfig.PISTACK_CENTROID_DIST_MAX])
write_hexagon_pdb(pll_max, "tmp/pll_max.pdb")

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
perp_max = hexagon3d(centroid=[0.0, 0.0, mastinxconfig.PISTACK_CENTROID_DIST_MAX],
                     y_theta=math.pi*0.5)
write_hexagon_pdb(perp_max, "tmp/perp_max.pdb")

# give a perpendicular example with closest bisecting atom near just
# under the maximum distance
perp_close_max = hexagon3d(centroid=[0.0, 0.0,
                                     mastinxconfig.PISTACK_CENTROID_DIST_MAX -
                                     benzene_bond_length],
                           y_theta=math.pi*0.5)
write_hexagon_pdb(perp_close_max, "tmp/perp_close_max.pdb")

# give a perpendicular example at maximum angular deviation
perp_close_max_dev = hexagon3d(centroid=[0.0, 0.0,
                                         mastinxconfig.PISTACK_CENTROID_DIST_MAX
                                         - benzene_bond_length],
                               y_theta=math.pi*0.5 +
                               np.radians(mastinxconfig.PISTACK_ANGLE_DEVIATION))

write_hexagon_pdb(perp_close_max_dev, "tmp/perp_close_max_dev.pdb")

# give a perpendicular example for different twist categories
perp_close_max_twist = hexagon3d(centroid=[0.0, 0.0,
                                           mastinxconfig.PISTACK_CENTROID_DIST_MAX -
                                           benzene_bond_length],
                                 y_theta=math.pi*0.5,
                                 z_theta=math.pi*0.5)
write_hexagon_pdb(perp_close_max_twist, "tmp/perp_close_max_twist.pdb")


# example of close but not projected onto the other ring, parallel
pll_displaced = hexagon3d(centroid=[benzene_bond_length,
                                     0.0, 5.0])
write_hexagon_pdb(pll_displaced, "tmp/pll_displaced.pdb")

pll_off_center = hexagon3d(centroid=[benzene_bond_length + mastinxconfig.PISTACK_OFFSET_MAX,
                                     0.0, 5.0])
write_hexagon_pdb(pll_off_center, "tmp/pll_off_center.pdb")
