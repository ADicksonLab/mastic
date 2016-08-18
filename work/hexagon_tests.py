import numpy as np
import math

def rotation_matrix(axis, theta):
    """
    Return the rotation matrix associated with counterclockwise rotation about
    the given axis by theta radians.
    """
    axis = np.asarray(axis)
    theta = np.asarray(theta)
    axis = axis/math.sqrt(np.dot(axis, axis))
    a = math.cos(theta/2.0)
    b, c, d = -axis*math.sin(theta/2.0)
    aa, bb, cc, dd = a*a, b*b, c*c, d*d
    bc, ad, ac, ab, bd, cd = b*c, a*d, a*c, a*b, b*d, c*d
    return np.array([[aa+bb-cc-dd, 2*(bc+ad), 2*(bd-ac)],
                     [2*(bc-ad), aa+cc-bb-dd, 2*(cd+ab)],
                     [2*(bd+ac), 2*(cd-ab), aa+dd-bb-cc]])

axes = np.array([[1.0, 0.0, 0.0],
                 [0.0, 1.0, 0.0],
                 [0.0, 0.0, 1.0]])

origin = np.array([0.0, 0.0, 0.0])
hex_theta = (1.0/3.0)*math.pi
hex_rot_matrices = [rotation_matrix(axes[2], hex_theta*(i+1.0)) for i in range(6)]
unit_point = np.array([1.0, 1.0, 0.0])

unit_hexagon = np.array([np.dot(hex_rot_matrix, unit_point) for hex_rot_matrix in hex_rot_matrices])

def hexagon3d(centroid=[0.0,0.0,0.0], x_theta=0.0, y_theta=0.0, z_theta=0.0, radius=1.0):
    # read in centroid
    centroid = np.asarray(centroid)
    # make the first point that will be rotated
    point = origin + np.array([radius, radius, 0.0])
    # make the points of the hexagon
    hexagon = np.array([np.dot(rot_mat, point) for rot_mat in hex_rot_matrices])
    # translate centroid position
    hexagon = hexagon + centroid
    # rotations
    x_mat = rotation_matrix(axes[0], x_theta)
    hexagon = [np.dot(x_mat, hex_point) for hex_point in hexagon]
    y_mat = rotation_matrix(axes[1], y_theta)
    hexagon = [np.dot(y_mat, hex_point) for hex_point in hexagon]
    z_mat = rotation_matrix(axes[2], z_theta)
    hexagon = [np.dot(z_mat, hex_point) for hex_point in hexagon]

    return hexagon

def pdb_row(i, x, y, z, el='C', color=0.0):
    pdb_row = "ATOM    {0:>3}  {4}   UNK C   1    {1:>8.3f}{2:>8.3f}{3:>8.3f}  1.00 {5:>5.2f}           C  \n".format(i, x, y, z, el, color)
    return pdb_row


end_row = "END                                                                             "
def hexagon_pdb_lines(hexagon, centroid=None):
    pdb_rows = []
    for i, coord in enumerate(hexagon):
        pdb_rows.append(pdb_row(i, coord[0], coord[1], coord[2], color=float(i)+1))
    if centroid:
        pdb_rows.append(pdb_row(len(hexagon)+1, centroid[0], centroid[1], centroid[2], el='H', color=10.0))
    pdb_rows.append(end_row)
    return pdb_rows

def write_hexagon_pdb(hexagon, file_path, centroid=None):
    pdb_lines = hexagon_pdb_lines(hexagon, centroid=centroid)
    with open(file_path, 'w') as wf:
        wf.writelines(pdb_lines)


if __name__ == "__main__":
    with open("tmp/origin.pdb", 'w') as wf:
        wf.write(pdb_row(0, 0.0, 0.0, 0.0, el='S'))
        wf.write(end_row)
    write_hexagon_pdb(unit_hexagon, "tmp/unit_hex.pdb")
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
