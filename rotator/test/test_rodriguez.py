from rotator import *
import numpy as np


def test_rodriguez():
    mol1 = read_geom("water.xyz")
    geom1 = gen_geom(mol1)
    mat90x1 = s_rot_matrix(degree=90, axis="x")
    mat90x2 = g_rot_matrix(degree=90, axis=[1, 0, 0])
    assert np.allclose(mat90x1, mat90x2)
    geom1_90x = np.dot(mat90x1, geom1)
    geom2_90x = np.dot(mat90x2, geom1)
    assert np.allclose(geom1_90x, geom2_90x)
