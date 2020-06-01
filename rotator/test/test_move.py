from rotator import *
import numpy as np


def test_move():
    mol1 = read_geom("water.xyz")
    geom1 = gen_geom(mol1)
    geom2 = g_displace(geom1, vec=np.asarray([1, 1, 1]))
    geom3 = g_displace(geom2, vec=np.asarray([-1, -1, -1]))
    assert np.allclose(geom1, geom3)
