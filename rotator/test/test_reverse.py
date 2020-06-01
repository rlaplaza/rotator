from rotator import *
import numpy as np


def test_reverse():
    mol1 = read_geom("water_opt.fchk")
    coords1 = gen_geom(mol1, verb_lvl=3)
    mol2 = put_geom(mol1, coords1, verb_lvl=3)
    coords2 = gen_geom(mol2, verb_lvl=3)
    assert np.allclose(coords1, coords2)
    rotmat = g_rot_matrix(verb_lvl=3)
    coords3 = np.dot(coords1, rotmat)
    mol2 = put_geom(mol1, coords3, verb_lvl=3)
    coords2 = gen_geom(mol2, verb_lvl=3)
    assert np.allclose(coords1, coords2)
