from rotator import *
import numpy as np


def test_fchk_molden():
    mol1 = read_geom("water_opt.fchk")
    mol2 = read_geom("water_opt.molden")
    coords1 = gen_geom(mol1, verb_lvl=3)
    coords2 = gen_geom(mol2, verb_lvl=3)
    assert np.allclose(coords1, coords2)
