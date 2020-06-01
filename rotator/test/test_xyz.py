from rotator import *
import numpy as np


def test_xyz():
    rotatexyz("water.xyz", degree="90", axis="x", filename2="water_90x.xyz")
    rotatexyz("water.xyz", degree="90", axis="y", filename2="water_90y.xyz")
    rotatexyz("water.xyz", degree="90", axis="z", filename2="water_90z.xyz")
    rotatexyz("water.xyz", degree="180", axis="x", filename2="water_180x.xyz")
    rotatexyz("water.xyz", degree="180", axis="y", filename2="water_180y.xyz")
    rotatexyz("water.xyz", degree="180", axis="z", filename2="water_180z.xyz")
    rotatexyz("water.xyz", degree="270", axis="x", filename2="water_270x.xyz")
    rotatexyz("water.xyz", degree="270", axis="y", filename2="water_270y.xyz")
    rotatexyz("water.xyz", degree="270", axis="z", filename2="water_270z.xyz")
    # Lets test consistency in the x axis
    mol1 = read_geom("water_90x.xyz")
    geom1 = gen_geom(mol1)
    mol2 = read_geom("water_180x.xyz")
    geom2 = gen_geom(mol2)
    mol3 = read_geom("water_270x.xyz")
    geom3 = gen_geom(mol3)
    xmat90 = s_rot_matrix(degree="90")
    geom1_90x = np.dot(xmat90, geom1)
    geom1_180x = np.dot(xmat90, geom1_90x)
    # Rotating 90 on 90 or 180 is the same
    assert np.allclose(geom1_90x, geom2)
    # Rotating 90 on 90 or 180 is the same
    assert np.allclose(geom1_180x, geom3)
    # Lets test consistency in the y axis
    mol1 = read_geom("water_90y.xyz")
    geom1 = gen_geom(mol1)
    mol2 = read_geom("water_180y.xyz")
    geom2 = gen_geom(mol2)
    mol3 = read_geom("water_270y.xyz")
    geom3 = gen_geom(mol3)
    ymat90 = s_rot_matrix(degree="90", axis="y")
    geom1_90y = np.dot(ymat90, geom1)
    geom1_180y = np.dot(ymat90, geom1_90y)
    # Rotating 90 on 90 or 180 is the same
    assert np.allclose(geom1_90y, geom2)
    # Rotating 90 on 90 or 180 is the same
    assert np.allclose(geom1_180y, geom3)
    # Lets test consistency in the z axis
    mol1 = read_geom("water_90z.xyz")
    geom1 = gen_geom(mol1)
    mol2 = read_geom("water_180z.xyz")
    geom2 = gen_geom(mol2)
    mol3 = read_geom("water_270z.xyz")
    geom3 = gen_geom(mol3)
    zmat90 = s_rot_matrix(degree="90", axis="z")
    geom1_90z = np.dot(zmat90, geom1)
    geom1_180z = np.dot(zmat90, geom1_90z)
    # Rotating 90 on 90 or 180 is the same
    assert np.allclose(geom1_90z, geom2)
    # Rotating 90 on 90 or 180 is the same
    assert np.allclose(geom1_180z, geom3)
