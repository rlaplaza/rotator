from rotator import *
import numpy as np


def test_flower():
   mol1 = read_geom("water.xyz")
   natoms = 4 * mol1.atnums.size  # IOData is very handy!
   geom1 = gen_geom(mol1)
   print(geom1.T)
   mat = s_rot_matrix(degree="90", axis="x")
   geom2 = np.dot(mat, geom1)
   print(geom2.T)
   geom3 = np.dot(mat, geom2)
   print(geom3.T)
   geom4 = np.dot(mat, geom3)
   print(geom4.T)
   geom2 = g_displace(geom2, vec=np.asarray([0, 0.5, 0.5]))
   print(geom2.T)
   geom3 = s_displace(geom3, axis="z", norm=1.0)
   print(geom3.T)
   geom4 = g_displace(geom4, vec=np.asarray([0, -0.5, 0.5]))
   print(geom4.T)
   merge = (geom1.T, geom2.T, geom3.T, geom4.T)
   geomflower = np.concatenate(merge, axis=0)
   # print(geomflower)
   f = open("flower.xyz", "w")
   f.write("" + str(natoms) + "\n")
   f.write("A beautiful flower of waters\n")
   for i in range(natoms):
       f.write("C")
       for j in range(0, 3):
           a = np.format_float_positional(geomflower[i, j], precision=4)
           f.write("   " + str(a))
       f.write("\n")
   f.close()
