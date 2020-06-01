import numpy as np
from numpy import linalg as la
from iodata import load_one, dump_one, IOData
import math
import os
import pprint
from rotator import *

pp = pprint.PrettyPrinter(indent=4)


class writererror(Exception):
    """ Exception class for errors in the reader module.
    """

    pass


def write_geom(mol, filename: str, verb_lvl=0):
    """Very simple wrapper for iodata dump_one.

    Parameters
    ----------
    filename
        A string that contains the path to an output file.
    verb
        Verbosity level integer flag.
    mol
        An IOdata molecule object.

    """
    dump_one(mol, filename)


def rotatexyz(filename1: str, degree="0", axis="x", verb_lvl=0, filename2="output.xyz"):
    """Read an xyz file, rotate it some degrees around some axis and write it.
    
     Parameters
     ----------
     filename1
        A string that contains the path to an input xyz file.
     filename2
        A string that contains the path to the output xyz file. By default will be called out.xyz 
     degree
        Angle of the rotation in degrees.
     axis
        Axis of the rotation. Can be a string x/y/z to use those axis or a vector.
     verb_lvl
        Verbosity level integer flag.
   """
    mol = read_geom(filename1, verb=verb_lvl)
    geom = gen_geom(mol, verb_lvl)
    if isinstance(axis, str):
        mat = s_rot_matrix(degree, axis, verb_lvl)
    else:
        mat = g_rot_matrix(degree, axis, verb_lvl)
    newgeom = np.dot(mat, geom)
    newmol = put_geom(mol, newgeom, verb_lvl)
    write_geom(newmol, filename2, verb_lvl)
