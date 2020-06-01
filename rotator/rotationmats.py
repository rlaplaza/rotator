import numpy as np
from numpy import linalg as la
from iodata import load_one, dump_one, IOData
import math
import os
import pprint

pp = pprint.PrettyPrinter(indent=4)


class rotatorerror(Exception):
    """ Exception class for errors in the rotator module.
        Usually means something is very weird.
    """

    pass


def g_rot_matrix(degree="0.0", axis=np.asarray([1, 1, 1]), verb_lvl=0):
    """
       Return the rotation matrix associated with counterclockwise rotation about
       the given axis by theta radians.
    Parameters
    ----------
    axis
        Axis of the rotation. Array or list.
    degree
        Angle of the rotation in degrees.
    verb_lvl
        Verbosity level integer flag.

    Returns
    -------
    rot 
        Rotation matrix to be used.
    """
    try:
        theta = degree * (np.pi / 180)
    except:
        degree = float(degree)
        theta = degree * (np.pi / 180)
    axis = np.asarray(axis)
    axis = axis / math.sqrt(np.dot(axis, axis))
    a = math.cos(theta / 2.0)
    b, c, d = -axis * math.sin(theta / 2.0)
    aa, bb, cc, dd = a * a, b * b, c * c, d * d
    bc, ad, ac, ab, bd, cd = b * c, a * d, a * c, a * b, b * d, c * d
    rot = np.array(
        [
            [aa + bb - cc - dd, 2 * (bc + ad), 2 * (bd - ac)],
            [2 * (bc - ad), aa + cc - bb - dd, 2 * (cd + ab)],
            [2 * (bd + ac), 2 * (cd - ab), aa + dd - bb - cc],
        ]
    )
    if verb_lvl > 1:
        print("Rotation matrix generated.")
    if verb_lvl > 2:
        pp.pprint(rot)
    return rot


def s_rot_matrix(degree="0.0", axis="x", verb_lvl=0):
    """
       Return the rotation matrix associated with counterclockwise rotation about
       the given axis by theta radians.
    Parameters
    ----------
    axis
        Axis of the rotation, string x, y or z. 
    theta
        Angle of the rotation in degrees.
    verb_lvl
        Verbosity level integer flag.

    Returns
    -------
    rot 
        Rotation matrix to be used.
    """
    if not isinstance(axis, str):
        raise rotatorerror("This function takes axis=x/y/z only.")
    try:
        theta = degree * (np.pi / 180)
        c, s = math.cos(theta), math.sin(theta)
    except:
        degree = float(degree)
        theta = degree * (np.pi / 180)
        c, s = math.cos(theta), math.sin(theta)
    if axis == "x":
        rot = np.array([[1.0, 0, 0], [0, c, -s], [0, s, c]])
    elif axis == "y":
        rot = np.array([[c, 0, s], [0, 1.0, 0], [-s, 0, c]])
    elif axis == "z":
        rot = np.array([[c, -s, 0], [s, c, 0], [0, 0, 1.0]])
    else:
        raise rotatorerror("This function takes axis=x/y/z only.")
    if verb_lvl > 1:
        print("Rotation matrix generated.")
    if verb_lvl > 2:
        pp.pprint(rot)
    return rot


def g_displace(coordmat, vec=np.asarray([1, 1, 1]), verb_lvl=0, norm=None):
    """Displace the geometry matrix following a displacement vector. 
     Alternatively, it can take the direction from any vector and renormalize to a norm.
    
     Parameters
     ----------
    coordmat
        Geometry matrix.
    norm
        The norm in angstrom. Optional, can be given in in the vector.
    vec
        The vector of the direction of the displacement. Array.
    verb_lvl
        Verbosity level integer flag.

    Returns
    -------
    coordmat 
        Displaced geometry matrix.
   """
    vec = np.asarray(vec)
    if norm is not None:
        try:
            norm = float(norm)
        except:
            raise rotatorerror("This function needs a real or integer norm.")
        vec = vec / (np.la.norm(vec) + 1e-16)
        vec = vec * norm
    vec.shape = (3,1)
    if verb_lvl > 2:
        pp.pprint(vec)
    coordmat += vec
    return coordmat


def s_displace(coordmat, axis="x", norm=1, verb_lvl=0):
    """Displace the geometry matrix following a displacement vector x/y/z 
       using a norm in angstroms.
    
     Parameters
     ----------
    coordmat
        Geometry matrix.
    norm
        The norm in angstrom. Optional, can be given in in the vector.
    vec
        The vector of the direction of the displacement, string x, y or z. 
    verb_lvl
        Verbosity level integer flag.

    Returns
    -------
    coordmat 
        Displaced geometry matrix.
   """
    if axis == "x":
        vec = np.array([1, 0, 0])
    elif axis == "y":
        vec = np.array([0, 1, 0])
    elif axis == "z":
        vec = np.array([0, 0, 1])
    else:
        raise rotatorerror("This function takes axis=x/y/z only.")
    try:
        norm = float(norm)
    except:
        raise rotatorerror("This function needs a real or integer norm.")
    vec = vec * norm
    vec.shape = (3,1)
    if verb_lvl > 2:
        pp.pprint(vec)
    coordmat += vec
    return coordmat


def gen_geom(mol, verb_lvl=0):
    """
       Return the geometry matrix from the molecule object.
    Parameters
    ----------
    mol
        IOData molecule object.
    verb_lvl
        Verbosity level integer flag.

    Returns
    -------
    coordmat
        Geometry matrix.
    """
    if not isinstance(mol, IOData):
        raise rotatorerror("Something other than an IOData mol object passed.")
    coords = []
    coordmat = np.empty(shape=[3, 3])
    # coordmat = mol.atcoords.T*0.52917721092  # that is it, this works
    coords = np.asarray(
        [
            [mol.atnums[i], mol.atcoords[i] * 0.52917721092]
            for i in range(mol.atnums.size)
        ]
    )
    for i in range(
        0, 3
    ):  # Its perfectly possible to simply transpose mol.atcoords, this is for hookability
        for j in range(mol.atnums.size):
            coordmat[i, j] = coords[j, 1][i]
    if verb_lvl > 1:
        print("Geometry matrix generated.")
    if verb_lvl > 2:
        pp.pprint(coordmat)
        pp.pprint(mol.atcoords)
    return coordmat


def put_geom(mol, coordmat, verb_lvl=0):
    """
       Put a new geometry matrix into the molecule object.
    Parameters
    ----------
    mol
        IOData molecule object.
    coordmat
        Geometry matrix.
    verb_lvl
        Verbosity level integer flag.

    Returns
    -------
    mol
        IOData molecule object with new geometry.

    """
    if not isinstance(mol, IOData):
        raise rotatorerror("Something other than an IOData mol object passed.")
    for i in range(mol.atnums.size):
        for j in range(
            0, 3
        ):  # Its perfectly possible to simply coordmat to mol.atcoords, this is for hookability
            mol.atcoords[i][j] = coordmat[j, i] * 1.8897259886
    if verb_lvl > 1:
        print("Geometry matrix updated.")
    if verb_lvl > 2:
        pp.pprint(coordmat)
        pp.pprint(mol.atcoords)
    return mol
