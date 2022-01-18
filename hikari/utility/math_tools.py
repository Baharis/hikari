"""
This file contains constants and functions containing purely mathematical
information / methods used in the package.
"""

import numpy as np
import random


def angle2rad(value):
    """
    Interpret a unit of given angle value and return this value in radians.
    Values in range between -3.15 and 3.15 are interpreted as given in radians.
    Values outside this range are interpreted as given in degrees.

    The function can be used to dynamically interpret and accept
    both radian and degrees and return the value in angles. Should you not
    need this functionality, please use numpy's rad2deg and deg2rad instead.

    :Example:

    >>> angle2rad(np.pi)
    3.141592653589793
    >>> angle2rad(180)
    3.141592653589793
    >>> angle2rad(2*np.pi)
    0.1096622711232151

    Please mind that in the last case the value given in radians was treated
    as given in degrees and wrongly recalculated to radians again.

    :param value: Angle value given in angles or radians
    :type value: float
    :return: Angle value expressed in radians
    :rtype: float
    """
    return value if -3.15 < value < 3.15 else np.deg2rad(value)


def fibonacci_sphere(samples=1, seed=1337):
    """
    Return a 3D cartesian coordinates of *samples* number of points
    evenly distributed on a surface of a unit sphere at (0, 0, 0).

    The algorithm utilised in this function gives a uniform distribution on the
    sphere. Contrary to an uniform distribution on coordinates, this algorithm
    does not favour points close to poles. For more details see `this article
    <http://extremelearning.com.au/evenly-distributing-points-on-a-sphere/>`_.

    The function is written so that it does always return the same set of points
    in order to ensure reproducibility of other methods. A list of other points
    can be obtained by changing (or even randomising) the value of a *seed*.

    :Example:

    >>> fibonacci_sphere(4)
    array([[ 0.62069, -0.75, -0.22857],
           [-0.44395, -0.25,  0.86047],
           [ 0.41275,  0.25,  0.87587],
           [-0.61208,  0.75, -0.25072]])
    >>> fibonacci_sphere(4)
    array([[ 0.62069, -0.75, -0.22857],
           [-0.44395, -0.25,  0.86047],
           [ 0.41275,  0.25,  0.87587],
           [-0.61208,  0.75, -0.25072]])
    >>> fibonacci_sphere(4, seed=420)
    array([[ 0.64040, -0.75,  0.16550],
           [-0.85489, -0.25,  0.45460],
           [ 0.32329,  0.25, -0.91268],
           [ 0.25831,  0.75,  0.60892]])

    :param samples: Number of points to be generated.
    :type samples: int
    :param seed: A seed value used once to slightly randomise point position.
    :type seed: any
    :return: (x,y,z) arrays representing points' cartesian coordinates.
    :rtype: np.array
    """
    random.seed(seed)
    rnd = random.random() * samples
    offset = 2. / samples
    increment = np.pi * (3. - np.sqrt(5.))
    i = np.arange(0, samples)
    y = ((i * offset) - 1) + (offset / 2)
    r = np.sqrt(1 - pow(y, 2))
    phi = ((i + rnd) % samples) * increment
    x = np.cos(phi) * r
    z = np.sin(phi) * r
    return np.vstack([x, y, z]).T


def euler_rodrigues_matrix(a, b, c, d):
    aa, bb, cc, dd = a * a, b * b, c * c, d * d
    bc, ad, ac, ab, bd, cd = b * c, a * d, a * c, a * b, b * d, c * d
    return np.array([[aa + bb - cc - dd, 2 * (bc + ad),     2 * (bd - ac)],
                     [2 * (bc - ad),     aa + cc - bb - dd, 2 * (cd + ab)],
                     [2 * (bd + ac),     2 * (cd - ab),     aa + dd - bb - cc]])


def rotation_from(from_, to):
    """
    Return matrix associated with rotation of vector `from_` onto `to`.
    :param from_: A 3D vector which will have been rotated onto `to`.
    :type from_: np.ndarray
    :param to: A 3D vector onto which `from` will have been rotated.
    :type to: np.ndarray
    :return: Matrix associated with rotation of `from` onto `to.
    :rtype: np.ndarray
    """
    if from_.size != 3:
        raise IndexError(f'Size of `from_` should be 3, but is {from_.size}')
    if to.size != 3:
        raise IndexError(f'Size of `to` should be 3, but is {to.size}')
    f = from_ / np.linalg.norm(from_)
    t = to / np.linalg.norm(to)
    axis = np.cross(f, t)
    if np.isclose(sum(axis), 0.):
        if np.allclose(f, t):
            return np.eye(3)
        else:
            e1, e2 = np.array([1, 0, 0]), np.array([0, 1, 0])
            e = e2 if np.isclose(sum(np.cross(e1, f)), 0.) else e1
            return rotation_from(f, to=e) @ rotation_from(e, to=t)
    else:
        axis = axis / np.linalg.norm(axis)
        angle = np.arcsin(np.linalg.norm(np.cross(f, t)))
        return euler_rodrigues_matrix(np.cos(angle / 2.),
                                      *(-axis * np.sin(angle / 2.)))


def rotation_around(axis, by):
    """
    Return matrix associated with counterclockwise rotation about `axis` through
    `by` radians. See Euler-Rodrigues form.,https://stackoverflow.com/q/6802577/
    :param axis: A 3D vector about which rotation is performed
    :type axis: np.ndarray
    :param by: Angle of rotation in radians
    :type by: float
    :return: Matrix associated with rotation around `axis` through `by`.
    :rtype: np.ndarray
    """
    if axis.size != 3:
        raise IndexError(f'Size of `axis` should be 3, but is {axis.size}')
    axis = axis / np.linalg.norm(axis)
    return euler_rodrigues_matrix(np.cos(by / 2.0), *(-axis * np.sin(by / 2.0)))


# TODO later check scipy.spatial.transform.Rotation.from_rotvec as substitute
