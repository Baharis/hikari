"""
This file contains constants and functions containing purely mathematical
information / methods used in the package.
"""

import numpy as np
import random
import uncertainties
from typing import Union


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
    return value if -3.15 < value < 3.15 else value * np.pi / 180


def cart2sph(x, y, z):
    """
    Convert Cartesian coordinates x, y, z
    to conventional spherical coordinates r, p, a

    :param x: Cartesian coordinate or vector x
    :type x: float or np.ndarray
    :param y: Cartesian coordinate or vector y
    :type y: float or np.ndarray
    :param z: Cartesian coordinates or vector z
    :type z: float or np.ndarray
    :return: Spherical coordinates: radius, polar angle, and azimuth angle
    :rtype: np.ndarray
    """
    r = (x ** 2 + y ** 2 + z ** 2) ** 0.5
    p = np.arccos(z / r)
    a = np.arctan2(y, x)
    return np.array([r, p, a])


def sph2cart(r, p, a):
    """
    Convert conventional spherical coordinates r, p, a
    to Cartesian coordinates x, y, z

    :param r: Spherical coordinate or vector radius
    :type r: float or np.ndarray
    :param p: Spherical coordinate or vector polar angle
    :type p: float or np.ndarray
    :param a: Spherical coordinate or vector azimuth angle
    :type a: float or np.ndarray
    :return: Cartesian coordinates: x, y, and z
    :rtype: np.ndarray
    """
    x = r * np.cos(a) * np.sin(p)
    y = r * np.sin(a) * np.sin(p)
    z = r * np.cos(p)
    return np.array([x, y, z])


def det3x3(matrix):
    """
    Calculate a determinant of a 3x3 matrix. Should be usually substituted
    by `numpy.linalg.det`, but is indispensable for matrices with uncertainties.

    :param matrix: 3x3 array/matrix which allows for 2d indexing
    :type matrix: numpy.ndarray or uncertainties.unumpy.matrix
    :return: Determinant of the matrix
    :rtype: int or float or uncertainties.core.Variable
    """
    m11, m12, m13 = matrix[0, 0], matrix[1, 0], matrix[2, 0]
    m21, m22, m23 = matrix[0, 1], matrix[1, 1], matrix[2, 1]
    m31, m32, m33 = matrix[0, 2], matrix[1, 2], matrix[2, 2]
    return m11 * m22 * m33 + m12 * m23 * m31 + m13 * m21 * m32 \
        - m13 * m22 * m31 - m12 * m21 * m33 - m11 * m23 * m32


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


def euler_rodrigues_matrix(
        a: Union[int, float, uncertainties.UFloat],
        b: Union[int, float, uncertainties.UFloat],
        c: Union[int, float, uncertainties.UFloat],
        d: Union[int, float, uncertainties.UFloat]) -> np.ndarray:
    """
    Return a rotation matrix based on a Euler-Rodrigues parametrisation. For
    details, see https://en.wikipedia.org/wiki/Euler%E2%80%93Rodrigues_formula.

    :param a: Euler-Rodrigues parameter *a*
    :param b: Euler-Rodrigues parameter *b*
    :param c: Euler-Rodrigues parameter *c*
    :param d: Euler-Rodrigues parameter *d*
    :return: 3x3 matrix describing rotation
    """
    aa, bb, cc, dd = a * a, b * b, c * c, d * d
    bc, ad, ac, ab, bd, cd = b * c, a * d, a * c, a * b, b * d, c * d
    return np.array([[aa + bb - cc - dd, 2 * (bc + ad), 2 * (bd - ac)],
                     [2 * (bc - ad), aa + cc - bb - dd, 2 * (cd + ab)],
                     [2 * (bd + ac), 2 * (cd - ab), aa + dd - bb - cc]])


def rotation_from(from_, to):
    """
    Return matrix associated with rotation of vector `from_` onto `to`.

    :param from_: A 3D vector which will have been rotated onto `to`.
    :type from_: np.ndarray
    :param to: A 3D vector onto which `from` will have been rotated.
    :type to: np.ndarray
    :return: Matrix associated with rotation of `from` onto `to`.
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
            e = np.cross(e2, f) if np.isclose(sum(np.cross(e1, f)), 0.) \
                else np.cross(e1, f)
            return rotation_around(e, by=np.pi)
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


def weighted_quantile(values, quantiles, weights=None):
    """
    A function to approximate quantiles based on many weighted points. Quantiles
    are derived exactly using method 7 of H&F from an unweighted distribution,
    otherwise they are approximated by comparing two opposing interpolations.

    Let A: {1, 2, 2, 2, 3} and B: {1, 2(x3), 3} be seemingly equivalent sets
    where (x3) denotes triple weight of the "2". Let q[n] denote quantiles.
    For set A, q[1/4](A) = q[3/4](A) = 2. However, for set B we need to
    extrapolate and read quantile from PDF, obtained via linear interpolation.
    As a result, q[1/4](B) = 4/3 < q[3/4](B) = 2 since f(0) = 1 and f(3/4) = 2
    at f(1/4) and the results are always underestimated.

    To combat this effect, we can repeat the process for -B: {-3, -2(x3), -1}.
    Here, q[1/4](-B) = -8/3 and q[3/4](-B) = -2. As a result, the quantile
    q[n](A) is described a bit better using formula (q[n](B) - q[1-n](-B)) / 2.

    :param values: data to be evaluated
    :type values: tuple or list or np.ndarray
    :param quantiles: iterable containing desired quartiles between 0 and 1
    :type quantiles: tuple or list or np.ndarray
    :param weights: weights of data in the same order
    :type weights: tuple or list or np.ndarray
    :return: array containing approximated quartiles
    :rtype: np.ndarray
    """

    values = np.array(values)
    quantiles = np.array(quantiles)
    w = np.ones_like(values) if weights is None else np.array(weights)
    if np.any(quantiles < 0) or np.any(quantiles > 1):
        raise ValueError(f'quantiles ({quantiles}) should be in range 0 to 1')

    sorter = np.argsort(values)
    values = values[sorter]
    w = w[sorter]

    weighted_cumsum_0to1 = (np.cumsum(w) - w[0]) / np.sum(w[1:])
    weighted_cumsum_1to0 = (np.cumsum(np.flip(w)) - w[-1]) / np.sum(w[:-1])
    results_0to1 = np.interp(quantiles, weighted_cumsum_0to1, values)
    results_1to0 = np.interp(1-quantiles, weighted_cumsum_1to0, np.flip(values))
    return (results_0to1 + results_1to0) / 2.0
