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

    The algorithm utilised in this function gives a distribution which is
    even on the surface of a sphere. Contrary to an even distribution of two
    angle values in spherical representation, this algorithm does not favour
    points in a proximity to sphere poles. For more details please refer to
    a fantastic article by Martin Roberts `here
    <http://extremelearning.com.au/evenly-distributing-points-on-a-sphere/>`_.

    The function is written so that it does always return the same set of points
    in order to ensure reproducibility of other methods. A list of other points
    can be obtained by changing (or even randomising) the value of a *seed*.

    :Example:

    >>> fibonacci_sphere(4)
    [(0.62069, -0.75, -0.22857),
     (-0.44395, -0.25, 0.86047),
     (0.41275, 0.25, 0.87587),
     (-0.61208, 0.75, -0.25072)]
    >>> fibonacci_sphere(4)
    [(0.62069, -0.75, -0.22857),
     (-0.44395, -0.25, 0.86047),
     (0.41275, 0.25, 0.87586),
     (-0.61208, 0.75, -0.25072)]
    >>> fibonacci_sphere(4, seed=420)
    [(0.64040, -0.75, 0.16550),
     (-0.85489, -0.25, 0.45460),
     (0.32329, 0.25, -0.91268),
     (0.25831, 0.75, 0.60892)]

    :param samples: Number of points to be generated.
    :type samples: int
    :param seed: A seed value used once to slightly randomise point position.
    :type seed: any
    :return: A list of 3-element tuple containing points' cartesian coordinates.
    :rtype: list
    """
    random.seed(seed)
    rnd = random.random() * samples
    points = []
    offset = 2. / samples
    increment = np.pi * (3. - np.sqrt(5.))
    for i in range(samples):
        y = ((i * offset) - 1) + (offset / 2)
        r = np.sqrt(1 - pow(y, 2))
        phi = ((i + rnd) % samples) * increment
        x = np.cos(phi) * r
        z = np.sin(phi) * r
        points.append((x, y, z))
    return points


def is2n(integer):
    """
    Check whether supplied *integer* is divisible by 2.
    For each divisible element, return True.
    For each indivisible element, return False.
    Equivalent to `*int_or_iterable* % 2 == 0`.

    :example:

    >>> is2n(6)
    True
    >>> is2n(-8)
    True
    >>> is2n(7)
    False

    :param integer: Value to be checked for divisibility.
    :type integer: int
    :return: True if remainder of *integer* / 2 == 0; False otherwise.
    :rtype: bool
    """
    return integer % 2 == 0


def is3n(integer):
    """
    Check whether supplied *integer* is divisible by 3.
    For each divisible element, return True.
    For each indivisible element, return False.
    Equivalent to `*int_or_iterable* % 3 == 0`.

    :example:

    >>> is3n(6)
    True
    >>> is3n(-8)
    False
    >>> is3n(7)
    False

    :param integer: Value to be checked for divisibility.
    :type integer: int
    :return: True if remainder of *integer* / 3 == 0; False otherwise.
    :rtype: bool
    """
    return integer % 3 == 0


def is4n(integer):
    """
    Check whether supplied *integer* is divisible by 4.
    For each divisible element, return True.
    For each indivisible element, return False.
    Equivalent to `*int_or_iterable* % 4 == 0`.

    :example:

    >>> is4n(6)
    False
    >>> is4n(-8)
    True
    >>> is4n(7)
    False

    :param integer: Value to be checked for divisibility.
    :type integer: int
    :return: True if remainder of *integer* / 4 == 0; False otherwise.
    :rtype: bool
    """
    return integer % 4 == 0


def is6n(integer):
    """
    Check whether supplied *integer* is divisible by 6.
    For each divisible element, return True.
    For each indivisible element, return False.
    Equivalent to `*int_or_iterable* % 6 == 0`.

    :example:

    >>> is6n(6)
    True
    >>> is6n(-8)
    False
    >>> is6n(7)
    False

    :param integer: Value to be checked for divisibility.
    :type integer: int
    :return: True if remainder of *integer* / 6 == 0; False otherwise.
    :rtype: bool
    """
    return integer % 6 == 0
