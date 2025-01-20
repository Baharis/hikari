"""
This file contains class and definitions of 3-dimensional symmetry operations.
"""

import numpy as np
from fractions import Fraction
from enum import Enum
from typing import Sequence, Union


class Operation:
    """
    Class storing information about symmetry operations, with clear string
    representation and intuitive syntax for combining / utilising operations:
    "=" to compare two symmetry operations for logical equivalence,
    "*" to combine two symmetry operations or transform a vector,
    "** n" to apply symmetry operation n times,
    "% n" to restrict symmetry operation to n unit cells.
    Some of the functions may work incorrectly for rhombohedral unit cells #TODO
    """
    class Type(Enum):
        """Enumerator class storing information about type of operation"""
        rotoinversion = 4
        identity = 3
        reflection = 2
        rotation = 1
        inversion = 0
        rototranslation = -1
        transflection = -2
        translation = -3

    def __init__(self, transformation, translation=np.array([0, 0, 0])):
        self.tf = np.array(transformation)
        self.tl = np.array(translation)

    def __eq__(self, other):
        return np.array_equal(self.tf, other.tf) and \
               np.array_equal(self._tl24, other._tl24)

    def __mul__(self, other):
        assert isinstance(other, Operation)
        return self.__class__(self.tf @ other.tf, self.tf @ other.tl + self.tl)

    def __pow__(self, power, modulo=None):
        return self.__class__.from_matrix(np.linalg.matrix_power(self.matrix, power))

    # TODO removed __mod__ operator: replace all uses of % 1 with new `bounded()`
    # def __mod__(self, other):
    #     return self.__class__(self.tf, np.mod(self._tl24, 24 * other) / 24)

    def __repr__(self):
        return f'{self.__class__.__name__}(np.{repr(self.tf)}, np.{repr(self.tl)})'.\
            replace('\n', '').replace(' ', '')

    def __str__(self):
        origin = ','.join([str(Fraction(o).limit_denominator(12))
                           for o in self.origin])
        return self.name + ': ' + self.code + ' (' + origin + ')'

    def __hash__(self):
        return hash(repr(self))

    @classmethod
    def from_code(cls, code):
        """
        Create new symmetry operation using symmetry code "x',y',z'" containing
        transformation of each individual coordinate from x,y,z to x',y',z'.

        :param code: string representing new coordinates after operation
        :type code: str
        :return: Symmetry operation generated from given coordinate triplet code
        :rtype: Operation
        """
        coords = str(code).replace(';', ',').replace(' ', '').lower().split(',')
        tf = np.zeros((3, 3))
        tl = np.zeros(3)
        for i, coord in enumerate(coords):
            coord = coord if coord[0] in '+-' else '+' + coord
            tf[i, 0] = 1 if '+x' in coord else -1 if '-x' in coord else 0
            tf[i, 1] = 1 if '+y' in coord else -1 if '-y' in coord else 0
            tf[i, 2] = 1 if '+z' in coord else -1 if '-z' in coord else 0
            coord = coord.replace('+x', '').replace('+y', '').replace('+z', '')\
                .replace('-x', '').replace('-y', '').replace('-z', '')
            tl[i] = 0 if coord == '' else float(Fraction(coord))
        return cls(tf, tl)

    @classmethod
    def from_matrix(cls, matrix):
        """
        Create new symmetry operation using augmented 4x4 transformation matrix

        :param matrix: augmented 4x4 matrix
        :type matrix: np.ndarray
        :return: Symmetry operation generated based on augmented matrix
        :rtype: Operation
        """
        return cls(matrix[0:3, 0:3], matrix[0:3, 3])

    @classmethod
    def from_pair(cls, matrix, vector):
        """
        Create new symmetry operation using point transformation 3x3 matrix
        and 3-length translation vector. Alias for standard creation method.

        :param matrix: 3x3 point transformation matrix
        :type matrix: np.ndarray
        :param vector: 3-length translation vector
        :type vector: np.ndarray
        :return: Symmetry operation generated based on matrix - vector pair
        :rtype: Operation
        """
        return cls(matrix, vector)

    @staticmethod
    def _row_to_str(xyz: np.ndarray[int], r: float) -> str:
        """
        Convert xyz: 3-el. list and r - number to single element of code triplet

        :param xyz: 3-element list of coordinates, e.g.: [1,-1,0]
        :param r: translation applied to the row
        :return: string representing change of coordinate
        """
        s = ''
        for var, char in zip(xyz, 'xyz'):
            if np.isclose(abs(var), 0):
                continue
            s += '{:+f}'.format(var)[0]
            if np.isclose(abs(var), 1):
                s += char
                continue
            s += str(Fraction(abs(var)).limit_denominator(9)) + char
        if not np.isclose(abs(r), 0):
            s += '{:+f}'.format(float(r))[0]
            s += str(Fraction(abs(r)).limit_denominator(9))
        return s.strip('+')

    @staticmethod
    def _project(vector: np.ndarray, onto: np.ndarray) -> np.ndarray:
        """Return projection of np.ndarray "vector" to np.ndarray "onto" """
        return (np.dot(vector, onto) / np.sqrt(sum(onto ** 2)) ** 2) * onto

    @property
    def code(self) -> str:
        return ','.join([self._row_to_str(xyz, r) for xyz, r
                         in zip(self.tf, self.tl)])

    @property
    def tf(self) -> np.ndarray:
        return self._tf

    @tf.setter
    def tf(self, value: np.ndarray) -> None:
        self._tf = np.rint(value).astype(int)

    @property
    def tl(self) -> np.ndarray:
        return self._tl24 / 24

    @tl.setter
    def tl(self, value: np.ndarray) -> None:
        self._tl24 = np.rint(value * 24).astype(int)

    @property
    def matrix(self) -> np.ndarray:
        """
        :return: Augmented 4 x 4 transformation matrix with float-type values
        :rtype: np.ndarray
        """
        matrix = np.eye(4, dtype=float)
        matrix[0:3, 0:3] = self.tf
        matrix[0:3, 3] = self.tl
        return matrix

    @property
    def det(self) -> np.float64:
        """
        :return: determinant of 3x3 transformation part of operation's matrix
        :rtype: int
        """
        return np.linalg.det(self.tf)

    @property
    def typ(self) -> Type:
        """
        :return: crystallographic type of this symmetry operation (see `Type`)
        """
        _trans = self.translational
        if self.trace == 3:
            return self.Type.translation if _trans else self.Type.identity
        elif self.trace == -3:
            return self.Type.inversion
        elif len(self.invariants) == 0:
            return self.Type.rotoinversion
        elif self.det < 0:
            return self.Type.transflection if _trans else self.Type.reflection
        else:
            return self.Type.rototranslation if _trans else self.Type.rotation

    @property
    def name(self) -> str:
        """
        :return: short name of symmetry operation, e.g.: "m", "3" or "2_1"
        """
        _glide = self.glide
        _glide_dir = 'n' if np.linalg.norm(_glide) < 1e-8 else 'x'
        d = {(1, 0, 0): 'a', (0, 1, 0): 'b', (0, 0, 1): 'c',
             (0, 1, 1): 'A', (1, 0, 1): 'B', (1, 1, 0): 'C', (1, 1, 1): 'I',
             (1, 2, 2): 'R', (2, 1, 1): 'R', (2, 1, 2): 'R', (1, 2, 1): 'R'}
        for key, value in d.items():
            if np.allclose(_glide, self._project(_glide, np.array(key))):
                _glide_dir = value

        if self.typ is self.Type.identity:
            return '1'
        elif self.typ is self.Type.reflection:
            return 'm'
        elif self.typ is self.Type.rotation:
            return str(self.fold) + self.sense
        elif self.typ is self.Type.inversion:
            return '-1'
        elif self.typ is self.Type.rototranslation:
            return str(self.fold) + \
                   str(round(np.linalg.norm(_glide)*self.fold)) + self.sense
        elif self.typ is self.Type.transflection:
            return _glide_dir.replace('A', 'n').replace('B', 'n').\
                replace('C', 'n') if self.glide_fold == 2 else 'd'
        elif self.typ is self.Type.translation:
            return _glide_dir
        elif self.typ is self.Type.rotoinversion:
            return str(-self.fold) + self.sense
        else:
            return '?'

    @property
    def fold(self) -> int:
        """
        Number of times operation must be repeated to become identity, inversion
        or translation: n for n-fold axes, 2 for reflections, 1 for other (max 6)
        """
        op = BoundedOperation(self.tf if self.det > 0 else -self.tf)
        for f in (1, 2, 3, 4, 5, 6):
            if (op ** f).trace == 3:
                return f
        raise NotImplementedError('fold is not in range 1 to 6')

    @property
    def order(self) -> int:
        """
        Number of times operation has to be repeated to become
        a translation, e.g.: n for all n-fold axes, 2 for other (max 6)
        """
        op = self.bounded()
        for f in (1, 2, 3, 4, 5, 6):
            if pow(op, f, 1).typ is self.Type.identity:
                return f
        raise NotImplementedError('order is not in range 1 to 6')

    @property
    def glide(self) -> np.ndarray:
        """
        Part of the translation vector stemming from operations' glide
        """
        return (self.unbounded() ** 24)._tl24 / 576

    @property
    def glide_fold(self) -> int:
        """
        Number of types glide component of the operation must be repeated
        in order to contain only integer values, e.g.: 3 for "6_2", 4 for "d"
        """
        return max([24 // t for t in [*self._tl24, 24] if t != 0])

    @property
    def origin(self) -> np.ndarray:
        """Selected point that remains on the symmetry element after operation"""
        return np.linalg.pinv(np.eye(3) - self.tf) @ self.tl

    # TODO Likely incorrect, see ITC 5.2.1. Transformations and example for P4/n
    # TODO Fixed? Test rigorously.

    @property
    def reciprocal(self):
        """
        :return: relevant symmetry operation in reciprocal hkl space
        :rtype: Operation
        """
        return PointOperation(np.linalg.inv(self.tf).T)

    @property
    def trace(self):
        """
        :return: trace of 3x3 transformation part of operation's matrix
        :rtype:
        """
        return np.trace(self.tf)

    @property
    def translational(self):
        """
        :return: True if operation has any glide component, False otherwise
        :rtype: bool
        """
        return not np.allclose(self.glide, 0)

    @property
    def invariants(self):
        """
        :return: List of directions not affected by this symmetry operation
        :rtype: list[np.ndarray]
        """
        eigenvalues, eigenvectors = np.linalg.eig(self.tf)
        return eigenvectors.T[np.isclose(eigenvalues.real, 1)].real

    @property
    def orientation(self):
        """
        :return: Direction of symmetry element (if can be defined) else None
        :rtype: Union[np.ndarray, None]
        """
        if self.det < 0:
            o = Operation(-1 * self.tf).orientation
        elif self.typ in {self.Type.rotation, self.Type.rototranslation}:
            m = self.tf  # method: see https://math.stackexchange.com/q/3441262
            t = m + m.T - (m.trace() - 1) * np.eye(3)
            o = [_ for _ in t if not np.isclose(sum(abs(_)), 0)][0]
        else:
            return None
        if o is None:
            return o
        o = o if sum(o) > 0 else -o
        return o / np.sqrt(sum(o*o))

    @property
    def sense(self):
        """
        :return: "+" or "-", the "sense" of rotation, as given in ITC A, 11.1.2
        :rtype: str
        """
        unique = np.array([np.sqrt(2), np.e, np.pi])
        rotation = self.tf if self.det > 0 else -self.tf
        sign = np.dot(self.orientation, np.cross(unique, rotation @ unique))
        return '' if np.isclose(sign, 0) else '+' if sign > 0 else '-'

    @property
    def is_bounded(self) -> bool:
        """
        True if the operation is a coset representatives i.e. if all
        translation components lie in range -0.5 exclusive to 0.5 inclusive.
        """
        return all(self._tl24 > -12) and all(self._tl24 <= 12)

    def bounded(self) -> 'BoundedOperation':
        return self if isinstance(self, BoundedOperation) \
            else BoundedOperation(self.tf, self.tl)

    def unbounded(self) -> 'Operation':
        return Operation(self.tf, self.tl)

    def at(self, point: np.ndarray) -> 'Operation':
        """
        Transform operation as little as possible so that its symmetry element
        contains "point". To be used after "into" if used together.
        Based on "ITC 5.2.1. Transformations".

        :param point: Target coordinates of point which should lie in element
        :return: New symmetry operation which contains "point" in its element
        """
        shift = np.array(point) - self.origin
        for invariant in self.invariants:
            shift -= self._project(shift, onto=invariant)
        return Operation(np.eye(3), shift) * self * Operation(np.eye(3), -shift)

    def into(self, direction, hexagonal=False):
        """
        Rotate operation so that its orientation changes to "direction", while
        preserving fractional glide. To be used before respective "at" method.
        Will most likely not work for unimplemented rhombohedral unit cells.

        :param direction: Target orientation for element of symmetry operation
        :type direction: np.ndarray
        :param hexagonal: True if operation is defined in hexagonal coordinates
        :type hexagonal: bool
        :return: New symmetry operation whose orientation is "direction"
        :rtype: Operation
        """
        rebase = np.array(((1, -1/2, 0), (0, np.sqrt(3)/2, 0), (0, 0, 1))) \
            if hexagonal else np.eye(3)
        d = rebase @ np.array(direction)
        d /= np.linalg.norm(d)
        o = rebase @ self.orientation
        o /= np.linalg.norm(o)

        def are_parallel(_v, _w):
            return np.isclose(np.dot(_v, _w), 1)

        if o is None:
            return self
        elif are_parallel(d, o):
            return self
        elif are_parallel(-d, o):
            x, y = np.array((1, 0, 0)), np.array((0, 1, 0))
            temp_vector = y if are_parallel(d, x) or are_parallel(d, -x) else x
            return self.into(temp_vector, hexagonal).into(d, hexagonal)
        else:
            p = np.cross(o, d)
            v = np.array(((0, -p[2], p[1]), (p[2], 0, -p[0]), (-p[1], p[0], 0)))
            rot = np.eye(3) + v + (1 / (1 + np.dot(o, d))) * v @ v
            rot_h = np.linalg.inv(rebase) @ rot @ rebase

            if np.allclose(self.glide, 0):
                new_glide = np.array((0, 0, 0))
            else:
                new_glide = rot_h @ self.glide
                glide_els = [abs(g) for g in new_glide if not np.isclose(g, 0)]
                new_glide = (new_glide / min(glide_els)) / self.glide_fold

            return Operation(rot_h @ self.tf @ np.linalg.inv(rot_h),
                             2 * self.origin + new_glide)

    def transform(self, other):
        """
        Transform a column containing rows of coordinate points

        :param other: A vertical numpy array of coordinate triplets kept in rows
        :type other: np.ndarray
        :return: Same-shaped array of coordinate triplets transformed by self
        :rtype: np.ndarray
        """
        if other.shape[1] == 3:
            return (self.tf @ other.T).T + self.tl
        elif other.shape[1] == 4:
            return self.matrix @ other.T
        raise TypeError('Cannot transform "{}"'.format(other))

    def extincts(self, hkl: np.ndarray) -> np.ndarray:
        """
        Return boolean array with truth whenever reflection should be extinct

        :param hkl: An array containing one hkl or multiple hkls in columns
        :return: array of booleans, where extinct reflections are marked as True
        """
        cond1 = np.isclose(hkl, self.reciprocal.transform(hkl)).all(axis=1)
        cond2 = ~np.isclose(np.dot(hkl, self.glide) % 1, 0)
        return cond1 & cond2

    def distance_to_point(self, point: Sequence[Union[int, float]]) -> float:
        if self.typ in {self.Type.rotation, self.Type.rototranslation, self.Type.rotoinversion}:
            point = np.array(point) - self.origin
            projection = np.dot(point, self.orientation) * self.orientation
            orthogonal_vector = point - projection
            return float(np.linalg.norm(orthogonal_vector))
        return NotImplemented


class BoundedOperation(Operation):
    """
    A subclass of `Operation` where all three elements of the translation
    vector are bounded between 0 (inclusive) and 1 (exclusive).
    This class is suitable for handling coset representatives of space groups,
    since upon binding operations related by unit translation become equivalent.
    """

    @property
    def tl(self):
        return self._tl24 / 24

    @tl.setter
    def tl(self, value):
        """
        This setter automatically handles binding the translation component
        of the `Operation` into the [0; 1) range.
        This makes the `Operation` a "bounded" "coset representative".
        With this operation 0 casts onto 0, +/-0.5 onto 0.5, and 1 back onto 0.
        """
        self._tl24 = np.rint(value * 24).astype(int) % 24


class PointOperation(BoundedOperation):
    """A subclass of `BoundedOperation`, asserts translation vector = [0,0,0]"""
    @property
    def tl(self):
        return np.array([0, 0, 0], dtype=float)

    @tl.setter
    def tl(self, value):
        if np.any(value):
            raise ValueError('Translation in `PointOperation` must be [0,0,0]')
        self._tl24 = np.array([0, 0, 0], dtype=int)
