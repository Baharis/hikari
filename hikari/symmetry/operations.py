"""
This file contains class and definitions of 3-dimensional symmetry operations.
"""

import numpy as np
from fractions import Fraction
from enum import Enum


class SymmOp:
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
        self.tf = transformation
        self.tl = translation

    def __eq__(self, other):
        return np.array_equal(self.tf, other.tf) and \
               np.array_equal(self.__tl24, other.__tl24)

    def __mul__(self, other):
        assert isinstance(other, SymmOp)
        return SymmOp(self.tf @ other.tf, self.tf @ other.tl + self.tl)

    def __pow__(self, power, modulo=None):
        return SymmOp.from_matrix(np.linalg.matrix_power(self.matrix, power))

    def __mod__(self, other):
        return SymmOp(self.tf, np.mod(self.__tl24, 24 * other) / 24)

    def __repr__(self):
        return f'SymmOp(np.{repr(self.tf)}, np.{repr(self.tl)})'.\
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
        :rtype: SymmOp
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
            tl[i] = 0 if coord is '' else float(Fraction(coord))
        return cls(tf, tl)

    @classmethod
    def from_matrix(cls, matrix):
        """
        Create new symmetry operation using augmented 4x4 transformation matrix

        :param matrix: augmented 4x4 matrix
        :type matrix: np.ndarray
        :return: Symmetry operation generated based on augmented matrix
        :rtype: SymmOp
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
        :rtype: SymmOp
        """
        return cls(matrix, vector)

    @staticmethod
    def _row_to_str(xyz, r):
        """
        Convert xyz: 3-el. list and r - number to single element of code triplet

        :param xyz: 3-element list of coordinates, eg
        :type xyz: Union[list, np.ndarray]
        :param r: translation applied to the row
        :type r: Union[float, Fraction]
        :return: string representing change of coordinate
        :rtype: str
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
    def _project(vector, onto):
        """Return projection of np.ndarray "vector" to np.ndarray "onto" """
        return (np.dot(vector, onto) / np.sqrt(sum(onto ** 2)) ** 2) * onto

    @property
    def code(self):
        return ','.join([self._row_to_str(xyz, r) for xyz, r
                         in zip(self.tf, self.tl)])

    @property
    def tf(self):
        return self.__tf

    @tf.setter
    def tf(self, value):
        self.__tf = np.rint(value).astype(int)

    @property
    def tl(self):
        return self.__tl24 / 24

    @tl.setter
    def tl(self, value):
        self.__tl24 = np.rint(value * 24).astype(int)

    @property
    def matrix(self):
        """
        :return: Augmented 4 x 4 transformation matrix with float-type values
        :rtype: np.ndarray
        """
        matrix = np.eye(4, dtype=float)
        matrix[0:3, 0:3] = self.tf
        matrix[0:3, 3] = self.tl
        return matrix

    @property
    def det(self):
        """
        :return: determinant of 3x3 transformation part of operation's matrix
        :rtype: int
        """
        return np.linalg.det(self.tf)

    @property
    def typ(self):
        """
        :return: crystallographic type of this symmetry operation
        :rtype: SymmOp.Type
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
    def name(self):
        """
        :return: short name of symmetry operation, eg.: "m", "3" or "2_1"
        :rtype: str
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
        op = SymmOp(self.tf) if self.det > 0 else SymmOp(-self.tf)
        for f in (1, 2, 3, 4, 5, 6):
            if (op ** f).trace == 3:
                return f
        raise NotImplementedError('fold is not in range 1 to 6')

    @property
    def order(self) -> int:
        """
        Number of times operation has to be repeated to become
        a translation, eg.: n for all n-fold axes, 2 for other (max 6)
        """
        for f in (1, 2, 3, 4, 5, 6):
            if pow(self, f, 1).typ is self.Type.identity:
                return f
        raise NotImplementedError('order is not in range 1 to 6')

    @property
    def glide(self) -> np.ndarray:
        """
        Part of the translation vector stemming from operations' glide
        """
        return (self ** 24).__tl24 / 576 % 1

    @property
    def glide_fold(self) -> int:
        """
        Number of types glide component of the operation must be repeated
        in order to contain only integer values, eg.: 3 for "6_2", 4 for "d"
        """
        return max([24 // t for t in [*self.__tl24, 24] if t != 0])

    @property
    def origin(self) -> np.ndarray:
        """
        Selected point invariant to the symmetry operation
        """
        return (self.tl - self.glide) * 1 / 2

    @property
    def reciprocal(self):
        """
        :return: relevant symmetry operation in reciprocal hkl space
        :rtype: SymmOp
        """
        return SymmOp(np.linalg.inv(self.tf).T)

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
            o = SymmOp(-1 * self.tf).orientation
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

    def at(self, point):
        """
        Transform operation as little as possible so that its symmetry element
        contains "point". To be used after "into" if used together.

        :param point: Target coordinates of point which should lie in element
        :type point: np.array
        :return: New symmetry operation which contains "point" in its element
        :rtype: SymmOp
        """
        shift = np.array(point) - self.origin
        for invariant in self.invariants:
            shift -= self._project(shift, onto=invariant)
        return SymmOp(np.eye(3), shift) * self * SymmOp(np.eye(3), -shift)

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
        :rtype: SymmOp
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

            return SymmOp(rot_h @ self.tf @ np.linalg.inv(rot_h),
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

    def extincts(self, hkl):
        """
        Return boolean array with truth whenever reflection should be extinct

        :param hkl: An array containing one hkl or multiple hkls in columns
        :type hkl: np.ndarray
        :return: array of booleans, where extinct reflections are marked as True
        :rtype: np.array
        """
        cond1 = np.isclose(hkl, self.reciprocal.transform(hkl)).all(axis=1)
        cond2 = ~np.isclose(np.dot(hkl, self.glide) % 1, 0)
        return cond1 & cond2
