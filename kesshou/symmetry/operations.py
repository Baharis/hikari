"""
This file contains the definition of 3-dimensional point symmetry operations.
"""

import numpy as np
from fractions import Fraction
from enum import Enum

symm_ops = dict()
"""
Dictionary containing 3x3 matrices of all 3-dimensional point-symmetry 
generators. Please mind that this file does not contain all symmetry operations,
but only those that are necessary to generate the point groups.

The names of symmetry operations follow the international crystallographic
naming convention wherever possible. For trigonal and hexagonal system case,
the names are additionally prefixed with the letter 'h'. 
"""
symm_ops['1'] = np.array([[1, 0, 0], [0, 1, 0], [0, 0, 1]])
symm_ops['-1'] = np.array([[-1, 0, 0], [0, -1, 0], [0, 0, -1]])
symm_ops['2_x'] = np.array([[1, 0, 0], [0, -1, 0], [0, 0, -1]])
symm_ops['2_y'] = np.array([[-1, 0, 0], [0, 1, 0], [0, 0, -1]])
symm_ops['2_z'] = np.array([[-1, 0, 0], [0, -1, 0], [0, 0, 1]])
symm_ops['2_xy'] = np.array([[0, 1, 0], [1, 0, 0], [0, 0, -1]])
symm_ops['-3_z'] = np.array([[1, 1, 0], [-1, 0, 0], [0, 0, -1]])
symm_ops['3_xyz'] = np.array([[0, 0, 1], [1, 0, 0], [0, 1, 0]])
symm_ops['-3_xyz'] = np.array([[0, 0, -1], [-1, 0, 0], [0, -1, 0]])
symm_ops['4_z'] = np.array([[0, -1, 0], [1, 0, 0], [0, 0, 1]])
symm_ops['-4_z'] = np.array([[0, 1, 0], [-1, 0, 0], [0, 0, -1]])
symm_ops['m_x'] = np.array([[-1, 0, 0], [0, 1, 0], [0, 0, 1]])
symm_ops['m_y'] = np.array([[1, 0, 0], [0, -1, 0], [0, 0, 1]])
symm_ops['m_z'] = np.array([[1, 0, 0], [0, 1, 0], [0, 0, -1]])
symm_ops['m_xy'] = np.array([[0, -1, 0], [-1, 0, 0], [0, 0, 1]])
symm_ops['h2_x'] = np.array([[1, -1, 0], [0, -1, 0], [0, 0, -1]])
symm_ops['h2_x2y'] = np.array([[-1, 1, 0], [0, 1, 0], [0, 0, -1]])
symm_ops['h3_z'] = np.array([[0, -1, 0], [1, -1, 0], [0, 0, 1]])
symm_ops['h-3_z'] = np.array([[0, 1, 0], [-1, 1, 0], [0, 0, -1]])
symm_ops['h6_z'] = np.array([[1, -1, 0], [1, 0, 0], [0, 0, 1]])
symm_ops['h-6_z'] = np.array([[-1, 1, 0], [-1, 0, 0], [0, 0, -1]])
symm_ops['hm_x'] = np.array([[-1, 1, 0], [0, 1, 0], [0, 0, 1]])
symm_ops['hm_x2y'] = np.array([[1, -1, 0], [0, -1, 0], [0, 0, 1]])
symm_ops['hm_z'] = np.array([[1, 0, 0], [0, 1, 0], [0, 0, -1]])


class SymmOpType(Enum):
    """Enumerator class storing information about type of symmetry operation"""
    rotoinversion = 4
    identity = 3
    reflection = 2
    rotation = 1
    inversion = 0
    rototranslation = -1
    transflection = -2
    translation = -3


class SymmOp:
    """
    Class storing information about symmetry operations, with clear string
    representation and intuitive syntax for combining / utilising operations:
    "=" to compare two symmetry operations for logical equivalence,
    "*" to combine two symmetry operations or transform a vector,
    "** n" to apply symmetry operation n times,
    "% n" to restrict symmetry operation to n unit cells.
    """

    def __init__(self, transformation, translation=np.array([0, 0, 0])):
        self.transformation = np.array(transformation, dtype=float)
        self.translation = np.array(translation, dtype=float)

    def __eq__(self, other):
        return np.allclose(self.matrix, other.matrix)

    def __mul__(self, other):
        if isinstance(other, SymmOp):
            return SymmOp.from_matrix(self.matrix @ other.matrix)
        if isinstance(other, np.ndarray):
            if other.shape[0] == 3:  # works only for 1,could work for multiple?
                return self.transformation @ other + self.translation
            elif other.shape[0] == 4:
                return self.matrix @ other
        raise TypeError('Cannot multiply "{}" and "{}"'.format(self, other))

    def __pow__(self, power, modulo=None):
        return SymmOp.from_matrix(np.linalg.matrix_power(self.matrix, power))

    def __mod__(self, other):
        return SymmOp(self.transformation, np.mod(self.translation, other))

    def __str__(self):
        code = ','.join([self._row_to_str(xyz, r) for xyz, r
                         in zip(self.transformation, self.translation)])
        origin = ','.join([str(Fraction(o).limit_denominator(9))
                           for o in self.origin])
        return self.name + ': ' + code + ' (' + origin + ')'

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
    def matrix(self):
        """
        :return: Augmented 4 x 4 transformation matrix with float-type values
        :rtype: np.ndarray
        """
        matrix = np.eye(4, dtype=float)
        matrix[0:3, 0:3] = self.transformation
        matrix[0:3, 3] = self.translation
        return matrix

    @property
    def det(self):
        """
        :return: determinant of 3x3 transformation part of operation's matrix
        :rtype: int
        """
        return int(round(np.linalg.det(self.transformation)))

    @property
    def typ(self):
        """
        :return: crystallographic type of this symmetry operation
        :rtype: SymmOpType
        """
        _trans = self.translational
        if np.isclose(self.trace, 3):
            return SymmOpType.translation if _trans else SymmOpType.identity
        elif np.isclose(self.trace, -3):
            return SymmOpType.inversion
        elif len(self.invariants) == 0:
            return SymmOpType.rotoinversion
        elif self.det < 0:
            return SymmOpType.transflection if _trans else SymmOpType.reflection
        else:
            return SymmOpType.rototranslation if _trans else SymmOpType.rotation

    @property
    def name(self):
        """
        :return: short name of symmetry operation, eg.: "m", "3" or "2_1"
        :rtype: str
        """
        _g = self.glide
        _glide_dir = 'n' if np.linalg.norm(_g) < 1e-8 else 'x'
        d = {(1, 0, 0): 'a', (0, 1, 0): 'b', (0, 0, 1): 'c',
             (0, 1, 1): 'A', (1, 0, 1): 'B', (1, 1, 0): 'C', (1, 1, 1): 'I'}
        for key, value in d.items():
            if np.allclose(_g, self._project(_g, np.array(key))):
                _glide_dir = value

        if self.typ is SymmOpType.identity:
            return '1'
        elif self.typ is SymmOpType.reflection:
            return 'm'
        elif self.typ is SymmOpType.rotation:
            return str(self.fold)
        elif self.typ is SymmOpType.inversion:
            return '-1'
        elif self.typ is SymmOpType.rototranslation:
            return str(self.fold) + '_' + str(min(np.rint([g * self.glide_fold
                            for g in _g if not np.isclose(g, 0)]).astype(int)))
        elif self.typ is SymmOpType.transflection:
            return _glide_dir.replace('A', 'n').replace('B', 'n').\
                replace('C', 'n') if self.glide_fold is 2 else 'd'
        elif self.typ is SymmOpType.translation:
            return 't_' + _glide_dir
        elif self.typ is SymmOpType.rotoinversion:
            return str(-self.fold)
        else:
            return '?'

    @property
    def fold(self):
        """
        :return: number of times operation has to be repeated to become
        translation, eg.: n for all n-fold axes, 2 for reflections (max 6)
        :rtype: int
        """
        point_self = SymmOp(self.transformation)
        for f in (1, 2, 3, 4, 5, 6):
            if point_self ** f == SymmOp(np.eye(3)):
                return f

    @property
    def order(self):
        """
        :return: number of times operation has to be repeated to become
        identity: equal to fold for non-translational, inf for translational ops
        :rtype: Union[int, float]
        """
        for f in (1, 2, 3, 4, 5, 6):
            if self ** f == SymmOp(np.eye(3)):
                return f
        return float('Inf')

    @property
    def glide(self):
        """
        :return: part of the translation vector stemming from operations' glide
        :rtype: np.ndarray
        """
        return (self ** self.fold).translation / self.fold

    @property
    def glide_fold(self):
        """
        :return: number of types glide component of the operation must be
        repeated to contain only integer values, eg.: 3 for "6_2", 4 for "d"
        :rtype:
        """
        return max([Fraction(comp).limit_denominator(9).denominator
                    for comp in self.glide])

    @property
    def origin(self):
        """
        :return: selected point belonging to symmetry element of the operation
        :rtype: np.ndarray
        """
        return (self.translation - self.glide) * 1/2

    @property
    def reciprocal(self):
        """
        :return: relevant symmetry operation in reciprocal hkl space
        :rtype: SymmOp
        """
        return SymmOp(np.linalg.inv(self.transformation).T)

    @property
    def trace(self):
        """
        :return: trace of 3x3 transformation part of operation's matrix
        :rtype:
        """
        return np.trace(self.transformation)

    @property
    def translational(self):
        """
        :return: True if operation has any glide component, False otherwise
        :rtype: bool
        """
        return not np.allclose(self.glide, np.zeros_like(self.glide))

    @property
    def invariants(self):
        """
        :return: List of directions not affected by this symmetry operation
        :rtype: list[np.ndarray]
        """
        eigenvalues, eigenvectors = np.linalg.eig(self.transformation)
        return eigenvectors.T[np.isclose(eigenvalues.real, 1)].real

    @property
    def orientation(self):
        """
        :return: Direction of symmetry element (if can be defined) else None
        :rtype: Union[np.ndarray, None]
        """
        if self.typ in {SymmOpType.rotation, SymmOpType.rototranslation}:
            return self.invariants[0]
        elif self.typ in {SymmOpType.reflection, SymmOpType.transflection,
                          SymmOpType.rotoinversion}:
            return SymmOp(-1 * self.transformation).orientation
        else:
            return None

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
        return SymmOp(self.transformation, self.translation + 2 * shift)

    def into(self, direction, hexagonal=False):
        """
        Rotate operation so that its orientation changes to "direction", while
        preserving fractional glide. To be used before respective "at" method.
        :param direction: Target orientation for element of symmetry operation
        :type direction: np.ndarray
        :param hexagonal: True if operation is defined in hexagonal coordinates
        :type hexagonal: bool
        :return: New symmetry operation whose orientation is "direction"
        :rtype: SymmOp
        """
        d = np.array(direction) / np.linalg.norm(direction)
        o = self.orientation / np.linalg.norm(self.orientation)
        ref_change = np.array(((1, -1/2, 0), (0, np.sqrt(3)/2, 0), (0, 0, 1))) \
            if hexagonal else np.eye(3)
        d = ref_change @ d
        d /= np.linalg.norm(d)
        o = ref_change @ o / np.linalg.norm(o)
        d /= np.linalg.norm(o)

        def are_parallel(v, w):
            return np.isclose(np.dot(v, w), 1)

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
            rot_m = np.eye(3) + v + (1 / (1 + np.dot(o, d))) * v @ v
            new_origin = self.origin
            new_glide = np.linalg.inv(ref_change) @ rot_m @ ref_change @ self.glide
            new_glide_comps = [abs(comp) for comp in new_glide if not np.isclose(comp, 0)]
            new_glide = new_glide if len(new_glide_comps) == 0 else (new_glide / min(new_glide_comps)) / self.glide_fold
            return SymmOp(np.linalg.inv(ref_change) @ rot_m @ ref_change @ self.transformation @ np.linalg.inv(ref_change) @ np.linalg.inv(rot_m) @ ref_change,
                          2 * new_origin + new_glide)
        # TODO clean this bullcrap and implement hkl/xyz transformation

    def _transform_hkl(self, hkl):
        raise NotImplementedError
        #return self.reciprocal.transformation @ np.array(hkl)

    def _transform_xyz(self, xyz):
        raise NotImplementedError

    def _extincts_hkl(self, hkl):
        raise NotImplementedError
        # hkl = np.array(hkl)
        # cond1 = np.allclose(hkl, self._transform_hkl(hkl))
        # cond2 = np.isclose(np.dot(hkl, self.glide) % 1, 0)
        # return cond1 and cond2


if __name__ == '__main__':
    t = SymmOp.from_code('-x,y+1/2,-z')
    u = SymmOp.from_code('x,-y,z+1/2')

    o1 = SymmOp.from_code('x+1/2, y+1/2, z')
    o2 = SymmOp.from_code('x+1/2, y, z+1/2')
    o3 = SymmOp.from_code('x, y+1/2, z+1/2')
    o4 = SymmOp.from_code('1/4-x, 1/4+y, 1/4+z')
    o5 = SymmOp.from_code('1/4+x, 1/4-y, 1/4+z')
    o6 = SymmOp.from_code('-x, -y, z')
    # print('glide fold:', t.glide_fold)
    # print(t.invariants)
    # print(t.matrix4x4)
    # print(t.det)
    # print(t, t.glide)
    # print('----------')
    # print(t.typ)
    # print(t.orientation)

    def build(*gens):
        new_gens = list(gens)
        for g1 in gens:
            for g2 in gens:
                if not g1 * g2 % 1 in new_gens:
                    new_gens.append(g1 * g2 % 1)
        if new_gens == list(gens):
            print('\n'.join([str(g) for g in gens]))
        else:
            build(*new_gens)

    build(o1, o2, o3, o4, o5, o6)

    # print('t._transform_hkl((1, 0, 0))')
    # print(t._transform_hkl((0, 0, 1)))
    # print('-' * 30)
    #
    # print('t._extincts_hkl((0, 0, 1))')
    # print(t._extincts_hkl((0, 0, 1)))
    # print('-' * 30)
    #
    # print('t._extincts_hkl((0, 0, 3))')
    # print(t._extincts_hkl((0, 0, 3)))
    # print('-' * 30)
    #
    # print('t._extincts_hkl((0, 0, 4))')
    # print(t._extincts_hkl((0, 0, 4)))
    # print('-' * 30)
    #
    # print('t._extincts_hkl((0, 1, 4))')
    # print(t._extincts_hkl((0, 1, 4)))
    # print('-'*30)
    #
    # print(t * np.array((1.6, 1.8, 2.2)))
    # print(t**2 * np.array((1.6, 1.8, 2.2)))
    # print(t**3 * np.array((1.6, 1.8, 2.2)))
    # print(t**4 * np.array((1.6, 1.8, 2.2)))
    # print(t ** 5 * np.array((1.6, 1.8, 2.2)))
