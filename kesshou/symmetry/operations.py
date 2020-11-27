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
    rotoinversion = 4
    identity = 3
    reflection = 2
    rotation = 1
    inversion = 0
    rototranslation = -1
    transflection = -2
    translation = -3


class SymmOp:
    def __init__(self, transformation, translation=np.array([0, 0, 0])):
        self.transformation = np.array(transformation, dtype=float)
        self.translation = np.array(translation, dtype=float)

    def __eq__(self, other):
        return np.allclose(self.matrix4x4, other.matrix4x4)

    def __mul__(self, other):
        if isinstance(other, SymmOp):
            return SymmOp.from_matrix(self.matrix4x4 @ other.matrix4x4)
        if isinstance(other, np.ndarray):
            if other.shape[0] == 3:
                return self.matrix3x3 @ other + self.translation #this works only for 1x1, could work for multiple? see
            elif other.shape[0] == 4:
                return self.matrix4x4 @ other
        raise TypeError('Cannot multiply "{}" and "{}"'.format(self, other))

    def __pow__(self, power, modulo=None):
        return SymmOp.from_matrix(np.linalg.matrix_power(self.matrix4x4, power))

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
        matrix = np.array(matrix)
        if matrix.shape[0] == matrix.shape[1] == 3:
            return cls(matrix)
        if matrix.shape[0] == matrix.shape[1] == 4:
            return cls(matrix[0:3, 0:3], matrix[0:3, 3])

    @staticmethod
    def _row_to_str(xyz, r):
        s = ''
        for var, char in zip(xyz, 'xyz'):
            if np.isclose(abs(var), 0):
                continue
            s += '{:+f}'.format(var)[0]
            if np.isclose(abs(var), 1):
                s += char
                continue
            s += str(Fraction(abs(var)).limit_denominator(9)) + char
        if np.isclose(abs(r), 0):
            pass
        else:
            s += '{:+f}'.format(float(r))[0]
            s += str(Fraction(abs(r)).limit_denominator(9))
        return s.strip('+')

    @staticmethod
    def _project(vector, onto):
        return (np.dot(vector, onto) / np.sqrt(sum(onto ** 2)) ** 2) * onto

    @property
    def matrix3x3(self):
        return self.transformation

    @property
    def matrix4x4(self):
        matrix = np.eye(4, dtype=float)
        matrix[0:3, 0:3] = self.transformation
        matrix[0:3, 3] = self.translation
        return matrix

    @property
    def det(self):
        return int(round(np.linalg.det(self.transformation)))

    @property
    def typ(self):
        return SymmOpType((3 if self.trace == 3 else 0 if self.trace == -3 else
                           4 if len(self.invariants) == 0 else 2 if self.det < 0
                           else 1) * (1 if self.translational else -1))

    @property
    def name(self):
        _g = self.glide
        _glide_fold = ((Fraction(np.max(_g)).limit_denominator(9) % 1) *
                       self.fold).numerator
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
            return str(self.fold) + '_' + str(_glide_fold)
        elif self.typ is SymmOpType.transflection:
            return _glide_dir.replace('A', 'n').replace('B', 'n').\
                replace('C', 'n') if self.fold is 2 else 'd'
        elif self.typ is SymmOpType.translation:
            return 't_' + _glide_dir
        elif self.typ is SymmOpType.rotoinversion:
            return str(-self.fold)
        else:
            return '?'

    @property
    def fold(self):
        point_self = SymmOp(self.transformation)
        for f in (1, 2, 3, 4, 5, 6):
            if point_self ** f == SymmOp(np.eye(3)):
                return f

    @property
    def order(self):
        for f in (1, 2, 3, 4, 5, 6):
            if self ** f == SymmOp(np.eye(3)):
                return f
        return float('Inf')

    @property
    def glide(self):
        return (self ** self.fold).translation / self.fold

    @property
    def origin(self):
        return (self.translation - self.glide) * 1/2

    @property
    def reciprocal(self):
        return SymmOp(np.linalg.inv(self.transformation).T)

    @property
    def trace(self):
        return int(np.trace(self.transformation))

    @property
    def translational(self):
        return np.linalg.norm(self.glide) < 1e-8

    @property
    def _eigenpairs(self):
        return np.linalg.eig(self.matrix3x3)

    @property
    def invariants(self):
        eigenvalues, eigenvectors = np.linalg.eig(self.matrix3x3)
        return eigenvectors.T[np.isclose(eigenvalues.real, 1)].real

    @property
    def orientation(self):
        if self.typ in {SymmOpType.rotation, SymmOpType.rototranslation}:
            return self.invariants[0]
        elif self.typ in {SymmOpType.reflection, SymmOpType.transflection,
                          SymmOpType.rotoinversion}:
            return SymmOp(-1 * self.transformation).orientation
        else:
            return None

    def at(self, point):
        shift = np.array(point) - self.origin
        for invariant in self.invariants:
            shift -= self._project(shift, onto=invariant)
        return SymmOp(self.transformation, self.translation + 2 * shift)

    def into(self, direction, hexagonal=False):
        """Warning: works only if axis are orthogonal!"""
        d = np.array(direction) / np.linalg.norm(direction)
        o = self.orientation / np.linalg.norm(self.orientation)
        ref_change = np.array(((1, -1/2, 0), (0, np.sqrt(3)/2, 0), (0, 0, 1))) \
            if hexagonal else np.eye(3)
        d = ref_change @ d
        d /= np.linalg.norm(d)
        o = ref_change @ o / np.linalg.norm(o)
        d /= np.linalg.norm(o)
        print(d, 'new direction')
        print(o, 'old orientation')

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

            # print(o)
            # print(d)
            # print(rot_m, ":rot_m")
            # print(np.linalg.det(rot_m), ":det rot_m")
            # print(np.linalg.eig(rot_m)[0], ":eigval rot_m")
            # print(np.linalg.eig(rot_m)[1], ":eigvec rot_m")
            # print(SymmOp(rot_m).fold, ":fold of rot_m")

            #dubugging block - WORKS AS LONG AS ORIGINAL OPERATION IS IN ORTHORTOMBIC (eg 2||z)
            print('_' * 40)
            x = np.array((0.2, 0.4, 0.0))
            print(' mój punkt to:', x)
            print('zmieniam jego współrzędne na takie w układzie ortogonalnym')
            print(str(ref_change @ x))
            print('obracam go wstecz od d', str(d), 'do poprzedniego kierunku o:', str(o))
            print(str(np.linalg.inv(rot_m) @ ref_change @ x))
            print('obracam go wokół osi wokół poprzedniego kierunku o:', o)
            print(str(self.transformation @ np.linalg.inv(rot_m) @ ref_change @ x))
            print('obracam go od o', str(o), 'do nowego kierunku d:', str(d))
            print(str(rot_m @ self.transformation @ np.linalg.inv(rot_m) @ ref_change @ x))
            print('zmieniam współrzędne spowrotem na heksagonalne; mój punkt to:')
            print(str(np.linalg.inv(ref_change) @ rot_m @ self.transformation @ np.linalg.inv(rot_m) @ ref_change @ x))
            print('_'*40)

            return SymmOp(np.linalg.inv(ref_change) @ rot_m @ ref_change @ self.transformation @ np.linalg.inv(ref_change) @ np.linalg.inv(rot_m) @ ref_change,
                          np.linalg.inv(ref_change) @ rot_m @ self.translation)

    def _transform_hkl(self, hkl):
        return self.reciprocal.transformation @ np.array(hkl)

    def _transform_xyz(self, xyz):
        pass

    def _extincts_hkl(self, hkl):
        hkl = np.array(hkl)
        cond1 = np.allclose(hkl, self._transform_hkl(hkl))
        cond2 = np.isclose(np.dot(hkl, self.glide) % 1, 0)
        return cond1 and cond2


if __name__ == '__main__':
    #t = SymmOp.from_code('x-y, x, z+2/3')
    t = SymmOp.from_code('-x,-y,z').into2((1, 1, 0), hexagonal=True).into2((1, 0, 0), hexagonal=True)
    #t = SymmOp.from_code('y,x,-z').into((0, 0, 1), hexagonal=True)
    print(t.invariants)
    print(t.matrix4x4)
    print(t.det)
    print(t, t.glide)
    print('----------')
    print(t.typ)
    print(t.orientation)

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
