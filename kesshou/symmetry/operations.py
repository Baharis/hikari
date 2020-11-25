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
        return ','.join([self._row_to_str(xyz, r) for xyz, r
                         in zip(self.transformation, self.translation)])

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
            s += '{:+d}'.format(int(var))[0]
            if np.isclose(abs(var), 1):
                s += char
                continue
            s += str(Fraction(abs(var)).limit_denominator(9))
        if np.isclose(abs(r), 0):
            pass
        else:
            s += '{:+f}'.format(float(r))[0]
            s += str(Fraction(abs(r)).limit_denominator(9))
        return s.strip('+')

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
        return int(np.linalg.det(self.transformation))

    @property
    def typ(self):
        return SymmOpType((3+self.trace)/2 * (1 if self.translational else -1))

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

    def at(self, pt):
        shift = np.array(pt) - self.origin
        return SymmOp(self.transformation, self.glide + 2 * shift)

    # botched approach to try to conside that you cannot move operation
    # in invariant directions, works only for glide vectrors
    #
    # def at2(self, pt):
    #     shift = np.array(pt) - self.origin
    #     if np.linalg.norm(self.glide) > 0:
    #         glide_dir = self.glide / np.linalg.norm(self.glide)
    #         shift_perp_to_glide = shift - glide_dir * np.dot(glide_dir, shift)
    #         return SymmOp(self.transformation,
    #                       self.glide + 2 * shift_perp_to_glide)
    #     else:
    #         return SymmOp(self.transformation,
    #                       self.glide + 2 * shift)

    def _transform_hkl(self, hkl):
        return self.reciprocal.transformation @ np.array(hkl)

    def _extincts_hkl(self, hkl):
        hkl = np.array(hkl)
        cond1 = np.allclose(hkl, self.reciprocal.transformation @ hkl)
        cond2 = np.isclose(np.dot(hkl, self.glide) % 1, 0)
        return cond1, cond2


if __name__ == '__main__':
    t = SymmOp.from_code('-x, -y, z').at((1, 1, 0))
    print(t)
    print(t.matrix4x4)
    print(t.glide)
    print(t.origin)
    print(t.typ)


    print('t')
    print(t)
    print('-'*30)

    print('t.matrix4x4')
    print(t.matrix4x4)
    print('-'*30)

    print('t.order')
    print(t.order)
    print('-' * 30)

    print('t.fold')
    print(t.fold)
    print('-' * 30)

    print('t.glide')
    print(t.glide)
    print('-'*30)

    print('t.origin')
    print(t.origin)
    print('-'*30)

    print('t.at((2, 0, 0))')
    print(t.at((2, 0, 0)))
    print('-' * 30)

    print('t.reciprocal')
    print(t.reciprocal)
    print('-' * 30)

    print('t._transform_hkl((1, 0, 0))')
    print(t._transform_hkl((0, 0, 1)))
    print('-' * 30)

    print('t._extincts_hkl((0, 0, 1))')
    print(t._extincts_hkl((0, 0, 1)))
    print('-' * 30)

    print('t._extincts_hkl((0, 0, 4))')
    print(t._extincts_hkl((0, 0, 4)))
    print('-' * 30)

    print('t._extincts_hkl((0, 1, 4))')
    print(t._extincts_hkl((0, 1, 4)))
    print('-'*30)

    print(t * np.array(((1.6, 1.8, 2.2), (1.6, 1.8, 2.2))).T)
    print(t**2 * np.array((1.6, 1.8, 2.2)))
    print(t**3 * np.array((1.6, 1.8, 2.2)))
    print(t**4 * np.array((1.6, 1.8, 2.2)))
    print(t ** 5 * np.array((1.6, 1.8, 2.2)))