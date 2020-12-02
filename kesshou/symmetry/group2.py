"""
This file contains class definition and necessary tools for constructing
and evaluating all symmetry groups.
"""
import numpy as np
from itertools import product as itertools_product
from enum import Enum
from kesshou.symmetry.operations import SymmOp


class CrystalSystem(Enum):
    triclinic = 0
    monoclinic = 1
    orthorhombic = 2
    trigonal = 3
    tetragonal = 4
    cubic = 5
    hexagonal = 6


class Group:
    """
    Base class containing information about symmetry groups.
    It stores information for point and space groups.

    In the future it is planned to expand this group using getter/setter
    operations to automatically refresh operations and symmetry equivalents
    whenever a list of generators is changed.
    """

    def __init__(self, *generators):
        """
        :param generators: List of operations necessary to construct whole group
        :type generators: List[SymmOp]
        """
        self.generators = generators

    @property
    def generators(self):
        return self.__generators

    @generators.setter
    def generators(self, new_generators):

        self.__generators = []
        for gen in new_generators:
            if gen not in self.__generators:
                self.__generators.append(gen)

        def _find_new_product(current_ops):
            if len(current_ops) > 200:
                raise ValueError('Generated group order exceeds size of 200')
            new_ops = []
            for op1, op2 in itertools_product(current_ops, current_ops):
                new_op = op1 * op2 % 1
                if new_op not in [*current_ops, *new_ops]:
                    new_ops.append(new_op)
            return _find_new_product([*current_ops, *new_ops]) \
                if len(new_ops) > 0 else current_ops

        self.__operations = _find_new_product(new_generators)

    @property
    def operations(self):
        return self.__operations

    @property
    def order(self):
        return len(self.__operations)

    @property
    def chiral_operations(self):
        """
        :return: A subgroup with only operations whose determinant is positive.
        :rtype: list
        """
        return Group(*[op for op in self.operations if op.det > 0])

    @property
    def is_centrosymmetric(self):
        """
        :return: True if group has centre of symmetry; False otherwise.
        :rtype: bool
        """
        return any(np.isclose(op.trace, -3) for op in self.operations)

    @property
    def is_chiral(self):
        """
        :return: True if determinant of all operations in group are positive.
        :rtype: bool
        """
        return all(op.det > 0 for op in self.operations)
        # TODO: differentiate between chiral and centrosymmetric?

    @property
    def is_polar(self):
        """
        :return: True if group can preserve polar properties.
        :rtype: bool
        """
        # return False if group has any inversion or rotoinversion
        if any(op.orientation is None for op in self.operations if op.det < 0):
            return False
        # return False if group has rotations around non-parallel axes
        axes_orientation = [op.orientation for op in self.operations
                            if op.det > 0 and op.orientation is not None]
        for or1, or2 in itertools_product(axes_orientation, axes_orientation):
            if not np.isclose(abs(np.dot(or1, or2)), 1):
                return False
        return True

    @property
    def system(self):
        folds = [op.fold for op in self.operations]
        orients = [op.orientation for op in self.operations]
        if 6 in folds:
            return CrystalSystem.hexagonal
        elif 3 in folds:
            orients_of_3 = len({o for f, o in zip(folds, orients) if f == 3})
            return CrystalSystem.cubic if orients_of_3 > 1 \
                else CrystalSystem.trigonal
        elif 4 in folds:
            return CrystalSystem.tetragonal
        elif 2 in folds:
            orients_of_2 = len({o for f, o in zip(folds, orients) if f == 2})
            return CrystalSystem.orthorhombic if orients_of_2 > 1 \
                else CrystalSystem.monoclinic
        else:
            return CrystalSystem.triclinic

    def lauefy(self):
        """
        :return: New PointGroup with centre of symmetry added to generators.
        :rtype: Group
        """
        inv = SymmOp.from_code('-x, -y, -z')
        return Group(*self.operations, inv)

    def reciprocate(self):
        new_generators = [op.reciprocal for op in self.generators]
        return Group(*new_generators)


if __name__ == '__main__':
    o1 = SymmOp.from_code('-x, -y, z')
    o2 = SymmOp.from_code('-x, -y, -z')
    o3 = SymmOp.from_code('z, x, y')
    o4 = SymmOp.from_code('x-y, x, z+1/6')
    g = Group(o4)
    # h = Group(o4).lauefy()
    [print(g) for g in g.operations]
    # print(len([g for g in g.operations]))
    # print(g.is_centrosymmetric)
    # print(g.is_chiral)
    # print(g.is_polar)
    # print('-'*40)
    # [print(h) for h in h.operations]
    # print(len([h for h in h.operations]))
    # print(h.is_centrosymmetric)
    # print(h.is_chiral)
    # print(h.is_polar)