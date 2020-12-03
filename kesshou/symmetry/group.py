"""
This file contains class definition and necessary tools for constructing
and evaluating all symmetry groups.
"""
import numpy as np
from itertools import product as itertools_product
from enum import Enum
from kesshou.symmetry.operations import SymmOp


class Group:
    """
    Base class containing information about symmetry groups.
    It stores information for point and space groups and, among others,
    allows for iteration over its elements from `kesshou.symmetry.SymmOp`.
    """

    class CrystalSystem(Enum):
        triclinic = 0
        monoclinic = 1
        orthorhombic = 2
        trigonal = 3
        tetragonal = 4
        cubic = 5
        hexagonal = 6

    def __init__(self, *generators):
        """
        :param generators: List of operations necessary to construct whole group
        :type generators: List[SymmOp]
        """
        self.generators = generators

    def __iter__(self):
        return iter(self.operations)

    def __str__(self):
        s = 'A centrosymmetric ' if self.is_centrosymmetric else 'A '
        s += 'enantiomorphic ' if self.is_enantiomorphic else ''
        s += 'polar ' if self.is_polar else ''
        return s + 'group of order {}.'.format(self.order)

    @property
    def generators(self):
        return self.__generators

    @generators.setter
    def generators(self, new_generators):

        generator_list = []
        for gen in new_generators:
            if gen % 1 not in generator_list:
                generator_list.append(gen % 1)

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

        self.__generators = generator_list
        self.__operations = _find_new_product(generator_list)

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
    def is_enantiomorphic(self):
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
        """
        :return: Predicted crystal system associated with this group
        :rtype: self.CrystalSystem()
        """
        folds = [op.fold for op in self.operations]
        orients = [op.orientation for op in self.operations]
        if 6 in folds:
            return self.CrystalSystem.hexagonal
        elif 3 in folds:
            orients_of_3 = len({o for f, o in zip(folds, orients) if f == 3})
            return self.CrystalSystem.cubic if orients_of_3 > 1 \
                else self.CrystalSystem.trigonal
        elif 4 in folds:
            return self.CrystalSystem.tetragonal
        elif 2 in folds:
            orients_of_2 = len({o for f, o in zip(folds, orients) if f == 2})
            return self.CrystalSystem.orthorhombic if orients_of_2 > 1 \
                else self.CrystalSystem.monoclinic
        else:
            return self.CrystalSystem.triclinic

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
