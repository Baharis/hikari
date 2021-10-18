"""
This file contains class definition and necessary tools for constructing
and evaluating all symmetry groups.
"""
import numpy as np
from pathlib import Path
from itertools import product as itertools_product
from enum import Enum
from hikari.symmetry.operations import SymmOp
import json
import pickle

# TODO add class method or creation method "from operators", without generating
# TODO and hard-code all space groups in separate json file. (pickle?)


def unpack_group_dict_from_csv(filename):
    path = Path(__file__).parent.absolute().joinpath(filename)
    with open(path) as file:
        json_dict = json.load(file)
    group_dict = {}
    for json_key, json_group in json_dict.items():
        sg_name = json_group["H-M_short"]
        sg_number = json_group["number"]
        sg_gens = [SymmOp.from_code(c) for c in json_group["generators"]]
        sg_ops = [SymmOp.from_code(o) for o in json_group["operations"]]
        sg_object = Group.create_manually(generators=sg_gens, operators=sg_ops)
        group_dict[json_key] = sg_object
        group_dict[sg_name] = sg_object
        group_dict[sg_number] = sg_object
    return group_dict


def _unpack_group_dict_from_pickle(filename):
    path = Path(__file__).parent.absolute().joinpath(filename)
    return pickle.load(open(path, 'rb'))


class Group:
    """
    Base immutable class containing information about symmetry groups.
    It stores information for point and space groups and, among others,
    allows for iteration over its elements from `hikari.symmetry.SymmOp`.
    """

    class System(Enum):
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

        generator_list = []
        for gen in generators:
            if gen % 1 not in generator_list:
                generator_list.append(gen % 1)

        def _find_new_product(ops):
            if len(ops) > 200:
                raise ValueError('Generated group order exceeds size of 200')
            new_ops = list({o1 * o2 % 1 for o1, o2 in itertools_product(ops, ops)})
            new_ops = set(ops).union(new_ops)
            # TODO this part gets slow for large groups; require "1" instead?
            # print([str(o) for o in new_ops])
            return _find_new_product(new_ops) if len(new_ops) > len(ops) else ops

        self.__generators = tuple(generator_list)
        self.__operations = tuple(_find_new_product(generator_list))

    @classmethod
    def create_manually(cls, generators, operators):
        new_group = cls()
        new_group.__generators = generators
        new_group.__operations = operators
        return new_group

    def __iter__(self):
        return iter(self.operations)

    def __str__(self):
        s = 'A centrosymmetric ' if self.is_centrosymmetric else 'A '
        s += 'enantiomorphic ' if self.is_enantiomorphic else ''
        s += 'polar ' if self.is_polar else ''
        return s + 'group of order {}.'.format(self.order)

    def __hash__(self):
        return sum(hash(o) for o in self.operations)

    @property
    def generators(self):
        return self.__generators

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

        def _many_orients_in(_orients):
            return any([np.dot(_orients[0], o) < 0.99 for o in _orients[1:]])

        if 6 in folds:
            return self.System.hexagonal
        elif 3 in folds:
            orients_of_3 = [o for f, o in zip(folds, orients) if f == 3]
            return self.System.cubic if _many_orients_in(orients_of_3) \
                else self.System.trigonal
        elif 4 in folds:
            return self.System.tetragonal
        elif 2 in folds:
            orients_of_2 = [o for f, o in zip(folds, orients) if f == 2]
            return self.System.orthorhombic if _many_orients_in(orients_of_2) \
                else self.System.monoclinic
        else:
            return self.System.triclinic

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
