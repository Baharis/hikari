"""
This file contains class definition and necessary tools for constructing
and evaluating all symmetry groups.
"""
import numpy as np
from pathlib import Path
from itertools import product as itertools_product
from enum import Enum
from hikari.symmetry.operations import SymmOp
from hikari.utility.list_tools import find_best
import json
import pickle


def unpack_group_dictionary_from_csv(filename):
    """Development function used to get PG and SG from csv to pickle it later"""
    path = Path(__file__).parent.absolute().joinpath(filename)
    with open(path) as file:
        json_dict = json.load(file)
    group_dict = {}
    for json_key, json_group in json_dict.items():
        g_name = json_group["H-M_short"]
        g_number = json_group["number"]
        g_gens = [SymmOp.from_code(g) for g in json_group["generators"]]
        g_ops = [SymmOp.from_code(o) for o in json_group["operations"]]
        g = Group.create_manually(generators=g_gens, operations=g_ops)
        g.name = json_group["H-M_short"]
        g.number = abs(g_number)
        group_dict[json_key] = g
        if g_number > 0:
            group_dict[g_name] = g
            group_dict[g_number] = g
    return group_dict


def unpack_group_dictionary_from_pickle(filename):
    """Function used to unpack point and space group pickles into their dicts"""
    path = Path(__file__).parent.absolute().joinpath(filename)
    return pickle.load(open(path, 'rb'))


def pack_group_dictionary_to_pickle(group_dict, filename):
    """Development function used to obtain point and space group pickles"""
    pickle.dump(group_dict, open(filename, 'wb'), protocol=4)


class Group:
    """
    Base immutable class containing information about symmetry groups.
    It stores information for point and space groups and, among others,
    allows for iteration over its elements from `hikari.symmetry.SymmOp`.
    """

    class System(Enum):
        """Enumerator class with information about associated crystal system"""
        triclinic = 0
        monoclinic = 1
        orthorhombic = 2
        trigonal = 3
        tetragonal = 4
        cubic = 5
        hexagonal = 6

        @property
        def directions(self):
            _a = np.array((1, 0, 0))
            _b = np.array((0, 1, 0))
            _c = np.array((0, 0, 1))
            _ab = np.array((1 / np.sqrt(2), 1 / np.sqrt(2), 0))
            _abc = np.array((1 / np.sqrt(3), 1 / np.sqrt(3), 1 / np.sqrt(3)))
            return [(), (_b, ), (_a, _b, _c), (_c, _a, _ab),
                    (_c, _a, _ab), (_c, _abc, _ab), (_c, _a, _ab)][self.value]

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
            new = list({o1 * o2 % 1 for o1, o2 in itertools_product(ops, ops)})
            new = set(ops).union(new)
            return _find_new_product(new) if len(new) > len(ops) else ops

        self.__generators = tuple(generator_list)
        self.__operations = tuple(_find_new_product(generator_list))
        self.name = self.auto_generated_name
        self.number = 0

    @classmethod
    def create_manually(cls, generators, operations):
        """
        Generate group using already complete list of generators and operators.
        :param generators: A complete list of group generators
        :type generators: List[np.ndarray]
        :param operations: A complete list of group operations
        :type operations: List[np.ndarray]
        :return:
        :rtype:
        """
        new_group = cls()
        new_group.__generators = generators
        new_group.__operations = operations
        return new_group

    def __iter__(self):
        return iter(self.operations)

    def __str__(self):
        s = 'A centrosymmetric ' if self.is_centrosymmetric else 'A '
        s += 'Sohncke ' if self.is_sohncke else ''
        s += 'polar ' if self.is_polar else ''
        return s + 'group of order {}.'.format(self.order)

    def __hash__(self):
        return sum(hash(o) for o in self.operations)

    @property
    def auto_generated_name(self):
        """Name of the group generated automatically. Use only as approx."""
        # TODO: enantiomorphs like P41 / P43 not recognised
        # TODO: 'e' found always whenever 'a' and 'b' present
        # TODO: some mistakes occur in trigonal crystal system (see SG149+)
        tl = ([o.name for o in self.operations if o.typ is o.Type.translation])
        tl.append('H' if self.system is self.System.trigonal else 'P')
        name = find_best(tl, 'A+B+C=F>I>C>B>A>H>P')
        for d in self.system.directions:
            ops = [o.name.partition(':')[0] for o in self.operations
                   if o.orientation is not None and
                   np.isclose(abs(np.dot(np.abs(o.orientation), np.abs(d))), 1)]
            axis_rules = '6>61>62>63>64>65>-6>4>41>42>43>-4>-3>3>31>32>2>21'
            plane_rules = 'm>a+b=e>a+c=e>b+c=e>a>b>c>n>d'
            best_axis = find_best(ops, axis_rules)
            best_plane = find_best(ops, plane_rules)
            sep = '/' if len(best_axis) > 0 and len(best_plane) > 0 else ''
            name += ' ' + best_axis + sep + best_plane
        return name

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
    def is_centrosymmetric(self):
        """
        :return: True if group has centre of symmetry; False otherwise.
        :rtype: bool
        """
        return any(np.isclose(op.trace, -3) for op in self.operations)

    @property
    def is_enantiogenic(self):
        """
        :return: True if determinant of all operations in group are positive.
        :rtype: bool
        """
        return any(op.det < 0 for op in self.operations)

    @property
    def is_sohncke(self):
        """
        :return: True if determinant of all operations in group are positive.
        :rtype: bool
        """
        return all(op.det > 0 for op in self.operations)

    @property
    def is_achiral(self):      # TODO See dictionary.iucr.org/Chiral_space_group
        return NotImplemented  # TODO and dx.doi.org/10.1524/zkri.2006.221.1.1

    @property
    def is_chiral(self):       # TODO See dictionary.iucr.org/Chiral_space_group
        return NotImplemented  # TODO and dx.doi.org/10.1524/zkri.2006.221.1.1

    @property
    def is_symmorphic(self):
        return all(g.typ not in {g.Type.rototranslation, g.Type.transflection}
                   and sum(g.origin) == 0 for g in self.generators)

    @property
    def is_polar(self):
        if any(op.typ in {op.Type.rotoinversion, op.Type.inversion}
               for op in self.operations):
            return False
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

        def _is_many(_orients):
            return any([np.dot(_orients[0], o) < 0.99 for o in _orients[1:]])

        if 6 in folds:
            return self.System.hexagonal
        elif 3 in folds:
            orients_of_3 = [o for f, o in zip(folds, orients) if f == 3]
            return self.System.cubic if _is_many(orients_of_3) \
                else self.System.trigonal
        elif 4 in folds:
            return self.System.tetragonal
        elif 2 in folds:
            orients_of_2 = [o for f, o in zip(folds, orients) if f == 2]
            return self.System.orthorhombic if _is_many(orients_of_2) \
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

    def transform(self, m):
        """
        Transform the group using 4x4 matrix. For reference, see `bilbao
        resources <https://www.cryst.ehu.es/cryst/trmatrix.html>`_ or `IUCr
        pamphlet no. 22 <https://www.iucr.org/education/pamphlets/22>`_.

        :example:

        >>> import numpy
        >>> from hikari.symmetry import SG
        >>> m = numpy.array([(1,0,1,0),(0,1,0,0),(-1,0,0,0),(0,0,0,1)])
        >>> SG['P21/c'].transform(m).auto_generated_name
        P 21/n

        :param m: A 4x4 array containing information about new base and origin.
        :type m: np.ndarray
        :return: Group with new, transformed basis and origin.
        :rtype: Group
        """
        transformed_group = Group.create_manually(
            generators=[SymmOp.from_matrix(np.linalg.inv(m) @ g.matrix @ m)
                        for g in self.generators],
            operations=[SymmOp.from_matrix(np.linalg.inv(m) @ o.matrix @ m)
                    for o in self.operations])
        transformed_group.name = self.name + ' @ ' + str(m)
        transformed_group.number = -abs(self.number)
        return transformed_group
