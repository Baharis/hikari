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


PG = {
    # TRICLINIC
    '1': Group(SymmOp.from_code('x,y,z')),
    '-1': Group(SymmOp.from_code('-x,-y,-z')),
    # MONOCLINIC
    '2': Group(SymmOp.from_code('-x,y,-z')),
    'm': Group(SymmOp.from_code('x,-y,z')),
    '2/m': Group(SymmOp.from_code('-x,-y,-z'), SymmOp.from_code('x,-y,z')),
    # ORTHORHOMBIC
    '222': Group(SymmOp.from_code('x,-y,-z'), SymmOp.from_code('-x,y,-z')),
    'mm2': Group(SymmOp.from_code('-x,y,z'), SymmOp.from_code('-x,-y,z')),
    'mmm': Group(SymmOp.from_code('-x,y,z'), SymmOp.from_code('x,-y,z'),
                 SymmOp.from_code('x,y,-z')),
    # TETRAGONAL
    '4': Group(SymmOp.from_code('-y,x,z')),
    '-4': Group(SymmOp.from_code('y,-x,-z')),
    '4/m': Group(SymmOp.from_code('-y,x,z'), SymmOp.from_code('x,y,-z')),
    '422': Group(SymmOp.from_code('-y,x,z'), SymmOp.from_code('-x,y,-z')),
    '4mm': Group(SymmOp.from_code('-y,x,z'), SymmOp.from_code('x,-y,z')),
    '-42m': Group(SymmOp.from_code('y,-x,-z'), SymmOp.from_code('x,-y,-z')),
    '-4m2': Group(SymmOp.from_code('y,-x,-z'), SymmOp.from_code('-x,y,z')),
    '4/mmm': Group(SymmOp.from_code('-y,x,z'), SymmOp.from_code('x,y,-z'),
                   SymmOp.from_code('-x,y,z')),
    # TRIGONAL
    '3': Group(SymmOp.from_code('-y,x-y,z')),
    '-3': Group(SymmOp.from_code('y,-x+y,-z')),
    '321': Group(SymmOp.from_code('-y,x-y,z'), SymmOp.from_code('y,x,-z')),
    '312': Group(SymmOp.from_code('-y,x-y,z'), SymmOp.from_code('-y,-x,-z')),
    '3m1': Group(SymmOp.from_code('-y,x-y,z'), SymmOp.from_code('-y,-x,z')),
    '31m': Group(SymmOp.from_code('-y,x-y,z'), SymmOp.from_code('y,x,z')),
    '-3m1': Group(SymmOp.from_code('y,-x+y,-z'), SymmOp.from_code('-y,-x,z')),
    '-31m': Group(SymmOp.from_code('y,-x+y,-z'), SymmOp.from_code('y,x,z')),
    # HEXAGONAL
    '6': Group(SymmOp.from_code('x-y,x,z')),
    '-6': Group(SymmOp.from_code('-x+y,-x,-z')),
    '6/m': Group(SymmOp.from_code('x-y,x,z'), SymmOp.from_code('x,y,-z')),
    '622': Group(SymmOp.from_code('x-y,x,z'), SymmOp.from_code('x-y,-y,-z')),
    '6mm': Group(SymmOp.from_code('x-y,x,z'), SymmOp.from_code('-x+y,y,z')),
    '-6m2': Group(SymmOp.from_code('-x+y,-x,-z'), SymmOp.from_code('-x+y,y,z')),
    '-62m': Group(SymmOp.from_code('-x+y,-x,-z'), SymmOp.from_code('y,x,-z')),
    '6/mmm': Group(SymmOp.from_code('x-y,x,z'), SymmOp.from_code('x,y,-z'),
                   SymmOp.from_code('-x+y,y,z')),
    # CUBIC
    '23': Group(SymmOp.from_code('-x,-y,z'), SymmOp.from_code('z,x,y')),
    'm-3': Group(SymmOp.from_code('x,y,-z'), SymmOp.from_code('-z,-x,-y')),
    '432': Group(SymmOp.from_code('-y,x,z'), SymmOp.from_code('z,x,y')),
    '-43m': Group(SymmOp.from_code('y,-x,-z'), SymmOp.from_code('z,x,y')),
    'm-3m': Group(SymmOp.from_code('x,y,-z'), SymmOp.from_code('-z,-x,-y'),
                  SymmOp.from_code('-y,-x,z'))}
"""
Dictionary containing all known point groups written as :class:`Group`
along with alternative axis settings. The point groups in this dictionary
can be accessed using their short Hermann-Maugin notation, as presented below.

+-------+---------------+----------------+---------------+-----------------+
| No.   | CRYSTAL       | Hermann-Maugin | Schoenflies   | Can be accessed |
|       | SYSTEM        | notation       | notation      | using           |
+-------+---------------+----------------+---------------+-----------------+
| 1     | triclinic     | 1              | C1            | `PG['1']`       |
+-------+               +----------------+---------------+-----------------+
| 2     |               | -1             | Ci            | `PG['-1']`      |
+-------+---------------+----------------+---------------+-----------------+
| 3     | monoclinic    | 2              | C2            | `PG['2']`       |
+-------+               +----------------+---------------+-----------------+
| 4     |               | m              | Cs            | `PG['m']`       |
+-------+               +----------------+---------------+-----------------+
| 5     |               | 2/m            | C2h           | `PG['2/m']`     |
+-------+---------------+----------------+---------------+-----------------+
| 6     | orthorhombic  | 222            | D2            | `PG['222']`     |
+-------+               +----------------+---------------+-----------------+
| 7     |               | mm2            | C2v           | `PG['mm2']`     |
+-------+               +----------------+---------------+-----------------+
| 8     |               | mmm            | D2h           | `PG['mmm']`     |
+-------+---------------+----------------+---------------+-----------------+
| 9     | tetragonal    | 4              | C4            | `PG['4']`       |
+-------+               +----------------+---------------+-----------------+
| 10    |               | -4             | S4            | `PG['-4']`      |
+-------+               +----------------+---------------+-----------------+
| 11    |               | 4/m            | C4h           | `PG['4/m']`     |
+-------+               +----------------+---------------+-----------------+
| 12    |               | 422            | D4            | `PG['422']`     |
+-------+               +----------------+---------------+-----------------+
| 13    |               | 4mm            | C4v           | `PG['4mm']`     |
+-------+               +----------------+---------------+-----------------+
| 14    |               | -42m           | D2d           | `PG['-42m']`    |
+-------+               +----------------+---------------+-----------------+
| 14*   |               | -4m2           | D2d           | `PG['-4m2']`    |
+-------+               +----------------+---------------+-----------------+
| 15    |               | 4/mmm          | D4h           | `PG['4/mmm']`   |
+-------+---------------+----------------+---------------+-----------------+
| 16    | trigonal      | 3              | C3            | `PG['3']`       |
+-------+               +----------------+---------------+-----------------+
| 17    |               | -3             | C3i           | `PG['-3']`      |
+-------+               +----------------+---------------+-----------------+
| 18    |               | 32             | D3            | `PG['321']`     |
+-------+               +----------------+---------------+-----------------+
| 18*   |               | 32             | D3            | `PG['312']`     |
+-------+               +----------------+---------------+-----------------+
| 19    |               | 3m             | C3v           | `PG['3m1']`     |
+-------+               +----------------+---------------+-----------------+
| 19*   |               | 3m             | C3v           | `PG['31m']`     |
+-------+               +----------------+---------------+-----------------+
| 20    |               | -3m            | D3d           | `PG['-3m1']`    |
+-------+               +----------------+---------------+-----------------+
| 20*   |               | -3m            | D3d           | `PG['-31m']`    |
+-------+---------------+----------------+---------------+-----------------+
| 21    | hexagonal     | 6              | C6            | `PG['6']`       |
+-------+               +----------------+---------------+-----------------+
| 22    |               | -6             | C3h           | `PG['-6']`      |
+-------+               +----------------+---------------+-----------------+
| 23    |               | 6/m            | C6h           | `PG['6/m']`     |
+-------+               +----------------+---------------+-----------------+
| 24    |               | 622            | D6            | `PG['622']`     |
+-------+               +----------------+---------------+-----------------+
| 25    |               | 6mm            | C6v           | `PG['6mm']`     |
+-------+               +----------------+---------------+-----------------+
| 26    |               | -6m2           | D3h           | `PG['-6m2']`    |
+-------+               +----------------+---------------+-----------------+
| 26*   |               | -62m           | D3h           | `PG['-62m']`    |
+-------+               +----------------+---------------+-----------------+
| 27    |               | 6/mmm          | D6h           | `PG['6/mmm']`   |
+-------+---------------+----------------+---------------+-----------------+
| 28    | cubic         | 23             | T             | `PG['23']`      |
+-------+               +----------------+---------------+-----------------+
| 29    |               | m-3            | Th            | `PG['m-3']`     |
+-------+               +----------------+---------------+-----------------+
| 30    |               | 432            | O             | `PG['432']`     |
+-------+               +----------------+---------------+-----------------+
| 31    |               | -43m           | Td            | `PG['-43m']`    |
+-------+               +----------------+---------------+-----------------+
| 32    |               | m-3m           | Oh            | `PG['m-3m']`    |
+-------+---------------+----------------+---------------+-----------------+

Asterisk (*) denotes alternative choice of axes.
"""


if __name__ == '__main__':
    o1 = SymmOp.from_code('x+1/2, y+1/2, z')
    o2 = SymmOp.from_code('x+1/2, y, z+1/2')
    o3 = SymmOp.from_code('z, x, y')
    o4 = SymmOp.from_code('x-y, x, z-1/6')
    g = Group(o4)
    h = g.reciprocate()
    # h = Group(o4).lauefy()
    [print(op, op.orientation) for op in g]
    print(str(g))
    print('-' * 40)
    [print(op, op.orientation) for op in h]
    print(str(h))
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