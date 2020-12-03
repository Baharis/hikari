"""
This file contains class definition and necessary tools for point groups.
"""

import copy
import numpy as np
from .group_old import Group
from .operations import symm_ops as so


class PointGroup(Group):
    """
    Class containing information about three dimensional point groups. It
    inherits a lots of properties after :class:`kesshou.symmetry.group.Group`.
    """

    unique_point = np.array([np.pi / 100, 0.1 + np.pi / 100, 0.2 + np.pi / 100])
    """Unique point in the space used for the sake of drawing etc."""

    def __init__(self, generators):
        super().__init__(generators=generators)

    @property
    def is_centrosymmetric(self):
        """
        Check whether operations of the space group contain centre of symmetry.

        :return: True if group has centre of symmetry; False otherwise.
        :rtype: bool
        """
        return any(np.array_equal(op, so['-1']) for op in self.operations)

    def lauefy(self):
        """
        Return itself with a centre of symmetry added to the set of generators.

        :return: New PointGroup with operations of self and centre of symmetry.
        :rtype: PointGroup
        """
        pg_with_added_inversion = copy.deepcopy(self)
        pg_with_added_inversion.generators.append(so['-1'])
        pg_with_added_inversion.construct()
        return pg_with_added_inversion


def _initiate_point_group_dictionary():
    """
    Define a dictionary of point groups PG.

    :return: A dictionary of known point groups with some non-standard settings.
    :rtype: Dict[str, kesshou.symmetry.pointgroup.PointGroup]
    """
    # TRICLINIC
    pg1 = PointGroup(generators=[so['1']])
    pg_1 = PointGroup(generators=[so['-1']])
    # MONOCLINIC
    pg2 = PointGroup(generators=[so['2_y']])
    pgm = PointGroup(generators=[so['m_y']])
    pg2om = PointGroup(generators=[so['2_y'], so['m_y']])
    # ORTHORHOMBIC
    pg222 = PointGroup(generators=[so['2_x'], so['2_y'], so['2_z']])
    pgmm2 = PointGroup(generators=[so['m_x'], so['m_y'], so['2_z']])
    pgmmm = PointGroup(generators=[so['m_x'], so['m_y'], so['m_z']])
    # TETRAGONAL
    pg4 = PointGroup(generators=[so['4_z']])
    pg_4 = PointGroup(generators=[so['-4_z']])
    pg4om = PointGroup(generators=[so['4_z'], so['m_z']])
    pg422 = PointGroup(generators=[so['4_z'], so['2_x'], so['2_xy']])
    pg4mm = PointGroup(generators=[so['4_z'], so['m_x'], so['m_xy']])
    pg_42m = PointGroup(generators=[so['-4_z'], so['2_x'], so['m_xy']])
    pg_4m2 = PointGroup(generators=[so['-4_z'], so['m_x'], so['2_xy']])
    pg4ommm = PointGroup([so['4_z'], so['m_z'], so['m_x'], so['m_xy']])
    # TRIGONAL
    pg3 = PointGroup(generators=[so['h3_z']])
    pg_3 = PointGroup(generators=[so['h-3_z']])
    pg321 = PointGroup(generators=[so['h3_z'], so['h2_x2y']])
    pg312 = PointGroup(generators=[so['h3_z'], so['h2_x2y']])
    pg3m1 = PointGroup(generators=[so['h3_z'], so['hm_x']])
    pg31m = PointGroup(generators=[so['h3_z'], so['hm_x2y']])
    pg_3m1 = PointGroup(generators=[so['h-3_z'], so['hm_x']])
    pg_31m = PointGroup(generators=[so['h-3_z'], so['hm_x2y']])
    # HEXAGONAL
    pg6 = PointGroup(generators=[so['h6_z']])
    pg_6 = PointGroup(generators=[so['h-6_z']])
    pg6om = PointGroup(generators=[so['h6_z'], so['hm_z']])
    pg622 = PointGroup(generators=[so['h6_z'], so['h2_x'], so['h2_x2y']])
    pg6mm = PointGroup(generators=[so['h6_z'], so['hm_x'], so['hm_x2y']])
    pg_6m2 = PointGroup(generators=[so['h-6_z'], so['hm_x'], so['h2_x2y']])
    pg_62m = PointGroup(generators=[so['h-6_z'], so['h2_x'], so['hm_x2y']])
    pg6ommm = PointGroup([so['h6_z'], so['hm_z'], so['hm_x'], so['hm_x2y']])
    # CUBIC
    pg23 = PointGroup(generators=[so['2_z'], so['3_xyz']])
    pgm_3 = PointGroup(generators=[so['m_z'], so['-3_xyz']])
    pg432 = PointGroup(generators=[so['4_z'], so['3_xyz'], so['2_xy']])
    pg_43m = PointGroup(generators=[so['-4_z'], so['3_xyz'], so['m_xy']])
    pgm_3m = PointGroup(generators=[so['m_z'], so['-3_xyz'], so['m_xy']])
    pg_dict = {'1': pg1, '-1': pg_1,
               '2': pg2, 'm': pgm, '2/m': pg2om,
               '222': pg222, 'mm2': pgmm2, 'mmm': pgmmm,
               '4': pg4, '-4': pg_4, '4/m': pg4om, '422': pg422,
               '4mm': pg4mm, '-42m': pg_42m, '-4m2': pg_4m2, '4/mmm': pg4ommm,
               '3': pg3, '-3': pg_3, '321': pg321, '312': pg312,
               '3m1': pg3m1, '31m': pg31m, '-3m1': pg_3m1, '-31m': pg_31m,
               '6': pg6, '-6': pg_6, '6/m': pg6om, '622': pg622,
               '6mm': pg6mm, '-6m2': pg_6m2, '-62m': pg_62m, '6/mmm': pg6ommm,
               '23': pg23, 'm-3': pgm_3, '432': pg432,
               '-43m': pg_43m, 'm-3m': pgm_3m}
    return pg_dict


PG = _initiate_point_group_dictionary()
"""
Dictionary containing all known point groups written as :class:`PointGroup`
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
