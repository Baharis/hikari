import copy
import numpy as np
from .group import Group
from .symmetry_operations import symm_ops as so


class PointGroup(Group):
    """Basic Point Group class info holder"""

    unique_point = np.array([np.pi / 100, 0.1 + np.pi / 100, 0.2 + np.pi / 100])

    def __init__(self, generators):
        super().__init__(generators=generators)

    @property
    def is_centrosymmetric(self):
        return any(np.array_equal(op, so['-1']) for op in self.operations)

    def lauefy(self):
        pg_with_added_inversion = copy.deepcopy(self)
        pg_with_added_inversion.generators.append(so['-1'])
        pg_with_added_inversion.construct()
        return pg_with_added_inversion

def initiate_point_group_dictionary():
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


PG = initiate_point_group_dictionary()
