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
        self.generators.append(so['-1'])
        self.construct()


PG = dict()
# TRICLINIC
PG['1'] = PointGroup(generators=[so['1']])
PG['-1'] = PointGroup(generators=[so['-1']])
# MONOCLINIC
PG['2'] = PointGroup(generators=[so['2_y']])
PG['m'] = PointGroup(generators=[so['m_y']])
PG['2/m'] = PointGroup(generators=[so['2_y'], so['m_y']])
# ORTHORHOMBIC
PG['222'] = PointGroup(generators=[so['2_x'], so['2_y'], so['2_z']])
PG['mm2'] = PointGroup(generators=[so['m_x'], so['m_y'], so['2_z']])
PG['mmm'] = PointGroup(generators=[so['m_x'], so['m_y'], so['m_z']])
# TETRAGONAL
PG['4'] = PointGroup(generators=[so['4_z']])
PG['-4'] = PointGroup(generators=[so['-4_z']])
PG['4/m'] = PointGroup(generators=[so['4_z'], so['m_z']])
PG['422'] = PointGroup(generators=[so['4_z'], so['2_x'], so['2_xy']])
PG['4mm'] = PointGroup(generators=[so['4_z'], so['m_x'], so['m_xy']])
PG['-42m'] = PointGroup(generators=[so['-4_z'], so['2_x'], so['m_xy']])
PG['-4m2'] = PointGroup(generators=[so['-4_z'], so['m_x'], so['2_xy']])
PG['4/mmm'] = PointGroup([so['4_z'], so['m_z'], so['m_x'], so['m_xy']])
# TRIGONAL
PG['3'] = PointGroup(generators=[so['h3_z']])
PG['-3'] = PointGroup(generators=[so['h-3_z']])
PG['321'] = PointGroup(generators=[so['h3_z'], so['h2_x2y']])
PG['312'] = PointGroup(generators=[so['h3_z'], so['h2_x2y']])
PG['3m1'] = PointGroup(generators=[so['h3_z'], so['hm_x']])
PG['31m'] = PointGroup(generators=[so['h3_z'], so['hm_x2y']])
PG['-3m1'] = PointGroup(generators=[so['h-3_z'], so['hm_x']])
PG['-31m'] = PointGroup(generators=[so['h-3_z'], so['hm_x2y']])
# HEXAGONAL
PG['6'] = PointGroup(generators=[so['h6_z']])
PG['-6'] = PointGroup(generators=[so['h-6_z']])
PG['6/m'] = PointGroup(generators=[so['h6_z'], so['hm_z']])
PG['622'] = PointGroup(generators=[so['h6_z'], so['h2_x'], so['h2_x2y']])
PG['6mm'] = PointGroup(generators=[so['h6_z'], so['hm_x'], so['hm_x2y']])
PG['-6m2'] = PointGroup(generators=[so['h-6_z'], so['hm_x'], so['h2_x2y']])
PG['-62m'] = PointGroup(generators=[so['h-6_z'], so['h2_x'], so['hm_x2y']])
PG['6/mmm'] = PointGroup([so['h6_z'], so['hm_z'], so['hm_x'], so['hm_x2y']])
# CUBIC
PG['23'] = PointGroup(generators=[so['2_z'], so['3_xyz']])
PG['m-3'] = PointGroup(generators=[so['m_z'], so['-3_xyz']])
PG['432'] = PointGroup(generators=[so['4_z'], so['3_xyz'], so['2_xy']])
PG['-43m'] = PointGroup(generators=[so['-4_z'], so['3_xyz'], so['m_xy']])
PG['m-3m'] = PointGroup(generators=[so['m_z'], so['-3_xyz'], so['m_xy']])
