import numpy as np
import numpy.linalg as lin
from .group import Group
from .symm_ops import symm_ops


class PointGroup(Group):
    """Basic Point Group class info holder"""

    unique_point = np.array([0.0 + np.pi / 100,
                             0.1 + np.pi / 100,
                             0.2 + np.pi / 100])

    def __init__(self, generators):
        super().__init__(generators=generators)

    @property
    def hp_disc_symm_ops(self):
        hp_disc_operations = [op for op in self.operations
                              if lin.det(op) > 0 and np.trace(op) < 3]
        return hp_disc_operations

    @property
    def is_centrosymmetric(self):
        return any(np.array_equal(op, symm_ops['-1']) for op in self.operations)

    def lauefy(self):
        self.generators.append(symm_ops['-1'])
        self.construct()


# TRICLINIC
PG1 = PointGroup(generators=[symm_ops['1']])
PG_1 = PointGroup(generators=[symm_ops['-1']])
# MONOCLINIC
PG2 = PointGroup(generators=[symm_ops['2_y']])
PGm = PointGroup(generators=[symm_ops['m_y']])
PG2pm = PointGroup(generators=[symm_ops['2_y'], symm_ops['m_y']])
# ORTHORHOMBIC
PG222 = PointGroup(generators=[symm_ops['2_x'], symm_ops['2_y'], symm_ops['2_z']])
PGmm2 = PointGroup(generators=[symm_ops['m_x'], symm_ops['m_y'], symm_ops['2_z']])
PGmmm = PointGroup(generators=[symm_ops['m_x'], symm_ops['m_y'], symm_ops['m_z']])
# TETRAGONAL
PG4 = PointGroup(generators=[symm_ops['4_z']])
PG_4 = PointGroup(generators=[symm_ops['-4_z']])
PG4pm = PointGroup(generators=[symm_ops['4_z'], symm_ops['m_z']])
PG422 = PointGroup(generators=[symm_ops['4_z'], symm_ops['2_x'], symm_ops['2_xy']])
PG4mm = PointGroup(generators=[symm_ops['4_z'], symm_ops['m_x'], symm_ops['m_xy']])
PG_42m = PointGroup(generators=[symm_ops['-4_z'], symm_ops['2_x'], symm_ops['m_xy']])
PG_4m2 = PointGroup(generators=[symm_ops['-4_z'], symm_ops['m_x'], symm_ops['2_xy']])
PG4pmmm = PointGroup(generators=[symm_ops['4_z'], symm_ops['m_z'], symm_ops['m_x'], symm_ops['m_xy']])
# TRIGONAL
PG3 = PointGroup(generators=[symm_ops['h3_z']])
PG_3 = PointGroup(generators=[symm_ops['h-3_z']])
PG321 = PointGroup(generators=[symm_ops['h3_z'], symm_ops['h2_x2y']])
PG312 = PointGroup(generators=[symm_ops['h3_z'], symm_ops['h2_x2y']])
PG3m1 = PointGroup(generators=[symm_ops['h3_z'], symm_ops['hm_x']])
PG31m = PointGroup(generators=[symm_ops['h3_z'], symm_ops['hm_x2y']])
PG_3m1 = PointGroup(generators=[symm_ops['h-3_z'], symm_ops['hm_x']])
PG_31m = PointGroup(generators=[symm_ops['h-3_z'], symm_ops['hm_x2y']])
# HEXAGONAL
PG6 = PointGroup(generators=[symm_ops['h6_z']])
PG_6 = PointGroup(generators=[symm_ops['h-6_z']])
PG6pm = PointGroup(generators=[symm_ops['h6_z'], symm_ops['hm_z']])
PG622 = PointGroup(generators=[symm_ops['h6_z'], symm_ops['h2_x'], symm_ops['h2_x2y']])
PG6mm = PointGroup(generators=[symm_ops['h6_z'], symm_ops['hm_x'], symm_ops['hm_x2y']])
PG_6m2 = PointGroup(generators=[symm_ops['h-6_z'], symm_ops['hm_x'], symm_ops['h2_x2y']])
PG_62m = PointGroup(generators=[symm_ops['h-6_z'], symm_ops['h2_x'], symm_ops['hm_x2y']])
PG6pmmm = PointGroup(generators=[symm_ops['h6_z'], symm_ops['hm_z'], symm_ops['hm_x'], symm_ops['hm_x2y']])
# CUBIC
PG23 = PointGroup(generators=[symm_ops['2_z'], symm_ops['3_xyz']])
PGm_3 = PointGroup(generators=[symm_ops['m_z'], symm_ops['-3_xyz']])
PG432 = PointGroup(generators=[symm_ops['4_z'], symm_ops['3_xyz'], symm_ops['2_xy']])
PG_43m = PointGroup(generators=[symm_ops['-4_z'], symm_ops['3_xyz'], symm_ops['m_xy']])
PGm_3m = PointGroup(generators=[symm_ops['m_z'], symm_ops['-3_xyz'], symm_ops['m_xy']])


if __name__ == '__main__':
    print(PG3.hp_disc_symm_ops)
