import numpy as np
import numpy.linalg as lin
from itertools import product as itertools_product

symm_ops = dict()
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


class PointGroup:
    """Basic Point Group class info holder"""

    def __init__(self, generators):
        self.operations = list()
        self.generators = generators
        self.operations = self.generate_group()

    def generate_group(self):
        """generate whole point group based on its generators"""
        operations = self.generators

        def _is_in(op, ops):
            return any([np.array_equal(op, op1) for op1 in ops])

        # then perform recursive group table multiplication
        def _find_new_product(current_ops):
            new_op = None
            for op1, op2 in itertools_product(current_ops, current_ops):
                op = np.dot(op1, op2)
                if not _is_in(op, current_ops):
                    new_op = op
                    break
            if new_op is None:
                return current_ops
            else:
                return _find_new_product([*current_ops, new_op])

        operations = _find_new_product(operations)
        return operations

    @property
    def hp_disc_transforming_symm_ops(self):
        hp_disc_operations = [op for op in self.operations
                              if lin.det(op) > 0 and np.trace(op) < 3]
        return hp_disc_operations

    @property
    def is_centrosymmetric(self):
        return any(np.array_equal(op, symm_ops['-1']) for op in self.operations)

    @property
    def is_enantiomorphic(self):
        return all(lin.det(op) > 0 for op in self.operations)

    @property
    def is_polar(self):
        polar_prop = np.array([1, 3, 7])
        sum_of_props = sum(np.dot(op, polar_prop) for op in self.operations)
        # print([op for op in self.operations])
        # print([np.dot(op, polar_prop) for op in self.operations])
        return lin.norm(sum_of_props) > 0.1

    def lauefy(self):
        self.generators.append(symm_ops['-1'])
        self.operations = self.generate_group()


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
    print(PG3.hp_disc_transforming_symm_ops)
