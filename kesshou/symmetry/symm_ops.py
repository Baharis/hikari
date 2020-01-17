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


symm_ops4 = dict()
symm_ops4['21_y'] = np.array([[-1, 0, 0, 0],
                              [0, 1, 0, 1/2],
                              [0, 0, -1, 0],
                              [0, 0, 0, 1]])
symm_ops4['c_y'] = np.array([[1, 0, 0, 0],
                             [0, -1, 0, 0],
                             [0, 0, 1, 1/2],
                             [0, 0, 0, 1]])


class Extinction:
    """This class tests symmetry operations and returns extinction conditions"""
    def __init__(self):
        testers = dict()
        pass

    def generate_tester(self, id, extinct, tester):
        pass

    def test(self, equivalent_points):
        pass

    extinct_refs = dict()
    extinct_refs['h = 2n @ h00'] = [(1, 0, 0), (3, 0, 0), (5, 0, 0)]
    extinct_refs['k = 2n @ 0k0'] = [(0, 1, 0), (0, 3, 0), (0, 5, 0)]
    extinct_refs['l = 2n @ 00l'] = [(0, 0, 1), (0, 0, 3), (0, 0, 5)]
    extinct_refs['h + k = 2n @ hk0'] = [(3, 1, 0), (6, 4, 0), (7, 1, 0)]
    extinct_refs['h + l = 2n @ h0k'] = [(3, 1, 0), (6, 4, 0), (7, 1, 0)]
    extinct_refs['h + k = 2n @ 0kl'] = [(3, 1, 0), (6, 4, 0), (7, 1, 0)]
