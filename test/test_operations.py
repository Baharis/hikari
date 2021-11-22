import unittest
from hikari.symmetry import SymmOp
import numpy as np


class TestSymmOp(unittest.TestCase):
    sg40_a_tf = np.array([[1, 0, 0], [0, -1, 0], [0, 0, 1]])
    sg40_a_tl = np.array([1 / 2, 0, 0])
    sg40_a_m = np.array([[1, 0, 0, .5], [0, -1, 0, 0],
                         [0, 0, 1, 0], [0, 0, 0, 1]])
    sg40_a_code = '1/2+x,-y,z'
    sg40_21_code = '-x,-y+1/2,z+1/2'
    sg40_n_code = '-x+1/2,y-1/2,z+1/2'

    def test_from_pair(self):
        o = SymmOp.from_pair(self.sg40_a_tf, self.sg40_a_tl)
        self.assertTrue(np.allclose(self.sg40_a_tf, o.tf))
        self.assertTrue(np.allclose(self.sg40_a_tl, o.tl))

    def test_from_code(self):
        o = SymmOp.from_code(self.sg40_a_code)
        self.assertTrue(np.allclose(self.sg40_a_tf, o.tf))
        self.assertTrue(np.allclose(self.sg40_a_tl, o.tl))

    def test_from_matrix(self):
        o = SymmOp.from_matrix(self.sg40_a_m)
        self.assertTrue(np.allclose(self.sg40_a_tf, o.tf))
        self.assertTrue(np.allclose(self.sg40_a_tl, o.tl))

    def test_equal(self):
        o1 = SymmOp.from_code(self.sg40_a_code)
        o2 = SymmOp.from_matrix(self.sg40_a_m)
        self.assertEqual(o1, o2)

    def test_mul_and_mod(self):
        o1 = SymmOp.from_code(self.sg40_a_code)
        o2 = SymmOp.from_code(self.sg40_21_code)
        o3 = SymmOp.from_code(self.sg40_n_code)
        self.assertEqual(o1 * o2 % 1, o3 % 1)


if __name__ == '__main__':
    unittest.main()
