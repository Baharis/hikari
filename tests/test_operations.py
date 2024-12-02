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
    sg147_m3_code = 'x - y, x, -z'

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

    def test_maths(self):
        o1 = SymmOp.from_code(self.sg40_a_code)
        o2 = SymmOp.from_code(self.sg40_21_code)
        o3 = SymmOp.from_code(self.sg40_n_code)
        self.assertEqual(o1 * o2 % 1, o3 % 1)
        self.assertEqual(o1*o1, o1**2)
        self.assertEqual(o1*o1*o1*o1*o1, o1**5)

    def test_repr(self):
        o1 = SymmOp.from_code(self.sg40_a_code)
        o2 = SymmOp.from_code(self.sg40_21_code)
        o3 = SymmOp.from_code(self.sg40_n_code)
        self.assertEqual(o1, eval(repr(o1)))
        self.assertEqual(o2, eval(repr(o2)))
        self.assertEqual(o3, eval(repr(o3)))

    def test_det_and_trace(self):
        o1 = SymmOp.from_code(self.sg40_a_code)
        o2 = SymmOp.from_code(self.sg40_21_code)
        self.assertEqual(o1.det, -1)
        self.assertEqual(o2.det, 1)
        self.assertEqual(o1.trace, 1)
        self.assertEqual(o2.trace, -1)

    def test_fold_and_order(self):
        o1 = SymmOp.from_code(self.sg40_a_code)
        o2 = SymmOp.from_code(self.sg147_m3_code)
        self.assertEqual(o1.fold, 2)
        self.assertEqual(o2.fold, 3)
        self.assertEqual(o1.order, 2)
        self.assertEqual(o2.order, 6)

    def test_typ(self):
        o1 = SymmOp.from_code(self.sg40_a_code)
        o2 = SymmOp.from_code(self.sg40_21_code)
        o3 = SymmOp.from_code(self.sg147_m3_code)
        self.assertIs(o1.typ, o1.Type.transflection)
        self.assertIs(o2.typ, o2.Type.rototranslation)
        self.assertIs(o3.typ, o3.Type.rotoinversion)


if __name__ == '__main__':
    unittest.main()
