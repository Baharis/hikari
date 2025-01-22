import unittest

import numpy as np

from tests import NumpyTestCase
from hikari.symmetry import Operation, BoundedOperation, PointOperation


p100 = np.array([1, 0, 0])
p010 = np.array([0, 1, 0])
p110 = np.array([1, 1, 0])
p111 = np.array([1, 1, 1])

x = np.array([1, 0, 0])
y = np.array([0, 1, 0])
z = np.array([0, 0, 1])

sg40_a_tf = np.array([[1, 0, 0], [0, -1, 0], [0, 0, 1]])
sg40_a_tl = np.array([1 / 2, 0, 0])
sg40_a_m = np.array([[1, 0, 0, .5], [0, -1, 0, 0], [0, 0, 1, 0], [0, 0, 0, 1]])
sg40_a_code = 'x+1/2,-y,z'
sg40_21_code = '-x,-y+1/2,z+1/2'
sg40_n_code = '-x+1/2,y-1/2,z+1/2'
sg147_m3_code = 'x-y,x,-z'


class TestOperationsInit(unittest.TestCase):
    def test_from_pair(self):
        o = Operation.from_pair(sg40_a_tf, sg40_a_tl)
        self.assertTrue(np.allclose(sg40_a_tf, o.tf))
        self.assertTrue(np.allclose(sg40_a_tl, o.tl))

    def test_from_code(self):
        o = Operation.from_code(sg40_a_code)
        self.assertTrue(np.allclose(sg40_a_tf, o.tf))
        self.assertTrue(np.allclose(sg40_a_tl, o.tl))

    def test_from_matrix(self):
        o = Operation.from_matrix(sg40_a_m)
        self.assertTrue(np.allclose(sg40_a_tf, o.tf))
        self.assertTrue(np.allclose(sg40_a_tl, o.tl))


class TestSelectedOperations(NumpyTestCase):
    """Define space group #40 'Ama2' operations here to avoid code repetition"""

    @classmethod
    def setUpClass(cls):
        cls.sg40_a = Operation.from_code(sg40_a_code)
        cls.sg40_21 = Operation.from_code(sg40_21_code)
        cls.sg40_n = Operation.from_code(sg40_n_code)
        cls.sg147_m3 = Operation.from_code(sg147_m3_code)
        cls.operations = [cls.sg40_a, cls.sg40_21, cls.sg40_n, cls.sg147_m3]


class TestOperations(TestSelectedOperations):
    """Test selected aspects specific to `Operation`s"""

    def test_equals(self):
        self.assertEqual(self.sg40_a, self.sg40_a)
        self.assertNotEqual(self.sg40_a, self.sg40_21)
        self.assertNotEqual(self.sg40_a, 42)

    def test_product(self):
        self.assertEqual(self.sg40_a * self.sg40_21, self.sg40_n)

    def test_power(self):
        for o in self.operations:
            self.assertEqual(o*o, o**2)
            self.assertEqual(o*o*o*o*o, o**5)

    def test_inverse(self):
        identity = Operation.from_code('x,y,z')
        for o in self.operations:
            self.assertEqual(o * (o ** -1), identity)
            self.assertEqual((o ** -1) * o, identity)
            self.assertEqual(o * pow(o, -1), identity)

    def test_repr(self):
        for o in self.operations:
            self.assertEqual(o, eval(repr(o)))

    def test_str(self):
        self.assertEqual(str(self.sg40_a), 'a: x+1/2,-y,z (0,0,0)')
        self.assertEqual(str(self.sg40_21), '21: -x,-y+1/2,z+1/2 (0,1/4,0)')
        self.assertEqual(str(self.sg40_n), 'x: -x+1/2,y-1/2,z+1/2 (1/4,0,0)')

    def test_hash(self):
        for o in self.operations:
            self.assertIsInstance(hash(o), int)

    def test_code(self):
        self.assertEqual(self.sg40_a.code, sg40_a_code)
        self.assertEqual(self.sg40_21.code, sg40_21_code)
        self.assertEqual(self.sg40_n.code, sg40_n_code)
        self.assertEqual(self.sg147_m3.code, sg147_m3_code)

    def test_det(self):
        self.assertEqual(self.sg40_a.det, -1)
        self.assertEqual(self.sg40_21.det, 1)
        self.assertEqual(self.sg40_n.det, -1)
        self.assertEqual(self.sg147_m3.det, -1)

    def test_typ(self):
        self.assertIs(self.sg40_a.typ, Operation.Type.transflection)
        self.assertIs(self.sg40_21.typ, Operation.Type.rototranslation)
        self.assertIs(self.sg40_n.typ, Operation.Type.transflection)
        self.assertIs(self.sg147_m3.typ, Operation.Type.rotoinversion)

    def test_fold(self):
        self.assertEqual(self.sg40_a.fold, 2)
        self.assertEqual(self.sg40_21.fold, 2)
        self.assertEqual(self.sg40_n.fold, 2)
        self.assertEqual(self.sg147_m3.fold, 3)

    def test_order(self):
        self.assertEqual(self.sg40_a.order, 2)
        self.assertEqual(self.sg40_21.order, 2)
        self.assertEqual(self.sg40_n.order, 2)
        self.assertEqual(self.sg147_m3.order, 6)

    def test_glide(self):
        self.assertEqual(self.sg40_a.glide[0], 0.5)
        self.assertEqual(self.sg40_21.glide[2], 0.5)
        self.assertEqual(self.sg40_n.glide[1], -0.5)
        self.assertEqual(self.sg40_n.glide[2], 0.5)
        self.assertEqual(sum(self.sg147_m3.glide), 0.)

    def test_trace(self):
        self.assertEqual(self.sg40_a.trace, 1)
        self.assertEqual(self.sg40_21.trace, -1)
        self.assertEqual(self.sg40_n.trace, 1)
        self.assertEqual(self.sg147_m3.trace, 0)

    def test_orientation(self):
        self.assertEqual(abs(np.dot(self.sg40_a.orientation, y)), 1)
        self.assertEqual(abs(np.dot(self.sg40_21.orientation, z)), 1)
        self.assertEqual(abs(np.dot(self.sg40_n.orientation, x)), 1)
        self.assertEqual(abs(np.dot(self.sg147_m3.orientation, x)), 1)

    def test_sense(self):
        self.assertEqual(self.sg40_a.sense, '')
        self.assertEqual(self.sg40_21.sense, '')
        self.assertEqual(self.sg40_n.sense, '')
        self.assertEqual(self.sg147_m3.sense, '-')

    def test_at(self):
        self.assertEqual(self.sg40_a.at(p111).code, 'x+1/2,-y+2,z')
        self.assertEqual(self.sg40_21.at(p111).code, '-x+2,-y+2,z+1/2')
        self.assertEqual(self.sg40_n.at(p111).code, '-x+2,y-1/2,z+1/2')
        self.assertEqual(self.sg147_m3.at(p111).code, 'x-y+1,x,-z+2')

    def test_into(self):
        self.assertEqual(self.sg40_a.into(z).code, 'x+1/2,y,-z')
        self.assertEqual(self.sg40_21.into(y).code, '-x,y+1,-z')
        self.assertEqual(self.sg40_n.into(y).code, 'x+1,-y,z+1/2')
        self.assertEqual(self.sg147_m3.into(y).code, '-y,x+y,-z')

    def test_transform(self):
        pp = np.expand_dims(p111, 1).T
        self.assertAllClose(self.sg40_a.transform(pp), np.array([[1.5, -1., 1.]]))
        self.assertAllClose(self.sg40_21.transform(pp), np.array([[-1, -0.5, 1.5]]))
        self.assertAllClose(self.sg40_n.transform(pp), np.array([[-0.5, 0.5, 1.5]]))
        self.assertAllClose(self.sg147_m3.transform(pp), np.array([[0., 1., -1.]]))

    def test_extincts(self):
        r = np.array([[0, 0, 1], [0, 0, 2], [1, 0, 1], [1, -1, 1]])
        self.assertEqual(list(self.sg40_a.extincts(r)), [False, False, True, False])
        self.assertEqual(list(self.sg40_21.extincts(r)), [True, False, False, False])
        self.assertEqual(list(self.sg40_n.extincts(r)), [True, False, False, False])
        self.assertEqual(list(self.sg147_m3.extincts(r)), [False, False, False, False])

    def test_distance_to_point(self):
        self.assertEqual(self.sg40_a.distance_to_point(p111), 1.0)
        self.assertEqual(self.sg40_21.distance_to_point(p110), 1.25)
        self.assertEqual(self.sg40_n.distance_to_point(p111), 0.75)
        self.assertEqual(self.sg147_m3.distance_to_point(p111), np.sqrt(2))


class TestBoundedOperations(TestSelectedOperations):
    """Test selected aspects specific to `BoundedOperation`s"""

    def test_binding(self):
        self.assertEqual(self.sg40_21, self.sg40_21.bounded)
        self.assertNotEqual(self.sg40_n, self.sg40_n.bounded)
        self.assertEqual(self.sg40_n.bounded.code, '-x+1/2,y+1/2,z+1/2')

    def test_unbinding(self):
        for op in self.operations:
            self.assertEqual((op.bounded ** 12).tl[0], 0)
            self.assertEqual((op.bounded ** 12).tl[1], 0)
            self.assertEqual((op.bounded ** 12).tl[2], 0)
        self.assertEqual((self.sg40_a.unbounded ** 12).tl[0], 6)
        self.assertEqual((self.sg40_21.unbounded ** 12).tl[2], 6)
        self.assertEqual((self.sg40_n.unbounded ** 12).tl[1], -6)
        self.assertEqual((self.sg40_n.unbounded ** 12).tl[2], 6)
        self.assertEqual((self.sg147_m3.unbounded ** 12).tl[0], 0)

    def test_is_bounded(self):
        self.assertTrue(self.sg40_a.is_bounded)
        self.assertTrue(self.sg40_21.is_bounded)
        self.assertFalse(self.sg40_n.is_bounded)
        self.assertTrue(self.sg147_m3.is_bounded)




class TestPointOperations(TestSelectedOperations):
    """Test selected aspects specific to `BoundedOperation`s"""

    def test_point_restraint(self):
        for op in [self.sg40_a, self.sg40_21, self.sg40_n]:
            with self.assertRaises(ValueError):
                _ = PointOperation(op.tf, op.tl)
        _ = PointOperation(self.sg147_m3.tf, self.sg147_m3.tl)


class TestCatalogOperations(unittest.TestCase):
    """Run complete rigorous tests on all operations defined in `SG` catalog"""

    @classmethod
    def setUpClass(cls):
        from hikari.symmetry import SG
        cls.operations: list[BoundedOperation] \
            = list({o for s in SG.values() for o in s.operations})

    def test_catalog_names(self):
        for op in self.operations:
            self.assertFalse(op.name.startswith('?'))

    def test_catalog_fold(self):
        for op in self.operations:
            self.assertIn(op.fold, {1, 2, 3, 4, 5, 6})

    def test_catalog_order(self):
        for op in self.operations:
            self.assertIn(op.order, {1, 2, 3, 4, 5, 6})

    def test_catalog_at_origin(self):
        for op in self.operations:
            self.assertEqual(op, op.at(op.origin))

    def test_catalog_into_orientation(self):
        for op in self.operations:
            if (orientation := op.orientation) is not None:
                self.assertEqual(op, op.into(orientation))

    def test_catalog_distance_to_point(self):
        p_pi = np.array([np.pi, np.pi, np.pi])
        for op in self.operations:
            self.assertNotEqual(op.distance_to_point(p_pi), 0.0)


if __name__ == '__main__':
    unittest.main()
