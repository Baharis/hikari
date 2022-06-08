import unittest
import numpy as np
from hikari.utility import *


class TestChemTools(unittest.TestCase):
    def test_chemical_elements(self):
        self.assertEqual(str(chemical_elements[0]), 'H')
        self.assertEqual(str(chemical_elements[99]), 'Fm')
        self.assertEqual(chemical_elements.index('Sn'), 49)
        self.assertEqual(chemical_elements.index('Og'), 117)

    def test_split_atom_label(self):
        self.assertEqual(split_atom_label('C12A'), ('C', '12', 'A'))
        self.assertEqual(split_atom_label('Fe1'), ('Fe', '1', ''))
        self.assertEqual(split_atom_label('H1O1'), ('H', '1', 'O1'))


class TestCertainFloat(unittest.TestCase):
    def test_on_float(self):
        self.assertEqual(cfloat('12.34'), 12.34)

    def test_on_shorthand_ufloat(self):
        self.assertEqual(cfloat('12.34(567)'), 12.34)

    def test_on_plus_minus_ufloat(self):
        self.assertEqual(cfloat('12.34+/-5.67'), 12.34)

    def test_on_plus_minus_spaced_ufloat(self):
        self.assertEqual(cfloat('12.34 +/- 5.67'), 12.34)

    def test_on_factored_exponent_ufloat(self):
        self.assertEqual(cfloat('1.234(567)e1'), 12.34)

    def test_on_pretty_print_ufloat(self):
        self.assertEqual(cfloat(u'12.34±5.67'), 12.34)


class TestInterval(unittest.TestCase):
    def test_creation(self):
        a = Interval(4, 8)
        self.assertIsInstance(a, Interval)

    def test_comparison_operations(self):
        a = Interval(4, 8)
        self.assertEqual(a, Interval(4, 8))
        self.assertTrue(a > 3)
        self.assertFalse(a > 6.5)
        self.assertTrue(a < 8.5)
        self.assertFalse(a < 7.5)
        self.assertTrue(a >= 4)
        self.assertFalse(a >= 6.5)
        self.assertTrue(a <= 8.1)
        self.assertFalse(a <= 7.9)

    def test_unary_and_arithmetic(self):
        a = Interval(4, 8)
        self.assertEqual(+a, Interval(4, 8))
        self.assertEqual(-a, Interval(-8, -4))
        self.assertNotEqual(a, -a)
        self.assertEqual(a + 1, Interval(5, 9))
        self.assertEqual(a - 1, Interval(3, 7))
        self.assertEqual(a * 2, Interval(8, 16))
        self.assertAlmostEqual((a / 2).left, Interval(2, 4).left)
        self.assertAlmostEqual((a / 2).right, Interval(2, 4).right)

    def test_str_and_repr(self):
        a = Interval(4, 8)
        self.assertEqual(str(a), '[4, 8]')
        self.assertEqual(repr(a), 'Interval(4, 8)')

    def test_iter_and_length(self):
        self.assertEqual([_ for _ in Interval(4, 8)], [4, 8])
        self.assertEqual(len(Interval(4, 8)), 2)

    def test_getitem(self):
        a = Interval(4, 8)
        self.assertEqual(a[0], 4)
        self.assertEqual(a[1], 8)
        self.assertEqual(a[0], a.left)
        self.assertEqual(a[1], a.right)
        self.assertEqual(a[:], [4, 8])
        with self.assertRaises(IndexError):
            _ = a[2]

    def test_setitem(self):
        a = Interval(4, 8)
        a[0], a[1] = 3, 7
        self.assertEqual(a, Interval(3, 7))
        a[:] = 2, 6
        self.assertEqual(a, Interval(2, 6))
        a.left, a.right = 1, 5
        self.assertEqual(a, Interval(1, 5))
        with self.assertRaises(IndexError):
            a[2] = 0

    def test_contains(self):
        a = Interval(4, 8)
        c = Interval(6, 6)
        self.assertIn(5, a)
        self.assertNotIn(3, a)
        self.assertIn(c, a)
        self.assertNotIn(a, c)
        self.assertIn([5, 7], a)

    def test_min_and_max(self):
        a = Interval(4, 8)
        self.assertEqual(min(a), 4)
        self.assertEqual(max(a), 8)

    def test_arange(self):
        a = Interval(1, 3)
        self.assertTrue(np.array_equal(a.arange(), np.array([1, 2, 3])))
        self.assertTrue(np.array_equal(a.arange(step=2), np.array([1, 3])))
        self.assertTrue(np.allclose(a.arange(step=0.5), np.arange(1, 3.1, 0.5)))

    def test_comb_with(self):
        a = Interval(1, 3)
        b = Interval(4, 5)
        c = Interval(6, 6)
        self.assertTrue(np.array_equal(a.comb_with(b)[0],
                                       np.array([1, 2, 3] * 2)))
        self.assertTrue(np.allclose(a.comb_with(b, step=2)[1], [4, 4]))
        self.assertTrue(np.allclose(a.comb_with(b, step=0.9)[1],
                                    np.array([4.0] * 3 + [4.9] * 3)))
        self.assertTrue(np.array_equal(a.comb_with(b, c)[:2],
                                       a.comb_with(b, c, c, c, c, c)[:2]))

    def test_mesh_with(self):
        a = Interval(1, 3)
        b = Interval(4, 5)
        self.assertTrue(np.array_equal(a.mesh_with(b)[0],
                                       np.array([[1, 2, 3], [1, 2, 3]])))
        self.assertTrue(np.allclose(a.mesh_with(b, step=0.8)[1],
                                    np.array([[4., 4., 4.], [4.8, 4.8, 4.8]])))
        self.assertTrue(np.array_equal(a.mesh_with(b)[0], b.mesh_with(a)[1].T))


class TestListTools(unittest.TestCase):
    def test_find_best(self):
        self.assertIs(find_best('12ab?#*', 'c>b>a'), 'b')
        self.assertIs(find_best('12ab?#*', '1+2=3>2>1'), '3')
        self.assertIs(find_best('12ab?#*', 'C>B>A'), '')

    def test_cubespace(self):
        cs1 = cubespace(10)
        cs2 = cubespace(0, -10, num=10, include_start=True)
        self.assertEqual(cs1[5], -cs2[5])
        self.assertEqual(len(cubespace(-5, 0, include_start=False)), 10)
        self.assertAlmostEqual(cs1[6]**3 - cs1[5]**3, cs1[4]**3 - cs1[3]**3)

    def test_rescale_list_to_range(self):
        self.assertEqual(len(rescale_list_to_range(list(range(9)), (0, 1j))), 9)
        self.assertAlmostEqual(rescale_list_to_range([1, 5, 10], (7, 14))[1],
                               rescale_list_to_range([-1, -6, -10], (7, 14))[1])

    def test_rescale_list_to_other(self):
        self.assertEqual(rescale_list_to_other([1, 9], list(range(9))), [0, 8])
        self.assertEqual(''.join(rescale_list_to_other([1, 4], 'Test')), 'Tt')


class TestMathTools(unittest.TestCase):
    def test_angle2rad(self):
        self.assertEqual(angle2rad(-5), -0.08726646259971647)
        self.assertEqual(angle2rad(-2.0), -2.0)
        self.assertEqual(angle2rad(3.14), 3.14)
        self.assertAlmostEqual(angle2rad(4), 0.06981317007977318)
        self.assertEqual(angle2rad(180), 3.141592653589793)

    def test_cart2sph(self):
        s = np.sqrt(2)
        self.assertTrue(np.allclose(cart2sph(1, 0, 0), [1, np.pi/2, 0]))
        self.assertTrue(np.allclose(cart2sph(0, 1, 0), [1, np.pi/2, np.pi/2]))
        self.assertTrue(np.allclose(cart2sph(0, 0, 1), [1, 0, 0]))
        self.assertTrue(np.allclose(cart2sph(0, -1, 1), [s, np.pi/4, -np.pi/2]))
        self.assertTrue(np.allclose(cart2sph(1, 0, -1), [s, 3*np.pi/4, 0]))

    def test_sph2cart(self):
        s = np.sqrt(2)
        self.assertTrue(np.allclose(sph2cart(1, np.pi/2, 0), [1, 0, 0]))
        self.assertTrue(np.allclose(sph2cart(1, np.pi/2, np.pi/2), [0, 1, 0]))
        self.assertTrue(np.allclose(sph2cart(1, 0, 0), [0, 0, 1]))
        self.assertTrue(np.allclose(sph2cart(s, np.pi/4, -np.pi/2), [0, -1, 1]))
        self.assertTrue(np.allclose(sph2cart(s, 3*np.pi/4, 0), [1, 0, -1]))

    def test_det3x3(self):
        a = np.array([(1, 2, 3), (4, 5, 6), (7, 8, 9)])
        self.assertEqual(det3x3(np.eye(3)), 1)
        self.assertEqual(det3x3(a), 0)
        self.assertEqual(det3x3(a ** 2), -216)

    def test_fibonacci_sphere(self):
        self.assertEqual(len(fibonacci_sphere(100)), 100)
        self.assertEqual(len(fibonacci_sphere(4)[0]), 3)
        self.assertAlmostEqual(fibonacci_sphere(4, seed=1)[0][0], 0.1833764108)
        self.assertAlmostEqual(
            (fibonacci_sphere(9)-fibonacci_sphere(9, seed=1337)).max().max(), 0)
        self.assertNotAlmostEqual(
            (fibonacci_sphere(9)-fibonacci_sphere(9, seed=123)).max().max(), 0)
        self.assertTrue(abs(fibonacci_sphere(99999).sum().sum() < 1.0))

    def test_rotation_from(self):
        a, b, z = np.array([1, 2, 3]), np.array([1, 2]), np.array([0, 0, 1])
        len_a = np.sqrt(14)
        self.assertTrue(np.allclose(rotation_from(a, to=a), np.eye(3)))
        self.assertTrue(np.allclose(rotation_from(a, to=z) @ a, z * len_a))
        self.assertTrue(np.allclose(rotation_from(a, to=-a) @ a, -a))
        self.assertTrue(np.allclose(
            rotation_from(np.array([1, 0, 0]), to=np.array([0, 1, 0])),
            np.array([(0, -1, 0), (1, 0, 0), (0, 0, 1)])))
        with self.assertRaises(IndexError):
            _ = rotation_from(from_=a, to=b)
        with self.assertRaises(IndexError):
            _ = rotation_from(from_=b, to=a)

    def test_rotation_around(self):
        x = np.array([1, 0, 0])
        z = np.array([0, 0, 1])
        self.assertTrue(np.allclose(rotation_around(z, by=0), np.eye(3)))
        self.assertTrue(np.allclose(rotation_around(z, by=np.pi) @ x, -x))
        self.assertTrue(np.allclose(
            rotation_around(z, by=np.pi/2),
            np.array([(0, -1, 0), (1, 0, 0), (0, 0, 1)])))
        with self.assertRaises(IndexError):
            _ = rotation_around(np.array([1, 2, 3, 4]), by=0)

    def test_weighted_quantile(self):
        v = [1, 2, 3, 4, 5]
        o = [1, 1, 1, 1, 1]
        self.assertEqual(weighted_quantile(values=v, quantiles=[0.5])[0], 3)
        self.assertEqual(weighted_quantile(v, [1.0])[0], 5)
        self.assertEqual(weighted_quantile(v, [0.5], weights=o)[0], 3)
        self.assertGreater(weighted_quantile(v, [0.5], weights=v)[0], 3)
        self.assertLess(weighted_quantile(v, [0.5], weights=v)[0], 4)
        with self.assertRaises(ValueError):
            _ = weighted_quantile(values=v, quantiles=[-0.1, 0.5, 1.1])


class TestOsTools(unittest.TestCase):
    def test_make_abspath(self):
        self.assertIsInstance(make_abspath(), str)
        self.assertIsInstance(make_abspath('~', '..'), str)
        self.assertEqual(make_abspath('~', 'f.ext'), make_abspath('~/f.ext'))


class TestPalettes(unittest.TestCase):
    def test_gnuplot_palette(self):
        self.assertIsInstance(gnuplot_map_palette[''], str)
        self.assertIsInstance(gnuplot_map_palette['h'], str)
        self.assertIsInstance(gnuplot_map_palette['k'], str)
        self.assertIsInstance(gnuplot_map_palette['l'], str)
        self.assertIsInstance(gnuplot_map_palette['hk'], str)
        self.assertIsInstance(gnuplot_map_palette['kl'], str)
        self.assertIsInstance(gnuplot_map_palette['hl'], str)

    def test_mpl_palette(self):
        try:
            _ = iter(mpl_map_palette['h'])
            _ = iter(mpl_map_palette['k'])
            _ = iter(mpl_map_palette['l'])
            _ = iter(mpl_map_palette['hk'])
            _ = iter(mpl_map_palette['kl'])
            _ = iter(mpl_map_palette['hl'])
        except TypeError:
            self.fail('mpl_map_palette is not iterable, but it should')


if __name__ == '__main__':
    unittest.main()
