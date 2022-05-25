import numpy as np
import unittest
from hikari.utility import *


class TestChemTools(unittest.TestCase):
    def test_chemical_elements(self):
        self.assertEqual(str(chemical_elements[0]), 'H')
        self.assertEqual(str(chemical_elements[49]), 'Sn')
        self.assertEqual(str(chemical_elements[99]), 'Fm')


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
        self.assertEqual(angle2rad(3.14), 3.14)
        self.assertAlmostEqual(angle2rad(4), 0.06981317007977318)

    def test_fibonacci_sphere(self):
        self.assertEqual(len(fibonacci_sphere(100)), 100)
        self.assertEqual(len(fibonacci_sphere(4)[0]), 3)
        self.assertAlmostEqual(fibonacci_sphere(4, seed=1)[0][0], 0.1833764108)
        self.assertAlmostEqual(
            (fibonacci_sphere(9)-fibonacci_sphere(9, seed=1337)).max().max(), 0)
        self.assertNotAlmostEqual(
            (fibonacci_sphere(9)-fibonacci_sphere(9, seed=123)).max().max(), 0)
        self.assertTrue(abs(fibonacci_sphere(99999).sum().sum() < 1.0))


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


class TestInterval(unittest.TestCase):
    def test_creation_and_access(self):
        a = Interval(4, 8)
        self.assertIsInstance(a, Interval)
        self.assertEqual(a[0], a.left)
        self.assertEqual(a[1], a.right)
        with self.assertRaises(IndexError):
            _ = a[2]
        self.assertIn(5, a)
        self.assertIn(4, a)
        self.assertNotIn(3, a)

    def test_comparison_operations(self):
        a = Interval(4, 8)
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
        self.assertEqual(-a, Interval(-8, -4))
        self.assertNotEqual(a, -a)
        self.assertEqual(a + 1, Interval(5, 9))
        self.assertEqual(a - 1, Interval(3, 7))
        self.assertEqual(a * 2, Interval(8, 16))
        self.assertAlmostEqual((a / 2).left, Interval(2, 4).left)
        self.assertAlmostEqual((a / 2).right, Interval(2, 4).right)

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


if __name__ == '__main__':
    unittest.main()
