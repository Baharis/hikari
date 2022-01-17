import unittest
from hikari.utility import make_abspath, angle2rad, fibonacci_sphere,\
    find_best, cubespace, rescale_list_to_range, rescale_list_to_other


class TestOsTools(unittest.TestCase):
    def test_make_abspath(self):
        self.assertIsInstance(make_abspath(), str)
        self.assertIsInstance(make_abspath('~', '..'), str)
        self.assertEqual(make_abspath('~', 'f.ext'), make_abspath('~/f.ext'))


class TestMathTools(unittest.TestCase):
    def test_angle2rad(self):
        self.assertEqual(angle2rad(3.14), 3.14)
        self.assertAlmostEqual(angle2rad(4), 0.06981317007977318)

    def test_fibonacci_sphere(self):
        self.assertEqual(len(fibonacci_sphere(100)), 100)
        self.assertEqual(len(fibonacci_sphere(4)[0]), 3)
        self.assertAlmostEqual(fibonacci_sphere(4, seed=1)[0][0], 0.1833764108)
        self.assertNotEqual(fibonacci_sphere(9), fibonacci_sphere(9, seed=123))


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


if __name__ == '__main__':
    unittest.main()
