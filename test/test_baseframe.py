import numpy as np
import unittest
from hikari.dataframes import BaseFrame


rad60 = 1.0471975511965976
rad70 = 1.2217304763960306
rad80 = 1.3962634015954636
rad90 = 1.5707963267948966


class TestBaseFrame(unittest.TestCase):
    def test_init(self):
        b = BaseFrame()
        self.assertAlmostEqual(b.a_d, 1.0)
        self.assertAlmostEqual(b.b_d, 1.0)
        self.assertAlmostEqual(b.c_d, 1.0)
        self.assertAlmostEqual(b.al_d, rad90)
        self.assertAlmostEqual(b.be_d, rad90)
        self.assertAlmostEqual(b.ga_d, rad90)

    def test_edit_cell(self):
        b = BaseFrame()
        b.edit_cell(a=6, b=7, c=8, al=60, be=70, ga=80)
        self.assertAlmostEqual(b.a_d, 6.0)
        self.assertAlmostEqual(b.b_d, 7.0)
        self.assertAlmostEqual(b.c_d, 8.0)
        self.assertAlmostEqual(b.al_d, rad60)
        self.assertAlmostEqual(b.be_d, rad70)
        self.assertAlmostEqual(b.ga_d, rad80)

    def test_reciprocal_cell(self):
        b = BaseFrame()
        b.edit_cell(a=6, b=7, c=8, al=60, be=70, ga=80)
        self.assertAlmostEqual(b.a_r, 0.1773638940193107)
        self.assertAlmostEqual(b.b_r, 0.1649580865234412)
        self.assertAlmostEqual(b.c_r, 0.1512680839136182)
        self.assertAlmostEqual(b.al_r, 2.067032910195549)
        self.assertAlmostEqual(b.be_r, 1.874672323582712)
        self.assertAlmostEqual(b.ga_r, 1.574038054665492)

    def test_vectors(self):
        b = BaseFrame()
        b.edit_cell(a=6, b=7, c=8, al=60, be=70, ga=80)
        expected_a_v = np.array([6.00000000,  0.00000000,  0.00000000])
        expected_b_v = np.array([1.21553724,  6.89365427,  0.00000000])
        expected_c_v = np.array([2.73616115,  3.57924741,  6.61077984])
        expected_a_w = np.array([0.16666667, -0.02938783, -0.05307098])
        expected_b_w = np.array([0.00000000,  0.14506094, -0.07853975])
        expected_c_w = np.array([0.00000000,  0.00000000,  0.15126808])
        self.assertAlmostEqual(sum(b.a_v - expected_a_v), 0.0)
        self.assertAlmostEqual(sum(b.b_v - expected_b_v), 0.0)
        self.assertAlmostEqual(sum(b.c_v - expected_c_v), 0.0)
        self.assertAlmostEqual(sum(b.a_w - expected_a_w), 0.0)
        self.assertAlmostEqual(sum(b.b_w - expected_b_w), 0.0)
        self.assertAlmostEqual(sum(b.c_w - expected_c_w), 0.0)

    def test_perpendicular(self):
        b = BaseFrame()
        b.edit_cell(a=6, b=7, c=8, al=60, be=70, ga=80)

        def norm(v):
            return v / np.sqrt(v.dot(v))
        self.assertAlmostEqual(sum(norm(np.cross(b.a_v, b.b_v))-norm(b.c_w)), 0)
        self.assertAlmostEqual(sum(norm(np.cross(b.b_v, b.c_v))-norm(b.a_w)), 0)
        self.assertAlmostEqual(sum(norm(np.cross(b.c_v, b.a_v))-norm(b.b_w)), 0)
        self.assertAlmostEqual(sum(norm(np.cross(b.a_w, b.b_w))-norm(b.c_v)), 0)
        self.assertAlmostEqual(sum(norm(np.cross(b.b_w, b.c_w))-norm(b.a_v)), 0)
        self.assertAlmostEqual(sum(norm(np.cross(b.c_w, b.a_w))-norm(b.b_v)), 0)

    def test_matrices(self):
        b = BaseFrame()
        b.edit_cell(a=6, b=7, c=8, al=60, be=70, ga=80)
        v = np.array([1.6180339887, 2.7182818285, 3.1415926536])
        self.assertAlmostEqual(v @ b.G_d @ v, 1797.2495106787824)
        self.assertAlmostEqual(v @ b.G_r @ v, 0.2238961843294296)
        self.assertAlmostEqual(sum(sum(b.G_d @ b.A_r - b.A_d)), 0.0)
        self.assertAlmostEqual(sum(sum(b.G_r @ b.A_d - b.A_r)), 0.0)


if __name__ == '__main__':
    unittest.main()
