import tempfile, unittest
import numpy as np
from hikari.dataframes import BaseFrame, HklFrame
from hikari.resources import nacl_hkl

rad60 = 1.0471975511965976
rad70 = 1.2217304763960306
rad80 = 1.3962634015954636
rad90 = 1.5707963267948966


class TestBaseFrame(unittest.TestCase):
    b = BaseFrame()

    def test_init(self):
        self.assertAlmostEqual(BaseFrame().a_d, 1.)
        self.assertAlmostEqual(BaseFrame().b_d, 1.)
        self.assertAlmostEqual(BaseFrame().c_d, 1.)
        self.assertAlmostEqual(BaseFrame().al_d, rad90)
        self.assertAlmostEqual(BaseFrame().be_d, rad90)
        self.assertAlmostEqual(BaseFrame().ga_d, rad90)

    def test_edit_cell(self):
        self.b.edit_cell(a=6, b=7, c=8, al=60, be=70, ga=80)
        self.assertAlmostEqual(self.b.a_d, 6.)
        self.assertAlmostEqual(self.b.b_d, 7.)
        self.assertAlmostEqual(self.b.c_d, 8.)
        self.assertAlmostEqual(self.b.al_d, rad60)
        self.assertAlmostEqual(self.b.be_d, rad70)
        self.assertAlmostEqual(self.b.ga_d, rad80)

    def test_reciprocal_cell(self):
        self.assertAlmostEqual(self.b.a_r, 0.1773638940193107)
        self.assertAlmostEqual(self.b.b_r, 0.1649580865234412)
        self.assertAlmostEqual(self.b.c_r, 0.1512680839136182)
        self.assertAlmostEqual(self.b.al_r, 2.067032910195549)
        self.assertAlmostEqual(self.b.be_r, 1.874672323582712)
        self.assertAlmostEqual(self.b.ga_r, 1.574038054665492)

    def test_vectors(self):
        expected_a_v = np.array([6.00000000,  0.00000000,  0.00000000])
        expected_b_v = np.array([1.21553724,  6.89365427,  0.00000000])
        expected_c_v = np.array([2.73616115,  3.57924741,  6.61077984])
        expected_a_w = np.array([0.16666667, -0.02938783, -0.05307098])
        expected_b_w = np.array([0.00000000,  0.14506094, -0.07853975])
        expected_c_w = np.array([0.00000000,  0.00000000,  0.15126808])
        self.assertAlmostEqual(sum(self.b.a_v - expected_a_v), 0.)
        self.assertAlmostEqual(sum(self.b.b_v - expected_b_v), 0.)
        self.assertAlmostEqual(sum(self.b.c_v - expected_c_v), 0.)
        self.assertAlmostEqual(sum(self.b.a_w - expected_a_w), 0.)
        self.assertAlmostEqual(sum(self.b.b_w - expected_b_w), 0.)
        self.assertAlmostEqual(sum(self.b.c_w - expected_c_w), 0.)

    def test_perpendicular(self):
        def norm(v):
            return v / np.sqrt(v.dot(v))
        av, bv, cv = self.b.a_v, self.b.b_v, self.b.c_v
        aw, bw, cw = self.b.a_w, self.b.b_w, self.b.c_w
        self.assertAlmostEqual(sum(norm(np.cross(av, bv)) - norm(cw)), 0.)
        self.assertAlmostEqual(sum(norm(np.cross(bv, cv)) - norm(aw)), 0.)
        self.assertAlmostEqual(sum(norm(np.cross(cv, av)) - norm(bw)), 0.)
        self.assertAlmostEqual(sum(norm(np.cross(aw, bw)) - norm(cv)), 0.)
        self.assertAlmostEqual(sum(norm(np.cross(bw, cw)) - norm(av)), 0.)
        self.assertAlmostEqual(sum(norm(np.cross(cw, aw)) - norm(bv)), 0.)

    def test_matrices(self):
        v = np.array([1.6180339887, 2.7182818285, 3.1415926536])
        self.assertAlmostEqual(v @ self.b.G_d @ v, 1797.2495106787824)
        self.assertAlmostEqual(v @ self.b.G_r @ v, 0.2238961843294296)
        self.assertAlmostEqual(sum(sum(self.b.G_d@self.b.A_r - self.b.A_d)), 0.)
        self.assertAlmostEqual(sum(sum(self.b.G_r@self.b.A_d - self.b.A_r)), 0.)


class TestHklFrame(unittest.TestCase):
    h = HklFrame()
    temporary_hkl_file = tempfile.NamedTemporaryFile('w+')

    @classmethod
    def setUpClass(cls) -> None:
        cls.temporary_hkl_file.write(nacl_hkl)
        cls.h.read(cls.temporary_hkl_file.name, hkl_format='shelx_4')

    @classmethod
    def tearDownClass(cls) -> None:
        cls.temporary_hkl_file.close()

    def test_len(self):
        self.assertEqual(len(self.h), 8578)

    def test_add(self):
        h2 = self.h + self.h
        self.assertEqual(len(h2), 17156)
        self.assertIsInstance(h2, HklFrame)

    def test_str(self):
        str_h = str(self.h)
        self.assertIn(str_h[:300], 'h')
        self.assertIn(str_h[:300], '52001.50')
        self.assertIn(str_h[-300:], '0.0')

    def test_read(self):
        self.h.read(self.temporary_hkl_file.name, hkl_format='free_4')
        self.assertEqual(self.h.table.__len__(), 8578)

# test other functions, in particular trim, dac, merge etc




if __name__ == '__main__':
    unittest.main()
