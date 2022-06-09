import copy
import pathlib
import tempfile
import unittest

import numpy as np
import uncertainties

from hikari.dataframes import BaseFrame, CifBlock, CifFrame, HklFrame, \
    UBaseFrame
from hikari.dataframes.cif import CifValidator
from hikari.symmetry import PG

RAD60 = 1.0471975511965976
RAD70 = 1.2217304763960306
RAD80 = 1.3962634015954636
RAD90 = 1.5707963267948966

U1 = uncertainties.ufloat(1, 1)
U6 = uncertainties.ufloat(6.0, 0.1)
U7 = uncertainties.ufloat(7.0, 0.1)
U8 = uncertainties.ufloat(8.0, 0.1)
U60 = uncertainties.ufloat(60.0, 1.0)
U70 = uncertainties.ufloat(70.0, 1.0)
U80 = uncertainties.ufloat(80.0, 1.0)
U90 = uncertainties.ufloat(90.0, 1.0)


class TempFile:
    temp_dir = tempfile.TemporaryDirectory()

    def __init__(self, name, content):
        self.name = name
        self.path = str(pathlib.Path(self.temp_dir.name) / name)
        self.file = open(self.path, 'w+')
        self.file.write(content)


nacl_cif_path = str(pathlib.Path(__file__).parent.joinpath('NaCl.cif'))
nacl_fcf_path = str(pathlib.Path(__file__).parent.joinpath('NaCl.fcf'))
nacl_hkl_path = str(pathlib.Path(__file__).parent.joinpath('NaCl.hkl'))


class TestBaseFrame(unittest.TestCase):
    def setUp(self) -> None:
        self.b = BaseFrame()
        self.b.edit_cell(a=6, b=7, c=8, al=60, be=70, ga=80)

    def test_init(self):
        b = BaseFrame()
        self.assertAlmostEqual(b.a_d, 1.)
        self.assertAlmostEqual(b.b_d, 1.)
        self.assertAlmostEqual(b.c_d, 1.)
        self.assertAlmostEqual(b.al_d, RAD90)
        self.assertAlmostEqual(b.be_d, RAD90)
        self.assertAlmostEqual(b.ga_d, RAD90)

    def test_edit_cell(self):
        self.assertAlmostEqual(self.b.a_d, 6.)
        self.assertAlmostEqual(self.b.b_d, 7.)
        self.assertAlmostEqual(self.b.c_d, 8.)
        self.assertAlmostEqual(self.b.al_d, RAD60)
        self.assertAlmostEqual(self.b.be_d, RAD70)
        self.assertAlmostEqual(self.b.ga_d, RAD80)
        with self.assertRaises(KeyError):
            BaseFrame().edit_cell(alpha=0.0)

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

    def test_volumes(self):
        self.assertAlmostEqual(self.b.v_d, 273.4345841924758)
        self.assertAlmostEqual(self.b.v_r, 0.003657181855591761)


class TestCifBlockRead(unittest.TestCase):
    def test_read(self):
        CifBlock().read(path=nacl_cif_path, block='NaCl')


class TestCifBlockGeneral(unittest.TestCase):
    b = CifBlock()

    @classmethod
    def setUpClass(cls) -> None:
        super().setUpClass()
        cls.b.read(path=nacl_cif_path, block='NaCl')

    def test_access(self):
        self.assertIn('_cell_length_a', self.b)
        self.assertEqual(self.b['_cell_length_a'], '5.64109(5)')
        self.assertEqual(len(self.b['_diffrn_refln_index_h']), 8578)

    def test_assign(self):
        self.b['_cell_angle_beta'] = '120'
        self.assertEqual(self.b['_cell_angle_beta'], '120')

    def test_get_as_type_single(self):
        k = '_cell_angle_alpha'
        self.assertEqual(self.b.get_as_type(k, typ=str), '90')
        self.assertEqual(self.b.get_as_type(k, typ=int), 90)
        u_typ = uncertainties.ufloat_fromstr
        self.assertEqual(repr(self.b.get_as_type(k, typ=u_typ)), repr(U90))

    def test_get_as_type_list(self):
        k = '_atom_site_occupancy'
        self.assertEqual(self.b.get_as_type(k, typ=str), ['1', '1'])
        self.assertEqual(self.b.get_as_type(k, typ=int), [1, 1])
        u_typ = uncertainties.ufloat_fromstr
        self.assertEqual(repr(self.b.get_as_type(k, typ=u_typ)), repr([U1, U1]))

    def test_get_as_type_exceptions(self):
        with self.assertRaises(KeyError):
            _ = self.b.get_as_type('_nonexistent_key', str)
        self.b['_cell_angle_beta'] = U90
        with self.assertRaises(TypeError):
            _ = self.b.get_as_type('_cell_angle_beta', str)

    def test_get_as_type_defaults(self):
        self.assertEqual(
            self.b.get_as_type('_nonexistent_key_with_default', str, 'default'),
            'default'
        )


class TestCifBlockInterface(unittest.TestCase):
    b = CifBlock()

    def test_base_frame_fill_from_cif_block_robust(self):
        self.b.read(path=nacl_cif_path, block='NaCl')
        base = BaseFrame()
        base.fill_from_cif_block(self.b, fragile=False)
        self.assertAlmostEqual(base.v_d, 179.51018149594705)
        self.assertAlmostEqual(base.orientation[0, 0], -0.0617076000)

    def test_base_frame_fill_from_cif_block_fragile(self):
        self.b.read(path=nacl_cif_path, block='NaCl_olex2_C2/m')
        base = BaseFrame()
        with self.assertRaises(KeyError):
            base.fill_from_cif_block(self.b, fragile=True)


class TestCifValidator(unittest.TestCase):
    def test_validator(self):
        c = CifValidator()
        self.assertIn('_atom_site_label', c)
        self.assertIsInstance(c['_atom_site_label'], CifBlock)


class TestCifFrame(unittest.TestCase):
    c_cif = CifFrame()
    c_fcf = CifFrame()

    def test_read_cif_file(self):
        self.c_cif.read(nacl_cif_path)
        self.assertIn('NaCl', self.c_cif)
        self.assertIsInstance(self.c_cif['NaCl'], CifBlock)

    def test_read_fcf_file(self):
        self.c_fcf.read(nacl_fcf_path)
        self.assertIn('NaCl', self.c_fcf)
        self.assertIsInstance(self.c_fcf['NaCl'], CifBlock)


class TestHklFrame(unittest.TestCase):
    h1 = HklFrame()
    h2 = HklFrame()

    @classmethod
    def setUpClass(cls) -> None:
        cls.h1.read(nacl_hkl_path, hkl_format='shelx_4')
        cls.h1.edit_cell(a=5.64109, b=5.64109, c=5.64109)

    def setUp(self) -> None:
        self.h2 = copy.deepcopy(self.h1)

    def test_len(self):
        self.assertEqual(len(self.h1), 8578)

    def test_add(self):
        h3 = self.h1 + self.h2
        self.assertEqual(len(h3), 17156)
        self.assertIsInstance(h3, HklFrame)

    def test_str(self):
        str_h = str(self.h1)
        self.assertIn('h', str_h[:300])
        self.assertIn('52001.5', str_h[:300])
        self.assertIn('419.465', str_h[-300:])

    def test_read(self):
        self.h2.read(nacl_hkl_path, hkl_format='free_4')
        self.assertEqual(self.h2.table.__len__(), 8578)

    def test_la_and_r_lim(self):
        self.h2.la = 0.50
        self.assertAlmostEqual(self.h2.la, 0.50)
        self.assertAlmostEqual(self.h2.r_lim, 4.0)
        self.h2.la = 'AgKa'
        self.assertAlmostEqual(self.h2.la, 0.5608)
        self.assertAlmostEqual(self.h2.r_lim, 3.56633380884)

    def test_copy(self):
        h = self.h1.copy()
        self.assertEqual(str(h), str(self.h1))

    def test_find_equivalents(self):
        self.h2.find_equivalents()
        self.assertEqual(self.h2.table['equiv'].nunique(), 2861)
        self.h2.find_equivalents(point_group=PG['m-3m'])
        self.assertEqual(self.h2.table['equiv'].nunique(), 111)

    def test_place(self):
        self.h2.place()
        xyz = self.h2.table.loc[:, ['x', 'y', 'z']].to_numpy()
        expected_xyz_sum = np.array([882.09902696, 1188.0682634, 516.21229231])
        expected_xyz_mean = np.array([0.102832710, 0.1385017800, 0.0601786300])
        self.assertAlmostEqual(xyz.max(), 2.4817898668519733)
        self.assertAlmostEqual(xyz.sum(), 2586.379582669306)
        self.assertAlmostEqual(sum(xyz.sum(axis=0) - expected_xyz_sum), 0.0)
        self.assertAlmostEqual(sum(xyz.mean(axis=0) - expected_xyz_mean), 0.0)

    def test_trim(self):
        self.h2.place()
        self.h2.trim(limit=1.2)
        self.h2.find_equivalents(point_group=PG['m-3m'])
        self.assertEqual(len(self.h2), 1761)
        self.assertEqual(self.h2.table['equiv'].nunique(), 18)

    def test_dac_trim(self):
        self.h2.place()
        self.h2.dac_trim(opening_angle=30, vector=np.array([1, 0, 0]))
        self.h2.find_equivalents(point_group=PG['m-3m'])
        self.assertEqual(len(self.h2), 434)
        self.assertEqual(self.h2.table['equiv'].nunique(), 10)

    def test_merge(self):
        self.h2.merge(point_group=PG['2/m'])
        self.assertEqual(len(self.h2.table), 831)
        self.h2.merge(point_group=PG['4/m'])
        self.assertEqual(len(self.h2.table), 390)
        self.h2.merge(point_group=PG['m-3'])
        self.assertEqual(len(self.h2.table), 159)
        self.h2.merge(point_group=PG['m-3m'])
        self.assertEqual(len(self.h2.table), 111)


class TestUBaseFrame(unittest.TestCase):

    def setUp(self) -> None:
        self.b = UBaseFrame()
        self.b.edit_cell(a=U6, b=U7, c=U8, al=U60, be=U70, ga=U80)

    def test_init(self):
        b = UBaseFrame()
        self.assertAlmostEqual(b.a_d.n, 1)
        self.assertAlmostEqual(b.b_d.n, 1)
        self.assertAlmostEqual(b.c_d.n, 1)
        self.assertAlmostEqual(b.al_d.n, np.pi / 2)
        self.assertAlmostEqual(b.be_d.n, np.pi / 2)
        self.assertAlmostEqual(b.ga_d.n, np.pi / 2)

    def test_edit_cell(self):
        self.assertAlmostEqual(self.b.a_d.n, 6)
        self.assertAlmostEqual(self.b.b_d.n, 7)
        self.assertAlmostEqual(self.b.c_d.n, 8)
        self.assertAlmostEqual(self.b.al_d.n, RAD60)
        self.assertAlmostEqual(self.b.be_d.n, RAD70)
        self.assertAlmostEqual(self.b.ga_d.n, RAD80)
        with self.assertRaises(KeyError):
            BaseFrame().edit_cell(alpha=0.0)

    def test_reciprocal_cell(self):
        self.assertAlmostEqual(self.b.a_r.n, 0.1773638940193107)
        self.assertAlmostEqual(self.b.b_r.n, 0.1649580865234412)
        self.assertAlmostEqual(self.b.c_r.n, 0.1512680839136182)
        self.assertAlmostEqual(self.b.al_r.n, 2.067032910195549)
        self.assertAlmostEqual(self.b.be_r.n, 1.874672323582712)
        self.assertAlmostEqual(self.b.ga_r.n, 1.574038054665492)

    def test_vectors(self):
        expected_a_v = np.array([6.00000000, 0.00000000, 0.00000000])
        expected_b_v = np.array([1.21553724, 6.89365427, 0.00000000])
        expected_c_v = np.array([2.73616115, 3.57924741, 6.61077984])
        expected_a_w = np.array([0.16666667, -0.02938783, -0.05307098])
        expected_b_w = np.array([0.00000000, 0.14506094, -0.07853975])
        expected_c_w = np.array([0.00000000, 0.00000000, 0.15126808])
        self.assertAlmostEqual(sum(self.b.a_v - expected_a_v).n, 0.)
        self.assertAlmostEqual(sum(self.b.b_v - expected_b_v).n, 0.)
        self.assertAlmostEqual(sum(self.b.c_v - expected_c_v).n, 0.)
        self.assertAlmostEqual(sum(self.b.a_w - expected_a_w).n, 0.)
        self.assertAlmostEqual(sum(self.b.b_w - expected_b_w).n, 0.)
        self.assertAlmostEqual(sum(self.b.c_w - expected_c_w).n, 0.)

    def test_perpendicular(self):
        def norm(v):
            return v / (v.dot(v) ** (1/2))
        av, bv, cv = self.b.a_v, self.b.b_v, self.b.c_v
        aw, bw, cw = self.b.a_w, self.b.b_w, self.b.c_w
        self.assertAlmostEqual(sum(norm(np.cross(av, bv)) - norm(cw)).n, 0.)
        self.assertAlmostEqual(sum(norm(np.cross(bv, cv)) - norm(aw)).n, 0.)
        self.assertAlmostEqual(sum(norm(np.cross(cv, av)) - norm(bw)).n, 0.)
        self.assertAlmostEqual(sum(norm(np.cross(aw, bw)) - norm(cv)).n, 0.)
        self.assertAlmostEqual(sum(norm(np.cross(bw, cw)) - norm(av)).n, 0.)
        self.assertAlmostEqual(sum(norm(np.cross(cw, aw)) - norm(bv)).n, 0.)

    def test_matrices(self):
        v = np.array([1.6180339887, 2.7182818285, 3.1415926536])
        self.assertAlmostEqual((v @ self.b.G_d @ v).n, 1797.2495106787824)
        self.assertAlmostEqual((v @ self.b.G_r @ v).n, 0.2238961843294296)
        self.assertAlmostEqual(sum(sum(self.b.G_d @ self.b.A_r - self.b.A_d)).n,
                               0.)
        self.assertAlmostEqual(sum(sum(self.b.G_r @ self.b.A_d - self.b.A_r)).n,
                               0.)

    def test_volumes(self):
        self.assertAlmostEqual(self.b.v_d.n, 273.4345841924758)
        self.assertAlmostEqual(self.b.v_r.n, 0.003657181855591761)


if __name__ == '__main__':
    unittest.main()
