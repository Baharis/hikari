import io
import pathlib
import sys
import tempfile
import unittest
from unittest import mock

import numpy as np

from hikari.scripts import calculate_similarity_indices, potency_map, \
    completeness_statistics, dac_statistics, reformat_hkl, simulate_dac


nacl_cif_path = str(pathlib.Path(__file__).parent.joinpath('NaCl.cif'))
nacl_hkl_path = str(pathlib.Path(__file__).parent.joinpath('NaCl.hkl'))
nacl_commons = {'a': 5.64109, 'b': 5.64109, 'c': 5.64109,
                'al': 90, 'be': 90, 'ga': 90, 'input_path': nacl_hkl_path,
                'input_format': 'shelx_4', 'input_wavelength': 'MoKa'}


class TestCompareADPsScripts(unittest.TestCase):
    temp_dir = tempfile.TemporaryDirectory()
    output1_path = str(pathlib.Path(temp_dir.name) / 'similarity1.out')
    output2_path = str(pathlib.Path(temp_dir.name) / 'similarity2.out')

    def test_calculate_similarity_index_for_implicit_cif2_path(self):
        calculate_similarity_indices(cif1_path=nacl_cif_path,
                                     output_path=self.output1_path)
        calculate_similarity_indices(cif1_path=nacl_cif_path,
                                     cif2_path=nacl_cif_path,
                                     cif1_block='NaCl',
                                     cif2_block='NaCl_olex2_C2/m',
                                     output_path=self.output2_path)
        with open(self.output1_path, 'r') as out1:
            with open(self.output2_path, 'r') as out2:
                self.assertEqual(out1.readlines()[4:], out2.readlines()[4:])

    def test_calculate_similarity_index_for_implicit_cif_blocks(self):
        calculate_similarity_indices(cif1_path=nacl_cif_path,
                                     cif2_path=nacl_cif_path,
                                     output_path=self.output1_path)
        calculate_similarity_indices(cif1_path=nacl_cif_path,
                                     cif2_path=nacl_cif_path,
                                     cif1_block='NaCl',
                                     cif2_block='NaCl',
                                     output_path=self.output2_path)
        with open(self.output1_path, 'r') as out1:
            with open(self.output2_path, 'r') as out2:
                self.assertEqual(out1.readlines()[4:], out2.readlines()[4:])

    def test_calculate_similarity_index_for_identical(self):
        calculate_similarity_indices(cif1_path=nacl_cif_path,
                                     cif2_path=nacl_cif_path,
                                     output_path=self.output1_path)
        with open(self.output1_path, 'r') as out:
            last_line_contents = out.readlines()[-1].strip().split()
            self.assertAlmostEqual(float(last_line_contents[-1]), 0.0)
            self.assertAlmostEqual(float(last_line_contents[-3]), 0.0)

    def test_calculate_similarity_index_for_different(self):
        calculate_similarity_indices(cif1_path=nacl_cif_path,
                                     output_path=self.output1_path)
        with open(self.output1_path, 'r') as out:
            last_line_contents = out.readlines()[-1].strip().split()
            self.assertNotAlmostEqual(float(last_line_contents[-1]), 0.0)
            self.assertNotAlmostEqual(float(last_line_contents[-3]), 0.0)

    def test_calculate_similarity_index_for_negative_adps(self):
        calculate_similarity_indices(cif1_path=nacl_cif_path,
                                     cif1_block='NaCl',
                                     cif2_block='NaCl_negative_ADPs',
                                     output_path=self.output1_path)
        with open(self.output1_path, 'r') as out:
            last_line_contents = out.readlines()[-1].strip().split()
            self.assertAlmostEqual(float(last_line_contents[-3]), 50.0)


class TestHklScripts(unittest.TestCase):
    temp_dir = tempfile.TemporaryDirectory()
    hkl_path = str(pathlib.Path(temp_dir.name) / 'temp.hkl')

    @mock.patch('sys.stdout', new_callable=io.StringIO)
    def get_stdout(self, callable_, kwargs, mock_stdout):
        callable_(**kwargs)
        return mock_stdout.getvalue()

    @unittest.skipIf(sys.version_info[:3] == (3, 11, 0), 'Faulty init.tcl')
    def test_potency_map_simple(self):
        potency_map(a=10, b=10, c=10, al=90, be=90, ga=90, space_group='P1',
                    path=self.hkl_path, output_quality=2, wavelength='MoKa')

    @unittest.skipIf(sys.version_info[:3] == (3, 11, 0), 'Faulty init.tcl')
    def test_potency_map_with_focus(self):
        ori = np.array([[-0.08263, +0.00536, -0.03667],
                        [-0.01671, -0.03679, -0.03238],
                        [+0.06030, +0.02438, -0.04088]])
        potency_map(a=10, b=10, c=10, al=90, be=90, ga=90, space_group='Pm-3m',
                    path=self.hkl_path, output_quality=2, wavelength='MoKa',
                    orientation=ori)

    def test_completeness_statistics(self):
        kwargs = dict({'space_group': 'Fm-3m'}, **nacl_commons)
        stdout = self.get_stdout(completeness_statistics, kwargs)
        line = '(2.373, 2.468]    287      6       6  11.039363'
        self.assertIn(line, stdout)
        # TODO: for some reason, on ubuntu latest ONLY, one refl is not read?
        # TODO: check that out eventually - the last reflection not read

    def test_dac_statistics(self):
        kwargs = dict({'opening_angle': 35, 'orientation': np.eye(3),
                       'resolution': 1.2}, **nacl_commons)
        stdout = self.get_stdout(dac_statistics, kwargs)
        line = '0.556991        26        60       122  0.433333  0.213115'
        self.assertIn(line, stdout)

    def test_reformat_hkl(self):
        reformat_hkl(input_path=nacl_hkl_path, input_format='shelx_4',
                     output_path=self.hkl_path, output_format='shelx_40')
        with open(nacl_hkl_path, 'r') as f1:
            with open(self.hkl_path, 'r') as f2:
                self.assertEqual(len(f1.readlines()), len(f2.readlines()))

    def test_simulate_dac(self):
        kwargs = dict({'opening_angle': 35, 'orientation': np.eye(3),
                       'resolution': 1.2, 'output_path': self.hkl_path,
                       'output_format': 'shelx_4'}, **nacl_commons)
        simulate_dac(**kwargs)
        with open(self.hkl_path, 'r') as f:
            self.assertEqual(len(f.readlines()), 553)


if __name__ == '__main__':
    unittest.main()
