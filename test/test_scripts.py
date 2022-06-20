import pathlib
import tempfile
import unittest

from hikari.scripts import calculate_similarity_indices


nacl_cif_path = str(pathlib.Path(__file__).parent.joinpath('NaCl.cif'))


class TestCompareADPsScripts(unittest.TestCase):
    temp_dir = tempfile.TemporaryDirectory()
    output1_path = str(pathlib.Path(temp_dir.name) / 'similarity.out')
    output2_path = str(pathlib.Path(temp_dir.name) / 'similarity.out')

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
                self.assertEqual(out1.read(), out2.read())

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
                self.assertEqual(out1.read(), out2.read())

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


if __name__ == '__main__':
    unittest.main()
