import numpy as np
from uncertainties import ufloat_fromstr, unumpy

from hikari.dataframes import BaseFrame, CifFrame
from hikari.utility import make_abspath, cfloat, det3x3, chemical_elements


def compare_adps(cif1_path, cif2_path=None, cif1_block=None, cif2_block=None,
                 output_path=None, uncertainties=True, normalize=False):
    """
    Compare Anisotropic Displacement Parameters of all atoms in two cif blocks
    and return a Similarity Index (SI) for these pairs of atoms that exist
    and are defined anisotropic (via "_atom_site_aniso_U_**") in both files.

    This script requires two independent cif blocks to run,
    but the location of both blocks can be specified in multiple ways:

    - If only `cif1_path` is given, the first and second data block from
      cif1 file will be used in SI evaluation;
    - If `cif1_path` and `cif2_path` are given, the first data block from
      each cif file will be used instead;
    - If `cif1_path`, `cif2_path`, and `cif1_block` are given,
      `cif1_block` from each `cif1_path` and `cif2_path` cif files will be used.
    - If `cif1_path`, `cif2_path`, `cif1_block`, and `cif2_block` are given,
      `cif1_block` in `cif1_path` and `cif2_block` in `cif2_path` will be used.

    For a cif file 'glycine.cif' with two data blocks, '100K' and 'RT',
    the three following commands will all yield the same results.

    :Example:

    >>> compare_adps('glycine.cif')
    >>> compare_adps('glycine.cif', 'glycine.cif', '100K', '100K')

    Similarily, for two files 'X-rays.cif' and 'neutrons.cif', both with only
    one data block '100K', the two following commands will have the same effect:

    :Example:

    >>> compare_adps('X-rays.cif', 'neutrons.cif')
    >>> compare_adps('X-rays.cif', 'neutrons.cif', '100K', '100K')

    The behaviour of script can be further altered via other parameters.
    If `output_file` is set True, the results will be written there instead of
    console. If `uncertainties` is set True, standard deviations of all ADPs
    as well as the unit cell parameters from 1st cif file will be assumed
    uncorrelated and used to estimate the uncertainty of every SI determination.

    For more information about Similarity Index (SI) itself, please consult
    Whitten and Spackman, Acta Cryst B **62**, 875 (2006)
    https://doi.org/10.1107/S0108768106020787.

    :param cif1_path: Absolute or relative path to the first cif file.
    :type cif1_path: str
    :param cif2_path: Absolute or relative path to the second cif file.
        If not specified, it is assumed equal to `cif1_path`.
    :type cif2_path: str
    :param cif1_block: Name of the first data block used in SI determination.
        It points to the data block inside the file specified by `cif1_path`.
        If not specified, the first block found in said file will be used.
    :type cif1_block: str
    :param cif2_block: Name of the second data block used in SI determination.
        It points to the data block inside the file specified by `cif2_path`.
        If not specified, the first unused block in said file will be used.
    :type cif2_block: str
    :param output_path: Path where the output of the program should be written.
    :type output_path: str
    :param uncertainties: If True, propagate the standard deviations of
        individual ADPs' and cif1's unit cell to estimate SI's uncertainties.
    :type uncertainties: bool
    :param normalize: If True, equalize the volume of displacement ellipsoids
        by normalizing ADP matrices' determinants to a common determinant before
        calculating the SI. Please mind that this invalidates the uncertainties.
    :type normalize: bool
    """

    u_type = ufloat_fromstr if uncertainties else cfloat

    def read_cif_block(cif_path, cif_block):
        c = CifFrame()
        c.read(make_abspath(cif_path))
        return c[cif_block]
    cif_block_1 = read_cif_block(cif1_path, cif1_block)
    cif_block_2 = read_cif_block(cif2_path, cif2_block)

    b = BaseFrame()
    b.fill_from_cif_block(cif_block_1)

    def find_common_labels(block1, block2):
        label_list1 = block1['_atom_site_aniso_label']
        label_list2 = block2['_atom_site_aniso_label']
        return [label for label in label_list1 if label in label_list2]
    common_labels = find_common_labels(cif_block_1, cif_block_2)

    def make_adp_fractional_arrays(block):
        labels = block['_atom_site_aniso_label']
        u11s = block.get_as_type(key='_atom_site_aniso_U_11', typ=u_type)
        u22s = block.get_as_type(key='_atom_site_aniso_U_22', typ=u_type)
        u33s = block.get_as_type(key='_atom_site_aniso_U_33', typ=u_type)
        u12s = block.get_as_type(key='_atom_site_aniso_U_12', typ=u_type)
        u13s = block.get_as_type(key='_atom_site_aniso_U_13', typ=u_type)
        u23s = block.get_as_type(key='_atom_site_aniso_U_23', typ=u_type)
        return {k: np.array([[u11, u12, u13], [u12, u22, u23], [u13, u23, u33]])
                for k, u11, u22, u33, u12, u13, u23
                in zip(labels, u11s, u22s, u33s, u12s, u13s, u23s)}
    adp_frac_dict_1 = make_adp_fractional_arrays(cif_block_1)
    adp_frac_dict_2 = make_adp_fractional_arrays(cif_block_2)

    def calculate_similarity_index(adp_frac_1, adp_frac_2):
        def adp_frac2cart(adp_frac):
            n = np.diag([b.a_r, b.b_r, b.c_r])
            return b.A_d.T @ n @ adp_frac @ n @ b.A_d
        adp_cart_1 = adp_frac2cart(adp_frac_1)
        adp_cart_2 = adp_frac2cart(adp_frac_2)
        if uncertainties:
            adp_inv_1 = unumpy.matrix(adp_cart_1).I
            adp_inv_2 = unumpy.matrix(adp_cart_2).I
        else:
            adp_inv_1 = np.linalg.inv(adp_cart_1)
            adp_inv_2 = np.linalg.inv(adp_cart_2)
        r12_num = 2 ** (3 / 2) * det3x3(adp_inv_1 @ adp_inv_2) ** (1 / 4)
        r12_den = det3x3(adp_inv_1 + adp_inv_2) ** (1 / 2)
        return 100 * (1 - r12_num / r12_den)
    sis = []
    for k in common_labels:
        si = calculate_similarity_index(adp_frac_dict_1[k], adp_frac_dict_2[k])
        if 'H' in k:
            sis.append(si)
        print(f'{k:>4}: {si}')
    avg_si = sum(sis) / len(sis)
    print(avg_si)


if __name__ == '__main__':
    # u1 = [0.050, 0.021, 0.042, -0.006, -0.008, 0.005]
    # loop_
    # _atom_site_Cryst_ADP2_U_label
    # _atom_site_Cryst_ADP2_U_11
    # _atom_site_Cryst_ADP2_U_22
    # _atom_site_Cryst_ADP2_U_33
    # _atom_site_Cryst_ADP2_U_12
    # _atom_site_Cryst_ADP2_U_13
    # _atom_site_Cryst_ADP2_U_23
    # _atom_site_Cryst_ADP2_U_11_esu
    # _atom_site_Cryst_ADP2_U_22_esu
    # _atom_site_Cryst_ADP2_U_33_esu
    # _atom_site_Cryst_ADP2_U_12_esu
    # _atom_site_Cryst_ADP2_U_13_esu
    # _atom_site_Cryst_ADP2_U_23_esu
    #  H18  0.050514  0.019578  0.030408  -0.011912  -0.010335   0.006249
    #       0.007746  0.005374  0.006059   0.005030   0.005525   0.004431
    # u2 = [0.051, 0.02, 0.03, 0.006, -0.01, -0.012]  # a
    # u2 = [0.050514, 0.019578, 0.030408, 0.006249, -0.011912, -0.010335]  # b
    # print(compare_adp(9.135, 8.814, 21.397, 90, 93.010, 90, u1, u2))
    compare_adps('~/x/HiPHAR/anders_script/rfpirazB_100K_SXD.cif',
                 'rfpirazB_100K_SXD',
                 '~/x/HiPHAR/anders_script/RFpirazB_cplt100.fractional.cif1',
                 'RFpirazB_cplt100')
