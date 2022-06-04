import numpy as np
from uncertainties import ufloat_fromstr, unumpy

from hikari.dataframes import BaseFrame, CifFrame
from hikari.utility import make_abspath, cfloat, det3x3


def compare_adp(a, b, c, al, be, ga, adp1, adp2):
    # read the unit cell
    f = BaseFrame()
    f.edit_cell(a=a, b=b, c=c, al=al, be=be, ga=ga)

    # read matrices from input
    u1_frac = np.array([[adp1[0], adp1[5], adp1[4]],
                        [adp1[5], adp1[1], adp1[3]],
                        [adp1[4], adp1[3], adp1[2]]], dtype=float)
    u2_frac = np.array([[adp2[0], adp2[5], adp2[4]],
                        [adp2[5], adp2[1], adp2[3]],
                        [adp2[4], adp2[3], adp2[2]]], dtype=float)

    # orthogonalise matrices and calculate inverse of cartesian adps
    def adp_frac2cart(f_, u_):
        n = np.diag([f_.a_r, f_.b_r, f_.c_r])
        return f_.A_d.T @ n @ u_ @ n @ f_.A_d
    u1_cart_inv = np.linalg.inv(adp_frac2cart(f, u1_frac))
    u2_cart_inv = np.linalg.inv(adp_frac2cart(f, u2_frac))

    # calculate the similarity index itself
    r12_num = 2 ** (3 / 2) * np.linalg.det(u1_cart_inv @ u2_cart_inv) ** (1 / 4)
    r12_den = np.linalg.det(u1_cart_inv + u2_cart_inv) ** (1 / 2)
    return 100 * (1 - r12_num / r12_den)


def compare_adps(cif1_path, cif1_block, cif2_path, cif2_block,
                 deviations=True):

    u_type = ufloat_fromstr if deviations else cfloat

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
        if deviations:
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
