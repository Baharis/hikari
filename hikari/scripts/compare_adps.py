import numpy as np
from hikari.dataframes import BaseFrame, CifFrame
from hikari.utility import make_abspath


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
    print(u1_frac)
    print(u2_frac)

    # calculate the similarity index itself
    r12_num = 2 ** (3 / 2) * np.linalg.det(u1_cart_inv @ u2_cart_inv) ** (1 / 4)
    r12_den = np.linalg.det(u1_cart_inv + u2_cart_inv) ** (1 / 2)
    return 100 * (1 - r12_num / r12_den)


# def compare_adps(cif1_path, cif1_block, cif2_path, cif2_block):
#     c1 = CifFrame()
#     c1.read(make_abspath(cif1_path), datablock=cif1_block)
#     c2 = CifFrame()
#     c2.read(make_abspath(cif2_path), datablock=cif2_block)
#     print(c1.data['_atom_site_fract_x'])
#     print(c1.data['_refine_ls_weighting_details'])


if __name__ == '__main__':
    u1 = [0.050, 0.021, 0.042, -0.006, -0.008, 0.005]
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
    u2 = [0.051, 0.02, 0.03, 0.006, -0.01, -0.012]  # a
    u2 = [0.050514, 0.019578, 0.030408, 0.006249, -0.011912, -0.010335]  # b
    print(compare_adp(9.135, 8.814, 21.397, 90, 93.010, 90, u1, u2))
    # compare_adps('~/x/HiPHAR/anders_script/rfpirazB_100K_SXD.cif',
    #              'rfpirazB_100K_SXD',
    #              '~/x/HiPHAR/anders_script/RFpirazB_cplt100.fractional.cif1',
    #              'RFpirazB_cplt100')
