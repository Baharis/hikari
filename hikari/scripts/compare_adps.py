from typing import Dict, Any, Iterable, Tuple, List

import matplotlib.pyplot as plt
from matplotlib import animation
from mpl_toolkits.mplot3d import axes3d
import numpy as np
import numpy.linalg as lin
from uncertainties import ufloat, ufloat_fromstr, unumpy, UFloat

from hikari.dataframes import BaseFrame, CifFrame, UBaseFrame
from hikari.utility import make_abspath, cfloat, det3x3, chemical_elements, \
    rotation_around


def calculate_similarity_indices(cif1_path: str,
                                 cif2_path: str = None,
                                 cif1_block: str = None,
                                 cif2_block: str = None,
                                 output_path: str = None,
                                 normalize: bool = False,
                                 uncertainties: bool = True) -> None:
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

    For a single cif file 'glycine.cif' with two data blocks named '100K' and
    'RT', the three following commands will have the same effect:

    :Example:

    >>> calculate_similarity_indices('glycine.cif')
    >>> calculate_similarity_indices('glycine.cif', 'glycine.cif')
    >>> calculate_similarity_indices('glycine.cif', 'glycine.cif', '100K', 'RT')

    For two cif files 'X-rays.cif' and 'neutrons.cif' with data block named 'RT'
    and '100K' (one each), the two following commands will have the same effect:

    :Example:

    >>> calculate_similarity_indices('X-rays.cif', 'neutrons.cif')
    >>> calculate_similarity_indices('X-rays.cif', 'neutrons.cif', 'RT', '100K')

    The behaviour of script can be further altered via other parameters.
    If `output_file` is set True, the results will be written there instead of
    console. If `uncertainties` is set True, standard deviations of all ADPs
    as well as the unit cell parameters from 1st cif file will be assumed
    uncorrelated and used to estimate the uncertainty of every SI determination.

    For more information about Similarity Index (SI) itself, please consult
    Whitten and Spackman, Acta Cryst B **62**, 875 (2006)
    https://doi.org/10.1107/S0108768106020787.

    :param cif1_path: Absolute or relative path to the first cif file.
    :param cif2_path: Absolute or relative path to the second cif file.
        If not specified, it is assumed equal to `cif1_path`.
    :param cif1_block: Name of the first data block used in SI determination.
        It points to the data block inside the file specified by `cif1_path`.
        If not specified, the first block found in said file will be used.
    :param cif2_block: Name of the second data block used in SI determination.
        It points to the data block inside the file specified by `cif2_path`.
        If not specified, the first unused block in said file will be used.
    :param output_path: Path where the output of the program should be written.
    :param uncertainties: If True, propagate the standard deviations of
        individual ADPs' and cif1's unit cell to estimate SI's uncertainties.
    :param normalize: If True, equalize the volume of displacement ellipsoids
        by normalizing the determinants of ADP matrices expressed in cartesian
        coordinates. As a result, SI is a function of displacement "shape" only.
    """

    u_type = ufloat_fromstr if uncertainties else cfloat

    def initialise_output():
        if output_path is None:
            f_ = None
        else:
            f_ = open(make_abspath(output_path), 'w+')
        return f_
    f = initialise_output()

    def interpret_paths():
        cif1p = make_abspath(cif1_path)
        cif2p = cif1p if cif2_path is None else make_abspath(cif2_path)
        cif1b = 1 if cif1_block is None else cif1_block
        cif2b = 2 if cif1b is 1 and cif1p is cif2p \
            else 1 if cif2_block is None else cif2_block
        return cif1p, cif2p, cif1b, cif2b
    cif1_path, cif2_path, cif1_block, cif2_block = interpret_paths()

    print(f"# Hikari will calculate similarity indices for atoms "
          f"in the following blocks:", file=f)
    print(f"# - Block {cif1_block} of file {cif1_path}", file=f)
    print(f"# - Block {cif2_block} of file {cif2_path}", file=f)
    print(f"# with the following settings:", file=f)
    print(f"# - normalize={normalize}", file=f)
    print(f"# - uncertainties={uncertainties}", file=f)

    def read_cif_block(path, block):
        c = CifFrame()
        c.read(path)
        block = list(c.keys())[block - 1] if isinstance(block, int) else block
        return c[block]
    cif_block_1 = read_cif_block(cif1_path, cif1_block)
    cif_block_2 = read_cif_block(cif2_path, cif2_block)

    def find_common_labels(block1, block2):
        label_list1 = block1['_atom_site_aniso_label']
        label_list2 = block2['_atom_site_aniso_label']
        return [label for label in label_list1 if label in label_list2]
    common_labels = find_common_labels(cif_block_1, cif_block_2)
    print(f"# Number of matching atoms found: {len(common_labels)}", file=f)

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

    base_frame1 = UBaseFrame() if uncertainties else BaseFrame()
    base_frame2 = UBaseFrame() if uncertainties else BaseFrame()
    base_frame1.fill_from_cif_block(cif_block_1)
    base_frame2.fill_from_cif_block(cif_block_2)

    def calculate_similarity_index(adp_frac_1, adp_frac_2) -> float or UFloat:
        def adp_frac2cart(adp_frac, base_frame):
            n = np.diag([base_frame.a_r, base_frame.b_r, base_frame.c_r])
            return base_frame.A_d.T @ n @ adp_frac @ n @ base_frame.A_d
        try:
            adp_cart_1 = adp_frac2cart(adp_frac_1, base_frame1)
            adp_cart_2 = adp_frac2cart(adp_frac_2, base_frame2)
            if normalize:
                adp_cart_1 /= det3x3(adp_cart_1) ** (1 / 3)
                adp_cart_2 /= det3x3(adp_cart_2) ** (1 / 3)
            if uncertainties:
                adp_inv_1 = unumpy.matrix(adp_cart_1).I
                adp_inv_2 = unumpy.matrix(adp_cart_2).I
            else:
                adp_inv_1 = np.linalg.inv(adp_cart_1)
                adp_inv_2 = np.linalg.inv(adp_cart_2)
            r12_num = 2 ** (3 / 2) * det3x3(adp_inv_1 @ adp_inv_2) ** (1 / 4)
            r12_den = det3x3(adp_inv_1 + adp_inv_2) ** (1 / 2)
        except ValueError:  # if matrix cannot be inverted, return S = 1.0
            r12_num = ufloat(0.0, 0) if uncertainties else 0.0
            r12_den = ufloat(1.0, 0) if uncertainties else 1.0
        return 100 * (1 - r12_num / r12_den)

    def si_string(si: float or UFloat) -> str:
        return f'{si:8.4f}' if isinstance(si, float) \
            else f'{si.n:8.4f} +/- {si.s:8.4f}'

    def si_header():
        return '#  label     SI.n +/-     SI.s' \
            if uncertainties else '#  label       SI'
    print(si_header(), file=f)

    def calculate_and_print_similarity_indices():
        sis: Dict[Any, float or UFloat] = {}
        for k in common_labels:
            sis[k] = calculate_similarity_index(adp_frac_dict_1[k],
                                                adp_frac_dict_2[k])
            print(f'{k:>8} {si_string(sis[k])}', file=f)
        sis_all = [v for k, v in sis.items()]
        sis_h = [v for k, v in sis.items()
                 if k[0] is 'H' and k[:2] not in chemical_elements]
        if sis_all:
            avg_si_all = sum(sis_all) / len(sis_all)
            print(f'# avg(*) {si_string(avg_si_all)}', file=f)
        else:
            print(f'# No atoms with matching names found', file=f)
        if sis_h:
            avg_si_h = sum(sis_h) / len(sis_h)
            print(f'# avg(H) {si_string(avg_si_h)}', file=f)
    calculate_and_print_similarity_indices()

    def terminate_output():
        if output_path is not None:
            f.close()
    terminate_output()


def animate_similarity_index(u_diag: Iterable,
                             transformations: List[np.ndarray],
                             output_path: str) -> None:
    """
    Make a gif presenting the change of similarity index against some arbitrary
    series of `transformations`. Initial displacement matrix must be diagonal,
    values given in `u_diag`, evaluated against a regular unit cell with a = 1.

    :param u_diag: triplet of u matrix diagonal elements
    :param transformations: list of transformations to be plotted every 0.1 sec
    :param output_path: path where resulting gif file will be saved
    """
    max_radius = 1.0

    fig = plt.figure(figsize=(12, 12))  # Square figure

    # Coefficients in a0 x**2 + a1 y**2 + a2 z**2 = 1 & corresponding radii
    rx, ry, rz = np.sqrt(u_diag)
    u = np.diag(u_diag)

    # Set of all spherical angles:
    v = np.linspace(0, 2 * np.pi, 100)
    w = np.linspace(0, np.pi, 100)

    # Cartesian coordinates that correspond to the spherical angles:
    # (this is the equation of an ellipsoid):
    x = rx * np.outer(np.cos(v), np.sin(w))
    y = ry * np.outer(np.sin(v), np.sin(w))
    z = rz * np.outer(np.ones_like(v), np.cos(w))

    # set animation settings
    fps = 10
    steps = len(transformations)

    def calculate_similarity_index(adp_frac_1, adp_frac_2) -> float or UFloat:
        base_frame = BaseFrame()

        def adp_frac2cart(adp_frac):
            n = np.diag([base_frame.a_r, base_frame.b_r, base_frame.c_r])
            return base_frame.A_d.T @ n @ adp_frac @ n @ base_frame.A_d
        adp_cart_1 = adp_frac2cart(adp_frac_1)
        adp_cart_2 = adp_frac2cart(adp_frac_2)
        adp_inv_1 = np.linalg.inv(adp_cart_1)
        adp_inv_2 = np.linalg.inv(adp_cart_2)
        r12_num = 2 ** (3 / 2) * det3x3(adp_inv_1 @ adp_inv_2) ** (1 / 4)
        r12_den = det3x3(adp_inv_1 + adp_inv_2) ** (1 / 2)
        return 100 * (1 - r12_num / r12_den)

    def get_xyz(step: int) -> Tuple[np.ndarray, np.ndarray, np.ndarray]:
        """Get positions of x, y, z nodes after `step` steps of animation"""
        trans_matrix = transformations[step]
        x_list = x.reshape(10_000)
        y_list = y.reshape(10_000)
        z_list = z.reshape(10_000)
        xyz = np.vstack([x_list, y_list, z_list]).T
        xyz_trans = xyz @ trans_matrix
        x_trans = xyz_trans[:, 0].reshape(100, 100)
        y_trans = xyz_trans[:, 1].reshape(100, 100)
        z_trans = xyz_trans[:, 2].reshape(100, 100)
        return x_trans, y_trans, z_trans

    def get_u(step: int) -> np.ndarray:
        """Get transformed u matrix after `step` steps of animation"""
        trans_matrix = transformations[step]
        return trans_matrix.T @ u @ trans_matrix

    # prepare list of similarity indices
    si = [calculate_similarity_index(u, get_u(i)) for i in range(steps)]
    si.append(si[0])

    def get_frame(step):
        """Clear previous axes and produce new at `step` step of animation"""
        for a in fig.axes:
            a.clear()
        x_, y_, z_ = get_xyz(step)
        ax1 = fig.add_subplot(221, projection='3d')
        ax2 = fig.add_subplot(222)
        ax3 = fig.add_subplot(223, projection='3d')
        ax4 = fig.add_subplot(224)
        ax1.plot_surface(x, y, z, rstride=4, cstride=4, color='r')
        ax2.plot(range(steps+1), si, 'b')
        ax2.set_xlabel(step)
        ax2.set_ylabel('similarity index')
        ax2.axvline(x=step, color='b')
        ax2.plot(step, si[step], 'bo')
        ax2.set_xlim([0, steps])
        ax2.set_ylim([0, 10])
        ax3.plot_surface(x_, y_, z_, rstride=4, cstride=4, color='b')
        u_ = get_u(step)
        ax4_font = {'family': 'serif', 'weight': 'normal', 'size': 16}
        ax4.clear()
        ax4.axis('off')
        ax4.text(0.2, 0.5, s='U = [', fontdict=ax4_font, ha='right')
        ax4.text(0.4, 0.7, s=f'{u_[0, 0]:6.3f}', fontdict=ax4_font, ha='right')
        ax4.text(0.6, 0.7, s=f'{u_[1, 0]:6.3f}', fontdict=ax4_font, ha='right')
        ax4.text(0.8, 0.7, s=f'{u_[2, 0]:6.3f}', fontdict=ax4_font, ha='right')
        ax4.text(0.4, 0.5, s=f'{u_[0, 1]:6.3f}', fontdict=ax4_font, ha='right')
        ax4.text(0.6, 0.5, s=f'{u_[1, 1]:6.3f}', fontdict=ax4_font, ha='right')
        ax4.text(0.8, 0.5, s=f'{u_[2, 1]:6.3f}', fontdict=ax4_font, ha='right')
        ax4.text(0.4, 0.3, s=f'{u_[0, 2]:6.3f}', fontdict=ax4_font, ha='right')
        ax4.text(0.6, 0.3, s=f'{u_[1, 2]:6.3f}', fontdict=ax4_font, ha='right')
        ax4.text(0.8, 0.3, s=f'{u_[2, 2]:6.3f}', fontdict=ax4_font, ha='right')
        ax4.text(0.9, 0.5, s=']', fontdict=ax4_font, ha='right')
        for axis in 'xyz':
            getattr(ax1, f'set_{axis}lim')((-max_radius, max_radius))
            getattr(ax3, f'set_{axis}lim')((-max_radius, max_radius))

    fig.tight_layout()
    fig.subplots_adjust(left=0.02, bottom=0.02, right=0.98, top=0.98,
                        wspace=None, hspace=None)
    ani = animation.FuncAnimation(
        fig=fig,
        func=get_frame,
        frames=steps,
        interval=1000/fps,
    )
    ani_path = make_abspath(output_path)
    ani.save(ani_path, writer='imagemagick', fps=fps)


if __name__ == '__main__':
    # calculate_similarity_indices('~/_/si/1.cif',
    #                              '~/_/si/2.cif',
    #                              output_path='~/_/si/output.txt')
    t1 = [rotation_around(np.array([0, 0, 1]), by=np.deg2rad(10 * i))
          for i in range(36)]
    t2 = [np.eye(3) * np.sqrt((1.5 - np.cos(np.deg2rad(10 * i)) / 2))
          for i in range(36)]
    t3 = [np.array([(1, 0, 0), (0, 1, 0), (0, 0, np.sqrt(2 ** np.sin(np.deg2rad(10 * i))))])
          for i in range(36)]
    animate_similarity_index(u_diag=(2, 1, 1),
                             transformations=t1,
                             output_path='~/_/si/S_animation_211_zrot.gif')
