import os
from pathlib import Path
import shutil

import numpy as np
import pandas as pd
import seaborn as sns
from matplotlib import pyplot, colors, cm
from mpl_toolkits.mplot3d import art3d
from numpy import linalg as lin

from hikari.dataframes import HklFrame, LstFrame
from hikari.symmetry import SG, Group
from hikari.utility import make_abspath, mpl_map_palette, gnuplot_map_palette, \
    fibonacci_sphere, rotation_around
from hikari.resources import potency_map_template


def r1_map(a, b, c, al, be, ga,
           space_group='P1',
           fix_scale=False,
           opening_angle=35,
           job_directory='~',
           job_name='cplt_map',
           output_quality=3,
           resolution=1.2,
           wavelength='MoKa'):
    """ An elementary script to map R1 as a function of dac orientation """
    hkl_path = make_abspath(job_directory, job_name + '.hkl')
    res_path = make_abspath(job_directory, job_name + '.res')
    dat_path = make_abspath(job_directory, job_name + '.dat')
    gnu_path = make_abspath(job_directory, job_name + '.gnu')
    lst_path = make_abspath(job_directory, job_name + '.lst')
    png_path = make_abspath(job_directory, job_name + '.png')
    sg = SG[space_group]
    lg = sg.reciprocate()

    def _read_hkl_frame():
        """Make ball or axis of hkl which will be cut in further steps"""
        _f = HklFrame()
        _f.edit_cell(a=a, b=b, c=c, al=al, be=be, ga=ga)
        _f.la = wavelength
        _f.read(hkl_path=hkl_path, hkl_format='shelx_4')
        _f.trim(limit=min(_f.r_lim, resolution))
        _f.find_equivalents(point_group=lg)
        return _f

    p = _read_hkl_frame()
    total_reflections = p.table['equiv'].nunique()

    def _determine_theta_and_phi_limits():
        """Define range of coordinates where potency map will be calculated.
        v1, v2, v3 are normal vectors pointing in zenith direction z*,
        orthogonal direction (x) and direction perpendicular to them both."""
        _v1 = p.c_w / lin.norm(p.c_w)
        _v2 = p.a_v / lin.norm(p.a_v)
        _v3 = np.cross(_v1, _v2)

        if sg.system in {Group.System.triclinic, Group.System.monoclinic,
                         Group.System.orthorhombic, Group.System.tetragonal,
                         Group.System.cubic}:
            _th_limits = [0, 90]
            _ph_limits = [0, 90]
        elif sg.system in {Group.System.trigonal,
                           Group.System.hexagonal}:
            _th_limits = [0, 90]
            _ph_limits = [0, 120]
        else:
            raise ValueError('Unknown crystal system (trigonal not supported)')
        return _v1, _v2, _v3, _th_limits, _ph_limits

    v1, v2, v3, th_limits, ph_limits = _determine_theta_and_phi_limits()

    def _make_theta_and_phi_mesh():
        """Define a list of theta and phi values to be investigated."""
        if output_quality not in {1, 2, 3, 4, 5}:
            raise KeyError('output_quality should be 1, 2, 3, 4 or 5')
        angle_res = {1: 15, 2: 10, 3: 5, 4: 2, 5: 1}[output_quality]
        _th_range = np.arange(th_limits[0], th_limits[1] + 0.001, angle_res)
        _ph_range = np.arange(ph_limits[0], ph_limits[1] + 0.001, angle_res)
        _th_mesh, _ph_mesh = np.meshgrid(_th_range, _ph_range)
        return _th_range, _ph_range, _th_mesh, _ph_mesh

    th_range, ph_range, th_mesh, ph_mesh = _make_theta_and_phi_mesh()
    data_dict = {'th': [], 'ph': [], 'cplt': [], 'r1': []}

    def _angles_to_vector(theta, phi):
        """Find the vector by rotating v1 by theta and then phi, in degrees."""
        _sin_th = np.sin(np.deg2rad(theta))
        _sin_ph = np.sin(np.deg2rad(phi))
        _cos_th = np.cos(np.deg2rad(theta))
        _cos_ph = np.cos(np.deg2rad(phi))
        _v = v1 * _cos_th + v2 * _sin_th
        _v_parallel = np.dot(_v, v1) * v1
        _v_perpend = lin.norm(_v - _v_parallel) * (_cos_ph * v2 + _sin_ph * v3)
        return _v_parallel + _v_perpend

    def _calculate_completeness_mesh():
        """Calculate completeness for each individual pair of theta and phi."""
        _r1_mesh = np.zeros_like(th_mesh)
        lst = open(lst_path, 'w+')
        lst.write('#     th      ph    cplt      R1\n')
        vectors = np.vstack([_angles_to_vector(th, ph) for th in th_range
                                             for ph in ph_range])

        for i, th in enumerate(th_range):
            for j, ph in enumerate(ph_range):
                subdir = job_name + f'_th{int(th+.1)}_ph{int(ph+.1)}'
                hkl_path2 = make_abspath(job_directory, subdir, job_name+'.hkl')
                ins_path2 = make_abspath(job_directory, subdir, job_name+'.ins')
                lst_path2 = make_abspath(job_directory, subdir, job_name+'.lst')
                dir_path2 = make_abspath(job_directory, subdir)
                Path(dir_path2).mkdir()
                shutil.copy(res_path, ins_path2)
                q = p.copy()
                q.dac_trim(opening_angle, _angles_to_vector(th, ph))
                q.write(hkl_path2, hkl_format='shelx_4')
                count = q.table['equiv'].nunique()
                os.system(f'cd {dir_path2}; shelxl {job_name}')
                r1 = LstFrame().read_r1(lst_path2)
                cplt = count / total_reflections
                data_dict['th'].append(th)
                data_dict['ph'].append(ph)
                data_dict['cplt'].append(cplt)
                data_dict['r1'].append(r1)
                lst.write(f'{th:8.0f}{ph:8.0f}{cplt:8.5f}{r1:8.4}\n')
                _r1_mesh[j][i] = r1
            lst.write('\n')
        index_max = np.unravel_index(np.argmax(_r1_mesh), _r1_mesh.shape)
        best_th, best_ph = th_range[index_max[1]], ph_range[index_max[0]]
        index_min = np.unravel_index(np.argmin(_r1_mesh), _r1_mesh.shape)
        worst_th, worst_ph = th_range[index_min[1]], ph_range[index_min[0]]
        max_p = max(data_dict['cplt'])
        min_p = min(data_dict['cplt'])
        mean_p = np.mean(data_dict['cplt'])
        lst.write(f'# best cplt = {max_p} for th = {best_th}, ph = {best_ph}\n'
                  f'# worst cplt = {min_p} for th = {worst_th}, ph = {worst_ph}'
                  f'\n# mean cplt = {mean_p}\n')
        lst.close()
        np.savetxt(dat_path, _r1_mesh)
        return data_dict, _r1_mesh

    data_dict, r1_mesh = _calculate_completeness_mesh()

    def _plot_in_matplotlib():
        """Plot the completeness map in radial coordinates using matplotlib"""
        fig = pyplot.figure(figsize=(5, 3))
        ax = fig.add_subplot(111, projection='3d')
        ax.view_init(elev=90 - sum(th_limits) / 2, azim=sum(ph_limits) / 2)
        ax.dist = 10
        ax.set_axis_off()

        # surface in cartesian coordinates
        x = np.sin(np.deg2rad(th_mesh)) * np.cos(np.deg2rad(ph_mesh))
        y = np.sin(np.deg2rad(th_mesh)) * np.sin(np.deg2rad(ph_mesh))
        z = np.cos(np.deg2rad(th_mesh))

        # wireframe
        np.warnings.filterwarnings('ignore',
                                   category=np.VisibleDeprecationWarning)
        ax.plot_wireframe(x, y, z, colors='k', linewidth=0.25)  # <- mpl warning

        # color map
        my_heatmap_colors = mpl_map_palette['']
        my_colormap = colors.LinearSegmentedColormap.from_list(
            'heatmapEX', my_heatmap_colors, N=256)
        m = cm.ScalarMappable(cmap=my_colormap)
        m.set_array(r1_mesh)
        if fix_scale is True:
            m.set_clim(0, 1)
            norm = colors.Normalize(vmin=0, vmax=1)
        else:
            norm = colors.Normalize(vmin=min(data_dict['cplt']),
                                    vmax=max(data_dict['cplt']))
        pyplot.colorbar(m, fraction=0.046, pad=0.04)

        # direction lines
        _len = 1.25
        _x = p.a_w / lin.norm(p.a_w)
        _y = p.b_w / lin.norm(p.b_w)
        _z = p.c_w / lin.norm(p.c_w)
        ax.add_line(art3d.Line3D((_x[0], _len * _x[0]), (_x[1], _len * _x[1]),
                                 (_x[2], _len * _x[2]), color='r', linewidth=5))
        ax.add_line(art3d.Line3D((_y[0], _len * _y[0]), (_y[1], _len * _y[1]),
                                 (_y[2], _len * _y[2]), color='g', linewidth=5))
        ax.add_line(art3d.Line3D((_z[0], _len * _z[0]), (_z[1], _len * _z[1]),
                                 (_z[2], _len * _z[2]), color='b', linewidth=5))

        # color mesh for heatmap
        color_mesh = r1_mesh[:-1, :-1]
        for i in range(r1_mesh.shape[0] - 1):
            for j in range(r1_mesh.shape[1] - 1):
                color_mesh[i, j] = (r1_mesh[i + 1, j] + r1_mesh[i, j + 1] +
                                    r1_mesh[i, j] + r1_mesh[
                                        i + 1, j + 1]) / 4

        # heatmap surface
        for item in [fig, ax]:
            item.patch.set_visible(False)
        ax.plot_surface(x, y, z, rstride=1, cstride=1, cmap=my_colormap,
                        linewidth=0, antialiased=False,
                        facecolors=my_colormap(norm(color_mesh)))
        pyplot.savefig(png_path, dpi=600, format='png', bbox_inches=None)
    _plot_in_matplotlib()

    def _prepare_gnuplot_input():
        """Prepare input to completeness map in radial coordinates in gnuplot"""
        gnu = open(gnu_path, 'w+')
        gnu.write(potency_map_template.format(
            axis_x1=(p.a_w / lin.norm(p.a_w))[0],
            axis_x2=(p.a_w / lin.norm(p.a_w))[1],
            axis_x3=(p.a_w / lin.norm(p.a_w))[2],
            axis_y1=(p.b_w / lin.norm(p.b_w))[0],
            axis_y2=(p.b_w / lin.norm(p.b_w))[1],
            axis_y3=(p.b_w / lin.norm(p.b_w))[2],
            axis_z1=(p.c_w / lin.norm(p.c_w))[0],
            axis_z2=(p.c_w / lin.norm(p.c_w))[1],
            axis_z3=(p.c_w / lin.norm(p.c_w))[2],
            cplt_min=0 if fix_scale else min(data_dict['r1']),
            cplt_max=1 if fix_scale else max(data_dict['r1']),
            job_name=job_name,
            min_ph=min(ph_limits),
            max_ph=max(ph_limits),
            min_th=min(th_limits),
            max_th=max(th_limits),
            palette=gnuplot_map_palette['']))

    _prepare_gnuplot_input()
    try:
        from os import system, getcwd
        _path = make_abspath(job_directory)
        os.system('cd ' + _path + '; gnuplot ' + gnu_path)
    except OSError:
        pass


if __name__ == '__main__':
    r1_map(5.64109, 5.64109, 5.64109, 90, 90, 90, space_group='Fm-3m',
           job_directory='~/_/NaCl', job_name='NaCl',
           output_quality=5, wavelength='MoKa')
    pass
