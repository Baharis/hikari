import os
from pathlib import Path
import shutil

import numpy as np
from numpy import linalg as lin

from hikari.dataframes import HklFrame, LstFrame
from hikari.symmetry import SG, Group
from hikari.utility import make_abspath, \
    GnuplotAngularHeatmapArtist, MatplotlibAngularHeatmapArtist


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
        lst.write('#     th      ph      R1    cplt\n')
        vectors = np.vstack([_angles_to_vector(th, ph) for th in th_range
                                             for ph in ph_range])
        #TODO rewrite using spherical2cartesian as in potency_map

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
                lst.write(f'{th:8.0f}{ph:8.0f}{r1:8.5}{cplt:8.5f}\n')
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

    ma = MatplotlibAngularHeatmapArtist()
    ma.x_axis = p.a_w / lin.norm(p.a_w)
    ma.y_axis = p.b_w / lin.norm(p.b_w)
    ma.z_axis = p.c_w / lin.norm(p.c_w)
    ma.heat_limits = (0 if fix_scale else min(data_dict['r1']),
                      1 if fix_scale else max(data_dict['r1']))
    ma.polar_limits = (min(th_limits), max(th_limits))
    ma.azimuth_limits = (min(ph_limits), max(ph_limits))
    ma.plot(png_path)

    # plot potency map using external gnuplot
    ga = GnuplotAngularHeatmapArtist()
    ga.x_axis = p.a_w / lin.norm(p.a_w)
    ga.y_axis = p.b_w / lin.norm(p.b_w)
    ga.z_axis = p.c_w / lin.norm(p.c_w)
    ga.heat_limits = (0 if fix_scale else min(data_dict['r1']) * 100,
                      1 if fix_scale else max(data_dict['r1']) * 100)
    ga.polar_limits = (min(th_limits), max(th_limits))
    ga.azimuth_limits = (min(ph_limits), max(ph_limits))
    ga.plot(png_path)


if __name__ == '__main__':
    r1_map(5.64109, 5.64109, 5.64109, 90, 90, 90, space_group='Fm-3m',
           job_directory='~/_/NaCl', job_name='NaCl',
           output_quality=5, wavelength='MoKa')
    pass
