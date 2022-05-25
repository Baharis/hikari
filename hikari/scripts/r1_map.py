import os
from pathlib import Path
import shutil

import numpy as np
from numpy import linalg as lin

from hikari.dataframes import HklFrame, LstFrame
from hikari.symmetry import SG, Group
from hikari.utility import make_abspath, sph2cart, weighted_quantile, \
    GnuplotAngularHeatmapArtist, MatplotlibAngularHeatmapArtist, Interval


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
        Unit vectors v1, v2, v3 point in zenith z*, orthogonal x and product."""
        _v1 = p.c_w / lin.norm(p.c_w)
        _v2 = p.a_v / lin.norm(p.a_v)
        _v3 = np.cross(_v1, _v2)

        if sg.system in {Group.System.triclinic}:
            _th_limits = Interval(0, 180)
            _ph_limits = Interval(-45, 135)
        elif sg.system in {Group.System.monoclinic}:
            _th_limits = Interval(0, 180)
            _ph_limits = Interval(0, 90)
        elif sg.system in {Group.System.orthorhombic, Group.System.tetragonal,
                         Group.System.cubic}:
            _th_limits = Interval(0, 90)
            _ph_limits = Interval(0, 90)
        elif sg.system in {Group.System.trigonal, Group.System.hexagonal}:
            _th_limits = Interval(0, 90)
            _ph_limits = Interval(0, 120)
        else:
            raise ValueError('Unknown crystal system (trigonal not supported)')
        return _v1, _v2, _v3, _th_limits, _ph_limits
    v1, v2, v3, th_limits, ph_limits = _determine_theta_and_phi_limits()

    def _determine_angle_res():
        if output_quality not in {1, 2, 3, 4, 5}:
            raise KeyError('output_quality should be 1, 2, 3, 4 or 5')
        return {1: 15, 2: 10, 3: 5, 4: 2, 5: 1}[output_quality]
    angle_res = _determine_angle_res()

    th_range = th_limits.arange(step=angle_res)
    ph_range = ph_limits.arange(step=angle_res)
    th_mesh, ph_mesh = th_limits.mesh_with(ph_limits, step=angle_res)
    data_dict = {'th': [], 'ph': [], 'cplt': [], 'r1': [], 'weight': []}

    def orientation_weight(th, ph):
        """Calculate how much each point should contribute to distribution"""
        def sphere_cutout_area(th1, th2, ph_span):
            """Calculate sphere area in specified ph and th degree range.
            For exact math, see articles about spherical cap and sector."""
            return np.deg2rad(abs(ph_span)) * \
                   abs(np.cos(np.deg2rad(th1)) - np.cos(np.deg2rad(th2)))
        th_max = min(th + angle_res / 2.0, th_limits[1])
        th_min = max(th - angle_res / 2.0, th_limits[0])
        ph_max = min(ph + angle_res / 2.0, ph_limits[1])
        ph_min = max(ph - angle_res / 2.0, ph_limits[0])
        return sphere_cutout_area(th_min, th_max, ph_max-ph_min)

    def _calculate_completeness_mesh():
        """Calculate completeness for each individual pair of theta and phi."""
        _r1_mesh = np.zeros_like(th_mesh, dtype=float)
        lst = open(lst_path, 'w+')
        lst.write('#     th      ph      R1    cplt\n')
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
                q.dac_trim(opening_angle, sph2cart(1.0, np.deg2rad(th),
                                                   np.deg2rad(ph)))
                q.write(hkl_path2, hkl_format='shelx_4')
                count = q.table['equiv'].nunique()
                os.system(f'cd {dir_path2}; shelxl {job_name}')
                r1 = LstFrame().read_r1(lst_path2)
                cplt = count / total_reflections
                data_dict['th'].append(th)
                data_dict['ph'].append(ph)
                data_dict['cplt'].append(cplt)
                data_dict['r1'].append(r1)
                data_dict['weight'].append(orientation_weight(th, ph))
                lst.write(f'{th:8.0f}{ph:8.0f}{r1:8.5}{cplt:8.5f}\n')
                _r1_mesh[j][i] = r1
            lst.write('\n')
        index_max = np.unravel_index(np.argmax(_r1_mesh), _r1_mesh.shape)
        best_th, best_ph = th_range[index_max[1]], ph_range[index_max[0]]
        index_min = np.unravel_index(np.argmin(_r1_mesh), _r1_mesh.shape)
        worst_th, worst_ph = th_range[index_min[1]], ph_range[index_min[0]]
        q1, q2, q3 = weighted_quantile(values=data_dict['r1'],
                                       quantiles=[0.25, 0.50, 0.75],
                                       weights=data_dict['weight'])
        avg_r = np.average(data_dict['r1'], weights=data_dict['weight'])
        max_r = max(data_dict['r1'])
        min_r = min(data_dict['r1'])
        s = f'# descriptive statistics for r1:\n' \
            f'# max ={max_r:8.5f} at th ={best_th :6.1f} ph ={best_ph :6.1f}\n'\
            f'# min ={min_r:8.5f} at th ={worst_th:6.1f} ph ={worst_ph:6.1f}\n'\
            f'# q_1 ={q1   :8.5f}\n' \
            f'# q_2 ={q2   :8.5f}\n' \
            f'# q_3 ={q3   :8.5f}\n' \
            f'# avg ={avg_r:8.5f}\n'
        lst.write(s)
        lst.close()
        np.savetxt(dat_path, _r1_mesh)
        return data_dict

    data_dict = _calculate_completeness_mesh()

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
