"""
This sub-module contains scripts made to work mainly with .hkl files,
containing information from single crystal diffraction experiments.
Due to authors interest in the subjects, majority of them aim to analyse
or model data from high-pressure experiments,
which utilise DAC - diamond anvil cell.
"""

import pandas as pd
import seaborn as sns
from kesshou.dataframes import HklFrame
from kesshou.symmetry import PG, SG, Group
from kesshou.utility import cubespace, fibonacci_sphere, home_directory, \
    make_absolute_path, gnuplot_cplt_map_palette, mpl_cplt_map_palette, \
    deep_palette, pastel_palette, cplt_map_template
from matplotlib import cm, colors, pyplot
from mpl_toolkits.mplot3d import Axes3D
from mpl_toolkits.mplot3d import art3d
from math import erf, erfc
import numpy as np
import numpy.linalg as lin
from scipy.special import erfinv
from scipy.stats import norm
from scipy.optimize import minimize


def completeness_map(a, b, c, al, be, ga,
                     space_group=SG['P1'],
                     axis='',
                     fix_scale=False,
                     opening_angle=35,
                     output_directory=home_directory,
                     output_name='cplt_map',
                     output_quality=3,
                     resolution=1.2,
                     wavelength='MoKa'):
    """
    Calculate and draw a map of predicted completeness for a given crystal
    encased in diamond anvil cell (dac) with a given opening angle
    as a function of crystal orientation in a dac.

    The script accepts and takes into consideration unit cell dimensions,
    laue group and extinctions, which allows to predict the completeness of
    fully merged data for investigated crystal in any space group.
    The results are printed and visualised using various formats,
    including raw text files and ready-to-review .png file.
    A gnuplot input is also provided for convenience.

    Program presents the results using a 2-dimensional spherical heatmap.
    Information about predicted completeness is visualised on a unit sphere,
    where each point on a sphere is associated with a certain orientation
    of the crystal relative to dac:

    - Firstly, the map presents the results in a reciprocal space.
      Therefore red, green and blue lines / vectors represent
      crystallographic directions **X\***, **Y\*** and **Z\***, respectively.

    - Since the *distribution* of reflections in reciprocal space has always
      at least a centre of inversion, it is unnecessary to show the whole
      sphere, as at least half of the data would be redundant.
      For this reason the script shows only a certain part of the unit sphere,
      whereas the information about other directions
      can be retrieved using reciprocal space symmetry for a given laue group.

    - Each point on the sphere is associated with a unit vector,
      which shows a direction perpendicular to a disc of dac-accessible volume.
      In other words, each point on the surface can be associated with a vector
      in reciprocal space from the origin towards said point. This vector is
      perpendicular to the dac-accessible plane (volume),
      which has been originally used to limit the sphere of hkl data.

    As an example, let's assume a orthorhombic cell with *a* = *b* = *c* = 10
    and laue group "mmm". Running the script and generating completeness the map
    yields the lowest values close to **X\***, **Y\*** and **Z\*** vector,
    while the highest values are observed inbeetween those vectors.
    Placing the crystal on its [100] face inside the dac will cause the
    dac-accessible plane to be placed perpendicularly to (100) direction
    in reciprocal space. Since the values close to **X\*** ((100) direction)
    are low, such a placement will allow us to collect data with low coverage.
    On the other hand, placing the crystal on its [111] face will cause the
    dac-acessible plane to be placed perpendicularly to (100) direction
    in reciprocal space. As on the sphere the values tend to rise the further
    we are from of **X\***, **Y\*** and **Z\***, we expect this crystal
    orientation to warrant a high completeness of collected data.

    The completeness is calculated as a ratio between the number of unique
    reflections inside the dac-accessible space and the number of unique
    reflections inside a reference sphere of a radius equal to `resolution`.


    Completeness is visualised using an extended rainbow heatmap,
    which utilises a wide color range to emphasize even small differences.
    Dy default, the color scale is dynamic and adapts to span between
    minimum and maximum values of found completeness,
    but it can be fixed to span between 0 and 100% using `fix_scale` parameter.

    Since the orientation is described using spherical coordinates,
    the exact positions of individual points on the sphere
    are described using *theta* and *phi* angles instead of using
    crystallographic directions. The *theta* and *phi* angles follow the
    physical definition of spherical coordinate system - for a given
    vector **v** perpendicular to the dac plane:

    - *theta* is the azimuth angle,
      i.e. angle between vector **Z\*** and **v**.
      *Theta* equals zero when **v** is parallel to **Z\***,
      equals *pi/2* when **v** is perpendicular to **Z\*** and
      equals *pi* when **v** is anti-parallel to **Z\***.
      Therefore *theta* should take values between zero and *pi*.

    - *phi* is the rotational angle,
      i.e. degree of rotation of **v** around **Z\***.
      *Phi* equals zero when **v** is in **X\*Z\*** plane,
      equals *beta\** when **v** is in **Y\*Z\*** plane,
      and equals *pi* when **v** is back in **X\*Z\*** plane.
      Therefore *phi* should take values between zero and *2 pi*.

    Finally, the script does not treat reflections in special positions (hk0,
    00l) in any special manner and thus does not consider lack of coverage in
    any crystallographic plane or direction as an issue. This information should
    be additionally considered while choosing crystal orientation and data
    collection strategy based on the script output.

    :param a: Unit cell parameter *a* in Angstrom.
    :type a: float
    :param b: Unit cell parameter *b* in Angstrom.
    :type b: float
    :param c: Unit cell parameter *c* in Angstrom.
    :type c: float
    :param al: Unit cell parameter *alpha* in degrees.
    :type al: float
    :param be: Unit cell parameter *alpha* in degrees.
    :type be: float
    :param ga: Unit cell parameter *alpha* in degrees.
    :type ga: float
    :param space_group: Instance of :class:`kesshou.symmetry.Group`
        describing symmetry of the crystal
    :type space_group: kesshou.symmetry.Group
    :param axis: area to calculate completeness of. Accepts 'x', 'y', 'z', 'xy',
        'xz', 'yz' or '' for whole sphere.
    :type axis: string
    :param fix_scale: If true, the colour scheme will not adapt to
        be fixed to the range from 0 to 100%
    :type fix_scale:
    :param opening_angle:
    :type opening_angle:
    :param output_directory:
    :type output_directory:
    :param output_name:
    :type output_name:
    :param output_quality:
    :type output_quality:
    :param resolution: If given, additionally limit data resolution to given
        value. Please provide the resolution as a distance from the origin
        in reciprocal space (twice the resolution in reciprocal angstrom).
    :type resolution: float
    :param wavelength: Wavelength of radiation to be simulated.
    :type wavelength: float or str
    :return: None
    :rtype: None
    """
    dat_path = make_absolute_path(output_directory, output_name + '.dat')
    gnu_path = make_absolute_path(output_directory, output_name + '.gnu')
    lst_path = make_absolute_path(output_directory, output_name + '.lst')
    png_path = make_absolute_path(output_directory, output_name + '.png')
    axis = axis.lower()
    laue_group = space_group.reciprocate()

    def _make_hkl_frame(ax=axis):
        """Make ball or axis of hkl which will be cut in further steps"""
        _f = HklFrame()
        _f.edit_cell(a=a, b=b, c=c, al=al, be=be, ga=ga)
        _f.la = wavelength
        _f.fill(radius=min(_f.r_lim, resolution))
        if ax in {'x'}:
            _f.table = _f.table.loc[_f.table['k'].eq(0) & _f.table['l'].eq(0)]
        elif ax in {'y'}:
            _f.table = _f.table.loc[_f.table['h'].eq(0) & _f.table['l'].eq(0)]
        elif ax in {'z'}:
            _f.table = _f.table.loc[_f.table['h'].eq(0) & _f.table['k'].eq(0)]
        elif ax in {'xy'}:
            _f.table = _f.table.loc[_f.table['l'].eq(0)]
        elif ax in {'xz'}:
            _f.table = _f.table.loc[_f.table['k'].eq(0)]
        elif ax in {'yz'}:
            _f.table = _f.table.loc[_f.table['h'].eq(0)]
        if ax in {'x', 'y', 'z', 'xy', 'xz', 'yz'}:
            _f.transform([o.tf for o in laue_group.operations])
        _f.extinct(space_group)
        return _f

    p = _make_hkl_frame()
    p.find_equivalents(point_group=laue_group)
    total_reflections = p.table['equiv'].nunique()
    assert total_reflections > 0, "No non-extinct reflections in this region"

    def _determine_theta_and_phi_limits():
        """Define the spherical coordinate system based on given point group.
        v1, v2, v3 are normal vectors pointing in zenith direction z*,
        orthogonal direction (x) and direction orthogonal to them both."""
        _v1 = p.z_w
        _v2 = p.x_v
        _v3 = np.cross(_v1, _v2)

        if space_group.system is Group.CrystalSystem.triclinic:
            _th_limits = [0, 180]
            _ph_limits = [0, 180]
        elif space_group.system in {Group.CrystalSystem.monoclinic,
                                    Group.CrystalSystem.orthorhombic,
                                    Group.CrystalSystem.tetragonal,
                                    Group.CrystalSystem.cubic}:
            _th_limits = [0, 90]
            _ph_limits = [0, 90]
        elif space_group.system in {Group.CrystalSystem.trigonal,
                                    Group.CrystalSystem.hexagonal}:
            _th_limits = [0, 90]
            _ph_limits = [0, 120]
        else:
            raise ValueError('Provided group is not one of known Laue groups')
        return _v1, _v2, _v3, _th_limits, _ph_limits

    v1, v2, v3, th_limits, ph_limits = _determine_theta_and_phi_limits()

    def _make_theta_and_phi_mesh():
        """Define a list of theta and phi values to be investigated.
        Theta and phi are defined as in physics ISO convention.
        Theta is "elevation angle" away from primary (zenith) direction,
        towards secondary direction, and takes values from 0 to 180 degrees.
        Phi is "azimuth angle" which rotates perpendicularly to primary,
        from secondary to tertiary direction. Takes values from 0 to 360."""
        assert output_quality in {1, 2, 3, 4, 5}, 'Quality not 1, 2, 3, 4 or 5'
        angle_res = {1: 15, 2: 10, 3: 5, 4: 2, 5: 1}[output_quality]
        _th_range = np.arange(th_limits[0], th_limits[1] + 0.001, angle_res)
        _ph_range = np.arange(ph_limits[0], ph_limits[1] + 0.001, angle_res)
        _th_mesh, _ph_mesh = np.meshgrid(_th_range, _ph_range)
        return _th_range, _ph_range, _th_mesh, _ph_mesh

    th_range, ph_range, th_mesh, ph_mesh = _make_theta_and_phi_mesh()
    data_dict = {'th': [], 'ph': [], 'cplt': [], 'reflns': []}

    def _translate_angles_to_vector(theta, phi):
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
        _cplt_mesh = np.zeros_like(th_mesh)
        lst = open(lst_path, 'w+', buffering=1)
        lst.write('#     th      ph    cplt  reflns\n')
        for i, th in enumerate(th_range):
            for j, ph in enumerate(ph_range):
                v = _translate_angles_to_vector(theta=th, phi=ph)
                q = p.duplicate()
                q.dac(opening_angle=opening_angle, vector=v)
                hkl_len = q.table['equiv'].nunique()
                data_dict['th'].append(th)
                data_dict['ph'].append(ph)
                data_dict['cplt'].append(hkl_len / total_reflections)
                data_dict['reflns'].append(hkl_len)
                lst.write('{:8.0f}'.format(th))
                lst.write('{:8.0f}'.format(ph))
                lst.write('{:8.5f}'.format(hkl_len / total_reflections))
                lst.write('{:8d}'.format(hkl_len))
                lst.write('\n')
                _cplt_mesh[j][i] = hkl_len / total_reflections
            lst.write('\n')
        index_max = np.unravel_index(np.argmax(_cplt_mesh), _cplt_mesh.shape)
        best_th, best_ph = th_range[index_max[1]], ph_range[index_max[0]]
        index_min = np.unravel_index(np.argmin(_cplt_mesh), _cplt_mesh.shape)
        worst_th, worst_ph = th_range[index_min[1]], ph_range[index_min[0]]
        lst.write('# best cplt = {_max} for th = {_th}, ph = {_ph}\n'.format(
            _max=max(data_dict['cplt']), _th=best_th, _ph=best_ph))
        lst.write('# worst cplt = {_min} for th = {_th}, ph = {_ph}\n'.format(
            _min=min(data_dict['cplt']), _th=worst_th, _ph=worst_ph))
        lst.write('# mean cplt = {}\n'.format(np.mean(data_dict['cplt'])))
        lst.close()
        np.savetxt(dat_path, _cplt_mesh)
        return data_dict, _cplt_mesh

    data_dict, cplt_mesh = _calculate_completeness_mesh()

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
        ax.plot_wireframe(x, y, z, colors='k', linewidth=0.25)

        # color map
        my_heatmap_colors = mpl_cplt_map_palette[axis]
        my_colormap = colors.LinearSegmentedColormap.from_list(
            'heatmapEX', my_heatmap_colors, N=256)
        m = cm.ScalarMappable(cmap=my_colormap)
        m.set_array(cplt_mesh)
        if fix_scale is True:
            m.set_clim(0, 1)
            norm = colors.Normalize(vmin=0, vmax=1)
        else:
            norm = colors.Normalize(vmin=min(data_dict['cplt']),
                                    vmax=max(data_dict['cplt']))
        pyplot.colorbar(m, fraction=0.046, pad=0.04)

        # direction lines
        _len = 1.25
        _x, _y, _z = p.x_w, p.y_w, p.z_w
        ax.add_line(art3d.Line3D((_x[0], _len * _x[0]), (_x[1], _len * _x[1]),
                                 (_x[2], _len * _x[2]), color='r', linewidth=5))
        ax.add_line(art3d.Line3D((_y[0], _len * _y[0]), (_y[1], _len * _y[1]),
                                 (_y[2], _len * _y[2]), color='g', linewidth=5))
        ax.add_line(art3d.Line3D((_z[0], _len * _z[0]), (_z[1], _len * _z[1]),
                                 (_z[2], _len * _z[2]), color='b', linewidth=5))

        # color mesh for heatmap
        color_mesh = cplt_mesh[:-1, :-1]
        for i in range(cplt_mesh.shape[0] - 1):
            for j in range(cplt_mesh.shape[1] - 1):
                color_mesh[i, j] = (cplt_mesh[i + 1, j] + cplt_mesh[i, j + 1] +
                                    cplt_mesh[i, j] + cplt_mesh[
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
        gnu = open(gnu_path, 'w+', buffering=1)
        gnu.write(cplt_map_template.format(
            axis_x1=p.x_w[0],
            axis_x2=p.x_w[1],
            axis_x3=p.x_w[2],
            axis_y1=p.y_w[0],
            axis_y2=p.y_w[1],
            axis_y3=p.y_w[2],
            axis_z1=p.z_w[0],
            axis_z2=p.z_w[1],
            axis_z3=p.z_w[2],
            cplt_min=0 if fix_scale else min(data_dict['cplt']) * 100,
            cplt_max=100 if fix_scale else max(data_dict['cplt']) * 100,
            job_name=output_name,
            min_ph=min(ph_limits),
            max_ph=max(ph_limits),
            min_th=min(th_limits),
            max_th=max(th_limits),
            palette=gnuplot_cplt_map_palette[axis.lower()]))

    _prepare_gnuplot_input()
    try:
        from os import system, getcwd
        _path = make_absolute_path(output_directory)
        system('cd ' + _path + '; gnuplot ' + gnu_path)
    except OSError:
        pass


def completeness_statistics(a, b, c, al, be, ga,
                            point_group=PG['-1'],
                            input_path='shelx.hkl',
                            input_format=4,
                            input_wavelength='CuKa'):
    """
    For a given experimental .hkl file
    calculate basic completeness statistics in equal-volume resolution shells.
    This script directly calls method :meth:`kesshou.dataframes.make_stats`.

    :param a: Unit cell parameter *a* in Angstrom.
    :type a: float
    :param b: Unit cell parameter *b* in Angstrom.
    :type b: float
    :param c: Unit cell parameter *c* in Angstrom.
    :type c: float
    :param al: Unit cell parameter *alpha* in degrees.
    :type al: float
    :param be: Unit cell parameter *alpha* in degrees.
    :type be: float
    :param ga: Unit cell parameter *alpha* in degrees.
    :type ga: float
    :param point_group: Point group of the crystal,
        defined as an instance of :class:`kesshou.symmetry.Group`
    :type point_group: kesshou.symmetry.Group
    :param input_path: Path to the input .hkl file.
    :type input_path: str
    :param input_format: Format of the .hkl file. For reference see
        :meth:`kesshou.dataframes.HklFrame.interpret_hkl_format`.
    :type input_format: int or str or dict
    :param input_wavelength: Wavelength of radiation utilised in experiment.
    :type input_wavelength: float or str
    :return: None
    """
    p = HklFrame()
    p.edit_cell(a=a, b=b, c=c, al=al, be=be, ga=ga)
    p.la = input_wavelength
    p.read(input_path, input_format)
    p.stats(space_group=point_group)


def dac_point_group_statistics(a, b, c, al=90, be=90, ga=90,
                               space_group=SG['P1'],
                               output_path='output.txt',
                               opening_angle=35,
                               precision=1000,
                               resolution=None,
                               wavelength='MoKa'):
    """
    Calculate max, min, avg completeness for multiple crystal orientations
    for a diamond anvil cell.future it is planned to merge it with the
    :func:`dac_completeness_map` script.

    :param a: Unit cell parameter *a* in Angstrom.
    :type a: float
    :param b: Unit cell parameter *b* in Angstrom.
    :type b: float
    :param c: Unit cell parameter *c* in Angstrom.
    :type c: float
    :param al: Unit cell parameter *alpha* in degrees.
    :type al: float
    :param be: Unit cell parameter *alpha* in degrees.
    :type be: float
    :param ga: Unit cell parameter *alpha* in degrees.
    :type ga: float
    :param space_group: Space group of the crystal,
        defined as an instance of :class:`kesshou.symmetry.Group`
    :type space_group: kesshou.symmetry.Group
    :param output_path: Path of created file containing calculated data.
    :type output_path: str
    :param opening_angle:
    :type opening_angle:
    :param precision: Number of vectors for which
        completeness in the dac will be calculated.
    :type precision: int
    :param resolution: If given, additionally limit data resolution to given
        value. Please provide the resolution as a distance from the origin
        in reciprocal space (twice the resolution in reciprocal angstrom).
    :type resolution: float
    :param wavelength: Wavelength of radiation to be simulated.
    :type wavelength: float or str
    :return: None
    """
    def _make_reference_ball():
        hkl_frame = HklFrame()
        hkl_frame.edit_cell(a=a, b=b, c=c, al=al, be=be, ga=ga)
        hkl_frame.la = wavelength
        hkl_frame.fill(radius=hkl_frame.r_lim)
        hkl_frame.merge()
        if not(resolution is None):
            hkl_frame.trim(resolution)
        return hkl_frame
    p = _make_reference_ball()
    p.find_equivalents(point_group=space_group.reciprocate())
    p.extinct(space_group=space_group)
    total_reflections = p.table['equiv'].nunique()
    max_resolution = max(p.table['r'])
    out = open(output_path, 'w', buffering=1)
    out.write('total_reflections: ' + str(total_reflections) + '\n')
    out.write('maximum_r_in_reciprocal_coordinates: ' + str(max_resolution))
    assert precision >= 100, 'Precision should be at least 100'
    vectors = [(1, 0, 0), (1, 0, 0), (0, 0, 1), (0.57735, 0.57735, 0.57735)]\
              + fibonacci_sphere(samples=precision-4, seed=1337)
    reflections = list()
    for vector in vectors:
        q = p.duplicate()
        q.dac(opening_angle=opening_angle, vector=vector)
        reflections.append(q.table['equiv'].nunique())
        out.write('\n' + str(vector) + ': ' + str(q.table['equiv'].nunique()))
        # q.to_res(output_path + '.res', colored='h')  # temp #
    out.write('\nmax_reflections: ' + str(max(reflections)))
    out.write('\nmin_reflections: ' + str(min(reflections)))
    out.write('\navg_reflections: ' + str(sum(reflections) / len(reflections)))
    out.close()
    return [r / total_reflections for r in reflections]


def dac_completeness_vs_opening_angle(output_path='output.txt',
                                      precision=91,
                                      resolution=None,
                                      wavelength='MoKa',
                                      theta=None):
    """
    Calculate completeness in P-1 as a function of DAC opening angle,
    assuming certain resolution and wavelength used.
    :param output_path: Path of created file containing calculated data.
    :type output_path: str
    :param precision: Number of probed opening angles between 0 and 90 degrees.
    :type precision: int
    :param resolution: If given, additionally limit data resolution to given
        value. Please provide the resolution as a distance from the origin
        in reciprocal space (twice the resolution in reciprocal angstrom).
    :type resolution: float
    :param wavelength: Wavelength of radiation to be simulated.
    :type wavelength: float or str
    :param theta: Resolution in theta angle in degrees instead radius in A-1
    :type theta: float
    :return: None
    """

    def _make_reference_ball():
        hkl_frame = HklFrame()
        hkl_frame.la = wavelength
        res = hkl_frame.r_lim * np.sin(np.deg2rad(theta)) \
            if theta is not None else resolution
        side = 10 * precision**(1/3) / res  # adapt to prec&la
        hkl_frame.edit_cell(a=side, b=side, c=side, al=90, be=90, ga=90)
        if res is not None:
            hkl_frame.fill(radius=res)
        else:
            hkl_frame.fill(radius=hkl_frame.r_lim)
        return hkl_frame
    p = _make_reference_ball()
    total_reflns = len(p.table)
    angles = np.linspace(start=90, stop=0, num=precision)
    out = open(output_path, 'w', buffering=1)
    out.write('#oa      cplt\n')
    for a in angles:
        p.dac(opening_angle=a, vector=np.array((1, np.pi/4, np.e/5)))  # rand. v
        out.write(' {a:7.4f} {c:7.5f}\n'.format(a=a, c=len(p.table)/total_reflns))
    out.close()


def completeness_violin_plot(input_path='',
                             output_path='output.txt',
                             opening_angle=35,
                             precision=1000,
                             resolution=None,
                             wavelength='MoKa'):
    if input_path is '':
        dpgs = dac_point_group_statistics
        kwargs = {'a': 20, 'b': 20, 'c': 20, 'al': 90, 'be': 90,
                  'opening_angle': opening_angle, 'precision': precision,
                  'output_path': output_path, 'resolution': resolution,
                  'wavelength': wavelength}
        cplt_dict = dict()
        cplt_dict['-1'] = dpgs(space_group=SG['P-1'], **kwargs)
        cplt_dict['2/m'] = dpgs(space_group=SG['P2/m'], **kwargs)
        cplt_dict['mmm'] = dpgs(space_group=SG['Pmmm'], **kwargs)
        cplt_dict['4/m'] = dpgs(space_group=SG['P4/m'], **kwargs)
        cplt_dict['4/mmm'] = dpgs(space_group=SG['P4/mmm'], **kwargs)
        cplt_dict['-3'] = dpgs(ga=120, space_group=SG['P-3'], **kwargs)
        cplt_dict['-3m'] = dpgs(ga=120, space_group=SG['P-3m1'], **kwargs)
        cplt_dict['6/m'] = dpgs(ga=120, space_group=SG['P6/m'], **kwargs)
        cplt_dict['6/mmm'] = dpgs(ga=120, space_group=SG['P6/mmm'], **kwargs)
        cplt_dict['m-3'] = dpgs(space_group=SG['Pm-3'], **kwargs)
        cplt_dict['m-3m'] = dpgs(space_group=SG['Pm-3m'], **kwargs)
        cplt_frame = pd.DataFrame.from_dict(cplt_dict)
        cplt_frame = cplt_frame.mul(100)
        cplt_frame.to_csv(output_path)
    else:
        cplt_frame = pd.read_csv(input_path, index_col=0)

    # GENERAL
    # cplt_frame.drop(labels='-1', axis=1, inplace=True)
    melt_frame = cplt_frame.melt(var_name='LG', value_name='Completeness')
    sns.set_style("whitegrid")
    # WIDE PLOT
    #pyplot.figure(figsize=(15, 6))
    #sns.set(font_scale=1.75)
    #ax = sns.violinplot(data=cplt_frame, scale='width', cut=0,
    #                    palette=deep_palette, inner='quartile')
    #ax.axhline(27.9458181944, linewidth=5, alpha=0.5)
    # pyplot.ylim(25, 101)
    # SMALL PLOT
    pyplot.figure(figsize=(8, 6))
    sns.set(font_scale=1.5)
    ax = sns.violinplot(x=melt_frame['Completeness'],
                        y=melt_frame['LG'], scale='width', cut=0,
                        palette=deep_palette, inner='quartile')
    for l in ax.lines[1::3]:
        l.set_linestyle('--')
        l.set_linewidth(1.25)
    ax.set_xlabel('')
    ax.set_ylabel('')
    for y, key in enumerate(cplt_frame.keys()):
        ax.axvline(x=cplt_frame[key].median(), color=deep_palette[key],
                   ymin=0, ymax=21/22-y/11, linewidth=5, alpha=0, zorder=0.99)
    ax.set_yticklabels([r"$\overline{1}$", r"$2/m$", r"$mmm$", r"$4/m$",
                        r"$4/mmm$", r"$\overline{3}$", r"$\overline{3}m$",
                        r"$6/m$", r"$6/mmm$",
                        r"$m\overline{3}$", r"$m\overline{3}m$"])
    ax.set_xticklabels(['10', r'$\mathcal{P}$ = 20       ', '30', '40',
                        '50', '60', '70', '80', '90', '100%'])
    pyplot.xlim(19, 101) #19
    # GENERAL
    sns.despine(left=True, top=False)
    pyplot.show()
    # use figsize=(5.3, 6),xlim(54, 101) and appropriate labels for thin version


def dac_statistics(a, b, c, al, be, ga,
                   space_group=SG['P1'],
                   opening_angle=35.0,
                   orientation=((1, 0, 0), (0, 1, 0), (0, 0, 1)),
                   input_path='shelx.hkl',
                   input_format=4,
                   input_wavelength='CuKa',
                   resolution=None):
    """
    For a given experimental .hkl file calculate number of experimentally found
    and theoretically possible reflections up to a given resolution.

    :param a: Unit cell parameter *a* in Angstrom.
    :type a: float
    :param b: Unit cell parameter *b* in Angstrom.
    :type b: float
    :param c: Unit cell parameter *c* in Angstrom.
    :type c: float
    :param al: Unit cell parameter *alpha* in degrees.
    :type al: float
    :param be: Unit cell parameter *alpha* in degrees.
    :type be: float
    :param ga: Unit cell parameter *alpha* in degrees.
    :type ga: float
    :param point_group: Point group of the crystal,
        defined as an instance of :class:`kesshou.symmetry.Group`
    :type point_group: kesshou.symmetry.Group
    :param opening_angle: Value of single opening angle as defined in
        :meth:`kesshou.dataframes.HklFrame.dac`.
    :type opening_angle: float
    :param orientation: Crystal orientation as defined in
        :class:`kesshou.dataframes.BaseFrame`
    :type orientation: tuple or numpy.array
    :param input_path: Path to the input .hkl file.
    :type input_path: str
    :param input_format: Format of the .hkl file. For reference see
        :meth:`kesshou.dataframes.HklFrame.interpret_hkl_format`.
    :type input_format: int or str or dict
    :param input_wavelength: Wavelength of radiation utilised in experiment.
    :type input_wavelength: float or str
    :return: None
    """

    point_group = space_group.reciprocate()#.lauefy()

    p = HklFrame()
    p.edit_cell(a=a, b=b, c=c, al=al, be=be, ga=ga)
    p.la = input_wavelength
    p.read(input_path, input_format)
    p.orientation = np.array(orientation)
    #p.extinct(space_group=space_group)

    resolution = p.r_lim if resolution is None else resolution
    p.trim(limit=resolution)
    p.merge(point_group=point_group)

    q = p.duplicate()
    q.fill(radius=resolution)
    q.dac(opening_angle=opening_angle)

    # uncomment this part for 2nd crystal in different orientation
    # q2 = p.duplicate()
    # q2.fill(radius=resolution)
    # q2.orientation = np.array(((0.08251, 0.01162, 0.03675), (0.01682, 0.03807, -0.03236), (-0.06016, 0.01982, 0.04088)))
    # q2.dac(opening_angle=opening_angle)
    # q = q + q2

    q.extinct(space_group=space_group)
    q.merge(point_group=point_group)

    b = p.duplicate()
    b.fill(radius=resolution)
    b.extinct(space_group=space_group)
    b.merge(point_group=point_group)

    r_max = resolution #max(b.table['r'])
    print('radius    experimnt theoryDAC theorBall DAC-Cplt  Ball-Cplt ')
    print('range     unique    unique    unique    exp/DAC   exp/Ball  ')
    for rad in reversed(cubespace(0, r_max, 10, include_start=False)):
        p.trim(rad)
        q.trim(rad)
        b.trim(rad)
        p_eq = p.table['equiv'].nunique()
        q_eq = q.table['equiv'].nunique()
        b_eq = b.table['equiv'].nunique()
        print(' {:9f}'.format(rad) + ' {:9d}'.format(p_eq) +
              ' {:9d}'.format(q_eq) + ' {:9d}'.format(b_eq) +
              ' {:9f}'.format(p_eq / q_eq) + ' {:9f}'.format(p_eq / b_eq))


def completeness_statistics_around_axis(a, b, c, al, be, ga,
                                        space_group=SG['P1'],
                                        opening_angle=35.0,
                                        wavelength='MoKa',
                                        vector=(1, 0, 0),
                                        topple=5):
    """
    For a given simulated .hkl file calculate average completeness of data
    obtained by toppling the crystal by "topple" degrees from the axis.

    :param a: Unit cell parameter *a* in Angstrom.
    :type a: float
    :param b: Unit cell parameter *b* in Angstrom.
    :type b: float
    :param c: Unit cell parameter *c* in Angstrom.
    :type c: float
    :param al: Unit cell parameter *alpha* in degrees.
    :type al: float
    :param be: Unit cell parameter *alpha* in degrees.
    :type be: float
    :param ga: Unit cell parameter *alpha* in degrees.
    :type ga: float
    :param space_group: Space group of the crystal,
        defined as an instance of :class:`kesshou.symmetry.Group`
    :type space_group: kesshou.symmetry.Group
    :param wavelength: Wavelength of radiation utilised in experiment.
    :type wavelength: float or str
    :param opening_angle: Value of single opening angle as defined in
        :meth:`kesshou.dataframes.HklFrame.dac`.
    :type opening_angle: float
    :param vector: Direction around which completeness will be calculated.
    :type vector: tuple
    :param topple: Angle by which vector will be toppled in varous directions.
    :type topple: float
    :return: None
    """

    point_group = space_group.reciprocate()#.lauefy()

    p = HklFrame()
    p.edit_cell(a=a, b=b, c=c, al=al, be=be, ga=ga)
    p.la = 'MoKa'
    p.fill(radius=p.r_lim)
    p.place()
    p.extinct(space_group=space_group)
    p.find_equivalents(point_group=point_group)

    # generate perpendicular vector for rotation
    v = np.array(vector) / lin.norm(np.array(vector))
    x = np.array((1, 0, 0))

    def are_parallel(v1, v2):
        return np.allclose(v1 / lin.norm(v1), v2 / lin.norm(v2))
    temp = x if not(are_parallel(v, x)) else np.array((0, 1, 0))
    perp = np.cross(v, temp)

    def rotate(v, k, angle):
        c = np.cos(np.deg2rad(angle))
        s = np.sin(np.deg2rad(angle))
        return v * c + np.cross(k, v) * s + k * np.dot(k, v) * (1 - c)

    # generate 360 toppled vectors
    toppled = rotate(v, perp, topple)
    toppleds = [rotate(toppled, v, i) for i in range(360)]

    # calculate the completeness for toppleds
    #rads = list(reversed(cubespace(0, 1.2, 10, include_start=False)))
    rads = [2.00, 1/0.7, 1.20, 1.00, 1/1.5]
    cplt = [0] * len(rads)
    p.trim(max(rads))
    for t in toppleds:
        q = p.duplicate()
        q.dac(opening_angle=opening_angle, vector=t)
        for cplt_bin, rad in enumerate(rads, start=1):
            q.trim(rad)
            cplt[-cplt_bin] += q.table['equiv'].nunique()

    # divide lists of completeness by total and return from sums to individuals

    full = [0] * len(rads)
    for cplt_bin, rad in enumerate(rads, start=1):
        p.trim(rad)
        full[-cplt_bin] = p.table['equiv'].nunique()
    print(cplt)
    print(full)
    cplt = [c2 - c1 for c1, c2 in zip([0] + cplt[:-1], cplt)]
    full = [f2 - f1 for f1, f2 in zip([0] + full[:-1], full)]
    cplt = [c / (360 * f) for c, f in zip(cplt, full)]

    print('Resolution limits (in distance to reflection, twice reciprocal:')
    print(np.array(list(reversed(rads))))
    print('Average shell completeness with {}degree topple'.format(topple))
    print(np.array(cplt))


def simulate_dac(a, b, c, al, be, ga,
                 opening_angle=35,
                 orientation=((1, 0, 0), (0, 1, 0), (0, 0, 1)),
                 vector=None,
                 resolution=None,
                 input_path='shelx.hkl',
                 input_format=4,
                 input_wavelength='CuKa',
                 output_path='output.hkl',
                 output_format=4):
    """
    For a given experimental .hkl file simulate a lack of completeness
    caused by a presence of high-pressure diamond anvil cell.

    :param a: Unit cell parameter *a* in Angstrom.
    :type a: float
    :param b: Unit cell parameter *b* in Angstrom.
    :type b: float
    :param c: Unit cell parameter *c* in Angstrom.
    :type c: float
    :param al: Unit cell parameter *alpha* in degrees.
    :type al: float
    :param be: Unit cell parameter *alpha* in degrees.
    :type be: float
    :param ga: Unit cell parameter *alpha* in degrees.
    :type ga: float
    :param opening_angle: Value of single opening angle as defined in
        :meth:`kesshou.dataframes.HklFrame.dac`.
    :type opening_angle: float
    :param orientation: Crystal orientation as defined in
        :class:`kesshou.dataframes.BaseFrame`
    :type orientation: tuple or numpy.array
    :param vector: If given, overwrite orientation to provide information
        about crystal placement in dac,
        as defined in :meth:`kesshou.dataframes.HklFrame.dac`.
    :type vector: tuple
    :param resolution: If given, additionally limit data resolution to given
        value. Please provide the resolution as a distance from the origin
        in reciprocal space (twice the resolution in reciprocal angstrom).
    :type resolution: float
    :param input_path: Path to the input .hkl file.
    :type input_path: str
    :param input_format: Format of the input .hkl file. For reference see
        :meth:`kesshou.dataframes.HklFrame.interpret_hkl_format`.
    :type input_format: int or str or dict
    :param input_wavelength: Wavelength of radiation utilised in experiment.
    :type input_wavelength: float or str
    :param output_path: Path to the output .hkl file.
    :type output_path: str
    :param output_format: Format of the input .hkl file. For reference see
        :meth:`kesshou.dataframes.HklFrame.interpret_hkl_format`.
    :type output_format: int or str or dict
    :return: None
    """
    p = HklFrame()
    p.edit_cell(a=a, b=b, c=c, al=al, be=be, ga=ga)
    p.la = input_wavelength
    p.read(input_path, input_format)
    p.orientation = np.array(orientation)
    if not(resolution is None):
        p.trim(resolution)
    p.dac(opening_angle=opening_angle, vector=vector)
    p.write(hkl_path=output_path, hkl_format=output_format)


def baycon_plot(x_key='ze', y_key='si',
                a=10.0, b=10.0, c=10.0, al=90.0, be=90.0, ga=90.0,
                input_path='shelx.fcf',
                input_format='shelx_fcf',
                input_wavelength='MoKa',
                output_path='baycon.png'):
    """
    For a given .fcf file prepare a bayesian conditional probability plot
    between x_key and y_key.

    :param x_key: Parameter of HklFrame which will be placed on x axis
    :type x_key: str
    :param y_key: Parameter of HklFrame which will be placed on x axis
    :type y_key: str
    :param a: Unit cell parameter *a* in Angstrom.
    :type a: float
    :param b: Unit cell parameter *b* in Angstrom.
    :type b: float
    :param c: Unit cell parameter *c* in Angstrom.
    :type c: float
    :param al: Unit cell parameter *alpha* in degrees.
    :type al: float
    :param be: Unit cell parameter *alpha* in degrees.
    :type be: float
    :param ga: Unit cell parameter *alpha* in degrees.
    :type ga: float
    :param input_path: Path to the input .fcf file.
    :type input_path: str
    :param input_format: Format of the input .fcf file. For reference see
        :meth:`kesshou.dataframes.HklFrame.interpret_hkl_format`.
    :type input_format: int or str or dict
    :param input_wavelength: Wavelength of radiation utilised in experiment.
    :type input_wavelength: float or str
    :param output_path: Path to the output .png file.
    :type output_path: str
    """
    no_of_bins = 10
    p = HklFrame()
    p.edit_cell(a=a, b=b, c=c, al=al, be=be, ga=ga)
    p.la = input_wavelength
    p.read(input_path, input_format)
    p.place()
    p.calculate_fcf_statistics()
    x = p.table.loc[:, x_key].rank(pct=True).to_numpy()
    y = p.table.loc[:, y_key].rank(pct=True).to_numpy()
    bins = np.zeros(shape=(no_of_bins, no_of_bins))
    lims = [-1.e-8] + [(i + 1) / no_of_bins for i in range(no_of_bins)]
    for i in range(no_of_bins):
        for j in range(no_of_bins):
            bins[i, j] = ((lims[i] < x) & (x <= lims[i+1]) &
                          (lims[j] < y) & (y <= lims[j+1])).sum()
    n_avg = len(x) / no_of_bins ** 2
    chi2 = np.sum((bins - n_avg) ** 2 / n_avg)
    fig = pyplot.figure()
    ax = fig.add_subplot(111, aspect='equal')
    pyplot.xlim(0, 1)
    pyplot.ylim(0, 1)
    h = ax.hist2d(x, y, bins=no_of_bins, alpha=0.25, cmap=cm.get_cmap('PiYG'))
    cb = pyplot.colorbar(h[3], ax=ax)
    cb.set_label('Number of observations')
    ax.scatter(x=x, y=y, s=5.0, c='#000080', marker='.', alpha=0.75)
    pyplot.title('Bayesian CoNditional probability, chi2 = {:.2f}'.format(chi2))
    pyplot.xlabel('"' +  x_key + '" rank')
    pyplot.ylabel('"' +  y_key + '" rank')
    pyplot.tight_layout()
    pyplot.savefig(fname=output_path, dpi=300)


def observed_vs_calculated_plot(input_path='shelx.fcf',
                                input_format='shelx_fcf',
                                output_path='Io_vs_Ic.png'):
    p = HklFrame()
    p.read(input_path, input_format)
    icalc = p.table.loc[:, 'Ic'].to_numpy()
    iobs = p.table.loc[:, 'I'].to_numpy()
    i_min = min(np.min(icalc[icalc > 0]), np.min(iobs[iobs > 0]))
    i_max = max(np.max(icalc[icalc > 0]), np.max(iobs[iobs > 0]))
    fig = pyplot.figure()
    ax = fig.add_subplot(111) # , aspect='equal'
    ax.set_xscale('log')
    ax.set_yscale('log')
    ax.set_xlim([i_min, i_max])
    ax.set_ylim([i_min, i_max])
    ax.plot(np.linspace(0, i_max), np.linspace(0, i_max), '-k', lw=1, zorder=0)
    ax.scatter(x=icalc, y=iobs, s=5.0, c='r', marker='.', alpha=0.75, zorder=10)
    pyplot.title('Calculated vs observed intensities plot')
    pyplot.xlabel('I_cal')
    pyplot.ylabel('I_obs')
    pyplot.tight_layout()
    pyplot.savefig(fname=output_path, dpi=300)


def normal_probability_plot(input_path='shelx.fcf',
                            input_format='shelx_fcf',
                            output_path='Io_vs_Ic.png'):

    # scale factors
    a = 0.1000
    b = 0.0

    p = HklFrame()
    p.read(input_path, input_format)
    i_obs = p.table.loc[:, 'I'].to_numpy()
    i_calc = p.table.loc[:, 'Ic'].to_numpy()
    si = p.table.loc[:, 'si'].to_numpy()
    p = 1/3 * i_obs + 2/3 * i_calc
    si = np.sqrt(si ** 2 + (a * p) ** 2 + b * p)

    # expected delta m
    def delta_m(f1, f2, k, si1, si2):
        return np.sort((f1 - k * f2) / np.sqrt(si1 ** 2 + k **2 * si2 ** 2))

    def sum_of_delta_m_squared(k):
        return np.sum(delta_m(i_obs, i_calc, k, si, np.zeros_like(si)) ** 2)

    def scale_factor():
        return minimize(sum_of_delta_m_squared, x0=np.array([1.0])).x[0]

    experiment_delta_m = delta_m(f1=i_obs, f2=i_calc, k=scale_factor(),
                               si1=si, si2=np.zeros_like(si))
    experiment_delta_m = experiment_delta_m / np.std(experiment_delta_m)

    # simulated delta m
    uniform = (np.arange(len(experiment_delta_m)) + 0.5) / len(experiment_delta_m)
    simulated_delta_m = [erfinv(-1 + 2 * q) for q in uniform]

    # drawing the plot
    fig = pyplot.figure()
    ax = fig.add_subplot(111, aspect='equal')
    ax.set_xlim([-5, 5])
    ax.set_ylim([-5, 5])
    pyplot.hist(experiment_delta_m, bins=100, density=True)
    ax.scatter(experiment_delta_m, simulated_delta_m, s=5.0, c='r', marker='.', alpha=0.75, zorder=10)
    ax.plot(np.linspace(-3, 3), np.linspace(-3, 3), '-k', lw=1, zorder=0)
    pyplot.plot(6 * uniform - 3, norm.pdf(6 * uniform - 3))
    pyplot.title('npp')
    pyplot.xlabel('delta_m experiment')
    pyplot.ylabel('delta_m simulated')
    pyplot.tight_layout()
    pyplot.savefig(fname=output_path, dpi=300)


def fcf_descriptors(input_path='shelx.fcf', input_format='shelx_fcf'):
    # scale factors
    a = 0.1000
    b = 0.0

    p = HklFrame()
    p.read(input_path, input_format)
    i_obs = p.table.loc[:, 'I'].to_numpy()
    i_calc = p.table.loc[:, 'Ic'].to_numpy()
    si = p.table.loc[:, 'si'].to_numpy()
    p = 1/3 * i_obs + 2/3 * i_calc
    si_weighted = np.sqrt(si ** 2 + (a * p) ** 2 + b * p)
    ze = (i_obs - i_calc) / si_weighted
    f_calc = np.sqrt(np.abs(i_calc)) * np.sign(i_calc)
    f_obs = np.sqrt(np.abs(i_obs)) * np.sign(i_obs)
    one_over_sf = (2 * abs(i_obs) ** 0.5) / si

    r1 = np.sum(np.abs(f_obs - f_calc)) / np.sum(np.abs(f_obs))
    wr2 = np.sqrt(
        np.sum(np.abs(si_weighted * np.abs(i_obs - i_calc) ** 2)) /
        np.sum(np.abs(si_weighted * i_obs ** 2)))
    awr2 = np.sqrt(
        (np.mean((i_obs - i_calc) ** 2) / np.mean(si_weighted ** 2)) /
        np.mean((i_obs / si_weighted) ** 2))
    gof_if_alpha_equal_one = np.sqrt(np.mean(ze ** 2))
    agof_if_alpha_equal_one = np.sqrt(
        np.mean((i_obs - i_calc) ** 2) /
        np.mean(si_weighted ** 2))

    print('R1    = {:f}'.format(r1))
    print('wR2   = {:f}'.format(wr2))
    print('awR2  = {:f}'.format(awr2))
    print('GoF*  = {:f}'.format(gof_if_alpha_equal_one))
    print('aGoF* = {:f}'.format(agof_if_alpha_equal_one))


if __name__ == '__main__':
    # TODO make sure unique is defined correctly (no cutoff errors in currently used float)
    # from os import system
    # sg = {'P-1': 'P-1',
    #       'P2/m': 'P2om',
    #       'Pmmm': 'Pmmm',
    #       'P4/m': 'P4om',
    #       'P4/mmm': 'P4ommm',
    #       'Pm-3': 'Pm-3',
    #       'Pm-3m': 'Pm-3m',
    #       'P-3': 'P-3',
    #       'P-3m1': 'P-3m1',
    #       'P6/m': 'P6om',
    #       'P6/mmm': 'P6ommm'}
    # kwargs = {'a': 20, 'b': 20, 'c': 20, 'al': 90, 'be': 90,
    #           'fix_scale': True, 'opening_angle': 45,
    #           'output_quality': 5, 'wavelength': 0.42, 'resolution': 2.0}
    # for k, v in sg.items():
    #     name = 'CpltMap_la42oa45res50_{}'.format(v)
    #     ga = 120 if v in {'P-3', 'P-3m1', 'P6om', 'P6ommm'} else 90
    #     completeness_map(space_group=SG[k], ga=ga,
    #                      output_name=name, **kwargs)
    pass


