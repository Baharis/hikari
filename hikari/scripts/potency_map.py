import numpy as np
from matplotlib import pyplot, colors, cm
from mpl_toolkits.mplot3d import art3d
from numpy import linalg as lin

from hikari.dataframes import HklFrame
from hikari.symmetry import SG, Group
from hikari.utility import make_abspath, mpl_map_palette, gnuplot_map_palette
from hikari.resources import potency_map_template


def potency_map(a, b, c, al, be, ga,
                space_group=SG['P1'],
                axis='',
                fix_scale=False,
                opening_angle=35,
                output_directory='~',
                output_name='cplt_map',
                output_quality=3,
                resolution=1.2,
                wavelength='MoKa'):
    """
    Calculate and draw a potency map for a given crystal in diamond anvil cell
    (dac) with a given opening angle, as a function of crystal orientation. For
    details see `this paper <https://doi.org/10.1107/S2052252521009532>`_.

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
    :param space_group: Instance of :class:`hikari.symmetry.Group`
        describing symmetry of the crystal
    :type space_group: hikari.symmetry.Group
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
    dat_path = make_abspath(output_directory, output_name + '.dat')
    gnu_path = make_abspath(output_directory, output_name + '.gnu')
    lst_path = make_abspath(output_directory, output_name + '.lst')
    png_path = make_abspath(output_directory, output_name + '.png')
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

        if space_group.system is Group.System.triclinic:
            _th_limits = [0, 180]
            _ph_limits = [0, 180]
        elif space_group.system in {Group.System.monoclinic,
                                    Group.System.orthorhombic,
                                    Group.System.tetragonal,
                                    Group.System.cubic}:
            _th_limits = [0, 90]
            _ph_limits = [0, 90]
        elif space_group.system in {Group.System.trigonal,
                                    Group.System.hexagonal}:
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
        my_heatmap_colors = mpl_map_palette[axis]
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
        gnu.write(potency_map_template.format(
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
            palette=gnuplot_map_palette[axis.lower()]))

    _prepare_gnuplot_input()
    try:
        from os import system, getcwd
        _path = make_abspath(output_directory)
        system('cd ' + _path + '; gnuplot ' + gnu_path)
    except OSError:
        pass


if __name__ == '__main__':
    pass
