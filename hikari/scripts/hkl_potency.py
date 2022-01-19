import numpy as np
import pandas as pd
import seaborn as sns
from matplotlib import pyplot, colors, cm
from mpl_toolkits.mplot3d import art3d
from numpy import linalg as lin

from hikari.dataframes import HklFrame
from hikari.symmetry import SG, Group
from hikari.utility import make_abspath, mpl_map_palette, gnuplot_map_palette, \
    fibonacci_sphere, rotation_around
from hikari.resources import potency_map_template


def potency_map(a, b, c, al, be, ga,
                space_group='P1',
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
    (DAC) with a given opening angle, as a function of crystal orientation.
    For details see `this paper <https://doi.org/10.1107/S2052252521009532>`_.

    The script accepts unit cell & space group information, and predicts the
    completeness of fully merged data for investigated crystal in said group.
    Results are logged into text files and drawn with gnuplot or matplotlib,
    depending on settings.

    Potency is calculated and visualised on a unit-sphere orientation heatmap.
    Each point **p** is associated with a certain crystal orientation in DAC,
    such that vector **v** from **0** to **p** acts as a symmetry axis for the
    dac-accessible volume traced in the reciprocal space up to `resolution`.
    Red / green / blue vectors represent crystallographic directions **X\***,
    **Y\*** and **Z\***, respectively. Since the distribution has some inherent
    symmetry, only a representative part of sphere (usually an octant) is shown.

    As an example, let's assume a orthorhombic cell with *a* = *b* = *c* = 10
    and laue group "mmm". Running the script and generating completeness the map
    yields the lowest values close to **X\***, **Y\*** and **Z\*** vector,
    while the highest values are observed between those vectors.
    Placing the crystal on its [100] face inside the dac will cause the
    dac-accessible plane to be placed perpendicularly to (100) direction
    in reciprocal space. Since the values close to **X\*** ((100) direction)
    are low, such a placement will allow us to collect data with low coverage.
    On the other hand, placing the crystal on its [111] face will cause the
    dac-accessible plane to be placed perpendicularly to (100) direction
    in reciprocal space. As on the sphere the values tend to rise the further
    we are from of **X\***, **Y\*** and **Z\***, we expect this crystal
    orientation to warrant a high completeness of collected data.

    The potency is calculated as a ratio of the number of unique reflections
    inside the DAC-accessible space to the number of unique reflections
    inside a reference sphere of a radius equal to `resolution`.
    Visualisation is performed using an extended rainbow heatmap,
    which utilises a wide color range to emphasize even small differences.
    Dy default, the color scale is dynamic and adapts to the range of calculated
    potency, but it can be fixed to 0-100% range using `fix_scale=True`.

    Since the orientation is given in spherical coordinates, the exact positions
    of individual points is given using *theta* and *phi* angles instead of
    crystallographic coordinates. The *theta* and *phi* angles here follow the
    physical definition (ISO 80000-2:2019). For a given DAC-axis vector **v**:

    - *theta* is the azimuth angle found between vector **Z\*** and **v**.
      It can assume values between zero and *pi* (*pi/2* in one octant).

    - *phi* is the rotational angle denoting rotation of **v** around **Z\***.
      It can assume values between zero and *2 pi* (*pi/2* in one octant).

    Finally it must be noted that higher potency setting is not always better.
    For example, for opening angle below 45 degrees, most potent orientation in
    laue class mmm renders all 0kl, h0l, hk0 reflections inaccessible.
    Potency map should be consulted before a high-pressure experiment, but it
    should not be treated as an universal quality indicator of given set-up.

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
    :param space_group: Short Hermann-Mauguin name or index of space group.
        For details see table in hikari.symmetry.space_groups.
    :type space_group: str or int
    :param axis: domain to calculate potency in. Accepts 'x'/ 'y'/ 'z' for h00/
        0k0/ 00l, 'xy'/'xz'/'yz' for hk0/ h0l/ 0kl, or '' for all reflections.
    :type axis: string
    :param fix_scale: If true, the colour scheme will fix to 0 - 100% range.
    :type fix_scale: bool
    :param opening_angle: Value of single opening angle as defined in
        :meth:`hikari.dataframes.HklFrame.dac`.
    :type opening_angle: float
    :param output_directory: Path to directory where output should be saved.
    :type output_directory: str
    :param output_name: Base name for files created in `output_directory`.
    :type output_name: str
    :param output_quality: Density of individual orientations to be considered.
        Should be in range from 1 (every 15 degrees) to 5 (every 1 degree).
    :type output_quality: int
    :param resolution: Upper limit of reflection resolution, given as a distance
        from zero to node in reciprocal space (one over plane spacing) in A-1.
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
    sg = SG[space_group]
    lg = sg.reciprocate()

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
            _f.transform([o.tf for o in lg.operations])
        _f.extinct(sg)
        return _f

    p = _make_hkl_frame()
    p.find_equivalents(point_group=lg)
    total_reflections = p.table['equiv'].nunique()
    if total_reflections == 0:
        raise KeyError('Specified part of reciprocal space contains zero nodes')

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
    data_dict = {'th': [], 'ph': [], 'cplt': [], 'reflns': []}

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
        _cplt_mesh = np.zeros_like(th_mesh)
        lst = open(lst_path, 'w+')
        lst.write('#     th      ph    cplt  reflns\n')
        vectors = np.vstack([_angles_to_vector(th, ph) for th in th_range
                                             for ph in ph_range])
        uniques = p.dacs_count(opening_angle=opening_angle, vectors=vectors)
        for i, th in enumerate(th_range):
            for j, ph in enumerate(ph_range):
                count = uniques[i * len(ph_range) + j]
                potency = count / total_reflections
                data_dict['th'].append(th)
                data_dict['ph'].append(ph)
                data_dict['cplt'].append(potency)
                data_dict['reflns'].append(count)
                lst.write(f'{th:8.0f}{ph:8.0f}{potency:8.5f}{count:8d}\n')
                _cplt_mesh[j][i] = potency
            lst.write('\n')
        index_max = np.unravel_index(np.argmax(_cplt_mesh), _cplt_mesh.shape)
        best_th, best_ph = th_range[index_max[1]], ph_range[index_max[0]]
        index_min = np.unravel_index(np.argmin(_cplt_mesh), _cplt_mesh.shape)
        worst_th, worst_ph = th_range[index_min[1]], ph_range[index_min[0]]
        max_p = max(data_dict['cplt'])
        min_p = min(data_dict['cplt'])
        mean_p = np.mean(data_dict['cplt'])
        lst.write(f'# best cplt = {max_p} for th = {best_th}, ph = {best_ph}\n'
                  f'# worst cplt = {min_p} for th = {worst_th}, ph = {worst_ph}'
                  f'\n# mean cplt = {mean_p}\n')
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
        np.warnings.filterwarnings('ignore',
                                   category=np.VisibleDeprecationWarning)
        ax.plot_wireframe(x, y, z, colors='k', linewidth=0.25)  # <- mpl warning

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


def potency_vs_dac_opening_angle(output_path='~/output.txt',
                                 precision=91,
                                 resolution=1.2,
                                 wavelength='MoKa',
                                 theta=None):
    """
    Calculate completeness in P1 as a function of DAC opening angle,
    assuming certain resolution and wavelength used.

    :param output_path: Path of created file containing calculated data.
    :type output_path: str
    :param precision: Number of probed opening angles between 0 and 90 degrees.
    :type precision: int
    :param resolution: Maximum distance from the origin in reciprocal space to
        reflection (twice the resolution in reciprocal angstrom). Default 1.2.
    :type resolution: float
    :param wavelength: Wavelength of radiation to be simulated.
    :type wavelength: float or str
    :param theta: If given, use max theta in degrees (instead of radius in A-1).
    :type theta: float
    :return: None
    """

    def _make_reference_ball():
        hkl_frame = HklFrame()
        hkl_frame.la = wavelength
        res = hkl_frame.r_lim * np.sin(np.deg2rad(theta)) \
            if theta is not None else resolution
        side = 10 * precision**(1/3) / res  # adapt to res&la
        hkl_frame.edit_cell(a=side, b=side, c=side, al=90, be=90, ga=90)
        hkl_frame.fill(radius=res)
        hkl_frame.find_equivalents()
        return hkl_frame

    p = _make_reference_ball()
    v = fibonacci_sphere(10)
    total = len(p)
    angles = np.linspace(start=90, stop=0, num=precision)
    out = open(make_abspath(output_path), 'w')
    out.write('#oa      cplt\n')
    for a in angles:  # for one random vector v
        c = p.dacs_count(opening_angle=a, vectors=v)
        out.write(' {a:7.4f} {c:7.5f}\n'.format(a=a, c=np.mean(c)/total))
    out.close()


laue_space_groups = 'P-1', 'P2/m', 'Pmmm', 'P4/m', 'P4/mmm', 'P-3', 'P-3m1',\
                    'P6/m', 'P6/mmm', 'Pm-3', 'Pm-3m'
laue_class_names = r"$\overline{1}$", r"$2/m$", r"$mmm$",\
                   r"$4/m$", r"$4/mmm$", r"$\overline{3}$", r"$\overline{3}m$",\
                   r"$6/m$", r"$6/mmm$", r"$m\overline{3}$", r"$m\overline{3}m$"


def potency_violin_plot(job_name='violin',
                        directory='~',
                        opening_angle=35,
                        precision=1000,
                        space_groups=laue_space_groups,
                        labels=laue_class_names,
                        resolution=None,
                        wavelength='MoKa'):
    """
    Calculate potency distribution for selected space groups and multiple sample
    orientations in a DAC, log it, and (re)draw appropriate violin plot.

    :param job_name: Base name for created files, defaults to 'violin'. If a
        .csv file with this name already exists, redraw instead of calculating.
    :type job_name: str
    :param directory: Target directory for saving or reading, defaults to '~'.
    :type directory: str
    :param opening_angle: Value of single opening angle as defined in
        :meth:`hikari.dataframes.HklFrame.dac`. Defaults to 35.
    :type opening_angle: float
    :param precision: Number of orientations to investigate, defaults to 1000.
    :type precision: int
    :param space_groups: List of space groups to investigate, strings or ints
        (see :class:`hikari.symmetry.Group`). Defaults to 11 "Laue" groups.
    :type space_groups: Tuple[int] or Tuple[str]
    :param labels: List of tex-style labels to be used for logging and plotting.
    :type labels: Tuple[str]
    :param resolution: If given, additionally limit data resolution to given
        value. Please provide the resolution as a distance from the origin
        in reciprocal space (twice the resolution in reciprocal angstrom).
    :type resolution: float
    :param wavelength: Wavelength of radiation to be simulated.
    :type wavelength: float or str
    :return: None
    :rtype: None
    """

    csv_path = make_abspath(directory, job_name + '.csv')
    png_path = make_abspath(directory, job_name + '.png')
    log_path = make_abspath(directory, job_name + '.log')
    space_groups = [SG[sg] for sg in space_groups]

    def _is_tri_or_hexagonal(sg):
        return sg.system in {sg.System.hexagonal, sg.System.trigonal}

    def _make_ball(ga):
        hkl_frame = HklFrame()
        hkl_frame.edit_cell(a=20, b=20, c=20, al=90, be=90, ga=ga)
        hkl_frame.la = wavelength
        hkl_frame.fill(radius=resolution if resolution else hkl_frame.r_lim)
        return hkl_frame

    h = _make_ball(ga=120) if any(_is_tri_or_hexagonal(sg)
                                  for sg in space_groups) else None
    c = _make_ball(ga=90) if not all(_is_tri_or_hexagonal(sg)
                                     for sg in space_groups) else None
    vectors = np.vstack([np.array([[1, 0, 0], [0, 1, 0], [0, 0, 1],
                                   [0.57735, 0.57735, 0.57735]]),
                         fibonacci_sphere(samples=precision-4, seed=1337)])

    try:
        cplt_frame = pd.read_csv(csv_path, index_col=0)
    except FileNotFoundError:
        cplt_dict = dict()
        log = open(log_path, 'w', buffering=1)
        for label, sg in zip(labels, space_groups):
            p = h.copy() if _is_tri_or_hexagonal(sg) else c.copy()
            p.find_equivalents(point_group=sg.reciprocate())
            p.extinct(space_group=sg)
            total_reflections = p.table['equiv'].nunique()
            log.write('space_group: ' + str(label) + '\n')
            log.write('total_reflections: ' + str(total_reflections) + '\n')
            log.write('max_r_in_reciprocal: ' + str(max(p.table['r'])) + '\n')
            reflections = p.dacs_count(opening_angle, vectors)
            for v, r in zip(vectors, reflections):
                log.write(str(v)+': '+str(r)+'\n')
            log.write('max_reflections: ' + str(max(reflections)) + '\n')
            log.write('min_reflections: ' + str(min(reflections)) + '\n')
            log.write('avg_reflections: ' + str(np.mean(reflections)) + '\n\n')
            cplt_dict[label] = [r / total_reflections for r in reflections]

        cplt_frame = pd.DataFrame.from_dict(cplt_dict)
        cplt_frame = cplt_frame.mul(100)
        cplt_frame.to_csv(csv_path)
        pd.set_option('display.expand_frame_repr', False)
        log.write('Summary:\n' + str(cplt_frame.describe().round(2)) + '\n')
        log.close()

    palette_number = {Group.System.triclinic: 0,
                      Group.System.monoclinic: 1,
                      Group.System.orthorhombic: 2,
                      Group.System.tetragonal: 3,
                      Group.System.trigonal: 4,
                      Group.System.hexagonal: 5,
                      Group.System.cubic: 6}
    palette = {label: sns.color_palette("deep")[palette_number[sg.system]]
               for label, sg in zip(labels, space_groups)}

    melt_frame = cplt_frame.melt(var_name='SG', value_name='Completeness')
    sns.set_style("whitegrid")
    pyplot.figure(figsize=(8, 6))
    sns.set(font_scale=1.5)
    ax = sns.violinplot(x=melt_frame['Completeness'],
                        y=melt_frame['SG'], scale='width', cut=0,
                        palette=palette, inner='quartile')
    for line in ax.lines[1::3]:
        line.set_linestyle('--')
        line.set_linewidth(1.25)
    ax.set_xlabel('')
    ax.set_ylabel('')
    for y, key in enumerate(cplt_frame.keys()):
        ax.axvline(x=cplt_frame[key].median(), color=palette[key],
                   ymin=0, ymax=21/22-y/11, linewidth=5, alpha=0, zorder=0.99)
    ax.set_yticklabels(labels)
    ax.set_xticklabels(['10', r'$\mathcal{P}$ = 20       ', '30', '40',
                        '50', '60', '70', '80', '90', '100%'])
    pyplot.xlim(19, 101)
    sns.despine(left=True, top=False)
    pyplot.savefig(png_path, dpi=600, format='png', bbox_inches=None)


def dac_potency_around_axis(a, b, c, al, be, ga,
                            space_group='P1',
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
    :param be: Unit cell parameter *beta* in degrees.
    :type be: float
    :param ga: Unit cell parameter *gamma* in degrees.
    :type ga: float
    :param space_group: Short Hermann-Mauguin name or index of space group.
        For details see table in hikari.symmetry.space_groups.
    :type space_group: str or int
    :param wavelength: Wavelength of radiation utilised in experiment.
    :type wavelength: float or str
    :param opening_angle: Value of single opening angle as defined in
        :meth:`hikari.dataframes.HklFrame.dac`.
    :type opening_angle: float
    :param vector: Direction around which completeness will be calculated.
    :type vector: tuple
    :param topple: Angle by which vector will be toppled in varous directions.
    :type topple: float
    :return: None
    """

    sg = SG[space_group]
    pg = sg.reciprocate()  # .lauefy()

    p = HklFrame()
    p.edit_cell(a=a, b=b, c=c, al=al, be=be, ga=ga)
    p.la = wavelength
    p.fill(radius=p.r_lim)
    p.place()
    p.extinct(space_group=sg)
    p.find_equivalents(point_group=pg)

    # generate perpendicular vector for rotation
    v = np.array(vector) / lin.norm(np.array(vector))
    x = np.array((1, 0, 0))

    def are_parallel(v1, v2):
        return np.allclose(v1 / lin.norm(v1), v2 / lin.norm(v2))
    temp = x if not(are_parallel(v, x)) else np.array((0, 1, 0))
    perp = np.cross(v, temp)

    def rotate(_v, _k, angle):
        return _v @ rotation_around(_k, by=np.deg2rad(angle))

    # generate 360 toppled vectors
    toppled = rotate(v, perp, topple)
    toppleds = [rotate(toppled, v, i) for i in range(360)]

    # calculate the completeness for toppleds
    rads = [2.00, 1/0.7, 1.20, 1.00, 1/1.5]
    cplt = [0] * len(rads)
    p.trim(max(rads))
    for t in toppleds:
        q = p.copy()
        q.dac_trim(opening_angle=opening_angle, vector=t)
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


if __name__ == '__main__':
    # potency_map(10, 10, 10, 90, 90, 90, space_group='Pmmm',
    #             resolution=1.2, output_directory='~/_/', output_name='_')
    pass
