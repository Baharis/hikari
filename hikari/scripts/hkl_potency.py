import numpy as np
import pandas as pd
import seaborn as sns
from matplotlib import pyplot
from numpy import linalg as lin

from hikari.dataframes import HklFrame
from hikari.symmetry import SG, Group
from hikari.utility import make_abspath, fibonacci_sphere, rotation_around
from hikari.scripts.angular_explorer import angular_property_explorer_factory


def potency_map(a, b, c, al, be, ga,
                space_group='P1',
                axis='',
                fix_scale=False,
                histogram=True,
                opening_angle=35,
                orientation=None,
                path='~/sortav.lst',
                output_quality=3,
                resolution=1.2,
                wavelength='MoKa'):
    """
    Calculate and draw a potency map for a given crystal in diamond anvil cell
    (DAC) with a given opening angle, as a function of crystal orientation.
    For details see `this paper <https://doi.org/10.1107/S2052252521009532>`_.

    The script accepts unit cell & space group information, and predicts max
    completeness of fully merged data for investigated crystal in said group.
    Results are logged into text files and drawn with gnuplot or matplotlib.

    Potency is calculated and visualised on a unit-sphere orientation heatmap.
    Each point **p** is associated with a certain crystal orientation in DAC,
    such that vector **v** from **0** to **p** acts as a symmetry axis for the
    dac-accessible volume traced in the reciprocal space up to `resolution`.
    Red / green / blue vectors represent crystallographic directions **X\***,
    **Y\*** and **Z\***, respectively. Since the distribution has some inherent
    symmetry, only a representative part of sphere (usually an octant) is shown.

    As an example, let's assume a orthorhombic cell with *a* = *b* = *c* = 10
    and laue group "mmm". Running the script and generating the potency the map
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
        For details see :py:mod:`hikari.symmetry.space_groups`.
    :type space_group: str or int
    :param axis: domain to calculate potency in. Accepts 'x'/ 'y'/ 'z' for h00/
        0k0/ 00l, 'xy'/'xz'/'yz' for hk0/ h0l/ 0kl, or '' for all reflections.
    :type axis: string
    :param fix_scale: If true, the colour scheme will fix to 0 - 100% range.
    :type fix_scale: bool
    :param histogram: If true, potency distribution will be plotted as histogram
    :type histogram: bool
    :param opening_angle: Value of single opening angle as defined in
        :meth:`hikari.dataframes.HklFrame.dac`.
    :type opening_angle: float
    :param orientation: either a cif-style 3x3 matrix of crystal orientation or
        a 3-length array with a diamond-perpendicular face to be marked on a map
    :type orientation: np.ndarray
    :param path: Path to a file where script is to be run. Extension is ignored,
        but file name will and must be the same for all input and output files.
    :type path: str
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
    kwargs = locals()
    ape = angular_property_explorer_factory.create(prop='potency')
    ape.set_up(**kwargs)
    ape.explore()
    ape.write_hist_file()
    ape.draw_matplotlib_map()
    ape.draw_gnuplot_map()


def potency_vs_dac_opening_angle(output_path='~/output.txt',
                                 precision=91,
                                 resolution=1.2,
                                 wavelength='MoKa',
                                 theta=None):
    """
    Calculate potency in P1 space group as a function of DAC opening angle,
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
    v = fibonacci_sphere(100)
    total = len(p)
    angles = np.linspace(start=90, stop=0, num=precision)
    with open(make_abspath(output_path), 'w', buffering=1) as out:
        out.write('#     oa potency\n')
        for a in angles:  # for one random vector v
            c = p.dacs_count(opening_angle=a, vectors=v)
            out.write(f' {a:7.4f} {np.mean(c)/total:7.5f}\n')


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
                            topple_angle=5):
    """
    For a given crystal, opening angle, and wavelength calculate average potency
    obtained by toppling the crystal by "topple_angle"° from the "vector" axis.

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
        For details see :py:mod:`hikari.symmetry.space_groups`.
    :type space_group: str or int
    :param wavelength: Wavelength of radiation utilised in experiment.
    :type wavelength: float or str
    :param opening_angle: Value of single opening angle as defined in
        :meth:`hikari.dataframes.HklFrame.dac`.
    :type opening_angle: float
    :param vector: Direction from which a theoretical crystal will be toppled.
    :type vector: tuple
    :param topple_angle: Angle by which theoretical crystal will be toppled.
    :type topple_angle: float
    :return: None
    """

    sg = SG[space_group]
    pg = sg.reciprocate()  # .lauefy()  # uncomment if hkl and -h-k-l are equiv.

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
    toppled_vector = rotate(v, perp, topple_angle)
    toppled_vectors = [rotate(toppled_vector, v, i) for i in range(360)]

    # calculate the potency for toppleds
    rads = [2.00, 1/0.7, 1.20, 1.00, 1/1.5]
    cplt = [0] * len(rads)
    p.trim(max(rads))
    for t in toppled_vectors:
        q = p.copy()
        q.dac_trim(opening_angle=opening_angle, vector=t)
        for cplt_bin, rad in enumerate(rads, start=1):
            q.trim(rad)
            cplt[-cplt_bin] += q.table['equiv'].nunique()

    # divide lists of potency by total and return from sums to individuals
    full = [0] * len(rads)
    for cplt_bin, rad in enumerate(rads, start=1):
        p.trim(rad)
        full[-cplt_bin] = p.table['equiv'].nunique()
    print(cplt)
    print(full)
    cplt = [c2 - c1 for c1, c2 in zip([0] + cplt[:-1], cplt)]
    full = [f2 - f1 for f1, f2 in zip([0] + full[:-1], full)]
    cplt = [c / (360 * f) for c, f in zip(cplt, full)]

    print('Resolution limits in distance to reflection, twice sin(θ/λ):')
    print(np.array(list(reversed(rads))))
    print(f'Average shell potency with {topple_angle} degree topple')
    print(np.array(cplt))


if __name__ == '__main__':
    potency_map(a=10, b=10, c=10, al=90, be=90, ga=120, space_group='P21/c',
                path='~/_/potency.hkl', output_quality=5, wavelength='MoKa')
