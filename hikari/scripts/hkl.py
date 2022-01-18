"""
This sub-module contains scripts made to work mainly with .hkl files,
containing information from single crystal diffraction experiments.
Due to authors interest in the subjects, majority of them aim to analyse
or model data from high-pressure experiments,
which utilise DAC - diamond anvil cell.
"""

import pandas as pd
import seaborn as sns
from hikari.dataframes import HklFrame
from hikari.symmetry import Group, PG, SG
from hikari.utility import cubespace, fibonacci_sphere, make_abspath
from matplotlib import pyplot
import numpy as np
import numpy.linalg as lin


def completeness_statistics(a, b, c, al, be, ga,
                            space_group=SG['P-1'],
                            input_path='shelx.hkl',
                            input_format=4,
                            input_wavelength='CuKa'):
    """
    For a given experimental .hkl file calculate basic completeness statistics
    in ten resolution shells of equal-volume.

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
    :param space_group: Space group of the crystal.
    :type space_group: hikari.symmetry.Group
    :param input_path: Path to the input .hkl file.
    :type input_path: str
    :param input_format: Format of the .hkl file. For reference see
        :meth:`hikari.dataframes.HklFrame.interpret_hkl_format`.
    :type input_format: int or str or dict
    :param input_wavelength: Wavelength of radiation utilised in experiment.
    :type input_wavelength: float or str
    :return: None
    """
    p = HklFrame()
    p.edit_cell(a=a, b=b, c=c, al=al, be=be, ga=ga)
    p.la = input_wavelength
    p.read(input_path, input_format)
    p.stats(space_group=space_group)


def dac_completeness_vs_opening_angle(output_path='output.txt',
                                      precision=91,
                                      resolution=1.2,
                                      wavelength='MoKa',
                                      theta=None):
    """
    Calculate completeness in P-1 as a function of DAC opening angle,
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
        side = 10 * precision**(1/3) / res  # adapt to prec&la
        hkl_frame.edit_cell(a=side, b=side, c=side, al=90, be=90, ga=90)
        hkl_frame.fill(radius=res)
        return hkl_frame
    p = _make_reference_ball()
    total = len(p.table)
    angles = np.linspace(start=90, stop=0, num=precision)
    out = open(output_path, 'w', buffering=1)
    out.write('#oa      cplt\n')
    for a in angles:
        p.dac(opening_angle=a, vector=np.array((1, np.pi/4, np.e/5)))  # rand. v
        out.write(' {a:7.4f} {c:7.5f}\n'.format(a=a, c=len(p.table)/total))
    out.close()


laue_space_groups = SG['P-1'], SG['P2/m'], SG['Pmmm'], \
                    SG['P4/m'], SG['P4/mmm'], SG['P-3'], SG['P-3m1'], \
                    SG['P6/m'], SG['P6/mmm'], SG['Pm-3'], SG['Pm-3m']
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
    :param space_groups: List of space groups to investigat, given as instances
        of :class:`hikari.symmetry.Group`. Defaults to 11 "Laue" groups.
    :type space_groups: list[hikari.symmetry.Group]
    :param labels: List of git-style labels to be used for logging and plotting.
    :type labels: list[str]
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
            reflections = list()
            total_reflections = p.table['equiv'].nunique()
            log.write('space_group: ' + str(label) + '\n')
            log.write('total_reflections: ' + str(total_reflections) + '\n')
            log.write('max_r_in_reciprocal: ' + str(max(p.table['r'])) + '\n')
            for vector in vectors:
                q = p.copy()
                q.dac(opening_angle=opening_angle, vector=vector)
                reflections.append(q.table['equiv'].nunique())
                log.write(str(vector)+': '+str(q.table['equiv'].nunique())+'\n')
            log.write('max_reflections: ' + str(max(reflections)) + '\n')
            log.write('min_reflections: ' + str(min(reflections)) + '\n')
            log.write('avg_reflections: ' + str(np.mean(reflections)) + '\n')
            cplt_dict[label] = [r / total_reflections for r in reflections]

        cplt_frame = pd.DataFrame.from_dict(cplt_dict)
        cplt_frame = cplt_frame.mul(100)
        cplt_frame.to_csv(csv_path)
        pd.set_option('display.expand_frame_repr', False)
        log.write('\nSummary:\n' + str(cplt_frame.describe().round(2)) + '\n')
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
    :param space_group: Space group of the crystal,
        defined as an instance of :class:`hikari.symmetry.Group`
    :type space_group: hikari.symmetry.Group
    :param opening_angle: Value of single opening angle as defined in
        :meth:`hikari.dataframes.HklFrame.dac`.
    :type opening_angle: float
    :param orientation: Crystal orientation as defined in
        :class:`hikari.dataframes.BaseFrame`
    :type orientation: tuple or numpy.array
    :param input_path: Path to the input .hkl file.
    :type input_path: str
    :param input_format: Format of the .hkl file. For reference see
        :meth:`hikari.dataframes.HklFrame.interpret_hkl_format`.
    :type input_format: int or str or dict
    :param input_wavelength: Wavelength of radiation utilised in experiment.
    :type input_wavelength: float or str
    :param resolution: If given, calculate statistics only up to this value.
    Please provide it as a distance in rec. space (twice the resolution in A-1).
    :type resolution: float
    :return: None
    """

    point_group = space_group.reciprocate()

    p = HklFrame()
    p.edit_cell(a=a, b=b, c=c, al=al, be=be, ga=ga)
    p.la = input_wavelength
    p.read(input_path, input_format)
    p.orientation = np.array(orientation)

    resolution = p.r_lim if resolution is None else resolution
    p.trim(limit=resolution)
    p.merge(point_group=point_group)

    q = p.copy()
    q.fill(radius=resolution)
    q.dac(opening_angle=opening_angle)

    # uncomment and fill this if you have 2 crystals in different orientations
    # q2 = p.copy()
    # q2.fill(radius=resolution)
    # q2.orientation = np.array(())
    # q2.dac(opening_angle=opening_angle)
    # q = q + q2

    q.merge(point_group=point_group)
    q.extinct(space_group=space_group)

    b = p.copy()
    b.fill(radius=resolution)
    b.extinct(space_group=space_group)
    b.merge(point_group=point_group)

    r_max = resolution
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
        defined as an instance of :class:`hikari.symmetry.Group`
    :type space_group: hikari.symmetry.Group
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

    point_group = space_group.reciprocate()  # .lauefy()

    p = HklFrame()
    p.edit_cell(a=a, b=b, c=c, al=al, be=be, ga=ga)
    p.la = wavelength
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

    def rotate(_v, _k, angle):
        _c = np.cos(np.deg2rad(angle))
        _s = np.sin(np.deg2rad(angle))
        return _v * _c + np.cross(_k, _v) * _s + _k * np.dot(_k, _v) * (1 - _c)

    # generate 360 toppled vectors
    toppled = rotate(v, perp, topple)
    toppleds = [rotate(toppled, v, i) for i in range(360)]

    # calculate the completeness for toppleds
    rads = [2.00, 1/0.7, 1.20, 1.00, 1/1.5]
    cplt = [0] * len(rads)
    p.trim(max(rads))
    for t in toppleds:
        q = p.copy()
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
        :meth:`hikari.dataframes.HklFrame.dac`.
    :type opening_angle: float
    :param orientation: Crystal orientation as defined in
        :class:`hikari.dataframes.BaseFrame`
    :type orientation: tuple or numpy.array
    :param vector: If given, overwrite orientation to provide information
        about crystal placement in dac,
        as defined in :meth:`hikari.dataframes.HklFrame.dac`.
    :type vector: tuple
    :param resolution: If given, additionally limit data resolution to given
        value. Please provide the resolution as a distance from the origin
        in reciprocal space (twice the resolution in reciprocal angstrom).
    :type resolution: float
    :param input_path: Path to the input .hkl file.
    :type input_path: str
    :param input_format: Format of the input .hkl file. For reference see
        :meth:`hikari.dataframes.HklFrame.interpret_hkl_format`.
    :type input_format: int or str or dict
    :param input_wavelength: Wavelength of radiation utilised in experiment.
    :type input_wavelength: float or str
    :param output_path: Path to the output .hkl file.
    :type output_path: str
    :param output_format: Format of the input .hkl file. For reference see
        :meth:`hikari.dataframes.HklFrame.interpret_hkl_format`.
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


if __name__ == '__main__':
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
