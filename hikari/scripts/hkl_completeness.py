"""
This sub-module contains scripts made to work mainly with .hkl files,
containing information from single crystal diffraction experiments.
Due to authors interest in the subjects, majority of them aim to analyse
or model data from high-pressure experiments,
which utilise DAC - diamond anvil cell.
"""

from hikari.dataframes import HklFrame
from hikari.symmetry import SG
from hikari.utility import cubespace
import numpy as np


def completeness_statistics(a, b, c, al, be, ga,
                            space_group='P-1',
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
    :param space_group: Short Hermann-Mauguin name or index of space group.
        For details see table in hikari.symmetry.space_groups.
    :type space_group: str or int
    :param input_path: Path to the input .hkl file.
    :type input_path: str
    :param input_format: Format of the .hkl file. For reference see
        :meth:`hikari.dataframes.HklFrame.interpret_hkl_format`.
    :type input_format: int or str or dict
    :param input_wavelength: Wavelength of radiation utilised in experiment.
    :type input_wavelength: float or str
    :return: None
    :rtype: None
    """
    p = HklFrame()
    p.edit_cell(a=a, b=b, c=c, al=al, be=be, ga=ga)
    p.la = input_wavelength
    p.read(input_path, input_format)
    p.stats(space_group=SG[space_group])


def dac_statistics(a, b, c, al, be, ga,
                   space_group='P1',
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
    :param space_group: Short Hermann-Mauguin name or index of space group.
        For details see table in hikari.symmetry.space_groups.
    :type space_group: str or int
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

    sg = SG[space_group]
    pg = sg.reciprocate()

    p = HklFrame()
    p.edit_cell(a=a, b=b, c=c, al=al, be=be, ga=ga)
    p.la = input_wavelength
    p.read(input_path, input_format)
    p.orientation = np.array(orientation)

    resolution = p.r_lim if resolution is None else resolution
    p.trim(limit=resolution)
    p.merge(point_group=pg)

    q = p.copy()
    q.fill(radius=resolution)
    q.dac_trim(opening_angle=opening_angle)

    # uncomment and fill this if you have 2 crystals in different orientations
    # q2 = p.copy()
    # q2.fill(radius=resolution)
    # q2.orientation = np.array(())
    # q2.dac(opening_angle=opening_angle)
    # q = q + q2

    q.merge(point_group=pg)
    q.extinct(space_group=sg)

    b = p.copy()
    b.fill(radius=resolution)
    b.extinct(space_group=sg)
    b.merge(point_group=pg)

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
    p.dac_trim(opening_angle=opening_angle, vector=vector)
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
