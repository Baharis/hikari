from kesshou.dataframes.hkl import HklFrame
import numpy as np


def simulate_dac(a, b, c, al, be, ga,
                 opening_angle=35,
                 orientation=((1, 0, 0), (0, 1, 0), (0, 0, 1)),
                 vector=(0, 0, 1),
                 resolution=None,
                 input_path=None,
                 input_format=4,
                 wavelength='CuKa',
                 output_path='output.hkl',
                 output_format=4):
    """Simulate high-pressure conditions using supplied or simulated data.
    If vector is given, use it to define a direction perpendicular to dac
    surface instead of using crystal orientation matrix.
    If input_path is given, use real hkl data instead of simulated one."""
    p = HklFrame()
    p.crystal.edit_cell(a=a, b=b, c=c, al=al, be=be, ga=ga)
    p.edit_wavelength(wavelength)
    if not(input_path is None):
        p.read(input_path, input_format)
    else:
        p.make_ball(p.r_lim)

    p.crystal.orient_matrix = np.array(orientation)
    p.extinct('000')
    if not(resolution is None):
        p.trim(resolution)
    p.dac(opening_angle=opening_angle, vector=vector)
    p.write(hkl_path=output_path, hkl_format=output_format)


if __name__ == '__main__':
    simulate_dac(a=10.0, b=10.0, c=10.0, al=90.0, be=90.0, ga=90.0)
