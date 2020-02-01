from kesshou.dataframes.hkl import HklFrame
from kesshou.symmetry import PG
from kesshou.utility import cubespace
import numpy as np


def dac_statistics(a, b, c, al, be, ga,
                   point_group=PG['-1'],
                   opening_angle=35,
                   orientation=((1, 0, 0), (0, 1, 0), (0, 0, 1)),
                   input_path='shelx.hkl',
                   input_format=4,
                   input_wavelength='CuKa'):

    p = HklFrame()
    p.crystal.edit_cell(a=a, b=b, c=c, al=al, be=be, ga=ga)
    p.edit_wavelength(input_wavelength)
    p.read(input_path, input_format)
    p.crystal.orient_matrix = np.array(orientation)
    p.merge()
    p.extinct('000')
    p.dac(opening_angle=opening_angle)
    p.find_equivalents(point_group=point_group)
    p.make_stats(point_group=point_group)

    q = p.duplicate()
    q.extinct()
    q.make_ball(radius=q.r_lim)
    q.extinct('000')
    q.dac(opening_angle=opening_angle)
    q.find_equivalents(point_group=point_group)

    r_max = max(q.data['r'])
    print('radius    exp.      theory    exp.      theory    uniqueCplt')
    print('range     all       all       unique    unique    exp/theory')
    for rad in reversed(cubespace(0, r_max, 10, include_start=False)):
        p.trim(rad)
        q.trim(rad)
        print(' {max_rad:9f} {exp_all:9d} {the_all:9d}' +
              ' {exp_uni:9d} {the_uni:9d} {cplt:9f}'.format(
                  max_rad=rad,
                  exp_all=len(p),
                  the_all=len(q),
                  exp_uni=p.data['equiv'].nunique(),
                  the_uni=q.data['equiv'].nunique(),
                  cplt=p.data['equiv'].nunique() / q.data['equiv'].nunique()))


if __name__ == '__main__':
    dac_statistics(a=10.0, b=10.0, c=10.0, al=90.0, be=90.0, ga=90.0)
