from kesshou.dataframes.hkl import HklFrame
from kesshou.symmetry.pointgroup import *
from kesshou.utility import cubespace
import copy


def dac_statistics(a, b, c, al, be, ga,
                   point_group=PG_1,
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
    p.extinct('000')
    p.dac(opening_angle=opening_angle)
    p.make_stats(point_group=point_group)

    q = copy.deepcopy(p)
    q.extinct()
    q.make_ball(radius=q.r_lim)
    q.dac(opening_angle=opening_angle)
    q.extinct('000')

    r = copy.deepcopy(q)
    r.resymmetrify(point_group.chiral_operations)

    r_max = max(q.data['r'])
    print('limiting_radius experiment theory_no_symm theory_w/_symm')
    for rad in reversed(cubespace(0, r_max, 10, include_start=False)):
        p.trim(rad)
        q.trim(rad)
        r.trim(rad)
        print(rad, len(p), len(r), len(q))


if __name__ == '__main__':
    dac_statistics(a=10.0, b=10.0, c=10.0, al=90.0, be=90.0, ga=90.0)
