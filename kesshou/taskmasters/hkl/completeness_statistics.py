from kesshou.dataframes.hkl import HklFrame
from kesshou.symmetry.pointgroup import *


def completeness_statistics(a, b, c, al, be, ga,
                            point_group=PG_1,
                            input_path='shelx.hkl',
                            input_format=4,
                            input_wavelength='CuKa'):

    p = HklFrame()
    p.crystal.edit_cell(a=a, b=b, c=c, al=al, be=be, ga=ga)
    p.edit_wavelength(input_wavelength)
    p.read(input_path, input_format)
    p.extinct('000')
    p.make_stats(point_group=point_group)


if __name__ == '__name__':
    completeness_statistics(a=10.0, b=10.0, c=10.0, al=90.0, be=90.0, ga=90.0)
