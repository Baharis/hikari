import numpy as np
from PIL import Image
from kesshou.utility import make_absolute_path


def merge_axis_completeness_maps(directory, r, g, b, out):
    image1 = Image.open(make_absolute_path(directory, r))
    image2 = Image.open(make_absolute_path(directory, g))
    image3 = Image.open(make_absolute_path(directory, b))
    o = make_absolute_path(directory, out)
    r_channel = image1.getchannel(0)
    g_channel = image2.getchannel(1)
    b_channel = image3.getchannel(2)
    try:
        a_channel = image1.getchannel(3)
        Image.merge('RGBA', (r_channel, g_channel, b_channel, a_channel))\
            .save(o, 'PNG')
    except ValueError:
        Image.merge('RGB', (r_channel, g_channel, b_channel)).save(o, 'PNG')


if __name__ == '__main__':
    start = 'CpltMap_Mo35a'
    group = '_Pm-3'
    end = '_gnu.png'
    merge_axis_completeness_maps(
        directory='/home/dtchon/x/HP/DAC_completeness/cplt_vs_orientation/spherical_Mo35a/',
        r=start+'X'+group+end,
        g=start+'Y'+group+end,
        b=start+'Z'+group+end,
        out=start+group+end
    )
