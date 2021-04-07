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
    merge_axis_completeness_maps(
        directory='/home/dtchon/x/_/',
        r='cplt_map_tests_x.png',
        g='cplt_map_tests_y.png',
        b='cplt_map_tests_z.png',
        out='cplt_map_tests_axes.png'
    )
