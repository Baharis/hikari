"""
This sub-module contains all scripts aiming to help the work with
CrysAlisPRO software. Contrary to many other scripts in the package,
this sub-module is independent, as it does not use information from
any dataframes, but due to its usefulness has been supplied as well.
"""

import numpy as np


def crysalis_edge_mask(frame_width=2048,
                       frame_height=2048,
                       opening_radius=1000,
                       output_path='edge_mask.mac',
                       resolution=100):
    """
    Generate a .mac file containing a circular detector mask
    using "dc rejectrect" commands.
    This script has been prepared in oder to exclude any unnecessary
    empty space in the corners of the frames.
    Such a case might occur if utilised detector is circular
    and a significant area of the frame is empty.

    Please mind that the resolution should not be higher than default 100,
    as such file can be incorrectly imported by CrysAlis.
    The script might also fail to work if `opening_radius` is larger than
    half the `frame_width` or `frame_height`.

    :param frame_width: A width of a single frame in pixels.
    :type frame_width: int
    :param frame_height: A height of a single frame in pixels.
    :type frame_height: int
    :param opening_radius: Radius of the circle on the frame
        containing an actual image instead of shadow in pixels.
    :type opening_radius: int
    :param output_path: A path to created .mac file.
    :type output_path: str
    :param resolution: A number of rectangles used to exclude the area
        outside the observed circle. While this number should be kept high,
        you might want to lower it if you are using other "dc rejectrect"
        commands and the .mac file is not loading correctly.
    :type resolution: int
    :return: None
    """

    angles = np.arange(0, 2 * np.pi, 2 * np.pi / resolution)
    file = open(output_path, 'w')
    for a in angles:
        circle_x = int(frame_width / 2 + np.sin(a) * opening_radius)
        circle_y = int(frame_height / 2 + np.cos(a) * opening_radius)
        if 0 <= a < np.pi / 2:
            x_pos = circle_x
            y_pos = circle_y
            x_span = frame_width
            y_span = frame_height
        elif np.pi / 2 <= a < np.pi:
            x_pos = circle_x
            y_pos = 0
            x_span = frame_width
            y_span = circle_y
        elif np.pi <= a < 3 * np.pi / 2:
            x_pos = 0
            y_pos = 0
            x_span = circle_x
            y_span = circle_y
        else:
            x_pos = 0
            y_pos = circle_y
            x_span = circle_x
            y_span = frame_height
        file.write('dc rejectrect {} {} {} {}\n'.format(
            x_pos, y_pos, x_span, y_span))
    file.close()


if __name__ == '__name__':
    crysalis_edge_mask(frame_width=2048,
                       frame_height=2048,
                       opening_radius=1000,
                       output_path='edge_mask.mac',
                       resolution=100)
