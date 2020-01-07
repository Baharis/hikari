# ~~~~~~~~~~~~~~~~~~~~~ IMPORT STATEMENTS - DO NOT CHANGE ~~~~~~~~~~~~~~~~~~~~ #
import numpy as np
# ~~~~~~~~~~~~~~~~~~~~ VARIABLES - CHANGE ONLY VALUES HERE ~~~~~~~~~~~~~~~~~~~ #

# Input details
opening_radius = 1000
frame_width = 2048
frame_height = 2048
resolution = 100
output_file_path = '/home/dtchon/edge_mask.mac'

# ~~~~~~~~~~~~~~~~~~~~~~~~ SCRIPT CODE - DO NOT CHANGE ~~~~~~~~~~~~~~~~~~~~~~~ #

angles = np.arange(0, 2*np.pi, 2*np.pi/resolution)
file = open(output_file_path, 'w')

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
    line = 'dc rejectrect {} {} {} {}\n'.format(x_pos, y_pos, x_span, y_span)
    file.write(line)

file.close()




# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ END OF FILE ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ #
