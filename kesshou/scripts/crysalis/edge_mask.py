import numpy as np


def edge_mask(frame_width=2048,
              frame_height=2048,
              opening_radius=1000,
              output_path='/home/dtchon/edge_mask.mac',
              resolution=100):

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
    edge_mask()
