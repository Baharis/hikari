gnuplot_map_palette = dict()
gnuplot_map_palette[''] = 'set palette model HSV defined (0 0.667 1 0,' + \
                          ' 20 0.667 1 1, 80 0.0 1 1, 100 0.0 0 1)'
gnuplot_map_palette['h'] = 'set palette model RGB defined (0 1 1 1, 1 1 0 0)'
gnuplot_map_palette['x'] = 'set palette model RGB defined (0 0 0 0, 1 1 0 0)'
gnuplot_map_palette['k'] = 'set palette model RGB defined (0 1 1 1, 1 0 1 0)'
gnuplot_map_palette['y'] = 'set palette model RGB defined (0 0 0 0, 1 0 1 0)'
gnuplot_map_palette['l'] = 'set palette model RGB defined (0 1 1 1, 1 0 0 1)'
gnuplot_map_palette['z'] = 'set palette model RGB defined (0 0 0 0, 1 0 0 1)'
gnuplot_map_palette['hk'] = gnuplot_map_palette['l']
gnuplot_map_palette['xy'] = gnuplot_map_palette['z']
gnuplot_map_palette['kl'] = gnuplot_map_palette['h']
gnuplot_map_palette['yz'] = gnuplot_map_palette['x']
gnuplot_map_palette['hl'] = gnuplot_map_palette['k']
gnuplot_map_palette['xz'] = gnuplot_map_palette['y']

mpl_map_palette = dict()
mpl_map_palette[''] = [(0.0, 0.0, 0.0),  # Black
                       (0.0, 0.0, 0.5),  # Indigo
                       (0.0, 0.0, 1.0),  # Blue
                       (0.0, 0.5, 1.0),  # Aqua
                       (0.0, 1.0, 1.0),  # Turquoise
                       (0.5, 1.0, 0.5),  # Lime
                       (1.0, 1.0, 0.0),  # Yellow
                       (1.0, 0.5, 0.0),  # Orange
                       (1.0, 0.0, 0.0),  # Red
                       (1.0, 0.5, 0.5),  # Light red
                       (1.0, 1.0, 1.0)]  # White
mpl_map_palette['h'] = [(1.0, 1.0, 1.0), (1.0, 0.0, 0.0)]
mpl_map_palette['x'] = [(0.0, 0.0, 0.0), (1.0, 0.0, 0.0)]
mpl_map_palette['k'] = [(1.0, 1.0, 1.0), (0.0, 1.0, 0.0)]
mpl_map_palette['y'] = [(0.0, 0.0, 0.0), (0.0, 1.0, 0.0)]
mpl_map_palette['l'] = [(1.0, 1.0, 1.0), (0.0, 0.0, 1.0)]
mpl_map_palette['z'] = [(0.0, 0.0, 0.0), (0.0, 0.0, 1.0)]
mpl_map_palette['hk'] = mpl_map_palette['l']
mpl_map_palette['xy'] = mpl_map_palette['z']
mpl_map_palette['kl'] = mpl_map_palette['h']
mpl_map_palette['yz'] = mpl_map_palette['x']
mpl_map_palette['hl'] = mpl_map_palette['k']
mpl_map_palette['xz'] = mpl_map_palette['y']
