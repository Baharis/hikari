gnuplot_cplt_map_palette = dict()
gnuplot_cplt_map_palette[''] = 'set palette model HSV defined (0 0.667 1 0, 5 0.667 1 1, 15 0.0 1 1, 20 0.0 0 1)' # , 25 0.0 0 0.5) to add gray on the end
gnuplot_cplt_map_palette['h'] = 'set palette model RGB defined (0 1 1 1, 1 1 0 0)'
gnuplot_cplt_map_palette['x'] = 'set palette model RGB defined (0 1 1 1, 1 1 0 0)'
gnuplot_cplt_map_palette['k'] = 'set palette model RGB defined (0 1 1 1, 1 0 1 0)'
gnuplot_cplt_map_palette['y'] = 'set palette model RGB defined (0 1 1 1, 1 0 1 0)'
gnuplot_cplt_map_palette['l'] = 'set palette model RGB defined (0 1 1 1, 1 0 0 1)'
gnuplot_cplt_map_palette['z'] = 'set palette model RGB defined (0 1 1 1, 1 0 0 1)'

mpl_cplt_map_palette = dict()
mpl_cplt_map_palette[''] = [(0.0, 0.0, 0.0),  # Black
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
mpl_cplt_map_palette['h'] = [(1.0, 1.0, 1.0), (1.0, 0.0, 0.0)]
mpl_cplt_map_palette['x'] = [(0.0, 0.0, 0.0), (1.0, 0.0, 0.0)]
mpl_cplt_map_palette['k'] = [(1.0, 1.0, 1.0), (0.0, 1.0, 0.0)]
mpl_cplt_map_palette['y'] = [(0.0, 0.0, 0.0), (0.0, 1.0, 0.0)]
mpl_cplt_map_palette['l'] = [(1.0, 1.0, 1.0), (0.0, 0.0, 1.0)]
mpl_cplt_map_palette['z'] = [(0.0, 0.0, 0.0), (0.0, 0.0, 1.0)]
