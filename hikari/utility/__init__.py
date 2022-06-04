"""
The module containing relatively simple functions,
which are not necessarily strictly connected
to the actual contents of the whole library,
but do significantly improve clarity of the code present in other modules.
"""

from .certain_float import cfloat
from .chem_tools import chemical_elements, split_atom_label
from .interval import Interval
from .math_tools import angle2rad, cart2sph, sph2cart, det3x3, \
    fibonacci_sphere, rotation_from, rotation_around, weighted_quantile
from .list_tools import cubespace, find_best, rescale_list_to_range,\
    rescale_list_to_other
from .os_tools import make_abspath
from .palettes import gnuplot_map_palette, mpl_map_palette
from .artists import artist_factory
