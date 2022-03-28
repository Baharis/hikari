"""
The module containing relatively simple functions,
which are not necessarily strictly connected
to the actual contents of the whole library,
but do significantly improve clarity of the code present in other modules.
"""

from .chem_tools import chemical_elements, split_atom_label
from .math_tools import angle2rad, fibonacci_sphere, \
    rotation_from, rotation_around
from .list_tools import cubespace, find_best, rescale_list_to_range,\
    rescale_list_to_other
from .os_tools import make_abspath
from .palettes import gnuplot_map_palette, mpl_map_palette
