"""
The module containing relatively simple functions,
which are not necessarily strictly connected
to the actual contents of the whole library,
but do significantly improve clarity of the code present in other modules.
"""

from .chem_tools import chemical_elements
from .math_tools import angle2rad, fibonacci_sphere, is2n, is3n, is4n, is6n
from .list_tools import cubespace, rescale_list_to_range, rescale_list_to_other
from .os_tools import make_absolute_path, home_directory
from .palettes import gnuplot_cplt_map_palette, mpl_cplt_map_palette

from pathlib import Path
cplt_map_template = open(Path(__file__).parent/'cplt_map_template.gnu').read()
