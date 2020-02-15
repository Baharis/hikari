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
