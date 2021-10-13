"""
This module contains all major script of Kesshou.
A script is a high-level hikari function, which utilises
a lot of low-level objects such as hikari :mod:`hikari.dataframe`
or :mod:`hikari.symmetry` in order to provide more sophisticated functionality.

Since HklFrame is the most developed :mod:`hikari.dataframe` of Kesshou,
majority of existing scripts revolve about modification of the .hkl file.
It should be mentioned that writing a hikari script does require elemental
knowledge of other objects within the package,
but should not be a problem after a short study of code.

Users are cordially invited to propose their own scripts or script ideas.
"""
from .hkl import completeness_map
from .hkl import completeness_statistics
from .hkl import dac_point_group_statistics
from .hkl import dac_statistics
from .hkl import simulate_dac
