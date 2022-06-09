"""
This module contains all major scripts of hikari.
A script is a high-level hikari function, which utilises
a lot of low-level objects such as hikari :mod:`hikari.dataframe`
or :mod:`hikari.symmetry` in order to provide more sophisticated functionality.

Since HklFrame is the most developed :mod:`hikari.dataframe` of hikari,
majority of existing scripts revolve about modification of the .hkl file.
It should be mentioned that writing a hikari script does require elemental
knowledge of other objects within the package,
but should not be a problem after a short study of code.

Users are cordially invited to propose their own scripts or script ideas.
"""
from .hkl_potency import potency_map, potency_vs_dac_opening_angle, \
    potency_violin_plot, dac_potency_around_axis
from .hkl_completeness import completeness_statistics, dac_statistics, \
    simulate_dac, reformat_hkl
from .compare_adps import calculate_similarity_indices
from .r1_map import r1_map
