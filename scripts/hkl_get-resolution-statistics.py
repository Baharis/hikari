# ~~~~~~~~~~~~~~~~~~~~~ IMPORT STATEMENTS - DO NOT CHANGE ~~~~~~~~~~~~~~~~~~~~ #
from kesshou.dataframes.hkl import HklFrame
from kesshou.symmetry.pointgroup import *

# ~~~~~~~~~~~~~~~~~~~~ VARIABLES - CHANGE ONLY VALUES HERE ~~~~~~~~~~~~~~~~~~~ #
# Unit Cell (in Angstrom in degrees)
unit_cell_a = 13.5495
unit_cell_b = 11.4153
unit_cell_c = 21.4555
unit_cell_al = 90.0
unit_cell_be = 90.0
unit_cell_ga = 90.0

# Crystal's point group
point_group = PGmmm

# Input hkl file details
input_hkl_path = '/home/dtchon/git/kesshou/test_data/' \
                 'AM_c1_v1_abs_m_resymmetrified_mmm.hkl'
input_hkl_format = 4
input_hkl_wavelength = 'AgKa'


# ~~~~~~~~~~~~~~~~~~~~~~~~ SCRIPT CODE - DO NOT CHANGE ~~~~~~~~~~~~~~~~~~~~~~~ #

# Load desired hkl file
p = HklFrame()
p.crystal.edit_cell(a=unit_cell_a, b=unit_cell_b, c=unit_cell_c,
                    al=unit_cell_al, be=unit_cell_be, ga=unit_cell_ga)
p.edit_wavelength(input_hkl_wavelength)
p.read(input_hkl_path, input_hkl_format)
p.drop_zero()

# Assign reflections their positions and make statistics
p._place()
p.make_stats(point_group=point_group)

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ END OF FILE ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ #
