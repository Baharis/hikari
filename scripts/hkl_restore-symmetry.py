# ~~~~~~~~~~~~~~~~~~~~~ IMPORT STATEMENTS - DO NOT CHANGE ~~~~~~~~~~~~~~~~~~~~ #
from kesshou.dataframes.hkl import HklFrame
from kesshou.symmetry.pointgroup import *

# ~~~~~~~~~~~~~~~~~~~~ VARIABLES - CHANGE ONLY VALUES HERE ~~~~~~~~~~~~~~~~~~~ #

# Point group used to transform the data
point_group = PG2pm
merge_at_the_end = False

# Input details
input_hkl_path = '/home/dtchon/git/kesshou/test_data/AM_c1_v1_abs_m.hkl'
input_hkl_format = 4

# Output details
output_hkl_path = '/home/dtchon/git/kesshou/test_data/' \
                  'AM_c1_v1_abs_m_resymmetrified_-1.hkl'
output_hkl_format = 4

# ~~~~~~~~~~~~~~~~~~~~~~~~ SCRIPT CODE - DO NOT CHANGE ~~~~~~~~~~~~~~~~~~~~~~~ #

# Prepare HklFrame objects
p = HklFrame()
p.read(input_hkl_path, input_hkl_format)

# Apply all symmetry elements to the hkl
p.resymmetrify(operations=point_group.operations)
if merge_at_the_end is True:
    p.merge()

# Write the output hkl file
p.write(output_hkl_path, output_hkl_format)

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ END OF FILE ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ #
