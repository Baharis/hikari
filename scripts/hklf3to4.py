# ~~~~~~~~~~~~~~~~~~~~~ IMPORT STATEMENTS - DO NOT CHANGE ~~~~~~~~~~~~~~~~~~~~ #
from kesshou.dataframes.hkl import HklFrame

# ~~~~~~~~~~~~~~~~~~~~ VARIABLES - CHANGE ONLY VALUES HERE ~~~~~~~~~~~~~~~~~~~ #

# Input details
input_hkl_path = '/home/dtchon/git/kesshou/test_data/AM_c1_v1_abs_m.hkl'

# ~~~~~~~~~~~~~~~~~~~~~~~~ SCRIPT CODE - DO NOT CHANGE ~~~~~~~~~~~~~~~~~~~~~~~ #

# Prepare HklFrame objects
p = HklFrame()
p.read(input_hkl_path, 3)
max_f = p.data['F'].max()
if max_f > 9999:
    p.rescale_f(9999 / max_f)
p.calculate_intensity_from_structure_factor()
p.write(input_hkl_path + 'f4', 4, columns_separator=False)

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ END OF FILE ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ #
