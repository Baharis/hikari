# ~~~~~~~~~~~~~~~~~~~~~ IMPORT STATEMENTS - DO NOT CHANGE ~~~~~~~~~~~~~~~~~~~~ #
from kesshou.dataframes.hkl import HklFrame
import numpy as np

# ~~~~~~~~~~~~~~~~~~~~ VARIABLES - CHANGE ONLY VALUES HERE ~~~~~~~~~~~~~~~~~~~ #

# Target completeness of the output file
# (the algorithm assumes that input is 100% complete merged hkl)
target_completeness = 0.50

# Input hkl file details
input_hkl_path = '/home/dtchon/git/kesshou/test_data/c1_p1_dt.hkl'
input_hkl_format = 4
input_hkl_wavelength = 'CuKa'

# Output hkl file details
output_hkl_path = '/home/dtchon/git/kesshou/test_data/output.hkl'
output_hkl_format = 4

# ~~~~~~~~~~~~~~~~~~~~~~~~ SCRIPT CODE - DO NOT CHANGE ~~~~~~~~~~~~~~~~~~~~~~~ #

p = HklFrame()
p.read(input_hkl_path, input_hkl_format)
p.thin_out(target_completeness=target_completeness)
p.write(output_hkl_path, output_hkl_format)

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ END OF FILE ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ #
