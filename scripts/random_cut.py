# ~~~~~~~~~~~~~~~~~~~~~ IMPORT STATEMENTS - DO NOT CHANGE ~~~~~~~~~~~~~~~~~~~~ #
from kesshou.dataframes.hkl import HklFrame
from copy import deepcopy

# ~~~~~~~~~~~~~~~~~~~~ VARIABLES - CHANGE ONLY VALUES HERE ~~~~~~~~~~~~~~~~~~~ #

# Target completeness of the output file
# (the algorithm assumes that input is 100% complete merged hkl)
target_completeness = [0.10]

# Input hkl file details
input_hkl_path = '/home/dtchon/x/HiPHAR/RFpirazB/hkl_preparation/full_merged/sortav.hkl'
input_hkl_format = 4

# Output hkl file details
output_name = 'RFpirazB'
output_directory = '/home/dtchon/x/HiPHAR/RFpirazB/hkl_preparation/'
output_hkl_format = 4

# ~~~~~~~~~~~~~~~~~~~~~~~~ SCRIPT CODE - DO NOT CHANGE ~~~~~~~~~~~~~~~~~~~~~~~ #

# Prepare HklFrame object
p = HklFrame()
p.read(input_hkl_path, input_hkl_format)

# Cut for subsequent target completeness values
try:
    _ = iter(target_completeness)
except TypeError:
    target_completeness = list(target_completeness)
for target in target_completeness:
    q = deepcopy(p)
    output_hkl_path = output_directory + output_name +\
                      '_cplt' + str(target)[0:6] + '.hkl'
    q.drop_zero()
    q.thin_out(target_completeness=target)
    q.write(output_hkl_path, output_hkl_format)

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ END OF FILE ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ #
