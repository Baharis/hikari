# ~~~~~~~~~~~~~~~~~~~~~ IMPORT STATEMENTS - DO NOT CHANGE ~~~~~~~~~~~~~~~~~~~~ #
from kesshou.dataframes.hkl import HklFrame

# ~~~~~~~~~~~~~~~~~~~~ VARIABLES - CHANGE ONLY VALUES HERE ~~~~~~~~~~~~~~~~~~~ #

# Input details
input_hkl_path = '/home/dtchon/x/HiPHAR/glycine/IAM_models/opt40/sortav.hkl'
input_hkl_format = 4

# Output details
output_hkl_path = '/home/dtchon/x/HiPHAR/glycine/HAR_models_ISO/opt40/gly_opt40_tonto.hkl'
output_hkl_format = 'tonto'


# ~~~~~~~~~~~~~~~~~~~~~~~~ SCRIPT CODE - DO NOT CHANGE ~~~~~~~~~~~~~~~~~~~~~~~ #

# Prepare HklFrame object
p = HklFrame()
p.read(input_hkl_path, input_hkl_format)
p.write(output_hkl_path, output_hkl_format)


# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ END OF FILE ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ #
