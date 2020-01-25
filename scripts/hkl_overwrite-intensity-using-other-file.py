# ~~~~~~~~~~~~~~~~~~~~~ IMPORT STATEMENTS - DO NOT CHANGE ~~~~~~~~~~~~~~~~~~~~ #
from kesshou.dataframes.hkl import HklFrame
from kesshou.symmetry.pointgroup import *

# ~~~~~~~~~~~~~~~~~~~~ VARIABLES - CHANGE ONLY VALUES HERE ~~~~~~~~~~~~~~~~~~~ #

# Prepare a list of symmetry operations characteristic for a given Laue group:
symmetry_operations = PG2pm.hp_disc_transforming_symm_ops

# Input_hkls
input_target_hkl = '/home/dtchon/x/HiPHAR/RFpirazB/hkl_preparation/ref_10geom_experimental/sortav_10kbar_10geom_experimental.hkl'
input_target_hkl_format = 4
input_overwriting_hkl = '/home/dtchon/x/HiPHAR/RFpirazB/hkl_preparation/full_merged/sortav.hkl'
input_overwriting_hkl_format = 4

# Output directory
output_overwritten_hkl = 'overwritten_refonto10.hkl'
output_directory = '/home/dtchon/x/HiPHAR/RFpirazB/hkl_preparation/ref_10geom_experimental/'
output_hkl_format = 4

# ~~~~~~~~~~~~~~~~~~~~~~~~ SCRIPT CODE - DO NOT CHANGE ~~~~~~~~~~~~~~~~~~~~~~~ #

# Prepare HklFrame objects
p, q = HklFrame(), HklFrame()
p.read(input_target_hkl, input_target_hkl_format)
q.read(input_overwriting_hkl, input_overwriting_hkl_format)
q.resymmetrify(symmetry_operations)
p.extinct('000')

# Perform the overwriting and prepare new file
p.overwrite_values(q, keys=['I', 'si'])
q.write(output_directory + output_overwritten_hkl, output_hkl_format)

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ END OF FILE ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ #
# might not work after last changes - to be checked