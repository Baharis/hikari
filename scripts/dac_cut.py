# ~~~~~~~~~~~~~~~~~~~~~ IMPORT STATEMENTS - DO NOT CHANGE ~~~~~~~~~~~~~~~~~~~~ #
from kesshou.dataframes.hkl import HklFrame
import numpy as np
from copy import deepcopy

# ~~~~~~~~~~~~~~~~~~~~ VARIABLES - CHANGE ONLY VALUES HERE ~~~~~~~~~~~~~~~~~~~ #
# Unit Cell (in Angstrom in degrees)
unit_cell_a = 9.133980
unit_cell_b = 8.809929
unit_cell_c = 21.418243
unit_cell_al = 90.0000
unit_cell_be = 93.0499
unit_cell_ga = 90.0000

# Glicyna
# unit_cell_a = 5.085678
# unit_cell_b = 11.803994
# unit_cell_c = 5.460576
# unit_cell_al = 90.0000
# unit_cell_be = 111.9796
# unit_cell_ga = 90.0000

# Crystal orientation matrix from .cif file
# UB_11 = 1.0
# UB_12 = 0.0
# UB_13 = 0.0
# UB_21 = 0.0
# UB_22 = 1.0
# UB_23 = 0.0
# UB_31 = 0.0
# UB_32 = 0.0
# UB_33 = 1.0
# OR perpendicular vector, if better suited
v1 = 1/np.sqrt(2)
v2 = 1/np.sqrt(2)
v3 = 0

# Opening angle in degrees
pressure_cell_oa = [35, 40]

# Input details
input_hkl_path = '/home/dtchon/x/HiPHAR/RFpirazB/hkl_preparation/full_unmerged/RFpirazB_full_unmerged.hkl'
input_hkl_format = 4
input_hkl_wavelength = 'AgKa'

# Output details
output_name = 'RFpirazB'
output_directory = '/home/dtchon/x/HiPHAR/RFpirazB/hkl_preparation/'
output_hkl_format = 4


# ~~~~~~~~~~~~~~~~~~~~~~~~ SCRIPT CODE - DO NOT CHANGE ~~~~~~~~~~~~~~~~~~~~~~~ #

# Prepare HklFrame object
p = HklFrame()
p.read(input_hkl_path, input_hkl_format)
p.crystal.edit_cell(a=unit_cell_a, b=unit_cell_b, c=unit_cell_c,
                    al=unit_cell_al, be=unit_cell_be, ga=unit_cell_ga)
p.edit_wavelength(input_hkl_wavelength)
# p.crystal.orient_matrix = np.array(((UB_11, UB_12, UB_13),
#                                     (UB_21, UB_22, UB_23),
#                                     (UB_31, UB_32, UB_33)))
p.drop_zero()
p.place()

# Prepare list of interesting projections
projections = (('h', 'k', 0), ('h', 0, 'l'), (0, 'k', 'l'),
               ('h', 'k', 1), ('h', 1, 'l'), (1, 'k', 'l'))


# Draw projections before dac operation
q = deepcopy(p)
q.reduce()
q.drop_zero()
for projection in projections:
    output_png_path = output_directory + output_name + '_full_' + \
                      str(projection[0]) + \
                      str(projection[1]) + \
                      str(projection[2])
    q.draw(colored='m', projection=projection, savepath=output_png_path)

# Cut and draw projections for subsequent DAC opening angles
try:
    _ = iter(pressure_cell_oa)
except TypeError:
    pressure_cell_oa = list(pressure_cell_oa)
for oa in pressure_cell_oa:
    q = deepcopy(p)
    vector = np.array((v1, v2, v3))
    q.dac(opening_angle=oa, vector=vector)
    output_hkl_path = output_directory + output_name +\
                      '_oa' + str(oa)[0:6] + '.hkl'
    q.write(hkl_path=output_hkl_path, hkl_format=output_hkl_format)
    q.reduce()
    q.drop_zero()
    for projection in projections:
        output_png_path = output_directory + output_name + \
                          '_oa' + str(oa)[0:6] + '_' + \
                          str(projection[0]) + \
                          str(projection[1]) + \
                          str(projection[2])
        q.draw(colored='m', projection=projection, savepath=output_png_path)


# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ END OF FILE ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ #
