# ~~~~~~~~~~~~~~~~~~~~~ IMPORT STATEMENTS - DO NOT CHANGE ~~~~~~~~~~~~~~~~~~~~ #
from kesshou.dataframes.hkl import HklFrame
import numpy as np
from copy import deepcopy

# ~~~~~~~~~~~~~~~~~~~~ VARIABLES - CHANGE ONLY VALUES HERE ~~~~~~~~~~~~~~~~~~~ #

# Unit Cell (in Angstrom in degrees)
unit_cell_a = 6.6139
unit_cell_b = 8.1151
unit_cell_c = 10.272
unit_cell_al = 109.439
unit_cell_be = 99.563
unit_cell_ga = 100.637

# Crystal orientation matrix from .cif file
UB_11 = -0.0483900000
UB_12 = 0.0237766000
UB_13 = -0.0253660000
UB_21 = 0.0588330000
UB_22 = 0.0310240000
UB_23 = -0.0129997000
UB_31 = -0.0119742000
UB_32 = -0.0528324000
UB_33 = -0.0431912000

# OR DAC perpendicular vector, if better suited
v1 = 1/np.sqrt(2)
v2 = 1/np.sqrt(2)
v3 = 0
use_vector_instead_of_orientation_matrix = False

# Opening angle in degrees
pressure_cell_oa = [40, ]

# Resolution limit as sin(th)/la, A-1
reslim = 0.8

# Input details
input_hkl_path = '/home/dtchon/x/HiPHAR/PA/hkl_preparation/0kbar_15geom_experimental/plate_2.hkl'
input_hkl_format = 5
input_hkl_wavelength = 0.4859 #'AgKa'

# Output details
output_name = '0kbar_15geom_experimental_x1'
output_directory = '/home/dtchon/x/HiPHAR/PA/hkl_preparation/0kbar_15geom_experimental/'
output_hkl_format = 4


# ~~~~~~~~~~~~~~~~~~~~~~~~ SCRIPT CODE - DO NOT CHANGE ~~~~~~~~~~~~~~~~~~~~~~~ #

# Prepare HklFrame object
p = HklFrame()
p.read(input_hkl_path, input_hkl_format)
p.crystal.edit_cell(a=unit_cell_a, b=unit_cell_b, c=unit_cell_c,
                    al=unit_cell_al, be=unit_cell_be, ga=unit_cell_ga)
p.edit_wavelength(input_hkl_wavelength)
if not use_vector_instead_of_orientation_matrix:
    p.crystal.orient_matrix = np.array(((UB_11, UB_12, UB_13),
                                        (UB_21, UB_22, UB_23),
                                        (UB_31, UB_32, UB_33)))
p.drop_zero()
p.place()

# Prepare list of interesting projections
projections = (('h', 'k', 0), ('h', 0, 'l'), (0, 'k', 'l'),
               ('h', 'k', 1), ('h', 1, 'l'), (1, 'k', 'l'))

# Trim data to resolution if necesary
p.trim(reslim)


# Draw projections before dac operation
q = deepcopy(p)
q.reduce()
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
    if use_vector_instead_of_orientation_matrix:
        vector = np.array((v1, v2, v3))
        q.dac(opening_angle=oa, vector=vector)
    else:
        q.dac(opening_angle=oa)
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
