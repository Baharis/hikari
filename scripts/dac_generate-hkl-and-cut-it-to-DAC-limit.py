# ~~~~~~~~~~~~~~~~~~~~~ IMPORT STATEMENTS - DO NOT CHANGE ~~~~~~~~~~~~~~~~~~~~ #
from kesshou.dataframes.hkl import HklFrame
import numpy as np
from copy import deepcopy

# ~~~~~~~~~~~~~~~~~~~~ VARIABLES - CHANGE ONLY VALUES HERE ~~~~~~~~~~~~~~~~~~~ #
# Unit Cell (in Angstrom in degrees)
unit_cell_a = 4.2815
unit_cell_b = 4.2815
unit_cell_c = 4.2815
unit_cell_al = 90.0000
unit_cell_be = 90.0000
unit_cell_ga = 90.0000

# Crystal orientation matrix from .cif file
UB_11 = 0.0025914000
UB_12 = -0.0132920000
UB_13 = -0.1299297000
UB_21 = 0.0567811000
UB_22 = 0.1172899000
UB_23 = -0.0107307000
UB_31 = 0.1174748000
UB_32 = -0.0564320000
UB_33 = 0.0081958000

# Opening angle in degrees
pressure_cell_oa = [37]

# Input details
input_hkl_wavelength = 'AgKa'

# Output details
output_name = 'Pm3m_AM'
output_directory = '/home/dtchon/x/HP/Pm3m_AM/'
output_hkl_format = 4

# ~~~~~~~~~~~~~~~~~~~~~~~~ SCRIPT CODE - DO NOT CHANGE ~~~~~~~~~~~~~~~~~~~~~~~ #

# Prepare HklFrame object
p = HklFrame()
p.edit_wavelength(input_hkl_wavelength)
p.crystal.edit_cell(a=unit_cell_a, b=unit_cell_b, c=unit_cell_c,
                    al=unit_cell_al, be=unit_cell_be, ga=unit_cell_ga)
p.crystal.orient_matrix = np.array(((UB_11, UB_12, UB_13),
                                    (UB_21, UB_22, UB_23),
                                    (UB_31, UB_32, UB_33)))
p.generate_ball(radius=2/p.meta['wavelength'])
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
