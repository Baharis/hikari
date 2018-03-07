# ~~~~~~~~~~~~~~~~~~~~~ IMPORT STATEMENTS - DO NOT CHANGE ~~~~~~~~~~~~~~~~~~~~ #
from kesshou.dataframes.hkl import HklFrame
import numpy as np

# ~~~~~~~~~~~~~~~~~~~~ VARIABLES - CHANGE ONLY VALUES HERE ~~~~~~~~~~~~~~~~~~~ #
# Unit Cell (in Angstrom in degrees)
unit_cell_a = 10.0
unit_cell_b = 10.0
unit_cell_c = 10.0
unit_cell_al = 90.0
unit_cell_be = 90.0
unit_cell_ga = 90.0

# Crystal orientation matrix from .cif file
UB_11 = -0.0117244000
UB_12 = -0.0419390000
UB_13 = -0.0270692000
UB_21 = -0.0571303000
UB_22 = -0.0081278000
UB_23 = 0.0181937000
UB_31 = -0.0264406000
UB_32 = 0.0361560000
UB_33 = -0.0272830000

# Opening angle in degrees
pressure_cell_oa = 35

# Input hkl file details
input_hkl_path = '/home/dtchon/git/kesshou/test_data/c1_p1_dt.hkl'
input_hkl_format = 4
input_hkl_wavelength = 'CuKa'


# ~~~~~~~~~~~~~~~~~~~~~~~~ SCRIPT CODE - DO NOT CHANGE ~~~~~~~~~~~~~~~~~~~~~~~ #


p = HklFrame()
p.crystal.edit_cell(a=unit_cell_a, b=unit_cell_b, c=unit_cell_c,
                    al=unit_cell_al, be=unit_cell_be, ga=unit_cell_ga)
p.edit_wavelength(input_hkl_wavelength)
p.read(input_hkl_path, input_hkl_format)
p.crystal.orient_matrix = np.array(((UB_11, UB_12, UB_13),
                                    (UB_21, UB_22, UB_23),
                                    (UB_31, UB_32, UB_33)))
p.drop_zero()
p.reduce()
p.place()

q = HklFrame()
p.crystal.edit_cell(a=unit_cell_a, b=unit_cell_b, c=unit_cell_c,
                    al=unit_cell_al, be=unit_cell_be, ga=unit_cell_ga)
q.edit_wavelength(input_hkl_wavelength)
q.generate_ball(radius=q.meta['wavelength'])
q.crystal.orient_matrix = np.array(((UB_11, UB_12, UB_13),
                                    (UB_21, UB_22, UB_23),
                                    (UB_31, UB_32, UB_33)))
q.drop_zero()
q.place()
print('oa', 'experiment', 'theory')
for angle in range(pressure_cell_oa, 0, -1):
    p.dac(opening_angle=angle)
    q.dac(opening_angle=angle)
    print(angle, len(p), len(q))

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ END OF FILE ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ #
