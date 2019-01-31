# ~~~~~~~~~~~~~~~~~~~~~ IMPORT STATEMENTS - DO NOT CHANGE ~~~~~~~~~~~~~~~~~~~~ #
from kesshou.dataframes.hkl import HklFrame
import numpy as np

# ~~~~~~~~~~~~~~~~~~~~ VARIABLES - CHANGE ONLY VALUES HERE ~~~~~~~~~~~~~~~~~~~ #
# Unit Cell (in Angstrom in degrees)
unit_cell_a = 13.5495
unit_cell_b = 11.4153
unit_cell_c = 21.4555
unit_cell_al = 90.0
unit_cell_be = 90.0
unit_cell_ga = 90.0

# Crystal orientation matrix from .cif file
UB_11 = -0.0110676000
UB_12 = -0.0348772000
UB_13 = -0.0169147000
UB_21 = 0.0281042000
UB_22 = -0.0302304000
UB_23 = 0.0102881000
UB_31 = -0.0281232000
UB_32 = -0.0165421000
UB_33 = 0.0169609000

# Opening angle in degrees
pressure_cell_oa = 45

# Input hkl file details
input_hkl_path = '/home/dtchon/git/kesshou/test_data/' \
                 'AM_c1_v1_abs_m_resymmetrified_1.hkl'
input_hkl_format = 4
input_hkl_wavelength = 'AgKa'


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
q.crystal.edit_cell(a=unit_cell_a, b=unit_cell_b, c=unit_cell_c,
                    al=unit_cell_al, be=unit_cell_be, ga=unit_cell_ga)
q.edit_wavelength(input_hkl_wavelength)
q.generate_ball(radius=2/q.meta['wavelength'])
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
