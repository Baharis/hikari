# ~~~~~~~~~~~~~~~~~~~~~ IMPORT STATEMENTS - DO NOT CHANGE ~~~~~~~~~~~~~~~~~~~~ #
from kesshou.dataframes.hkl import HklFrame
import numpy as np
import copy
from kesshou.symmetry.pointgroup import PGm_3m

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
UB_33 =  0.0081958000

# Prepare a list of symmetry operations characteristic for a given Laue group:
# remember half of them is unnecessary to retrieve the shape (disc has -1 symm)
symmetry_operations = PGm_3m.hp_disc_transforming_symm_ops

# Opening angle in degrees
pressure_cell_oa = 37

# Input hkl file details
input_hkl_path = '/home/dtchon/x/HP/Pm3m_AM/Pm3m_AM_oa37.hkl'
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
p.dac(opening_angle=pressure_cell_oa)


q = HklFrame()
q.crystal.edit_cell(a=unit_cell_a, b=unit_cell_b, c=unit_cell_c,
                    al=unit_cell_al, be=unit_cell_be, ga=unit_cell_ga)
q.edit_wavelength(input_hkl_wavelength)
#q.generate_ball(radius=2/q.meta['wavelength'])
q.generate_ball(radius=max(p.data['r']))
q.crystal.orient_matrix = np.array(((UB_11, UB_12, UB_13),
                                    (UB_21, UB_22, UB_23),
                                    (UB_31, UB_32, UB_33)))
q.drop_zero()
q.place()


# print('oa', 'experiment', 'theory')
# for angle in range(pressure_cell_oa, 0, -1):
#     p.dac(opening_angle=angle)
#     q.dac(opening_angle=angle)
#     r = copy.deepcopy(p).resymmetrify()
#     print(angle, len(p), len(r), len(q))

for resolution in np.arange(0.5, 1.5, 0.05):
    rad = 1/resolution
    p.trim(rad)
    q.trim(rad)
    r = copy.deepcopy(p)
    r.resymmetrify(symmetry_operations)
    print(resolution, len(p), len(r), len(q))

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ END OF FILE ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ #
