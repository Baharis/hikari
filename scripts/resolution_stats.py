# ~~~~~~~~~~~~~~~~~~~~~ IMPORT STATEMENTS - DO NOT CHANGE ~~~~~~~~~~~~~~~~~~~~ #
from kesshou.dataframes.hkl import HklFrame
import numpy as np
import copy

# ~~~~~~~~~~~~~~~~~~~~ VARIABLES - CHANGE ONLY VALUES HERE ~~~~~~~~~~~~~~~~~~~ #
# Unit Cell (in Angstrom in degrees)
unit_cell_a = 13.5495
unit_cell_b = 11.4153
unit_cell_c = 21.4555
unit_cell_al = 90.0
unit_cell_be = 90.0
unit_cell_ga = 90.0

# Input hkl file details
input_hkl_path = '/home/dtchon/git/kesshou/test_data/' \
                 'AM_c1_v1_abs_m_resymmetrified_mmm.hkl'
input_hkl_format = 4
input_hkl_wavelength = 'AgKa'


# ~~~~~~~~~~~~~~~~~~~~~~~~ SCRIPT CODE - DO NOT CHANGE ~~~~~~~~~~~~~~~~~~~~~~~ #

# at this point, only completeness(resolution) is calculated
# check eof for the previous code

# Load desired hkl file
p = HklFrame()
p.crystal.edit_cell(a=unit_cell_a, b=unit_cell_b, c=unit_cell_c,
                    al=unit_cell_al, be=unit_cell_be, ga=unit_cell_ga)
p.edit_wavelength(input_hkl_wavelength)
#q = copy.deepcopy(p)
p.read(input_hkl_path, input_hkl_format)
p.drop_zero()
p.reduce()
p.place()

# generate a big ball for comparison
#q.generate_ball(radius=2.0/q.meta['wavelength'])
#q.drop_zero()
#q.place()

# find an optimal range limit
max_rad = 2.0/p.meta['wavelength'] + 0.1 - (2.0/p.meta['wavelength']) % 0.1
max_rad = 2.0 #if you have desired limit, type it

print('res[A^-1]', 'experiment', 'theory')
for radius in np.arange(max_rad, 0, -0.05):
    p.trim(limit=radius)
    #q.trim(limit=radius)
    print(str(radius/2), len(p))#, len(q))

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ END OF FILE ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ #
#     def calculate_statistics(self, radius):
#         # TODO forgot to consider radius :)
#         pos_reflections = self.count_points_in_sphere(radius)
#         q = copy.deepcopy(self)
#         q.drop_zero()
#         q.data.sort_values(['r'], inplace=True)
#         all_reflections = 0
#         for index, row in q.data.iterrows():
#             if row['r'] <= radius:
#                 all_reflections += row['m']
#         q.reduce()
#         ind_reflections = 0
#         obs_reflections = 0
#         for index, row in q.data.iterrows():
#             if row['r'] <= radius:
#                 ind_reflections += 1
#                 try:
#                     if row['I']/row['si'] > 2:
#                         obs_reflections += row['m']
#                 except KeyError:
#                     if row['F']/row['si'] > 2:
#                         obs_reflections += row['m']
#         del(copy)
#         print('RADIUS: ' + str(radius))
#         print('all: ' + str(all_reflections))
#         print('independent: ' + str(ind_reflections))
#         print('observed: ' + str(obs_reflections))
#         print('possible: ' + str(pos_reflections))
#         print('completeness: ' + str(ind_reflections/pos_reflections))
#         print('redundancy: ' + str(all_reflections/ind_reflections) + '\n')