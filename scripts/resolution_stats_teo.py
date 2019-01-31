# ~~~~~~~~~~~~~~~~~~~~~ IMPORT STATEMENTS - DO NOT CHANGE ~~~~~~~~~~~~~~~~~~~~ #
from kesshou.dataframes.hkl import HklFrame
import numpy as np
import copy

# ~~~~~~~~~~~~~~~~~~~~ VARIABLES - CHANGE ONLY VALUES HERE ~~~~~~~~~~~~~~~~~~~ #
# Unit Cell (in Angstrom in degrees)
unit_cell_a = 10.0
unit_cell_b = 10.0
unit_cell_c = 10.0
unit_cell_al = 90.0
unit_cell_be = 90.0
unit_cell_ga = 90.0

# Crystal orientation matrix from .cif file
vector = [1.0, 1.0, 1.0]

# Opening angle in degrees
pressure_cell_oa = 40
hkl_wavelength = 0.30

# Prepare a list of symmetry operations characteristic for a given Laue group:
# remember half of them is unnecessary to retrieve the shape (disc has -1 symm)
# symmetry_operations = []  # LAUE GROUP -1:
# symmetry_operations = [  # LAUE GROUP 2/m:
#    [[1, 0, 0], [0, -1, 0], [0, 0, 1]]          # m_y plane
# ]
symmetry_operations = [  # LAUE GROUP mmm:
   [[-1, 0, 0], [0, 1, 0], [0, 0, 1]],         # m_x plane
   [[1, 0, 0], [0, -1, 0], [0, 0, 1]],         # m_y plane
   [[1, 0, 0], [0, 1, 0], [0, 0, -1]]         # m_z plane
]

# ~~~~~~~~~~~~~~~~~~~~~~~~ SCRIPT CODE - DO NOT CHANGE ~~~~~~~~~~~~~~~~~~~~~~~ #

# generate reference ball
p = HklFrame()
p.crystal.edit_cell(a=unit_cell_a, b=unit_cell_b, c=unit_cell_c,
                    al=unit_cell_al, be=unit_cell_be, ga=unit_cell_ga)
p.edit_wavelength(hkl_wavelength)
p.generate_ball(radius=1./p.meta['wavelength'])
p.drop_zero()
p.reduce()
p.place()

# find an optimal range limit
max_rad = 1.0/p.meta['wavelength'] + 0.1 - (1.0/p.meta['wavelength']) % 0.1

print('res[A^-1]', 'theory', 'dac', 'dac_ressym')
for radius in np.arange(max_rad, 0, -0.05):
    p.trim(limit=radius)
    reflections_ball = len(p)

    # generate not-ressymetrified hkl list
    q = copy.deepcopy(p)
    q.dac(opening_angle=pressure_cell_oa, vector=vector)
    reflections_dac = len(q)

    # ressymetrify
    q.resymmetrify(operations=symmetry_operations)
    reflections_ressym = len(q)

    print(str(radius/2), str(reflections_ball),
          str(reflections_dac), str(reflections_ressym))

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ END OF FILE ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ #
