# ~~~~~~~~~~~~~~~~~~~~~ IMPORT STATEMENTS - DO NOT CHANGE ~~~~~~~~~~~~~~~~~~~~ #
from kesshou.dataframes.hkl import HklFrame
import copy
import numpy as np
import random
from kesshou.symmetry.pointgroup import PGmmm

# ~~~~~~~~~~~~~~~~~~~~ VARIABLES - CHANGE ONLY VALUES HERE ~~~~~~~~~~~~~~~~~~~ #
# This script verifies the values presented in paper of Casati et. al.
# "Exploring charge density analysis [...]", Table 1

# Unit Cell (in Angstrom in degrees)
unit_cell_a = 10.0
unit_cell_b = 10.0
unit_cell_c = 10.0
unit_cell_al = 90.0
unit_cell_be = 90.0
unit_cell_ga = 90.0

# Random seed for selecting random orientation vectors
random.seed(1337)

# Prepare a list of symmetry operations characteristic for a given Laue group:
# remember half of them is unnecessary to retrieve the shape (disc has -1 symm)
point_group = PGmmm

# Multistats precision, how many random vectors should be checked
precision = 1000

# Checked parameters
pressure_cell_oa = 55
hkl_wavelength = 0.45
resolution = 0.40  # checked 0.50 or 0.83

# Output file details
output_txt_path = '/home/dtchon/git/kesshou/test_data/' \
                  'orthorhombic_oa55_la45_res40.txt'

# ~~~~~~~~~~~~~~~~~~~~~~~~ SCRIPT CODE - DO NOT CHANGE ~~~~~~~~~~~~~~~~~~~~~~~ #


# Generate a list of uniformly distributed vectors
# using Fibbonaci sphere alghoritm https://stackoverflow.com/questions/9600801/
def fibonacci_sphere(samples=1, randomize=True):
    rnd = random.random() * samples if randomize else 1.
    points = []
    offset = 2./samples
    increment = np.pi * (3. - np.sqrt(5.))
    for i in range(samples):
        y = ((i * offset) - 1) + (offset / 2)
        r = np.sqrt(1 - pow(y, 2))
        phi = ((i + rnd) % samples) * increment
        x = np.cos(phi) * r
        z = np.sin(phi) * r
        points.append([x, y, z])
    return points


# Open a file for writing the output
output_file = open(output_txt_path, 'w')

# Generate and measure reference ball
q = HklFrame()
q.crystal.edit_cell(a=unit_cell_a, b=unit_cell_b, c=unit_cell_c,
                    al=unit_cell_al, be=unit_cell_be, ga=unit_cell_ga)
q.edit_wavelength(hkl_wavelength)
q.make_ball(radius=1. / resolution)
q.reduce()
q.drop_zero()
total_reflections = len(q)
output_file.write('total_reflections: ' + str(total_reflections))

# Prepare a template sphere for dac-cutting
oa_in_radians = pressure_cell_oa * np.pi / 180.
observed_radius = min((2.001*np.sin(oa_in_radians)/q.meta['wavelength'],
                       1./resolution))
q.make_ball(radius=observed_radius)
q.place()
q.drop_zero()

# Generate a list of points on sphere and for each...
points = fibonacci_sphere(samples=precision)
for point in points:
    p = copy.deepcopy(q)
    p.dac(opening_angle=pressure_cell_oa, vector=point)
    p.resymmetrify(point_group.hp_disc_transforming_symm_ops)
    counted_reflections = len(p)
    output_file.write('\n' + str(point) + ' ' + str(counted_reflections))
    del p

output_file.close()

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ END OF FILE ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ #
