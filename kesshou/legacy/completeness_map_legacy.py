# ~~~~~~~~~~~~~~~~~~~~~ IMPORT STATEMENTS - DO NOT CHANGE ~~~~~~~~~~~~~~~~~~~~ #
from kesshou.dataframes.hkl import HklFrame
import numpy as np
import copy
from datetime import datetime as dt
from kesshou.symmetry.pointgroup_old import PG4pmmm

# ~~~~~~~~~~~~~~~~~~~~ VARIABLES - CHANGE ONLY VALUES HERE ~~~~~~~~~~~~~~~~~~~ #
# Unit Cell (in Angstrom in degrees)
unit_cell_a = 10.0
unit_cell_b = 10.0
unit_cell_c = 10.0
unit_cell_al = 90.0
unit_cell_be = 90.0
unit_cell_ga = 90.0

# Vector orientation matrix resolution
# fractional coordinates in reciprocal space
primary_orientation = [0, 0, 1]
secondary_orientation = [1, 0, 0]

phi_from = 1    # phi is rotation of DAC normal vector around primary vector,
                # where 0 corresponds to primary-secondary plane in recip. space
phi_to = 45    # "from" and "to" describe investigated range of phi
phi_step = 2   # "step" describes degrees step b/ phi values;approximates to int

psi_from = 1    # psi is rotation away from initial vector towards primary
                # vector on the primary-secondary plane
psi_to = 90    # "from" and "to" describe investigated range of phi
psi_step = 2   # "step" describes degrees step b/ phi values;approximates to int


# Experimental details in angstrom and degrees
pressure_cell_oa = 35
hkl_wavelength = 0.71037
resolution_cutoff = 0.8

# Prepare a list of symmetry operations characteristic for a given Laue group:
# remember half of them is unnecessary to retrieve the shape (disc has -1 symm)
point_group = PG4pmmm


# Output directory
output_name = 'laue4pmmm_oa35_laMoKa_sparse'
output_directory = '/home/dtchon/git/kesshou/test_data/'

# ~~~~~~~~~~~~~~~~~~~~~~~~ SCRIPT CODE - DO NOT CHANGE ~~~~~~~~~~~~~~~~~~~~~~~ #

time_start = dt.now()
print('START TIME: ', str(time_start))
# prepare output file
log_path = output_directory+output_name+'.log'
log = open(log_path, 'a+', buffering=1)
dat_path = output_directory+output_name+'.dat'
dat = open(dat_path, 'a+', buffering=1)
log.write(20 * '-' + ' B A S I C   I N F O ' + 20 * '-' + '\n')
log.write('KESSHOU ORIENTATION-DEPENDENT COMPLETENESS CALCULATIONS FOR DAC\n')
log.write('> start time: ' + str(time_start) + '\n')
log.write('> log file path: ' + str(log_path) + '\n')
log.write('> dat file path: ' + str(dat_path) + '\n')
log.write('> pressure cell opening angle: ' + str(pressure_cell_oa) + '\n')
log.write('> unit cell a: ' + str(unit_cell_a) + '\n')
log.write('> unit cell b: ' + str(unit_cell_b) + '\n')
log.write('> unit cell c: ' + str(unit_cell_c) + '\n')
log.write('> unit cell al: ' + str(unit_cell_al) + '\n')
log.write('> unit cell be: ' + str(unit_cell_be) + '\n')
log.write('> unit cell ga: ' + str(unit_cell_ga) + '\n')
log.write('> primary orientation (hkl): ' + str(primary_orientation) + '\n')
log.write('> secondary orientation (hkl): ' + str(secondary_orientation) + '\n')
log.write('> symmetry operations: ' + '\n')
for operation in point_group.chiral_operations:
    log.write('> ' + str(operation) + '\n')

# generate reference ball
p = HklFrame()
p.crystal.edit_cell(a=unit_cell_a, b=unit_cell_b, c=unit_cell_c,
                    al=unit_cell_al, be=unit_cell_be, ga=unit_cell_ga)
p.edit_wavelength(hkl_wavelength)
ball_radius = 2./p.meta['wavelength']
p.fill(radius=ball_radius)
ball_trim_range = 1./resolution_cutoff
p.drop_zero()
p.merge()
p.place()
p.trim(limit=ball_trim_range)
total_reflections = p.table.shape[0]
log.write(20 * '-' + ' R A D I A T I O N ' + 20 * '-' + '\n')
log.write('> radiation wavelength (A): ' + str(p.meta['wavelength']) + '\n')
log.write('< reference ball full radius (1/A): ' + str(ball_radius) + '\n')
log.write('> resolution cutoff (A): ' + str(resolution_cutoff) + '\n')
log.write('> reference ball cutoff (1/A): ' + str(ball_trim_range) + '\n')
log.write('< reference ball reflections: ' + str(total_reflections) + '\n')

# define primary and secondary vectors in reciprocal space
primary_v = np.array(primary_orientation[0] * p.crystal.a_w +
                     primary_orientation[1] * p.crystal.b_w +
                     primary_orientation[2] * p.crystal.c_w)
secondary_v = np.array(secondary_orientation[0] * p.crystal.a_w +
                       secondary_orientation[1] * p.crystal.b_w +
                       secondary_orientation[2] * p.crystal.c_w)

# normalise them and correct secondary vector to be perpendicular to primary
primary_v /= np.linalg.norm(primary_v)
secondary_v /= np.linalg.norm(secondary_v)
secondary_v = secondary_v - np.dot(primary_v, secondary_v) * primary_v
secondary_v /= np.linalg.norm(secondary_v)
tertiary_v = np.cross(primary_v, secondary_v)
log.write(20 * '-' + ' O R I E N T A T I O N ' + 20 * '-' + '\n')
log.write('< primary vector (cartesian*):' + str(primary_v) + '\n')
log.write('< secondary vector (cartesian*):' + str(secondary_v) + '\n')
log.write('< tertiary vector (cartesian*):' + str(tertiary_v) + '\n')

# define orientation points to be investigated
phi_range = np.arange(phi_from, phi_to, phi_step)
psi_range = np.arange(psi_from, psi_to, psi_step)
phi_psi_matrix = np.zeros(shape=(len(phi_range)+1, len(psi_range)+1))
phi_psi_matrix[1:, 0] = phi_range
phi_psi_matrix[0, 1:] = psi_range
orientations_to_calculate = len(phi_range) * len(psi_range)
log.write(20 * '-' + ' R E S O L U T I O N ' + 20 * '-' + '\n')
log.write('> phi low limit:' + str(phi_from) + '\n')
log.write('> phi high limit:' + str(phi_to) + '\n')
log.write('> phi step size:' + str(phi_step) + '\n')
log.write('< phi value list:' + str(phi_range) + '\n')
log.write('> psi low limit:' + str(psi_from) + '\n')
log.write('> psi high limit:' + str(psi_to) + '\n')
log.write('> psi step size:' + str(psi_step) + '\n')
log.write('< psi value list:' + str(psi_range) + '\n')
log.write('< no of points to evaluate:' + str(orientations_to_calculate) + '\n')

# calculate completeness for each investigated phi/psi pair
log.write(20 * '-' + ' R E F L E C T I O N S ' + 20 * '-' + '\n')
log.write('| WRITING DAT FILE' + '\n')
dat.write('#    phi     psi    cplt  reflns' + '\n')
time_loop = dt.now()
orientation_no = 0
for phi_index, phi in enumerate(phi_range):
    for psi_index, psi in enumerate(psi_range):
        phi_r, psi_r = np.deg2rad(phi), np.deg2rad(psi)
        orientation_no += 1
        # rotate the DAC orientation vector psi around primary-secondary plane
        DAC_v = primary_v * np.cos(psi_r) + secondary_v * np.sin(psi_r)
        # rotate the DAC orientation vector phi around the primary vector
        DAC_v_parallel = np.dot(DAC_v, primary_v) * primary_v
        DAC_v_perpend = DAC_v - DAC_v_parallel
        DAC_v_perpend = np.linalg.norm(DAC_v_perpend) * \
                        (np.cos(phi_r) * secondary_v +
                         np.sin(phi_r) * tertiary_v)
        DAC_v = DAC_v_parallel + DAC_v_perpend
        # trim a copy of ball and calculate completeness
        q = copy.deepcopy(p)
        q.dac(opening_angle_in_radians=pressure_cell_oa, vector=DAC_v)
        q.resymmetrify(operations=point_group.chiral_operations)
        phi_psi_matrix[phi_index+1, psi_index+1] = q.table.shape[0]
        dat.write('{:8.0f}'.format(phi))
        dat.write('{:8.0f}'.format(psi))
        dat.write('{:8.5f}'.format(q.table.shape[0] / total_reflections))
        dat.write('{:8d}'.format(q.table.shape[0]))
        dat.write('\n')
        log.write('| ORIENTATION ' + str(orientation_no) + '/' +
              str(orientations_to_calculate) + '; TIME PASSED: ' +
              str(dt.now()-time_start) + '; EST END: ' +
              str(time_loop + (dt.now()-time_loop) *
                  orientations_to_calculate / orientation_no) + '\n')
    dat.write('\n')

# just for safety print the table to the log
time_end = dt.now()
log.write(20 * '-' + ' FINAL MATRIX ' + 20 * '-' + '\n')
log.write('> phi psi matrix:' + '\n')
log.write(np.array2string(phi_psi_matrix, separator=' ', max_line_width=999999))
log.write('\n> end time: ' + str(time_end) + '\n')
log.close()
dat.close()

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ END OF FILE ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ #
