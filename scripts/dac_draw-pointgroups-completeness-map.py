# ~~~~~~~~~~~~~~~~~~~~~ IMPORT STATEMENTS - DO NOT CHANGE ~~~~~~~~~~~~~~~~~~~~ #
import copy
from kesshou.dataframes.hkl import HklFrame
from kesshou.symmetry.pointgroup import *
from mpl_toolkits.mplot3d import Axes3D
from mpl_toolkits.mplot3d import art3d
from matplotlib import cm, colors
import matplotlib.pyplot as plt
import numpy as np

# ~~~~~~~~~~~~~~~~~~~~ VARIABLES - CHANGE ONLY VALUES HERE ~~~~~~~~~~~~~~~~~~~ #
# Unit Cell (in Angstrom in degrees)
unit_cell_a = 2.456
unit_cell_b = 2.456
unit_cell_c = 6.694
unit_cell_al = 90.0
unit_cell_be = 90.0
unit_cell_ga = 120.0

# Laue group of the data
# select one from:  PG_1, PG2pm, PGmmm, PG4pm, PG4pmmm,
#                   PG_3, PG_3m1, PG6pm, PG6pmmm, PGm_3, PGm_3m]
laue_group = PG6pmmm

# List of diffraction pattern extinctions in a format "hk0: h+k=2n"
extinctions = ['hhl: l=2n',
               '00l: l=2n']

# Experimental details in angstrom and degrees
pressure_cell_oa = 35
hkl_wavelength = 'AgKa'
resolution_cutoff = 0.83

# Output details
output_name = 'graphite_oa35_AgKa_res83'
output_directory = '/home/dtchon/_/'
output_quality = 3
# integer from 1 (low) to 5 (high). Each next quality step increases time ~4x.

# ~~~~~~~~~~~~~~~~~~~~~~~~ SCRIPT CODE - DO NOT CHANGE ~~~~~~~~~~~~~~~~~~~~~~~ #

# prepare output files
dat_path = output_directory + output_name + '.dat'
gnu_path = output_directory + output_name + '.gnu'
lst_path = output_directory + output_name + '.lst'
png_path = output_directory + output_name + '.png'
png2_path = output_directory + output_name + '_gnu.png'
lst = open(lst_path, 'w+', buffering=1)

# generate reference ball
p = HklFrame()
p.crystal.edit_cell(a=unit_cell_a, b=unit_cell_b, c=unit_cell_c,
                    al=unit_cell_al, be=unit_cell_be, ga=unit_cell_ga)
p.edit_wavelength(hkl_wavelength)
p.make_ball(radius=min(p.r_lim, 1/resolution_cutoff))
p.extinct('000')
for extinction in extinctions:
    p.extinct(*extinction.split(':', 1), laue_group=laue_group)

# access and store number of all available reflections
total_reflections = len(p)

# define the spherical coordinate system based on given point group
# v1 is vector pointing in zenith direction (z* or y* for monoclinic)
# v2 is vector placed in orthogonal plane (x or z for monoclinic)
# TRICLINIC
if laue_group in {PG_1}:
    v1 = p.crystal.c_w
    v2 = p.crystal.a_v
    th_limits = [0, 180]
    ph_limits = [0, 180]
# MONOCLINIC
elif laue_group in {PG2pm}:
    v1 = p.crystal.c_w
    v2 = p.crystal.a_v
    th_limits = [0, 180]
    ph_limits = [0, 90]
# ORTHORHOMBIC / TETRAGONAL / CUBIC
elif laue_group in {PGmmm, PG4pm, PG4pmmm, PGm_3, PGm_3m}:
    v1 = p.crystal.c_w
    v2 = p.crystal.a_v
    # th_limits = [0, 180]
    # ph_limits = [0, 360]
    th_limits = [0, 90]
    ph_limits = [0, 90]
# TRIGONAL / CUBIC
elif laue_group in {PG_3, PG_3m1, PG6pm, PG6pmmm}:
    v1 = p.crystal.c_w
    v2 = p.crystal.a_v
    th_limits = [0, 90]
    ph_limits = [0, 120]
else:
    v1 = v2 = None
    raise ValueError('Provided group is not one of known Laue groups')

# define or redefine the _star vectors for the sake of drawing
# if laue_group in {PG2pm}:
#     x_star = p.crystal.c_w
#     y_star = p.crystal.a_w
#     z_star = p.crystal.b_w

# normalise the vectors and define the reference system
v1 /= np.linalg.norm(v1)
v2 /= np.linalg.norm(v2)
v3 = np.cross(v1, v2)

# define all orientation points (theta, phi) to be investigated
# phi and theta as in physics ISO convention (theta in 0 to 180,phi in 0 to 360)
# theta is rotation away from primary towards secondary vector (elevation angle)
# phi is rotation of DAC normal vector around primary vector (azimuth angle)
assert output_quality in {1, 2, 3, 4, 5}
angle_resolution = {1: 15, 2: 10, 3: 5, 4: 2, 5: 1}[output_quality]
th_range = np.arange(th_limits[0], th_limits[1] + 0.001, angle_resolution)
ph_range = np.arange(ph_limits[0], ph_limits[1] + 0.001, angle_resolution)
th_mesh, ph_mesh = np.meshgrid(th_range, ph_range)
th_list = np.concatenate(th_mesh)
ph_list = np.concatenate(ph_mesh)
cplt_mesh = np.zeros_like(th_mesh)
data_dict = {'th': [], 'ph': [], 'cplt': [], 'reflns': []}

# calculate completeness for each investigated phi / theta pair
lst.write('#     th     psi    cplt  reflns\n')
for i, th in enumerate(th_range):
    for j, ph in enumerate(ph_range):
        # calculate phi and theta in radians
        ph_r, th_r = np.deg2rad(ph), np.deg2rad(th)
        # find the DAC orientation vector around primary-secondary plane
        v = v1 * np.cos(th_r) + v2 * np.sin(th_r)
        # rotate the DAC orientation vector in the secondary / tertiary plane
        v_parallel = np.dot(v, v1) * v1
        v_perpend = v - v_parallel
        v_perpend = lin.norm(v_perpend) * (np.cos(ph_r) * v2 + np.sin(ph_r) * v3)
        v = v_parallel + v_perpend
        # trim a copy of ball and calculate completeness
        q = copy.deepcopy(p)
        q.dac(opening_angle=pressure_cell_oa, vector=v)
        q.resymmetrify(operations=laue_group.hp_disc_transforming_symm_ops,
                       reduce=True)
        # save and write the data
        data_dict['th'].append(th)
        lst.write('{:8.0f}'.format(th))
        data_dict['ph'].append(ph)
        lst.write('{:8.0f}'.format(ph))
        data_dict['cplt'].append(len(q) / total_reflections)
        lst.write('{:8.5f}'.format(len(q) / total_reflections))
        cplt_mesh[j][i] = len(q) / total_reflections
        data_dict['reflns'].append(len(q))
        lst.write('{:8d}'.format(len(q)))
        lst.write('\n')
    lst.write('\n')
lst.write('\n')
lst.close()

# PLOTTING IN MATPLOTLIB
# prepare the plot
fig = plt.figure()
fig.set_size_inches(5, 3)
ax = fig.add_subplot(111, projection='3d')
ax.view_init(elev=90-sum(th_limits)/2, azim=sum(ph_limits)/2)
ax.dist = 10
ax.set_axis_off()
norm = colors.Normalize(vmin=0, vmax=1, clip=False)

# prepare a list of cartesian coordinates
x = np.sin(np.deg2rad(th_mesh)) * np.cos(np.deg2rad(ph_mesh))
y = np.sin(np.deg2rad(th_mesh)) * np.sin(np.deg2rad(ph_mesh))
z = np.cos(np.deg2rad(th_mesh))

# add wireframe
wf_lw = min(angle_resolution * 0.1, 1)
ax.plot_wireframe(1.001 * x, 1.001 * y, 1.001 * z, colors='k', linewidth=0.25)

# add the colorbar
heatmapEX_colors = [(0.0, 0.0, 0.0),    # Black
                    (0.0, 0.0, 0.5),    # Indigo
                    (0.0, 0.0, 1.0),    # Blue
                    (0.0, 0.5, 1.0),    # Aqua
                    (0.0, 1.0, 1.0),    # Turquoise
                    (0.5, 1.0, 0.5),    # Lime
                    (1.0, 1.0, 0.0),    # Yellow
                    (1.0, 0.5, 0.0),    # Orange
                    (1.0, 0.0, 0.0),    # Red
                    (1.0, 0.5, 0.5),    # Light red
                    (1.0, 1.0, 1.0)]    # White
heatmapEX_colormap = colors.LinearSegmentedColormap.from_list(
    'heatmapEX', heatmapEX_colors, N=256)
m = cm.ScalarMappable(cmap=heatmapEX_colormap)
m.set_array(cplt_mesh)
plt.colorbar(m, fraction=0.046, pad=0.04)

# add rgb x*, y*, z* axis lines
line_len = 1.25
line_x = p.crystal.a_w / lin.norm(p.crystal.a_w)
line_y = p.crystal.b_w / lin.norm(p.crystal.b_w)
line_z = p.crystal.c_w / lin.norm(p.crystal.c_w)
line_x = art3d.Line3D((1 * line_x[0], line_len * line_x[0]),
                      (1 * line_x[1], line_len * line_x[1]),
                      (1 * line_x[2], line_len * line_x[2]),
                      color='r', linewidth=5)
line_y = art3d.Line3D((1 * line_y[0], line_len * line_y[0]),
                      (1 * line_y[1], line_len * line_y[1]),
                      (1 * line_y[2], line_len * line_y[2]),
                      color='g', linewidth=5)
line_z = art3d.Line3D((1 * line_z[0], line_len * line_z[0]),
                      (1 * line_z[1], line_len * line_z[1]),
                      (1 * line_z[2], line_len * line_z[2]),
                      color='b', linewidth=5)
ax.add_line(line_x)
ax.add_line(line_y)
ax.add_line(line_z)

# plot the colormap
for item in [fig, ax]:
    item.patch.set_visible(False)
plt.tight_layout()
ax.plot_surface(x, y, z, rstride=1, cstride=1, cmap=heatmapEX_colormap, linewidth=0,
                antialiased=False, facecolors=heatmapEX_colormap(cplt_mesh),
                vmin=min(data_dict['cplt']), vmax=max(data_dict['cplt']))
plt.savefig(png_path, dpi=600, format='png', bbox_inches=None)

# PLOTTING IN GNUPLOT
# prepare the gnuplot file which will draw the results from csv
gnu = open(gnu_path, 'w+', buffering=1)
gnu.write("""
# prepare output
reset
set encoding utf8
set terminal pngcairo size 800,600 enhanced font 'Sans,16' solid
set output '{gnu_output_name}'

# color definitions
set border lw 1.5
set style line 1 lt 1 lc rgb "#000000" lw 2     # x direction ring
set style line 2 lt 1 lc rgb "#000000" lw 2     # y direction ring
set style line 3 lt 1 lc rgb "#000000" lw 2     # z direction ring
set style arrow 1 lt 1 lc rgb "#ff0000" lw 5     # x arrow
set style arrow 2 lt 1 lc rgb "#008000" lw 5     # y arrow
set style arrow 3 lt 1 lc rgb "#0000ff" lw 5     # z arrow
#set style arrow 4 lt 1 lc rgb "#ffffff" lw 6     # arrow shade

# unset border and tics
unset key
unset border
set format x ''
set format y ''
set format z ''
set tics scale 0
set cbtics
set format cb "%.0f%%"
set colorbox user origin 0.825, 0.05 size 0.075, 0.9
set lmargin screen 0.04
set bmargin screen 0.00
set rmargin screen 0.76
set tmargin screen 0.90

# prepare mapping on sphere
set mapping spherical
set angles degrees
set xyplane at -1
set view {gnu_avg_th}, {gnu_avg_ph_plus90}

# prepare spherical coordinate system
set parametric
set isosamples 25
set urange[{gnu_min_ph}:{gnu_max_ph}]
set vrange[{gnu_min_th}:{gnu_max_th}]
#set urange[{gnu_min_ph_minus25}:{gnu_max_ph_plus25}]
#set vrange[{gnu_min_th_minus25}:{gnu_max_th_plus25}]
set cbrange[{gnu_cplt_min}:{gnu_cplt_max}]
set xrange[-1.1:1.1]
set yrange[-1.1:1.1]
set zrange[-1.1:1.1]
set palette maxcolors 50
set palette model HSV defined (0 0.66 1 0, 5 0.66 1 1, 15 0.0 1 1, 20 0.0 0 1) 
# , 25 0.0 0 0.5) to add gray on the end

# define axes and their labels 
#set label "x*" at 1.5,0.1,0 center front textcolor rgb "#ff0000"
set arrow from {gnu_arrow_x1},{gnu_arrow_x2},{gnu_arrow_x3} to {gnu_arrow_x4},{gnu_arrow_x5},{gnu_arrow_x6} as 1 front
#set label "y*" at 0.1,1.5,0 center front textcolor rgb "#008000"
set arrow from {gnu_arrow_y1},{gnu_arrow_y2},{gnu_arrow_y3} to {gnu_arrow_y4},{gnu_arrow_y5},{gnu_arrow_y6} as 2 front
#set label "z*" at -0.07,0.07,1.47 center front textcolor rgb "#0000ff"
set arrow from {gnu_arrow_z1},{gnu_arrow_z2},{gnu_arrow_z3} to {gnu_arrow_z4},{gnu_arrow_z5},{gnu_arrow_z6} as 3 front

# draw everything
# here for splot psi is redefined to match gnuplot reference frame and cplt
# is multiplied to be expressed in percents
r = 1.001
splot '{gnu_input_name}' using 2:(90-$1):(1):($3*100) with pm3d,\\
      r*cos(v)*cos(0),r*cos(v)*sin(0),r*sin(v) with line ls 1, \\
      r*cos(v)*cos({gnu_max_ph}),r*cos(v)*sin({gnu_max_ph}),r*sin(v) with line ls 2, \\
      r*cos(0)*cos(u),r*cos(0)*sin(u),r*sin(0) with line ls 3
      """.format(gnu_arrow_x1=p.crystal.a_w[0] / lin.norm(p.crystal.a_w),
                 gnu_arrow_x2=p.crystal.a_w[1] / lin.norm(p.crystal.a_w),
                 gnu_arrow_x3=p.crystal.a_w[2] / lin.norm(p.crystal.a_w),
                 gnu_arrow_x4=p.crystal.a_w[0] / lin.norm(p.crystal.a_w) * 1.2,
                 gnu_arrow_x5=p.crystal.a_w[1] / lin.norm(p.crystal.a_w) * 1.2,
                 gnu_arrow_x6=p.crystal.a_w[2] / lin.norm(p.crystal.a_w) * 1.2,
                 gnu_arrow_y1=p.crystal.b_w[0] / lin.norm(p.crystal.b_w),
                 gnu_arrow_y2=p.crystal.b_w[1] / lin.norm(p.crystal.b_w),
                 gnu_arrow_y3=p.crystal.b_w[2] / lin.norm(p.crystal.b_w),
                 gnu_arrow_y4=p.crystal.b_w[0] / lin.norm(p.crystal.b_w) * 1.2,
                 gnu_arrow_y5=p.crystal.b_w[1] / lin.norm(p.crystal.b_w) * 1.2,
                 gnu_arrow_y6=p.crystal.b_w[2] / lin.norm(p.crystal.b_w) * 1.2,
                 gnu_arrow_z1=p.crystal.c_w[0] / lin.norm(p.crystal.c_w),
                 gnu_arrow_z2=p.crystal.c_w[1] / lin.norm(p.crystal.c_w),
                 gnu_arrow_z3=p.crystal.c_w[2] / lin.norm(p.crystal.c_w),
                 gnu_arrow_z4=p.crystal.c_w[0] / lin.norm(p.crystal.c_w) * 1.2,
                 gnu_arrow_z5=p.crystal.c_w[1] / lin.norm(p.crystal.c_w) * 1.2,
                 gnu_arrow_z6=p.crystal.c_w[2] / lin.norm(p.crystal.c_w) * 1.2,
                 gnu_avg_th=max(th_limits) / len(th_limits),
                 gnu_avg_ph_plus90=max(ph_limits) / len(ph_limits) + 90,
                 gnu_cplt_min=min(data_dict['cplt']) * 100,
                 gnu_cplt_max=max(data_dict['cplt']) * 100,
                 gnu_input_name=lst_path,
                 gnu_min_ph=min(ph_range),
                 gnu_max_ph=max(ph_range),
                 gnu_min_th=min(th_range),
                 gnu_max_th=max(th_range),
                 gnu_min_ph_minus25=min(ph_range) - 25,
                 gnu_max_ph_plus25=max(ph_range) + 25,
                 gnu_min_th_minus25=min(th_range) - 25,
                 gnu_max_th_plus25=max(th_range) + 25,
                 gnu_output_name=png2_path))







# TODO there is still something wrong with the scale on gnuplot graphic (too large values)
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ END OF FILE ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ #
