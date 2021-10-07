# prepare output
reset
set encoding utf8
set terminal pngcairo size 800,600 enhanced font 'Sans,16' solid
set output '{job_name}_gnu.png'

# load variables from python
min_ph = {min_ph}
max_ph = {max_ph}
min_th = {min_th}
max_th = {max_th}
X1 = {axis_x1}
X2 = {axis_x2}
X3 = {axis_x3}
Y1 = {axis_y1}
Y2 = {axis_y2}
Y3 = {axis_y3}
Z1 = {axis_z1}
Z2 = {axis_z2}
Z3 = {axis_z3}


# color definitions
set border lw 1.5
set style line 1  lt 1 lc rgb "#000000" lw 2     # map edge
set style arrow 1 lt 1 lc rgb "#ff0000" lw 5     # x arrow
set style arrow 2 lt 1 lc rgb "#008000" lw 5     # y arrow
set style arrow 3 lt 1 lc rgb "#0000ff" lw 5     # z arrow

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
set view (max_th+min_th)/2, (max_ph+min_ph)/2+90

# prepare spherical coordinate system
set parametric
set isosamples 25
set urange[min_ph:max_ph]
set vrange[min_th:max_th]
set cbrange[{cplt_min}:{cplt_max}]
set xrange[-1.1:1.1]
set yrange[-1.1:1.1]
set zrange[-1.1:1.1]
set palette maxcolors 40
{palette}

# define axes
set arrow from X1, X2, X3 to 1.2*X1, 1.2*X2, 1.2*X3 as 1 front
set arrow from Y1, Y2, Y3 to 1.2*Y1, 1.2*Y2, 1.2*Y3 as 2 front
set arrow from Z1, Z2, Z3 to 1.2*Z1, 1.2*Z2, 1.2*Z3 as 3 front

# draw everything, minding different coordinate system of splot
r = 1.001
splot '{job_name}.lst' using 2:(90-$1):(1):($3*100) with pm3d, \
      r*cos(v)*cos(0),r*cos(v)*sin(0),r*sin(v) with line ls 1, \
      r*cos(v)*cos(max_ph),r*cos(v)*sin(max_ph),r*sin(v) with line ls 1, \
      r*cos(0)*cos(u),r*cos(0)*sin(u),r*sin(0) with line ls 1