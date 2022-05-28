# prepare output
reset
set encoding utf8
set terminal pngcairo size 800,600 enhanced font 'Sans,16' solid
set output '{job_name}_gnu.png'

# load variables from python
make_histogram = {histogram}
min_heat = 100. * {cplt_min}
max_heat = 100. * {cplt_max}
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
set format cb "%.0f%%"
set colorbox user origin 0.825, 0.05 size 0.075, 0.90
set lmargin screen 0.00
set bmargin screen 0.00
set rmargin screen 0.72
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
set cbrange[min_heat:max_heat]
set xrange[-1.1:1.1]
set yrange[-1.1:1.1]
set zrange[-1.1:1.1]
set palette maxcolors 32
{palette}

# define axes
set arrow from X1, X2, X3 to 1.2*X1, 1.2*X2, 1.2*X3 as 1 front
set arrow from Y1, Y2, Y3 to 1.2*Y1, 1.2*Y2, 1.2*Y3 as 2 front
set arrow from Z1, Z2, Z3 to 1.2*Z1, 1.2*Z2, 1.2*Z3 as 3 front
set label at 1.2*X1, 1.2*X2, 1.2*X3 '(100)' offset screen 0.01, 0 nopoint front
set label at 1.2*Y1, 1.2*Y2, 1.2*Y3 '(010)' offset screen 0.01, 0 nopoint front
set label at 1.2*Z1, 1.2*Z2, 1.2*Z3 '(001)' offset screen 0.01, 0 nopoint front

set style line 10 pt 12 ps 3 lw 3 lc rgb '#000000'
{focus_string}

set multiplot

# draw everything, minding different coordinate system of splot
r = 1.001
splot '{job_name}.lst' using 2:(90-$1):(1):($3*100) with pm3d, \
      r*cos(min_ph)*sin(v),r*sin(min_ph)*sin(v),r*cos(v) with line ls 1, \
      r*cos(max_ph)*sin(v),r*sin(max_ph)*sin(v),r*cos(v) with line ls 1, \
      r*cos(u)*sin(max_th),r*sin(u)*sin(max_th),r*cos(max_th) with line ls 1

# if asked, add histogram based on .his file
if (make_histogram == 1) {{
unset parametric
stats '{job_name}.his' u 3 prefix 'heat' nooutput
set lmargin screen 0.825 - 0.15 / heat_max
set bmargin screen 0.050
set rmargin screen 0.825
set tmargin screen 0.950
unset label
unset arrow
set xrange [1:0]
set yrange [min_heat:max_heat]
set style fill solid border lt 0
plot '{job_name}.his' using (0):1:(0):3:($1*100):($2*100):(($1+$2)*50.0) with boxxy lc pal notitle
}}

unset multiplot