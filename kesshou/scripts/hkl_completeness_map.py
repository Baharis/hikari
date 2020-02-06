from kesshou.dataframes import HklFrame
from kesshou.symmetry import PG
from kesshou.utility import make_absolute_path, home_directory
from matplotlib import cm, colors, pyplot
from mpl_toolkits.mplot3d import Axes3D
from mpl_toolkits.mplot3d import art3d
import numpy as np
import numpy.linalg as lin


def completeness_map(a, b, c, al, be, ga,
                     laue_group=PG['1'],
                     extinctions=tuple(),
                     fix_scale=False,
                     legacy_cplt=False,
                     opening_angle=35,
                     output_directory=home_directory,
                     output_name='cplt_map',
                     output_quality=3,
                     resolution=0.83,
                     wavelength='MoKa'):

    dat_path = make_absolute_path(output_directory, output_name + '.dat')
    gnu_path = make_absolute_path(output_directory, output_name + '.gnu')
    lst_path = make_absolute_path(output_directory, output_name + '.lst')
    png_path = make_absolute_path(output_directory, output_name + '.png')
    png2_path = make_absolute_path(output_directory, output_name + '_gnu.png')

    def _make_reference_ball():
        """Make ball of hkl which will be cut in further steps"""
        _hkl_frame = HklFrame()
        _hkl_frame.edit_cell(a=a, b=b, c=c, al=al, be=be, ga=ga)
        _hkl_frame.edit_wavelength(wavelength)
        _hkl_frame.make_ball(radius=min(_hkl_frame.r_lim, 1/resolution))
        _hkl_frame.extinct('000')
        for extinction in extinctions:
            _hkl_frame.extinct(extinction, point_group=laue_group)
        return _hkl_frame
    p = _make_reference_ball()
    p.find_equivalents(point_group=laue_group)
    total_reflections = len(p) if legacy_cplt else p.data['equiv'].nunique()

    def _determine_theta_and_phi_limits():
        """Define the spherical coordinate system based on given point group.
        v1, v2, v3 are normal vectors pointing in zenith direction z*,
        orthogonal direction (x) and direction orthogonal to them both."""
        _v1 = p.z_w
        _v2 = p.x_v
        _v3 = np.cross(_v1, _v2)
        if laue_group in {PG['-1']}:
            _th_limits = [0, 180]
            _ph_limits = [0, 180]
        # MONOCLINIC
        elif laue_group in {PG['2/m']}:
            _th_limits = [0, 180]
            _ph_limits = [0, 90]
        # ORTHORHOMBIC / TETRAGONAL / CUBIC
        elif laue_group in {PG['mmm'], PG['4/m'], PG['4/mmm'],
                            PG['m-3'], PG['m-3m']}:
            _th_limits = [0, 90]
            _ph_limits = [0, 90]
        # TRIGONAL / CUBIC
        elif laue_group in {PG['-3'], PG['-3m1'], PG['-31m'],
                            PG['6/m'], PG['6/mmm']}:
            _th_limits = [0, 90]
            _ph_limits = [0, 120]
        else:
            raise ValueError('Provided group is not one of known Laue groups')
        return _v1, _v2, _v3, _th_limits, _ph_limits
    v1, v2, v3, th_limits, ph_limits = _determine_theta_and_phi_limits()

    def _make_theta_and_phi_mesh():
        """Define a list of theta and phi values to be investigated.
        Theta and phi are defined as in physics ISO convention.
        Theta is "elevation angle" away from primary (zenith) direction,
        towards secondary direction, and takes values from 0 to 180 degrees.
        Phi is "azimuth angle" which rotates perpendicularly to primary,
        from secondary to tertiary direction. Takes values from 0 to 360."""
        assert output_quality in {1, 2, 3, 4, 5}, 'Quality not 1, 2, 3, 4 or 5'
        angle_res = {1: 15, 2: 10, 3: 5, 4: 2, 5: 1}[output_quality]
        _th_range = np.arange(th_limits[0], th_limits[1] + 0.001, angle_res)
        _ph_range = np.arange(ph_limits[0], ph_limits[1] + 0.001, angle_res)
        _th_mesh, _ph_mesh = np.meshgrid(_th_range, _ph_range)
        return _th_range, _ph_range, _th_mesh, _ph_mesh
    th_range, ph_range, th_mesh, ph_mesh = _make_theta_and_phi_mesh()
    data_dict = {'th': [], 'ph': [], 'cplt': [], 'reflns': []}

    def _translate_angles_to_vector(theta, phi):
        """Find the vector by rotating v1 by theta and then phi, in degrees."""
        _sin_th = np.sin(np.deg2rad(theta))
        _sin_ph = np.sin(np.deg2rad(phi))
        _cos_th = np.cos(np.deg2rad(theta))
        _cos_ph = np.cos(np.deg2rad(phi))
        _v = v1 * _cos_th + v2 * _sin_th
        _v_parallel = np.dot(_v, v1) * v1
        _v_perpend = lin.norm(_v - _v_parallel) * (_cos_ph * v2 + _sin_ph * v3)
        return _v_parallel + _v_perpend

    def _calculate_completeness_mesh():
        """Calculate completeness for each individual pair of theta and phi."""
        _cplt_mesh = np.zeros_like(th_mesh)
        lst = open(lst_path, 'w+', buffering=1)
        lst.write('#     th      ph    cplt  reflns\n')
        for i, th in enumerate(th_range):
            for j, ph in enumerate(ph_range):
                v = _translate_angles_to_vector(theta=th, phi=ph)
                q = p.duplicate()
                q.dac(opening_angle=opening_angle, vector=v)
                if legacy_cplt:
                    q.transform(operations=laue_group.chiral_operations)
                    q.merge()
                    hkl_len = len(q)
                else:
                    hkl_len = q.data['equiv'].nunique()
                data_dict['th'].append(th)
                data_dict['ph'].append(ph)
                data_dict['cplt'].append(hkl_len / total_reflections)
                data_dict['reflns'].append(hkl_len)
                lst.write('{:8.0f}'.format(th))
                lst.write('{:8.0f}'.format(ph))
                lst.write('{:8.5f}'.format(hkl_len / total_reflections))
                lst.write('{:8d}'.format(hkl_len))
                lst.write('\n')
                _cplt_mesh[j][i] = hkl_len / total_reflections
            lst.write('\n')
        index_max = np.unravel_index(np.argmax(_cplt_mesh), _cplt_mesh.shape)
        best_th, best_ph = th_range[index_max[1]], ph_range[index_max[0]]
        index_min = np.unravel_index(np.argmin(_cplt_mesh), _cplt_mesh.shape)
        worst_th, worst_ph = th_range[index_min[1]], ph_range[index_min[0]]
        lst.write('# best cplt = {_max} for th = {_th}, ph = {_ph}\n'.format(
            _max=max(data_dict['cplt']), _th=best_th, _ph=best_ph))
        lst.write('# worst cplt = {_min} for th = {_th}, ph = {_ph}\n'.format(
            _min=min(data_dict['cplt']), _th=worst_th, _ph=worst_ph))
        lst.write('# mean cplt = {}\n'.format(np.mean(data_dict['cplt'])))
        lst.close()
        np.savetxt(dat_path, _cplt_mesh)
        return data_dict, _cplt_mesh
    data_dict, cplt_mesh = _calculate_completeness_mesh()

    def _plot_in_matplotlib():
        """Plot the completeness map in radial coordinates using matplotlib"""
        fig = pyplot.figure(figsize=(5, 3))
        ax = fig.add_subplot(111, projection='3d')
        ax.view_init(elev=90-sum(th_limits)/2, azim=sum(ph_limits)/2)
        ax.dist = 10
        ax.set_axis_off()

        # surface in cartesian coordinates
        x = np.sin(np.deg2rad(th_mesh)) * np.cos(np.deg2rad(ph_mesh))
        y = np.sin(np.deg2rad(th_mesh)) * np.sin(np.deg2rad(ph_mesh))
        z = np.cos(np.deg2rad(th_mesh))

        # wireframe
        ax.plot_wireframe(x, y, z, colors='k', linewidth=0.25)

        # color map
        my_heatmap_colors = [(0.0, 0.0, 0.0),    # Black
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
        my_colormap = colors.LinearSegmentedColormap.from_list(
            'heatmapEX', my_heatmap_colors, N=256)
        m = cm.ScalarMappable(cmap=my_colormap)
        m.set_array(cplt_mesh)
        if fix_scale is True:
            m.set_clim(0, 1)
            norm = colors.Normalize(vmin=0, vmax=1)
        else:
            norm = colors.Normalize(vmin=min(data_dict['cplt']),
                                    vmax=max(data_dict['cplt']))
        pyplot.colorbar(m, fraction=0.046, pad=0.04)

        # direction lines
        _len = 1.25
        _x, _y, _z = p.x_w, p.y_w, p.z_w
        ax.add_line(art3d.Line3D((_x[0], _len * _x[0]), (_x[1], _len * _x[1]),
                                 (_x[2], _len * _x[2]), color='r', linewidth=5))
        ax.add_line(art3d.Line3D((_y[0], _len * _y[0]), (_y[1], _len * _y[1]),
                                 (_y[2], _len * _y[2]), color='g', linewidth=5))
        ax.add_line(art3d.Line3D((_z[0], _len * _z[0]), (_z[1], _len * _z[1]),
                                 (_z[2], _len * _z[2]), color='b', linewidth=5))

        # color mesh for heatmap
        color_mesh = cplt_mesh[:-1, :-1]
        for i in range(cplt_mesh.shape[0]-1):
            for j in range(cplt_mesh.shape[1]-1):
                color_mesh[i, j] = (cplt_mesh[i + 1, j] + cplt_mesh[i, j + 1] +
                                    cplt_mesh[i, j] + cplt_mesh[i + 1, j + 1])/4

        # heatmap surface
        for item in [fig, ax]:
            item.patch.set_visible(False)
        ax.plot_surface(x, y, z, rstride=1, cstride=1, cmap=my_colormap,
                        linewidth=0, antialiased=False,
                        facecolors=my_colormap(norm(color_mesh)))
        pyplot.savefig(png_path, dpi=600, format='png', bbox_inches=None)
    _plot_in_matplotlib()

    def _prepare_gnuplot_input():
        """Prepare input to completeness map in radial coordinates in gnuplot"""
        gnu = open(gnu_path, 'w+', buffering=1)
        arrow_len = 1.2
        gnu_cplt_max = 100 if fix_scale else max(data_dict['cplt']) * 100
        gnu_cplt_min = 0 if fix_scale else min(data_dict['cplt']) * 100
        gnu.write("""# gnuplot file generated by Kesshou package
            
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
            """.format(
            gnu_arrow_x1=p.x_w[0],
            gnu_arrow_x2=p.x_w[1],
            gnu_arrow_x3=p.x_w[2],
            gnu_arrow_x4=p.x_w[0] * arrow_len,
            gnu_arrow_x5=p.x_w[1] * arrow_len,
            gnu_arrow_x6=p.x_w[2] * arrow_len,
            gnu_arrow_y1=p.y_w[0],
            gnu_arrow_y2=p.y_w[1],
            gnu_arrow_y3=p.y_w[2],
            gnu_arrow_y4=p.y_w[0] * arrow_len,
            gnu_arrow_y5=p.y_w[1] * arrow_len,
            gnu_arrow_y6=p.y_w[2] * arrow_len,
            gnu_arrow_z1=p.z_w[0],
            gnu_arrow_z2=p.z_w[1],
            gnu_arrow_z3=p.z_w[2],
            gnu_arrow_z4=p.z_w[0] * arrow_len,
            gnu_arrow_z5=p.z_w[1] * arrow_len,
            gnu_arrow_z6=p.z_w[2] * arrow_len,
            gnu_avg_th=max(th_limits) / len(th_limits),
            gnu_avg_ph_plus90=max(ph_limits) / len(ph_limits) + 90,
            gnu_cplt_min=gnu_cplt_min,
            gnu_cplt_max=gnu_cplt_max,
            gnu_input_name=lst_path,
            gnu_min_ph=min(ph_limits),
            gnu_max_ph=max(ph_limits),
            gnu_min_th=min(th_limits),
            gnu_max_th=max(th_limits),
            gnu_min_ph_minus25=min(ph_limits) - 25,
            gnu_max_ph_plus25=max(ph_limits) + 25,
            gnu_min_th_minus25=min(th_limits) - 25,
            gnu_max_th_plus25=max(th_limits) + 25,
            gnu_output_name=png2_path))
    _prepare_gnuplot_input()


if __name__ == '__main__':
    completeness_map(a=10.0, b=10.0, c=10.0, al=90.0, be=90.0, ga=90.0,
                     extinctions=('h00: h=2n', '0k0: k=2n', '00l: l=2n'),
                     laue_group=PG['m-3m'], output_quality=3, fix_scale=True)
