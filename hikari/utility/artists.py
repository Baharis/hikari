import abc
from pathlib import Path
from matplotlib import pyplot, colors, cm
from mpl_toolkits.mplot3d import art3d
import numpy as np

from hikari.utility import gnuplot_map_palette, make_abspath, \
    mpl_map_palette, sph2cart
from hikari.resources import gnuplot_angular_heatmap_template


class ArtistError(Exception):
    """Exception raised when a problem with plotting in hikari occurs."""

    def __init__(self, message):
        super().__init__(message)


class ArtistFactory:
    """A factory method for creating artists."""
    def __init__(self):
        self._artists = {}

    def register(self, name, artist):
        self._artists[name] = artist

    def create(self, name, **kwargs):
        artist = self._artists.get(name)
        if not artist:
            raise ValueError(f'Artist called "{name}" has not been registered!')
        return artist(**kwargs)


class Artist:
    """Base class used for plotting matplotlib and gnuplot plots"""

    @staticmethod
    def _assert_is_iterable(iterable, length=0):
        try:
            _ = len(iterable)
            if not len(iterable) == length and length is not 0:
                raise TypeError()
        except TypeError:
            raise ArtistError(f'object {iterable} should be an iterable'
                              f' of length {length}' if length else '')

    @abc.abstractmethod
    def plot(self, path):
        pass


class AngularHeatmapArtist(Artist, abc.ABC):
    """Base class used for plotting angular heatmaps"""

    HEAT_EXTENSION = '.lst'
    HISTOGRAM_EXTENSION = '.his'

    def __init__(self):
        self.histogram = False
        self._x_axis = ()
        self._y_axis = ()
        self._z_axis = ()
        self._focus = ()
        self._heat_limits = ()
        self._polar_limits = ()
        self._azimuth_limits = ()

    @property
    def x_axis(self):
        return self._x_axis

    @x_axis.setter
    def x_axis(self, a):
        self._assert_is_iterable(a, 3)
        self._x_axis = tuple(a)

    @property
    def y_axis(self):
        return self._y_axis

    @y_axis.setter
    def y_axis(self, a):
        self._assert_is_iterable(a, 3)
        self._y_axis = tuple(a)

    @property
    def z_axis(self):
        return self._z_axis

    @z_axis.setter
    def z_axis(self, a):
        self._assert_is_iterable(a, 3)
        self._z_axis = tuple(a)

    @property
    def heat_limits(self):
        return self._heat_limits

    @heat_limits.setter
    def heat_limits(self, lims):
        self._assert_is_iterable(lims, 2)
        self._heat_limits = tuple(lims)

    @property
    def focus(self):
        return self._focus

    @focus.setter
    def focus(self, coords):
        self._assert_is_iterable(coords)
        [self._assert_is_iterable(c, 3) for c in coords]
        self._focus = tuple(coords)

    @property
    def polar_limits(self):
        return self._polar_limits

    @polar_limits.setter
    def polar_limits(self, lims):
        self._assert_is_iterable(lims, 2)
        self._polar_limits = tuple(lims)

    @property
    def azimuth_limits(self):
        return self._azimuth_limits

    @azimuth_limits.setter
    def azimuth_limits(self, lims):
        self._assert_is_iterable(lims, 2)
        self._azimuth_limits = tuple(lims)


class GnuplotArtist(abc.ABC):
    """Base class used for plotting gnuplot plots"""

    GNUPLOT_EXTENSION = '.gnu'

    def __init__(self):
        super().__init__()
        self._heat_palette = gnuplot_map_palette['']

    @property
    def heat_palette(self):
        return self._heat_palette

    @heat_palette.setter
    def heat_palette(self, name):
        self._heat_palette = gnuplot_map_palette[name]


class MatplotlibArtist(abc.ABC):
    """Base class used for plotting matplotlib plots"""

    def __init__(self):
        super().__init__()
        self._heat_palette = colors.LinearSegmentedColormap.from_list(
            'heatmap', mpl_map_palette[''], N=256)

    @property
    def heat_palette(self):
        return self._heat_palette

    @heat_palette.setter
    def heat_palette(self, name):
        self._heat_palette = colors.LinearSegmentedColormap.from_list(
            'heatmap', mpl_map_palette[name], N=256)


class GnuplotAngularHeatmapArtist(GnuplotArtist, AngularHeatmapArtist):
    """Base class used for plotting gnuplot angular heatmaps"""

    template = gnuplot_angular_heatmap_template

    @property
    def focus_string(self):
        label = "set label at {}, {}, {} '' point ls 10 front"
        return '\n'.join([label.format(*f) for f in self.focus])

    def plot(self, path):
        png_path = Path(make_abspath(path))
        directory, stem, ext = png_path.parent, png_path.stem, png_path.suffix
        gnu_name = png_path.stem + self.GNUPLOT_EXTENSION
        gnu_path = Path().joinpath(directory, gnu_name)
        s = gnuplot_angular_heatmap_template.format(
            axis_x1=self.x_axis[0], axis_x2=self.x_axis[1],
            axis_x3=self.x_axis[2], axis_y1=self.y_axis[0],
            axis_y2=self.y_axis[1], axis_y3=self.y_axis[2],
            axis_z1=self.z_axis[0], axis_z2=self.z_axis[1],
            axis_z3=self.z_axis[2],
            cplt_min=self.heat_limits[0],
            cplt_max=self.heat_limits[1],
            histogram=int(self.histogram),
            focus_string=self.focus_string,
            job_name=stem,
            min_ph=self.azimuth_limits[0],
            max_ph=self.azimuth_limits[1],
            min_th=self.polar_limits[0],
            max_th=self.polar_limits[1],
            palette=self.heat_palette)
        with open(gnu_path, 'w+') as f:
            f.write(s)
        try:
            from os import system, getcwd
            system('cd ' + str(directory) + '; gnuplot ' + gnu_name)
        except OSError:
            raise ArtistError(f'OSError passed: Cannot plot {gnu_name}')


class MatplotlibAngularHeatmapArtist(MatplotlibArtist, AngularHeatmapArtist):
    """Base class used for plotting matplotlib angular heatmaps"""

    MESH_EXTENSION = '.dat'

    def plot(self, path):
        # OS and I/O operations
        png_path = Path(make_abspath(path))
        directory, stem, ext = png_path.parent, png_path.stem, png_path.suffix
        mesh_name = png_path.stem + self.MESH_EXTENSION
        mesh_path = Path().joinpath(directory, mesh_name)
        heat_mesh = np.loadtxt(str(mesh_path))

        # set-up the plot
        fig = pyplot.figure(figsize=(8, 6))
        ax = fig.add_subplot(111, projection='3d')
        pyplot.rcParams.update({'font.size': 16})
        ax.view_init(elev=90 - sum(self.polar_limits) / 2,
                     azim=sum(self.azimuth_limits) / 2)
        ax.dist = 6.5
        ax.plot([1], [1], [1], 'w')
        ax.plot([-1], [-1], [-1], 'w')
        ax.set_axis_off()

        # prepare surface in cartesian coordinates
        polar_range = np.linspace(start=self.polar_limits[0],
                                  stop=self.polar_limits[1],
                                  num=heat_mesh.shape[1])  # 1=angle_res
        azimuth_range = np.linspace(start=self.azimuth_limits[0],
                                    stop=self.azimuth_limits[1],
                                    num=heat_mesh.shape[0])
        polar_mesh, azimuth_mesh = np.meshgrid(polar_range, azimuth_range)
        x_mesh, y_mesh, z_mesh = sph2cart(r=np.ones_like(polar_mesh),
                                          p=np.deg2rad(polar_mesh),
                                          a=np.deg2rad(azimuth_mesh))
        np.warnings.filterwarnings('ignore',  # mpl uses depreciated numpy
                                   category=np.VisibleDeprecationWarning)
        ax.plot_wireframe(x_mesh, y_mesh, z_mesh, colors='k', linewidth=0.25)

        # color map declarations
        m = cm.ScalarMappable(cmap=self.heat_palette)
        m.set_array(heat_mesh)
        m.set_clim(*self.heat_limits)
        pyplot.colorbar(m, fraction=0.06, pad=0.0, shrink=0.9)
        norm = colors.Normalize(*self.heat_limits)

        # draw (100), (010), (010) axes
        len_ = 1.25
        xa, xb = self.x_axis, [_ * len_ for _ in self.x_axis]
        ya, yb = self.y_axis, [_ * len_ for _ in self.y_axis]
        za, zb = self.z_axis, [_ * len_ for _ in self.z_axis]
        ax.add_line(art3d.Line3D((xa[0], xb[0]), (xa[1], xb[1]), (xa[2], xb[2]),
                                 color='r', linewidth=5, zorder=9))
        ax.text(*xb, '(100)', zorder=10, color='k')
        ax.add_line(art3d.Line3D((ya[0], yb[0]), (ya[1], yb[1]), (ya[2], yb[2]),
                                 color='g', linewidth=5, zorder=9))
        ax.text(*yb, '(010)', zorder=10, color='k')
        ax.add_line(art3d.Line3D((za[0], zb[0]), (za[1], zb[1]), (za[2], zb[2]),
                                 color='b', linewidth=5, zorder=9))
        ax.text(*zb, '(001)', zorder=10, color='k')

        # draw focus points
        if self.focus:
            xf, yf, zf = np.array(self.focus).T
            ax.plot(xf, yf, zf, linestyle='none', marker='D', markersize=10,
                    markerfacecolor='none', markeredgewidth=2,
                    markeredgecolor='k', zorder=8)

        # prepare smaller heat mesh for polygon centers and plot the heatmap
        face_heat_mesh = (heat_mesh[1:, 1:] + heat_mesh[1:, :-1] +
                          heat_mesh[:-1, 1:] + heat_mesh[:-1, :-1]) / 4
        color_mesh = self.heat_palette(norm(face_heat_mesh))
        for item in [fig, ax]:
            item.patch.set_visible(False)
        ax.plot_surface(x_mesh, y_mesh, z_mesh, rstride=1, cstride=1,
                        cmap=self.heat_palette, linewidth=0,
                        antialiased=False, facecolors=color_mesh)
        pyplot.subplots_adjust(left=0.0, bottom=0.0, right=0.95, top=1.0)
        pyplot.savefig(png_path, dpi=100, format='png', bbox_inches=None)


artist_factory = ArtistFactory()
artist_factory.register(name='gnuplot_angular_heatmap_artist',
                        artist=GnuplotAngularHeatmapArtist)
artist_factory.register(name='matplotlib_angular_heatmap_artist',
                        artist=MatplotlibAngularHeatmapArtist)
