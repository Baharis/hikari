import abc
from pathlib import Path
from .palettes import gnuplot_map_palette
from hikari.resources import gnuplot_angular_heatmap_template


class ArtistError(Exception):
    """Exception raised when a problem with plotting in hikari occurs."""
    def __init__(self, message):
        super().__init__(message)


class Artist:
    @staticmethod
    def _assert_is_iterable(iterable, length=0):
        try:
            _ = len(iterable)
            if not len(iterable) == length and length is not 0:
                raise TypeError()
        except TypeError:
            raise ArtistError(f'object {iterable} should be an iterable' + \
                              f'of length {length}' if length else '')

    @abc.abstractmethod
    def plot(self, path):
        pass


class AngularHeatmapArtist(Artist, abc.ABC):
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
    GNUPLOT_EXTENSION = '.gnu'

    def __init__(self):
        super().__init__()
        self._heat_palette = ''

    @property
    def heat_palette(self):
        return self._heat_palette

    @heat_palette.setter
    def heat_palette(self, name):
        self._heat_palette = gnuplot_map_palette[name]


class MatplotlibArtist(abc.ABC):
    pass


class GnuplotAngularHeatmapArtist(GnuplotArtist, AngularHeatmapArtist):
    template = gnuplot_angular_heatmap_template

    @property
    def focus_string(self):
        label = "set label at {}, {}, {} '' point ls 10 front"
        return '\n'.join([label.format(*f) for f in self.focus])

    def plot(self, path):
        png_path = Path(path)
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
    pass




