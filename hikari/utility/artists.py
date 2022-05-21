import abc
from pathlib import Path
from .palettes import gnuplot_map_palette
from hikari.resources import potency_map_template


class ArtistError(Exception):
    """Exception raised when a problem with plotting in hikari occurs."""
    def __init__(self, message):
        super().__init__(message)


class Artist:
    HEAT_EXTENSION = '.lst'

    @staticmethod
    def _assert_iterable_length(iterable, length):
        try:
            if not len(iterable) == length:
                raise TypeError()
        except TypeError:
            raise ArtistError(f'object {iterable} should be'
                                      f'an iterable of length {length}')

    @abc.abstractmethod
    def plot(self, path):
        pass


class AngularHeatmapArtist(Artist, abc.ABC):
    def __init__(self):
        self._x_axis = ()
        self._y_axis = ()
        self._z_axis = ()
        self._heat_limits = ()
        self._polar_limits = ()
        self._azimuth_limits = ()

    @property
    def x_axis(self):
        return self._x_axis

    @x_axis.setter
    def x_axis(self, a):
        self._assert_iterable_length(a, 3)
        self._x_axis = tuple(a)

    @property
    def y_axis(self):
        return self._y_axis

    @y_axis.setter
    def y_axis(self, a):
        self._assert_iterable_length(a, 3)
        self._y_axis = tuple(a)

    @property
    def z_axis(self):
        return self._x_axis

    @z_axis.setter
    def z_axis(self, a):
        self._assert_iterable_length(a, 3)
        self._z_axis = tuple(a)

    @property
    def heat_limits(self):
        return self._heat_limits

    @heat_limits.setter
    def heat_limits(self, lims):
        self._assert_iterable_length(lims, 2)
        self._heat_limits = tuple(lims)

    @property
    def polar_limits(self):
        return self._polar_limits

    @polar_limits.setter
    def polar_limits(self, lims):
        self._assert_iterable_length(lims, 2)
        self._polar_limits = tuple(lims)

    @property
    def azimuth_limits(self):
        return self._azimuth_limits

    @azimuth_limits.setter
    def azimuth_limits(self, lims):
        self._assert_iterable_length(lims, 2)
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
    template = potency_map_template

    def plot(self, path):
        p = Path(path)
        directory, stem, ext = p.parent, p.stem, p.suffix
        gnu_file_name = p.stem + self.GNUPLOT_EXTENSION
        s = potency_map_template.format(
            axis_x1=self.x_axis[0], axis_x2=self.x_axis[1],
            axis_x3=self.x_axis[2], axis_y1=self.y_axis[0],
            axis_y2=self.y_axis[1], axis_y3=self.y_axis[2],
            axis_z1=self.z_axis[0], axis_z2=self.z_axis[1],
            axis_z3=self.z_axis[2],
            cplt_min=self.heat_limits[0],
            cplt_max=self.heat_limits[1],
            job_name=stem,
            min_ph=self.azimuth_limits[0],
            max_ph=self.azimuth_limits[1],
            min_th=self.polar_limits[0],
            max_th=self.polar_limits[1],
            palette=self.heat_palette)
        with open(gnu_file_name, 'w+') as f:
            f.write(s)
        try:
            from os import system, getcwd
            system('cd ' + str(directory) + '; gnuplot ' + gnu_file_name)
        except OSError:
            raise ArtistError(f'OSError passed: Cannot plot {gnu_file_name}')


class MatplotlibAngularHeatmapArtist(MatplotlibArtist, AngularHeatmapArtist):
    pass




