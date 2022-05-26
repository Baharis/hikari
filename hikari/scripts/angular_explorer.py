"""This file contains tools for making property maps visualised on sphere"""

import abc
from pathlib import Path

import numpy as np
from numpy import linalg as lin

from hikari.dataframes import HklFrame
from hikari.symmetry import SG, Group
from hikari.utility import make_abspath, weighted_quantile, \
    fibonacci_sphere, rotation_around, sph2cart, cart2sph, Interval
from hikari.utility import GnuplotAngularHeatmapArtist, \
    MatplotlibAngularHeatmapArtist


class AngularPropertyExplorerFactory:
    """A factory method for creating angular property explorers."""
    def __init__(self):
        self._explorers = {}

    def register_explorer(self, prop, explorer):
        self._explorers[prop] = explorer

    def create(self, prop, **kwargs):
        explorer = self._explorers.get(prop)
        if not explorer:
            raise ValueError(f'Explorer for {prop} has not been registered!')
        return explorer(**kwargs)


class AngularPropertyExplorer:
    """
    An abstract base class for objects handling analysing parameters as
    a function of crystal orientation in a Diamond Anvil Cell.
    In order to generate a map of desired property,
    the following methods must be executed in order:

    - set_path
    - set_sample
    - set_output
    - explore
    """

    hkl_is_read_not_generated = True
    property_name = 'UNDEFINED'
    property_theoretical_minimum = 0
    property_theoretical_maximum = 1

    GNUPLOT_INPUT_EXTENSION = '.gnu'
    GNUPLOT_OUTPUT_EXTENSION = '.pnG'
    HISTOGRAM_EXTENSION = '.his'
    HKL_EXTENSION = '.hkl'
    LISTING_EXTENSION = '.lst'
    MATPLOTLIB_EXTENSION = '.png'
    MESH_EXTENSION = '.dat'

    POLAR_LIMIT_DICT = {Group.System.triclinic: Interval(0, 180),
                        Group.System.monoclinic: Interval(0, 180),
                        Group.System.orthorhombic: Interval(0, 90),
                        Group.System.tetragonal: Interval(0, 90),
                        Group.System.trigonal: Interval(0, 90),
                        Group.System.hexagonal: Interval(0, 90),
                        Group.System.cubic: Interval(0, 90)}

    AZIMUTH_LIMIT_DICT = {Group.System.triclinic: Interval(-45, 135),
                          Group.System.monoclinic: Interval(0, 90),
                          Group.System.orthorhombic: Interval(0, 90),
                          Group.System.tetragonal: Interval(0, 90),
                          Group.System.trigonal: Interval(0, 120),
                          Group.System.hexagonal: Interval(0, 120),
                          Group.System.cubic: Interval(0, 90)}

    def __init__(self):
        self._path = ''
        self.sg = Group()
        self.axis = ''
        self.hkl_frame = HklFrame()
        self.oa = 0.0
        self.orientation = None
        self.resolution = 1.2
        self.axis = ''
        self.fix_scale = False
        self.histogram = False
        self.output_quality = 1
        self.data_dict = {'th': [], 'ph': [], 'cplt': [], 'r1': [], 'weight': []}

    def set_experimental(self, opening_angle, orientation, resolution):
        self.oa = opening_angle
        self.orientation = orientation
        self.resolution = resolution

    def set_hkl_frame(self, a, b, c, al, be, ga, space_group, wavelength, axis):
        self.sg = SG[space_group]
        self.hkl_frame.edit_cell(a=a, b=b, c=c, al=al, be=be, ga=ga)
        self.hkl_frame.la = wavelength
        self.axis = axis
        if self.hkl_is_read_not_generated:
            self._read_hkl_frame(hkl_format='shelx_4')
        else:
            self._make_hkl_frame()
        self.hkl_frame.find_equivalents(point_group=self.pg)
        total_reflections = self.hkl_frame.table['equiv'].nunique()
        if total_reflections == 0:
            raise KeyError('Specified part of reciprocal space has zero nodes')

    def set_options(self, path, fix_scale, histogram, output_quality):
        self.path = path
        self.fix_scale = fix_scale
        self.histogram = histogram
        self.output_quality = output_quality

    @abc.abstractmethod
    def explore(self):
        """
        Main interface method which, apart from the setters, is run within
        the script itself and generates data and figures using class methods
        """
        pass

    @property
    def prop_limits(self):
        if not self.fix_scale and self.data_dict[self.property_name]:
            lower_property_limit = min(self.data_dict[self.property_name])
            upper_property_limit = max(self.data_dict[self.property_name])
        else:
            lower_property_limit = self.property_theoretical_minimum
            upper_property_limit = self.property_theoretical_maximum
        return Interval(lower_property_limit, upper_property_limit)

    def _make_hkl_frame(self):
        """Make ball or axis of hkl which will be cut in further steps"""
        f = self.hkl_frame
        f.fill(radius=min(self.hkl_frame.r_lim, self.resolution))
        if self.axis in {'x'}:
            f.table = f.table.loc[f.table['k'].eq(0) & f.table['l'].eq(0)]
        elif self.axis in {'y'}:
            f.table = f.table.loc[f.table['h'].eq(0) & f.table['l'].eq(0)]
        elif self.axis in {'z'}:
            f.table = f.table.loc[f.table['h'].eq(0) & f.table['k'].eq(0)]
        elif self.axis in {'xy'}:
            f.table = f.table.loc[f.table['l'].eq(0)]
        elif self.axis in {'xz'}:
            f.table = f.table.loc[f.table['k'].eq(0)]
        elif self.axis in {'yz'}:
            f.table = f.table.loc[f.table['h'].eq(0)]
        if self.axis in {'x', 'y', 'z', 'xy', 'xz', 'yz'}:
            f.transform([o.tf for o in self.pg.operations])
        f.extinct(self.sg)
        return f

    def _read_hkl_frame(self, hkl_format):
        """Read reflections of hkl which will be cut in further steps"""
        f = self.hkl_frame
        f.read(hkl_path=self.path + self.HKL_EXTENSION, hkl_format=hkl_format)
        f.trim(limit=min(f.r_lim, self.resolution))
        f.find_equivalents(point_group=self.lg)
        return f

    @property
    def th_limits(self):
        """Interval range of polar angle where property will be calculated"""
        return self.POLAR_LIMIT_DICT[self.sg.system]

    @property
    def ph_limits(self):
        """Interval range of azimuth angle where property will be calculated"""
        return self.AZIMUTH_LIMIT_DICT[self.sg.system]

    def _draw_matplotlib_map(self):
        ma = MatplotlibAngularHeatmapArtist()
        ma.x_axis = self.hkl_frame.a_w / lin.norm(self.hkl_frame.a_w)
        ma.y_axis = self.hkl_frame.b_w / lin.norm(self.hkl_frame.b_w)
        ma.z_axis = self.hkl_frame.c_w / lin.norm(self.hkl_frame.c_w)
        ma.focus = self.focus
        ma.heat_limits = self.prop_limits
        ma.heat_palette = self.axis
        ma.polar_limits = self.th_limits
        ma.azimuth_limits = self.ph_limits
        ma.plot(self.path + self.MATPLOTLIB_EXTENSION)

    def _draw_gnuplot_map(self):
        ga = GnuplotAngularHeatmapArtist()
        ga.x_axis = self.hkl_frame.a_w / lin.norm(self.hkl_frame.a_w)
        ga.y_axis = self.hkl_frame.b_w / lin.norm(self.hkl_frame.b_w)
        ga.z_axis = self.hkl_frame.c_w / lin.norm(self.hkl_frame.c_w)
        ga.focus = self.focus
        ga.heat_limits = self.prop_limits
        ga.heat_palette = self.axis
        ga.histogram = self.histogram
        ga.polar_limits = self.th_limits
        ga.azimuth_limits = self.ph_limits
        ga.plot(self.path + self.GNUPLOT_OUTPUT_EXTENSION)

    @property
    def angle_res(self):
        if self.output_quality not in {1, 2, 3, 4, 5}:
            raise KeyError('output_quality should be 1, 2, 3, 4 or 5')
        return {1: 15, 2: 10, 3: 5, 4: 2, 5: 1}[self.output_quality]

    @property
    def path(self):
        """
        Provides path = directory + stem to the workspace. They are common to
        all input and output files, the extensions are specified as class var.
        """
        return self._path

    @path.setter
    def path(self, path):
        abs_path = Path(make_abspath(path))
        dir_, stem, ext = abs_path.parent, abs_path.stem, abs_path.suffix
        self._path = str(Path().joinpath(dir_, stem))

    @property
    def pg(self):
        return self.sg.reciprocate()

    @property
    def lg(self):
        return self.pg.lauefy()

    @property
    def focus(self):
        if self.orientation is None:
            return []
        else:
            _focus = []
            a = lin.inv(self.orientation) @ np.array((1, 0, 0))
            for op in self.lg.operations:
                v = self.hkl_frame.A_r.T @ op.tf @ a
                c = np.rad2deg(cart2sph(*v))
                if c[1] in self.th_limits and c[2] in self.ph_limits:
                    _focus.append(v / lin.norm(v))
        return _focus

    @property
    def th_range(self):
        return self.th_limits.arange(step=self.angle_res)

    @property
    def ph_range(self):
        return self.ph_limits.arange(step=self.angle_res)

    @property
    def th_mesh(self):
        return self.th_limits.mesh_with(self.ph_limits, step=self.angle_res)[0]

    @property
    def ph_mesh(self):
        return self.th_limits.mesh_with(self.ph_limits, step=self.angle_res)[1]

    def orientation_weight(self, th, ph):
        """Calculate how much each point should contribute to distribution"""
        def sphere_cutout_area(th1, th2, ph_span):
            """Calculate sphere area in specified ph and th degree range.
            For exact math, see articles about spherical cap and sector."""
            return np.deg2rad(abs(ph_span)) * \
                   abs(np.cos(np.deg2rad(th1)) - np.cos(np.deg2rad(th2)))
        th_max = min(th + self.angle_res / 2., self.th_limits[1])
        th_min = max(th - self.angle_res / 2., self.th_limits[0])
        ph_max = min(ph + self.angle_res / 2., self.ph_limits[1])
        ph_min = max(ph - self.angle_res / 2., self.ph_limits[0])
        return sphere_cutout_area(th_min, th_max, ph_max-ph_min)


class AngularPotencyExplorer(AngularPropertyExplorer):
    hkl_is_read_not_generated = False
    property_name = 'cplt'
    property_theoretical_minimum = 0
    property_theoretical_maximum = 1


class AngularR1Explorer(AngularPropertyExplorer):
    hkl_is_read_not_generated = True
    property_name = 'r1'
    property_theoretical_minimum = 0
    property_theoretical_maximum = 1

