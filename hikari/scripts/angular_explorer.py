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

    GNUPLOT_INPUT_EXTENSION = '.gnu'
    GNUPLOT_OUTPUT_EXTENSION = '.pnG'
    HISTOGRAM_EXTENSION = '.his'
    HKL_EXTENSION = '.hkl'
    LISTING_EXTENSION = '.lst'
    MATPLOTLIB_EXTENSION = '.png'
    MESH_EXTENSION = '.dat'

    HKL_IS_READ_NOT_GENERATED = True

    def __init__(self):
        self.gnu_path = ''
        self.gnu_png_path = ''
        self.his_path = ''
        self.hkl_path = ''
        self.lst_path = ''
        self.mpl_png_path = ''
        self.mesh_path = ''
        self.sg = Group()
        self.pg = Group()
        self.lg = Group()
        self.heat_limits = Interval(0, 1)
        self.hkl_frame = HklFrame()
        self.oa = 0.0
        self.orientation = None
        self.resolution = 1.2
        self.th_limits = Interval(0, 0)
        self.ph_limits = Interval(0, 0)
        self.ax = ''
        self.fix_scale = False
        self.histogram = False
        self.output_quality = 1

    def set_path(self, path):
        """
        Set the path to the workspace. The directory and stem must be common
        for all input and output files, the extensions are set automatically.
        :param path: A path to any file with correct stem in working directory.
        :type path: str
        """
        png_path = Path(make_abspath(path))
        dir_, stem, ext = png_path.parent, png_path.stem, png_path.suffix
        base_path = str(Path().joinpath(dir_, stem))
        self.gnu_path = base_path + self.GNUPLOT_INPUT_EXTENSION
        self.gnu_png_path = base_path + self.GNUPLOT_OUTPUT_EXTENSION
        self.his_path = base_path + self.HISTOGRAM_EXTENSION
        self.hkl_path = base_path + self.HKL_EXTENSION
        self.lst_path = base_path + self.LISTING_EXTENSION
        self.mpl_png_path = base_path + self.MATPLOTLIB_EXTENSION
        self.mesh_path = base_path + self.MESH_EXTENSION

    def set_experimental(self, opening_angle, orientation, resolution):
        self.oa = opening_angle
        self.orientation = orientation
        self.resolution = resolution

    def set_hkl_frame(self, a, b, c, al, be, ga, space_group, wavelength, axis):
        self.sg = SG[space_group]
        self.pg = self.sg.reciprocate()
        self.lg = self.pg.lauefy()
        self.hkl_frame.edit_cell(a=a, b=b, c=c, al=al, be=be, ga=ga)
        self.hkl_frame.la = wavelength
        if self.HKL_IS_READ_NOT_GENERATED:
            self._read_hkl_frame(hkl_format='shelx_4')
        else:
            self._make_hkl_frame(ax=axis)
        self.hkl_frame.find_equivalents(point_group=self.pg)
        total_reflections = self.hkl_frame.table['equiv'].nunique()
        if total_reflections == 0:
            raise KeyError('Specified part of reciprocal space has zero nodes')
        self._determine_theta_and_phi_limits()

    def set_options(self, ax, fix_scale, histogram, output_quality):
        self.ax = ''
        self.fix_scale = False
        self.histogram = False
        self.output_quality = output_quality

    def explore(self):
        th_range = self.th_limits.arange(step=self.angle_res)
        ph_range = self.ph_limits.arange(step=self.angle_res)
        th_mesh, ph_mesh = self.th_limits.mesh_with(self.ph_limits,
                                                    step=self.angle_res)
        data_dict = {'th': [], 'ph': [], 'cplt': [], 'r1': [], 'weight': []}

        self.heat_limits = Interval(
            0 if self.fix_scale else min(data_dict['heat']),
            1 if self.fix_scale else max(data_dict['heat']))

    def _make_hkl_frame(self, ax):
        """Make ball or axis of hkl which will be cut in further steps"""
        f = self.hkl_frame
        f.fill(radius=min(self.hkl_frame.r_lim, self.resolution))
        if ax in {'x'}:
            f.table = f.table.loc[f.table['k'].eq(0) & f.table['l'].eq(0)]
        elif ax in {'y'}:
            f.table = f.table.loc[f.table['h'].eq(0) & f.table['l'].eq(0)]
        elif ax in {'z'}:
            f.table = f.table.loc[f.table['h'].eq(0) & f.table['k'].eq(0)]
        elif ax in {'xy'}:
            f.table = f.table.loc[f.table['l'].eq(0)]
        elif ax in {'xz'}:
            f.table = f.table.loc[f.table['k'].eq(0)]
        elif ax in {'yz'}:
            f.table = f.table.loc[f.table['h'].eq(0)]
        if ax in {'x', 'y', 'z', 'xy', 'xz', 'yz'}:
            f.transform([o.tf for o in self.pg.operations])
        f.extinct(self.sg)
        return f

    def _read_hkl_frame(self, hkl_format):
        """Read reflections of hkl which will be cut in further steps"""
        f = self.hkl_frame
        f.read(hkl_path=self.hkl_path, hkl_format=hkl_format)
        f.trim(limit=min(f.r_lim, self.resolution))
        f.find_equivalents(point_group=self.lg)
        return f

    def _determine_theta_and_phi_limits(self):
        """Define range of coordinates where potency map will be calculated."""
        if self.sg.system in {Group.System.triclinic}:
            self.th_limits = Interval(0, 180)
            self.limits = Interval(-45, 135)
        elif self.sg.system in {Group.System.monoclinic}:
            self.th_limits = Interval(0, 180)
            self.ph_limits = Interval(0, 90)
        elif self.sg.system in {Group.System.orthorhombic,
                                Group.System.tetragonal, Group.System.cubic}:
            self.th_limits = Interval(0, 90)
            self.ph_limits = Interval(0, 90)
        elif self.sg.system in {Group.System.trigonal, Group.System.hexagonal}:
            self.th_limits = Interval(0, 90)
            self.ph_limits = Interval(0, 120)
        else:
            raise ValueError('Unknown crystal system (trigonal not supported)')

    def _draw_matplotlib_map(self):
        ma = MatplotlibAngularHeatmapArtist()
        ma.x_axis = self.hkl_frame.a_w / lin.norm(self.hkl_frame.a_w)
        ma.y_axis = self.hkl_frame.b_w / lin.norm(self.hkl_frame.b_w)
        ma.z_axis = self.hkl_frame.c_w / lin.norm(self.hkl_frame.c_w)
        ma.focus = self.focus
        ma.heat_limits = self.heat_limits
        ma.heat_palette = self.axis
        ma.polar_limits = self.th_limits
        ma.azimuth_limits = self.ph_limits
        ma.plot(self.mpl_png_path)

    def _draw_gnuplot_map(self):
        ga = GnuplotAngularHeatmapArtist()
        ga.x_axis = self.hkl_frame.a_w / lin.norm(self.hkl_frame.a_w)
        ga.y_axis = self.hkl_frame.b_w / lin.norm(self.hkl_frame.b_w)
        ga.z_axis = self.hkl_frame.c_w / lin.norm(self.hkl_frame.c_w)
        ga.focus = self.focus
        ga.heat_limits = self.heat_limits
        ga.heat_palette = self.axis
        ga.histogram = self.histogram
        ga.polar_limits = self.th_limits
        ga.azimuth_limits = self.ph_limits
        ga.plot(self.gnu_png_path)

    @property
    def angle_res(self):
        if self.output_quality not in {1, 2, 3, 4, 5}:
            raise KeyError('output_quality should be 1, 2, 3, 4 or 5')
        return {1: 15, 2: 10, 3: 5, 4: 2, 5: 1}[self.output_quality]

    @property
    def focus(self):
        if self.orientation is None:
            return []
        else:
            _focus = []
            for op in self.lg.operations:
                v = self.hkl_frame.A_r.T @ op.tf @ lin.inv(self.orientation) \
                    @ np.array((1, 0, 0))
                c = cart2sph(*v)
                if np.rad2deg(c[1]) in self.th_limits and \
                        np.rad2deg(c[2]) in self.ph_limits:
                    _focus.append(v / lin.norm(v))
        return _focus


class AngularPotencyExplorer(AngularPropertyExplorer):
    HKL_IS_READ_NOT_GENERATED = False


class AngularR1Explorer(AngularPropertyExplorer):
    HKL_IS_READ_NOT_GENERATED = True


