"""This file contains tools for making property maps visualised on sphere"""

import abc
import os
import shutil
from pathlib import Path

import numpy as np
from numpy import linalg as lin

from hikari.dataframes import HklFrame, LstFrame
from hikari.symmetry import SG, Group
from hikari.utility import make_abspath, weighted_quantile, sph2cart, cart2sph,\
    Interval, artist_factory


class AngularPropertyExplorerFactory:
    """A factory method for creating angular property explorers."""
    def __init__(self):
        self._explorers = {}

    def register(self, prop, explorer):
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

    hkl_is_read_not_generated: bool
    property_name: str
    property_theoretical_limits: Interval

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
        self._sg = Group()
        self._pg = Group()
        self._lg = Group()
        self._th_limits = Interval(0, 0)
        self._ph_limits = Interval(0, 0)
        self.axis = ''
        self.hkl_frame = HklFrame()
        self.opening_angle = 0.0
        self.orientation = None
        self.resolution = 0.0
        self.axis = ''
        self.fix_scale = False
        self.histogram = False
        self.output_quality = 1
        self.data_dict = {'th': [], 'ph': [], 'potency': [], 'reflns': [],
                          'R1': [], 'weight': []}

    def set_up(self, a, b, c, al, be, ga, space_group, wavelength, axis,
               opening_angle, orientation, resolution,
               path, fix_scale, histogram, output_quality):
        self.opening_angle = opening_angle
        self.orientation = None if orientation is None \
            else np.array(orientation)
        self.resolution = resolution
        self.path = path
        self.fix_scale = fix_scale
        self.histogram = histogram
        self.output_quality = output_quality
        self.sg = SG[space_group]
        self.hkl_frame.edit_cell(a=a, b=b, c=c, al=al, be=be, ga=ga)
        self.hkl_frame.la = wavelength
        self.axis = axis
        if self.hkl_is_read_not_generated:
            self._read_hkl_frame(hkl_format='shelx_4')
        else:
            self._make_hkl_frame()
        self.hkl_frame.find_equivalents(point_group=self.pg)
        total_unique = self.hkl_frame.table['equiv'].nunique()
        if total_unique == 0:
            raise KeyError('Specified part of reciprocal space has zero nodes')

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
            _prop_limits = Interval(min(self.data_dict[self.property_name]),
                                    max(self.data_dict[self.property_name]))
        else:
            _prop_limits = self.property_theoretical_limits
        return _prop_limits

    def write_hist_file(self):
        hist_bins, hist_edges = np.histogram(
            self.data_dict[self.property_name], density=True,
            weights=self.data_dict['weight'], bins=32, range=self.prop_limits)
        hist_bins = hist_bins / sum(hist_bins)
        with open(self.path + self.HISTOGRAM_EXTENSION, 'w+') as h:
            h.write('#   from      to   prob.\n')
            for _f, _t, _p in zip(hist_edges[:-1], hist_edges[1:], hist_bins):
                h.write(f'{_f:8.5f}{_t:8.5f}{_p:8.5f}\n')

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
        return self._th_limits

    @property
    def ph_limits(self):
        """Interval range of azimuth angle where property will be calculated"""
        return self._ph_limits

    def _draw_map(self, artist, extension):
        a = artist
        a.x_axis = self.hkl_frame.a_w / lin.norm(self.hkl_frame.a_w)
        a.y_axis = self.hkl_frame.b_w / lin.norm(self.hkl_frame.b_w)
        a.z_axis = self.hkl_frame.c_w / lin.norm(self.hkl_frame.c_w)
        a.focus = self.focus
        a.heat_limits = self.prop_limits
        a.heat_palette = self.axis
        a.histogram = self.histogram
        a.polar_limits = self.th_limits
        a.azimuth_limits = self.ph_limits
        a.plot(self.path + extension)

    def draw_matplotlib_map(self):
        a = artist_factory.create('matplotlib_angular_heatmap_artist')
        self._draw_map(artist=a, extension=self.MATPLOTLIB_EXTENSION)

    def draw_gnuplot_map(self):
        a = artist_factory.create('gnuplot_angular_heatmap_artist')
        self._draw_map(artist=a, extension=self.GNUPLOT_OUTPUT_EXTENSION)

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
    def sg(self):
        return self._sg

    @sg.setter
    def sg(self, value):
        _pg = value.reciprocate()
        _sg_system = value.system
        self._sg = value
        self._pg = _pg
        self._lg = _pg.lauefy()
        self._th_limits = self.POLAR_LIMIT_DICT[_sg_system]
        self._ph_limits = self.AZIMUTH_LIMIT_DICT[_sg_system]

    @property
    def pg(self):
        return self._pg

    @property
    def lg(self):
        return self._lg

    @property
    def focus(self):
        _focus = []
        if self.orientation is not None:
            if len(self.orientation.shape) == 1:
                a = self.orientation
            elif len(self.orientation.shape) == 2:
                a = lin.inv(self.orientation) @ np.array((1, 0, 0))
            else:
                raise ValueError(f'Unknown orientation: {self.orientation}')
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
    def th_comb(self):
        return self.th_limits.comb_with(self.ph_limits, step=self.angle_res)[0]

    @property
    def ph_comb(self):
        return self.th_limits.comb_with(self.ph_limits, step=self.angle_res)[1]

    @property
    def th_mesh(self):
        return self.th_limits.mesh_with(self.ph_limits, step=self.angle_res)[0]

    @property
    def ph_mesh(self):
        return self.th_limits.mesh_with(self.ph_limits, step=self.angle_res)[1]

    def orientation_weights(self, th, ph):
        """Calculate how much each point should contribute to distribution"""
        def sphere_cutout_area(th1, th2, ph_span):
            """Calculate sphere area in specified ph and th degree range.
            For exact math, see articles about spherical cap and sector."""
            return np.deg2rad(abs(ph_span)) * \
                   abs(np.cos(np.deg2rad(th1)) - np.cos(np.deg2rad(th2)))
        th_max = (th + self.angle_res / 2.).clip(max=self.th_limits[1])
        th_min = (th - self.angle_res / 2.).clip(min=self.th_limits[0])
        ph_max = (ph + self.angle_res / 2.).clip(max=self.ph_limits[1])
        ph_min = (ph - self.angle_res / 2.).clip(min=self.ph_limits[0])
        return sphere_cutout_area(th_min, th_max, ph_max-ph_min)

    @property
    def descriptive_statistics_string(self):
        min_p = min(self.data_dict[self.property_name])
        worst_index = self.data_dict[self.property_name].index(min_p)
        worst_th = self.data_dict['th'][worst_index]
        worst_ph = self.data_dict['ph'][worst_index]

        max_p = max(self.data_dict[self.property_name])
        best_index = self.data_dict[self.property_name].index(max_p)
        best_th = self.data_dict['th'][best_index]
        best_ph = self.data_dict['ph'][best_index]

        avg_p = np.average(self.data_dict[self.property_name],
                           weights=self.data_dict['weight'])
        q1, q2, q3 = weighted_quantile(self.data_dict[self.property_name],
                                       quantiles=[0.25, 0.50, 0.75],
                                       weights=self.data_dict['weight'])

        s = f'# descriptive statistics for {self.property_name}:\n' \
            f'# max ={max_p:8.5f} at th ={best_th :6.1f} ph ={best_ph :6.1f}\n'\
            f'# min ={min_p:8.5f} at th ={worst_th:6.1f} ph ={worst_ph:6.1f}\n'\
            f'# q_1 ={q1   :8.5f}\n' \
            f'# q_2 ={q2   :8.5f}\n' \
            f'# q_3 ={q3   :8.5f}\n' \
            f'# avg ={avg_p:8.5f}\n'
        return s


class AngularPotencyExplorer(AngularPropertyExplorer):
    hkl_is_read_not_generated = False
    property_name = 'potency'
    property_theoretical_limits = Interval(0, 1)

    def explore(self):
        dat_path = self.path + self.MESH_EXTENSION
        lst_path = self.path + self.LISTING_EXTENSION

        potency_mesh = np.zeros_like(self.th_mesh, dtype=float)
        lst = open(lst_path, 'w+')
        lst.write('#     th      ph potency  reflns\n')
        vectors = sph2cart(r=np.ones_like(self.th_comb),
                           p=np.deg2rad(self.th_comb),
                           a=np.deg2rad(self.ph_comb)).T
        weights = self.orientation_weights(th=self.th_comb, ph=self.ph_comb)
        uniques = self.hkl_frame.dacs_count(self.opening_angle, vectors=vectors)
        total_unique = self.hkl_frame.table['equiv'].nunique('')

        for i, th in enumerate(self.th_range):
            for j, ph in enumerate(self.ph_range):
                unique = uniques[j * len(self.th_range) + i]
                weight = weights[j * len(self.th_range) + i]
                potency = unique / total_unique
                self.data_dict['th'].append(th)
                self.data_dict['ph'].append(ph)
                self.data_dict['potency'].append(potency)
                self.data_dict['reflns'].append(unique)
                self.data_dict['weight'].append(weight)
                lst.write(f'{th:8.0f}{ph:8.0f}{potency:8.5f}{unique:8d}\n')
                potency_mesh[j][i] = potency
            lst.write('\n')

        # noinspection PyTypeChecker
        np.savetxt(dat_path, potency_mesh)
        lst.write(self.descriptive_statistics_string)
        lst.close()


class AngularR1Explorer(AngularPropertyExplorer):

    hkl_is_read_not_generated = True
    property_name = 'R1'
    property_theoretical_limits = Interval(0, 1)

    def explore(self):
        dat_path = self.path + self.MESH_EXTENSION
        lst_path = self.path + self.LISTING_EXTENSION
        job_name = Path(self.path).stem

        r1_mesh = np.zeros_like(self.th_mesh, dtype=float)
        lst = open(lst_path, 'w+')
        lst.write('#     th      ph      R1  reflns\n')

        weights = self.orientation_weights(th=self.th_comb, ph=self.ph_comb)
        total_unique = self.hkl_frame.table['equiv'].nunique()

        for i, th in enumerate(self.th_range):
            for j, ph in enumerate(self.ph_range):
                subdir = self.path + f'_th{int(th+.1)}_ph{int(ph+.1)}'
                hkl_path2 = make_abspath(subdir, job_name+'.hkl')
                ins_path2 = make_abspath(subdir, job_name+'.ins')
                lst_path2 = make_abspath(subdir, job_name+'.lst')
                dir_path2 = make_abspath(subdir)
                Path(dir_path2).mkdir()
                shutil.copy(self.path + '.res', ins_path2)
                q = self.hkl_frame.copy()
                q.dac_trim(self.opening_angle,
                           sph2cart(1.0, np.deg2rad(th), np.deg2rad(ph)))
                q.write(hkl_path2, hkl_format='shelx_4')
                unique = q.table['equiv'].nunique()
                weight = weights[j * len(self.th_range) + i]
                os.system(f'cd {dir_path2}; shelxl {job_name}')
                r1 = LstFrame().read_r1(lst_path2)
                potency = unique / total_unique
                self.data_dict['th'].append(th)
                self.data_dict['ph'].append(ph)
                self.data_dict['potency'].append(potency)
                self.data_dict['R1'].append(r1)
                self.data_dict['weight'].append(weight)
                lst.write(f'{th:8.0f}{ph:8.0f}{r1:8.5}{potency:8.5f}\n')
                r1_mesh[j][i] = r1
            lst.write('\n')

        # noinspection PyTypeChecker
        np.savetxt(dat_path, r1_mesh)
        lst.write(self.descriptive_statistics_string)
        lst.close()


angular_property_explorer_factory = AngularPropertyExplorerFactory()
angular_property_explorer_factory.register(
    prop='potency', explorer=AngularPotencyExplorer)
angular_property_explorer_factory.register(
    prop='r1', explorer=AngularR1Explorer)

