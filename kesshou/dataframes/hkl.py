from collections import OrderedDict
from kesshou.dataframes import BaseFrame
from kesshou.utility import cubespace, is2n, is3n, is4n, is6n
from kesshou.symmetry import PG
from pathlib import Path
from typing import Union
import copy
import json
import random
import struct
import sys
import numpy as np
import numpy.linalg as lin
import pandas as pd
import matplotlib.cm
import matplotlib.pyplot as plt


class HklKeys:
    """
    A helper class supporting HklFrame.
    Menages properties and presence of keys
    present in HklFrame's dataframe
    """
    # DEFINE ALL POSSIBLE HKL KEYS
    __h = {
        'default': 0,
        'description': 'Reciprocal lattice index h',
        'imperative': True,
        'dtype': 'int8',
        'reduce_behaviour': 'keep',
        'type': int
    }
    __k = {
        'default': 0,
        'description': 'Reciprocal lattice index k',
        'imperative': True,
        'dtype': 'int8',
        'reduce_behaviour': 'keep',
        'type': int
    }
    __l = {
        'default': 0,
        'description': 'Reciprocal lattice index l',
        'imperative': True,
        'dtype': 'int8',
        'reduce_behaviour': 'keep',
        'type': int
    }
    __F = {
        'default': 1.0,
        'description': 'Structure factor',
        'imperative': False,
        'dtype': 'float64',
        'reduce_behaviour': 'average',
        'type': float
    }
    __I = {
        'default': 1.0,
        'description': 'Intensity',
        'imperative': False,
        'dtype': 'float64',
        'reduce_behaviour': 'average',
        'type': float
    }
    __si = {
        'default': 0.0,
        'description': 'Structure factor/intensity uncertainty',
        'imperative': False,
        'dtype': 'float64',
        'reduce_behaviour': 'average',
        'type': float
    }
    __u = {
        'default': 0.0,
        'description': 'Structure factor/intensity to uncertainty ratio',
        'imperative': False,
        'dtype': 'float64',
        'reduce_behaviour': 'average',
        'type': float
    }
    __b = {
        'default': 1,
        'description': 'Batch / run number',
        'imperative': False,
        'dtype': 'int16',
        'reduce_behaviour': 'discard',
        'type': int
    }
    __c = {
        'default': 1,
        'description': 'crystal or twin number',
        'imperative': False,
        'dtype': 'int16',
        'reduce_behaviour': 'discard',
        'type': int
    }
    __m = {
        'default': 1,
        'description': 'multiplicity',
        'imperative': True,
        'dtype': 'int16',
        'reduce_behaviour': 'add',
        'type': int
    }
    __la = {
        'default': 0.0,
        'description': 'wavelength',
        'imperative': False,
        'dtype': 'float64',
        'reduce_behaviour': 'discard',
        'type': float
    }
    __ph = {
        'default': 0.0,
        'description': 'reflection phase in radians',
        'imperative': False,
        'dtype': 'float64',
        'reduce_behaviour': 'average',
        'type': float
    }
    __r = {
        'default': 0.0,
        'description': 'distance from 000. Divide by two for resol. in A^-1',
        'imperative': False,
        'dtype': 'float32',
        'reduce_behaviour': 'keep',
        'type': float
    }
    __t = {
        'default': 0.0,
        'description': 'absorption weighted path length in centimeters',
        'imperative': False,
        'dtype': 'float64',
        'reduce_behaviour': 'keep',
        'type': float
    }
    __u1 = {
        'default': 0.0,
        'description': 'Direction cosines of a vector (see XD manual for ref.)',
        'imperative': False,
        'dtype': 'float64',
        'reduce_behaviour': 'keep',
        'type': float
    }
    __u2 = {
        'default': 0.0,
        'description': 'Direction cosines of a vector (see XD manual for ref.)',
        'imperative': False,
        'dtype': 'float64',
        'reduce_behaviour': 'keep',
        'type': float
    }
    __u3 = {
        'default': 0.0,
        'description': 'Direction cosines of a vector (see XD manual for ref.)',
        'imperative': False,
        'dtype': 'float64',
        'reduce_behaviour': 'keep',
        'type': float
    }
    __v1 = {
        'default': 0.0,
        'description': 'Direction cosines of a vector (see XD manual for ref.)',
        'imperative': False,
        'dtype': 'float64',
        'reduce_behaviour': 'keep',
        'type': float
    }
    __v2 = {
        'default': 0.0,
        'description': 'Direction cosines of a vector (see XD manual for ref.)',
        'imperative': False,
        'dtype': 'float64',
        'reduce_behaviour': 'keep',
        'type': float
    }
    __v3 = {
        'default': 0.0,
        'description': 'Direction cosines of a vector (see XD manual for ref.)',
        'imperative': False,
        'dtype': 'float64',
        'reduce_behaviour': 'keep',
        'type': float
    }
    __x = {
        'default': 0.0,
        'description': 'reciprocal position vector x value',
        'imperative': False,
        'dtype': 'float32',
        'reduce_behaviour': 'keep',
        'type': float
    }
    __y = {
        'default': 0.0,
        'description': 'reciprocal position vector y value',
        'imperative': False,
        'dtype': 'float32',
        'reduce_behaviour': 'keep',
        'type': float
    }
    __z = {
        'default': 0.0,
        'description': 'reciprocal position vector z value',
        'imperative': False,
        'dtype': 'float32',
        'reduce_behaviour': 'keep',
        'type': float
    }
    __equiv = {
        'default': (0, 0, 0),
        'description': 'tuple with lexicographically first equivalent hkl',
        'imperative': False,
        'dtype': None,
        'reduce_behaviour': 'keep',
        'type': tuple
    }
    defined_keys = {'h', 'k', 'l', 'F', 'I', 'si', 'b', 'm', 'la', 'ph',
                    'u', 'r', 't', 'u1', 'u2', 'u3', 'v1', 'v2', 'v3',
                    'x', 'y', 'z', 'equiv'}

    def __init__(self):
        # DEFINE ALL KNOWN KEYS
        self.imperatives = set()
        for key in HklKeys.defined_keys:
            if self.get_property(key, 'imperative'):
                self.imperatives.add(key)
        self.all = set()
        self.add(self.imperatives)

    def add(self, keys):
        """Add keys from a keys list to the HklKeys handler"""
        for key in keys:
            self.all.add(key)
        self.__refresh()

    def set(self, keys):
        """Set keys from a keys list in the HklKeys handler"""
        self.all = set()
        self.add(self.imperatives)
        for key in keys:
            self.all.add(key)
        self.__refresh()

    def remove(self, keys):
        """Remove keys from a keys list in the HklKeys handler"""
        for key in keys:
            self.all.discard(key)
        self.add(self.imperatives)
        self.__refresh()

    def get_property(self, key, prop):
        """Get desired key property imported from defined dictionary"""
        return getattr(self, '_HklKeys__'+key)[prop]

    def __refresh(self):
        """refresh defined dictionaries based on all keys' set"""

        # SET 'OF_TYPE' DICTIONARY:
        types = set()
        for key in HklKeys.defined_keys:
            typ = self.get_property(key, 'type').__name__
            types.add(typ)
        self.of_type = dict()
        for typ in types:
            self.of_type[typ] = set()
        for key in self.all:
            typ = self.get_property(key, 'type').__name__
            self.of_type[typ].add(key)

        # SET 'REDUCE BEHAVIOUR' DICTIONARY
        behaviours = set()
        for key in HklKeys.defined_keys:
            behaviour = self.get_property(key, 'reduce_behaviour')
            behaviours.add(behaviour)
        self.reduce_behaviour = dict()
        for behaviour in behaviours:
            self.reduce_behaviour[behaviour] = set()
        for key in self.all:
            behaviour = self.get_property(key, 'reduce_behaviour')
            self.reduce_behaviour[behaviour].add(key)


class HklFrame(BaseFrame):
    """
    This is a master object, which utilises all other
    Hkl classes to import, store, manipulate and output
    information about single-crystal diffraction patterns.

    An empty container which stores the diffraction data
    using python's Pandas library and elementary crystal cell data.
    Demanding methods belonging to this class are vectorized,
    providing relatively satisfactory performance and high data capacity.
    HklFrame methods are designed to work in-place, so the work strategy
    is to create a new instance of HklFrame for each reflection dataset,
    manipulate it using methods, eg. :func:`merge` or :func:`trim`,
    and than :func:`duplicate` or :func:`write` if needed.

    The HklFrame always initiates empty and does not accept any arguments.
    It also does not return anything, though some magic methods, such as
    :func:`__len__` and :func:`__add__` are defined.
    """
    def __init__(self):
        """HklFrame constructor"""
        super().__init__()

        self.data = pd.DataFrame()
        """Pandas dataframe containing diffraction data information."""

        self.hkl_limit = 999
        """Highest absolute value of h, k or l index,
        which can be interpreted correctly (default 999)."""

        self.keys = HklKeys()
        """Object managing keys (column names) of :attr:`a`."""

        self.la = 0.71069
        """Wavelength of radiation used in experiment."""

    def __add__(self, other):
        """
        Add magic method. Adds contents of two `self.data` while preserving
        meta-information from the first `HklFrame` object.
        :param other: HklFrame to be added to data
        :type other: HklFrame
        :return: HklFrame with concatenated HklFrame.data pandas dataframes
        :rtype: HklFrame
        """
        _copied = self.duplicate()
        _copied.data = pd.concat([self.data, other.data], ignore_index=True)
        return _copied

    def __len__(self):
        """
        Len magic method, defined as a number of individual reflections
        :return: Number of rows (individual reflections) in `self.data`
        :rtype: int
        """
        return self.data.shape[0]

    def __str__(self):
        """
        Str magic method, provides human-readable representation of data
        :return: Human-readable representation of `self.data`
        :rtype: str
        """
        return self.data.__str__()

    @property
    def r_lim(self):
        """
        Calculate limiting sphere radius based on :ivar:`la`.

        :return: Value of limiting sphere radius in reciprocal Angstrom
        :rtype: float
        """
        return 2 / self.la

    def _condition(self, equation=''):
        """
        This method limits the reflection data based on truth of equation
        :param equation: Equation which must be met for data to be preserved
        :type equation: str
        :return: Dataframe containing only reflections which meet the equation
        :rtype: pd.DataFrame
        """
        equation = equation.lower().replace(' ', '').replace('_', '')
        df = self.data
        # no variables return dataframe with no rows
        if equation == '':
            return df.iloc[0:0]
        # one variable
        if equation == 'h=2n':
            return df.loc[is2n(df['h'])]
        if equation == 'k=2n':
            return df.loc[is2n(df['k'])]
        if equation == 'l=2n':
            return df.loc[is2n(df['l'])]
        if equation == 'h=3n':
            return df.loc[is3n(df['h'])]
        if equation == 'k=3n':
            return df.loc[is3n(df['k'])]
        if equation == 'l=3n':
            return df.loc[is3n(df['l'])]
        if equation == 'h=4n':
            return df.loc[is4n(df['h'])]
        if equation == 'k=4n':
            return df.loc[is4n(df['k'])]
        if equation == 'l=4n':
            return df.loc[is4n(df['l'])]
        if equation == 'h=6n':
            return df.loc[is6n(df['h'])]
        if equation == 'k=6n':
            return df.loc[is6n(df['k'])]
        if equation == 'l=6n':
            return df.loc[is6n(df['l'])]
        # sum of variables
        if equation in ('h+k=2n', 'k+h=2n'):
            return df.loc[is2n(df['h'] + df['k'])]
        if equation in ('h+l=2n', 'l+h=2n'):
            return df.loc[is2n(df['h'] + df['l'])]
        if equation in ('k+l=2n', 'l+k=2n'):
            return df.loc[is2n(df['k'] + df['l'])]
        if equation in ('2h+k=4n', 'k+2h=4n'):
            return df.loc[is4n(df['h'] + df['h'] + df['k'])]
        if equation in ('h+2k=4n', '2k+h=4n'):
            return df.loc[is4n(df['h'] + df['k'] + df['k'])]
        if equation in ('2h+l=4n', 'l+2h=4n'):
            return df.loc[is4n(df['h'] + df['h'] + df['l'])]
        if equation in ('h+2l=4n', '2l+h=4n'):
            return df.loc[is4n(df['h'] + df['l'] + df['l'])]
        if equation in ('2k+l=4n', 'l+2k=4n'):
            return df.loc[is4n(df['k'] + df['k'] + df['l'])]
        if equation in ('k+2l=4n', '2l+k=4n'):
            return df.loc[is4n(df['k'] + df['l'] + df['l'])]
        if equation in ('h+k+l=2n', 'h+l+k=2n', 'k+h+l=2n',
                        'k+l+h=2n', 'l+h+k=2n', 'l+k+h=2n'):
            return df.loc[is2n(df['h'] + df['k'] + df['l'])]
        # mixed sum and multiple variables
        if equation in ('h,k=2n,h+k=4n', 'h+k=4n,h,k=2n',
                        'k,h=2n,h+k=4n', 'h+k=4n,k,h=2n',
                        'h,k=2n,k+h=4n', 'k+h=4n,h,k=2n',
                        'k,h=2n,k+h=4n', 'k+h=4n,k,h=2n'):
            return df.loc[(is2n(df['h'])) & (is2n(df['k'])) &
                          (is4n(df['h'] + df['k']))]
        if equation in ('h,l=2n,h+l=4n', 'h+l=4n,h,l=2n',
                        'l,h=2n,h+l=4n', 'h+l=4n,l,h=2n',
                        'h,l=2n,l+h=4n', 'l+h=4n,h,l=2n',
                        'l,h=2n,l+h=4n', 'l+h=4n,l,h=2n'):
            return df.loc[(is2n(df['h'])) & (is2n(df['l'])) &
                          (is4n(df['h'] + df['l']))]
        if equation in ('k,l=2n,k+l=4n', 'k+l=4n,k,l=2n',
                        'l,k=2n,k+l=4n', 'k+l=4n,l,k=2n',
                        'k,l=2n,l+k=4n', 'l+k=4n,k,l=2n',
                        'l,k=2n,l+k=4n', 'l+k=4n,l,k=2n'):
            return df.loc[(is2n(df['k'])) & (is2n(df['l'])) &
                          (is4n(df['k'] + df['l']))]
        # multiple variables
        if equation in ('h,k=2n', 'k,h=2n'):
            return df.loc[(is2n(df['h'])) & (is2n(df['k']))]
        if equation in ('h,l=2n', 'l,h=2n'):
            return df.loc[(is2n(df['h'])) & (is2n(df['l']))]
        if equation in ('k,l=2n', 'l,k=2n'):
            return df.loc[(is2n(df['k'])) & (is2n(df['l']))]
        if equation in ('h,k,l=2n', 'h,l,k=2n', 'k,h,l=2n',
                        'k,l,h=2n', 'l,h,k=2n', 'l,k,h=2n', ):
            return df.loc[(is2n(df['h'])) & (is2n(df['k'])) & (is2n(df['l']))]
        # multiple sums of variables
        if equation in ('h+k,h+l,k+l=2n', 'k+h,h+l,k+l=2n', 'h+k,l+h,k+l=2n',
                        'h+k,h+l,l+k=2n', 'k+h,l+h,k+l=2n', 'k+h,l+h,l+k=2n'):
            return df.loc[(is2n(df['h'] + df['k'])) &
                          (is2n(df['h'] + df['l'])) &
                          (is2n(df['k'] + df['l']))]
        # raise exception if the equation is unknown
        raise ValueError('Unknown condition equation have been supplied')

    def dac(self, opening_angle=35.0, vector=None):
        """
        Cut from the dataframe all the reflections,
        which lie outside the accessible volume of diamond anvil cell.

        The diamond anvil cell (DAC) accessible volume is described
        using single *opening angle* (0 to 90 degrees, do not mistake with
        double opening angle which takes values from 0 to 180 degrees)
        and crystal orientation. The orientation information can be supplied
        either via specifying crystal orientation in
        :class:`kesshou.dataframes.BaseFrame`, in :attr:`orientation`
        or by providing a *vector*. The *vector* is perpendicular to
        the dac-accesible space traced by the tori.

        In order to see further details about the shape of dac-accessible space
        and orientation matrix / vector please refer to
        *Merrill & Bassett, Review of Scientific Instruments 45, 290 (1974)*
        and *Paciorek et al., Acta Cryst. A55, 543 (1999)*, respectively.

        :param opening_angle: DAC single opening angle in degrees
        :type opening_angle: float
        :param vector: Provides information about orientation of crystal
            relative to DAC. If None, current :attr:`orientation` is used instead.
        :type vector: Tuple[float, float, float]
        """
        # make the abbreviations for opening angle and prepare objects
        oa = np.radians([float(opening_angle)])[0]
        self.extinct('000')
        if not('x' in self.data.columns):
            self.place()
        # calculate normal vector "n" in reciprocal lattice
        if vector is None:
            l_v = np.array((1.0, 0.0, 0.0))             # vec. parallel to beam
            h = np.dot(lin.inv(self.orientation), l_v)  # calculate hkl vector
            n = h[0] * self.a_w + h[1] * self.b_w + h[2] * self.c_w
        else:                                           # calculate xyz* vector
            n = np.array(vector)
        n = 1 / lin.norm(n) * n                        # normalise the n vector
        # extinct reflections further than max radius
        self.trim(limit=self.r_lim * np.sin(oa))
        # remove reflections perpendicular to disc's normal vector
        xyz = self.data.loc[:, ('x', 'y', 'z')].to_numpy()
        self.data = self.data.loc[~np.isclose(xyz @ n, self.data['r'])]
        # calculate reflection's position in 2D disc reference system "m"
        # m1 / m2 is a coordinate parallel / perpendicular to vector n
        xyz = self.data.loc[:, ('x', 'y', 'z')].to_numpy()
        m1 = np.outer(xyz @ n, n)
        m2 = xyz - m1
        m1 = m1 @ n
        m2 = lin.norm(m2, axis=1)
        # find the middle of two tori, which trace the DAC-limiting shape
        t1 = 1 / 2 * self.r_lim * np.cos(oa)
        t2 = 1 / 2 * self.r_lim * np.sin(oa)
        # check if points lie in one of two tori making the DAC shape
        in_torus1 = (m1 - t1) ** 2 + (m2 - t2) ** 2 <= (self.r_lim / 2) ** 2
        in_torus2 = (m1 + t1) ** 2 + (m2 - t2) ** 2 <= (self.r_lim / 2) ** 2
        # leave only points which lie in both tori
        self.data = self.data[in_torus1 * in_torus2]
        self.data.reset_index(drop=True, inplace=True)

    def _domain(self, address='hkl'):
        """
        This method limits the reflection data to the ones living in address
        :param address: Address which the reflections must have to be preserved
        :type address: str
        :return: Dataframe containing only reflections with given address
        :rtype: pd.DataFrame
        """

        address = address.lower().replace(' ', '').replace('_', '')
        df = self.data
        # no zeroes in address
        if address == 'hkl':
            return df
        if address in ('hhl', 'kkl'):
            return df.loc[df['h'] == df['k']]
        if address in ('hkh', 'lkl'):
            return df.loc[df['h'] == df['l']]
        if address in ('hkk', 'hll'):
            return df.loc[df['k'] == df['l']]
        # one zero in address
        if address == 'hk0':
            return df.loc[df['l'] == 0]
        if address == 'h0l':
            return df.loc[df['k'] == 0]
        if address == '0kl':
            return df.loc[df['h'] == 0]
        if address in ('hh0', 'kk0'):
            return df.loc[(df['h'] == df['k']) & (df['l'] == 0)]
        if address in ('h0h', 'l0l'):
            return df.loc[(df['h'] == df['l']) & (df['k'] == 0)]
        if address in ('0kk', '0ll'):
            return df.loc[(df['k'] == df['l']) & (df['h'] == 0)]
        # two zeroes in address
        if address == 'h00':
            return df.loc[(df['k'] == 0) & (df['l'] == 0)]
        if address == '0k0':
            return df.loc[(df['h'] == 0) & (df['l'] == 0)]
        if address == '00l':
            return df.loc[(df['h'] == 0) & (df['k'] == 0)]
        # three zeroes in address
        if address == '000':
            return df.loc[(df['h'] == 0) & (df['k'] == 0) & (df['l'] == 0)]
        # raise exception if the address is unknown
        raise ValueError('Unknown domain address have been supplied')

    def draw(self, alpha=False, colored='b', dpi=600, legend=True,
             master_key='I', projection=('h', 'k', 0),
             savepath=False, scale=1.0, showfig=False):
        """
        Draw a cross-section of reciprocal lattice for given pattern

        alpha        (string)   Value to be visualised as alpha or False
        color        (string)   Int value represented as colour or False
        dpi          (integer)  Dots Per Inch, quality of saved graphics
        legend       (boolean)  Legend of used colors should be printed
        master_key   (string)   Value to be visualised as radius, 'I' or 'F'
        projection   (tuple)    desired cross-section, default ('h', 'k', 0)
        savepath     (string)   Path to the file to save the image or False
        scale        (float)    Scale factor for the reflection size
        showfig      (boolean)  Figure should be shown in matplotlib window
        uncertainty  (boolean)  master/sigma(master) drawn as transparency

        This object should be deleted or cleared for the sake of release version
        """
        # TODO clear of delete

        # SET NECESSARY PARAMETERS
        color_scheme = 'gist_rainbow'
        distance = 1
        axes = list()
        self.extinct('000')
        fig = plt.figure()

        # INTERPRET PROJECTION
        for value in projection:
            try:
                distance = int(value)
            except ValueError:
                axes.append(value)
        direction = ({'h', 'k', 'l'} - set(axes)).pop()
        title = str(projection).replace(', ', '').replace('\'', '')
        ax = fig.add_subplot(1, 1, 1)

        # PREPARE MINIMA AND MAXIMA FOR SCALING PURPOSES
        data_minima, data_maxima = dict(), dict()
        for key in list(self.keys.all):
            data_minima[key] = self.data[key].min()
            data_maxima[key] = self.data[key].max()

        # PREPARE COLOUR PALETTE
        color_map = matplotlib.cm.get_cmap(color_scheme)
        color_range = 1
        if colored:
            color_range += int(data_maxima[colored] - data_minima[colored])
        colors, _rgb, _alp = [], (1.0, 0.0, 0.0), (1,)
        [colors.append(color_map(i / color_range)) for i in range(color_range)]

        # PREPARE NECESSARY LISTS FOR MATPLOTLIB
        _x, _y, _size, _color, _edge = list(), list(), list(), list(), list()
        actual_index = -1
        for index, row in self.data.iterrows():
            if row[direction] != distance:
                continue
            actual_index += 1
            coordinates = {'h': 'x', 'k': 'y', 'l': 'z'}
            _x.append(row[coordinates[axes[0]]])
            _y.append(row[coordinates[axes[1]]])
            _size.append(scale ** 2 * np.log(abs(row[master_key])+1.0) ** 2)
            if colored:
                _rgb = colors[int(row[colored] - data_minima[colored])][:3]
            if alpha:
                _alp = ((row[alpha]/data_maxima[alpha])**0.25, )
            _color.append(_rgb + _alp)
            _edge.append('None') if row[master_key] > 0 else _edge.append('k')

        # CHECKING IF LIST IS NOT EMPTY
        if len(_color) == 0:
            print("Reflection set intended to be drawn is empty! Aborting")
            return

        # DRAW THE PLOT
        directions = {'h': 'a* [A^-1]', 'k': 'b* [A^-1]', 'l': 'c* [A^-1]'}
        ax.set_title(title)
        ax.set_xlabel(directions[axes[0]])
        ax.set_ylabel(directions[axes[1]])
        ax.scatter(0, 0, s=20, c='k', marker='x')
        ax.scatter(_x, _y, s=_size, c=_color, marker='.',
                   edgecolors=_edge, linewidth=[0.05 * s for s in _size])
        ax.axis('equal')

        # ADD LEGEND IF APPLICABLE
        if legend and colored:
            legend_colors, legend_batches = [], []
            color_skip = int(color_range / 25) + 1
            for i in range(0, color_range, color_skip):
                legend_colors.append(plt.Rectangle((0, 0), 1, 1, fc=colors[i]))
                legend_batches.append(str(i + data_minima[colored]))
            # if color_range < 10:
            ax.legend(legend_colors, legend_batches, loc=1, prop={'size': 7})
            # elif color_range < 25:
            #     ax.legend(legend_colors, legend_batches, loc=1, ncol=2,
            #               prop={'size': 6})
            # else:
            #     ax.legend(legend_colors, legend_batches, loc=1, ncol=3,
            #               prop={'size': 5})

        # SAVE OR SHOW FIGURE IF APPLICABLE
        fig = plt.gcf()
        if savepath:
            fig.savefig(savepath, bbox_inches=None, dpi=dpi)
        if showfig:
            plt.show()

    def duplicate(self):
        """
        Make and return an exact deep copy of this HklFrame.

        :return: A copy of this HklFrame.
        :rtype: HklFrame
        """
        return copy.deepcopy(self)

    def edit_wavelength(self, wavelength):
        """
        Define wavelength of radiation used in the diffraction experiment
        for the sake of other methods.

        Alternative to manually setting `HklFrame.la`,
        but can accept and interpret string values such as "MoKa".
        Interprets the definitions of the alpha and beta1 radiation of: "Ag",
        "Co", "Cr", "Cu", "Fe", "Mn", "Mo", "Ni", "Pd", "Rh", "Ti" and "Zn".
        The values of wavelengths have been imported from International Tables
        of Crystallography, Volume C, Table 4.2.4.1, 3rd Edition.

        :param wavelength: Wavelength of the radiation used in the experiment.
        :type wavelength: str or float
        """
        characteristic_radiations = {'agka': 0.5608, 'agkb': 0.4970,
                                     'coka': 1.7905, 'cokb': 1.6208,
                                     'crka': 2.2909, 'crkb': 2.0848,
                                     'cuka': 1.5418, 'cukb': 1.3922,
                                     'feka': 1.9373, 'fekb': 1.7565,
                                     'mnka': 2.1031, 'mnkb': 1.9102,
                                     'moka': 0.7107, 'mokb': 0.6323,
                                     'nika': 1.6591, 'nikb': 1.5001,
                                     'pdka': 0.5869, 'pdkb': 0.5205,
                                     'rhka': 0.6147, 'rhkb': 0.5456,
                                     'tika': 2.7496, 'tikb': 2.5138,
                                     'znka': 1.4364, 'znkb': 1.2952}
        try:
            self.la = characteristic_radiations[wavelength[:4].lower()]
        except TypeError:
            self.la = float(wavelength)

    def extinct(self, rule='hkl:', point_group=PG['1']):
        """
        Removes from dataframe all reflections which are in a specified
        *domain*, but do not meet the *condition*.
        The *rules* have a format "domain:condition*,
        where the domain specifies a part of reciprocal space which
        we plan to extinct, while the condition describes all the reflections,
        which should **not** be extinct despite belonging to domain.

        The rule should be written using a format present in
        International Tables of Crystallography,
        which is easily accessible using the web page of
        `Department of Chemistry, University College London
        <http://img.chem.ucl.ac.uk/sgp/large/sgp.htm>`_
        In the current release the *domains* are hard-coded and thus not all
        theoretically possible domains can be specified.
        However, all domains specified in the International Tables
        should be accesible.

        :param rule: A string with *domain* and *condition* separated by colon.
        :type rule: str
        :param point_group: Point group containing information about present
            symmetry operations, necessary to correctly apply the extinction.
        :type point_group: kesshou.symmetry.pointgroup.PointGroup
        """

        def _interpret_rule():
            try:
                _split_rule = [*rule.split(':', 1), '']
            except AttributeError:
                _split_rule = rule[0:2]
            _domain_string = _split_rule[0].strip(' ')
            _condition_string = _split_rule[1].strip(' ')
            return _domain_string, _condition_string
        domain_string, condition_string = _interpret_rule()
        dom = self._domain(domain_string)
        con = self._condition(condition_string)
        # obtain a part of domain which does NOT meet the condition
        dom = dom.loc[~dom.isin(con).all(1)][['h', 'k', 'l']]
        # make a set with all hkls which should be extinct
        extinct = set()
        for op in point_group.operations:
            extinct = extinct | {tuple(row) for row in dom.to_numpy() @ op}
        # remove all reflections whose hkls are in extinct set
        for h, k, l in list(extinct):
            self.data = self.data[~((self.data['h'] == h) &
                                    (self.data['k'] == k) &
                                    (self.data['l'] == l))]
        # reset the indices for other methods to use
        self.data.reset_index(drop=True, inplace=True)

    def find_equivalents(self, point_group=PG['1']):
        """
        Assign each reflection its symmetry equivalent with highest values
        of h, k, l indices and store it in "equiv" column in the dataframe.

        In order to provide an information about equivalence, a *point_group*
        must be provided (default PG1). Point groups and their notation can
        be found within :mod:`kesshou.symmetry` sub-package.

        :param point_group: Point Group used to determine symmetry equivalence
        :type point_group: kesshou.symmetry.PointGroup
        """
        hkl_lim = self.hkl_limit
        self.keys.add(('equiv',))
        self.data['equiv'] = [(-hkl_lim, -hkl_lim, -hkl_lim)] * len(self.data)
        _hkl_matrix = self.data.loc[:, ('h', 'k', 'l')].to_numpy()
        for op in point_group.operations:
            new_hkl = pd.Series(map(tuple, _hkl_matrix @ op[0:3, 0:3]))
            _to_update = self.data['equiv'] < new_hkl
            self.data.loc[_to_update, 'equiv'] = new_hkl.loc[_to_update]

    def from_dict(self, dictionary):
        """
        Construct the self.data using information stored in dictionary.
        The dictionary keys must be valid strings, see :class:`HklKeys` for
        a list of valid keys. The dictionary values must be iterable of equal
        size, preferably `numpy.ndarray`.

        :param dictionary: Dictionary with "key - iterable of values" pairs.
        :type dictionary: Dict[str, numpy.ndarray]
        """
        new_data = pd.DataFrame()
        self.keys.add(dictionary.keys())
        for key, value in dictionary.items():
            typ = self.keys.get_property(key, 'dtype')
            new_data[key] = pd.Series(value, dtype=typ, name=key)
        self.data = new_data
        self.place()

    @staticmethod
    def interpret_hkl_format(hkl_format):
        """
        Interpret hkl format and return format strings, column labels etc.
        necessary for other parsers: readers and writers.

        The hkl_format might be integer, string or ordered dictionary.
        Available integer-type input formats consist of:

        - 2: Standard hkl2 format containing h, k, l (4 digits each),
          square of structure factor, its uncertainty (8 digits each),
          batch number (4 digits) and wavelength (8 digits).

        - 3: Standard hkl3 format containing h, k, l (4 digits each),
          structure factor, its uncertainty (8 digits each)
          and batch number (4 digits).

        - 4: Standard hkl4 format containing h, k, l (4 digits each),
          square of structure factor, its uncertainty (8 digits each)
          and batch number (4 digits).

        - 40: Modified hkl4 format containing h, k, l (4 digits each),
          square of structure factor and its uncertainty (8 digits each).

        - 5: Standard hkl5 format containing h, k, l (4 digits each),
          square of structure factor, its uncertainty (8 digits each)
          and crystal number (4 digits).

        - 6: Standard hkl6 format containing h, k, l (4 digits each),
          square of structure factor, its uncertainty (8 digits each)
          and multiplicity (4 digits).

        Available string-type input formats consist of:

        - 'xd': hkl format accepted by program "XD", containing
          h, k, l, batch number (5 digits each), square of structure factor
          and its uncertainty (10 digits each).

        - 'tonto': hkl format accepted by program "tonto", containing
          square of structure factor and its uncertainty (8 digits each),
          as well as relevant prefix and suffix lines.

        - 'free': hkl space-separated free format containing
          h, k, l, square of structure factor, its uncertainty and batch number.

        Available ordered dictionary-type input should contains *key-value*
        pairs, where key is a short string name of accepted information
        (eg. "h", "I" or "si"; for full list please refer to :class:`HklKeys`)
        and value is a length of given field in the hkl file.

        This function returns a tuple containing five objects in this order:

        - column_labels:: a tuple of keys - labels of subsequently read data

        - format_string:: string used by other methods to read/write data

        - file_prefix:: None or string which appears on beginning of the file

        - file_suffix:: None or string which appears on end of the file

        - zero_line:: True or False; should the file contain the 0, 0, 0 line.

        :param hkl_format: Format of the hkl file to be read or written
        :type hkl_format: int or str or OrderedDict
        :return: A tuple of format characteristics of the read/written file
        :rtype: tuple
        """
        if hkl_format == 2:
            column_labels = ('h', 'k', 'l', 'I', 'si', 'b', 'la')
            format_string = '4s 4s 4s 8s 8s 4s 8s'
            file_prefix = False
            file_suffix = False
            zero_line = True
        elif hkl_format == 3:
            column_labels = ('h', 'k', 'l', 'F', 'si', 'b')
            format_string = '4s 4s 4s 8s 8s 4s'
            file_prefix = False
            file_suffix = False
            zero_line = True
        elif hkl_format == 4:
            column_labels = ('h', 'k', 'l', 'I', 'si', 'b')
            format_string = '4s 4s 4s 8s 8s 4s'
            file_prefix = False
            file_suffix = False
            zero_line = True
        elif hkl_format == 40:
            column_labels = ('h', 'k', 'l', 'I', 'si')
            format_string = '4s 4s 4s 8s 8s'
            file_prefix = False
            file_suffix = False
            zero_line = True
        elif hkl_format == 5:
            column_labels = ('h', 'k', 'l', 'I', 'si', 'c')
            format_string = '4s 4s 4s 8s 8s 4s'
            file_prefix = False
            file_suffix = False
            zero_line = True
        elif hkl_format == 6:
            column_labels = ('h', 'k', 'l', 'I', 'si', 'm')
            format_string = '4s 4s 4s 8s 8s 4s'
            file_prefix = False
            file_suffix = False
            zero_line = True
        elif hkl_format in ('xd', 'Xd', 'xD', 'XD'):
            column_labels, format_string = tuple(), 'XD'
            file_prefix = False
            file_suffix = False
            zero_line = True
        elif hkl_format in ('tonto', 'Tonto', 'TONTO', 'HAR', 'har', 'Har'):
            column_labels = ('h', 'k', 'l', 'I', 'si')
            format_string = '4s 4s 4s 8s 8s'
            file_prefix = 'reflection_data= {\n' \
                          'keys= { h= k= l= i_exp= i_sigma= }\n' \
                          'data= {'
            file_suffix = '}\n' \
                          '}\n' \
                          'REVERT'
            zero_line = False
        elif hkl_format in {'FREE', 'Free', 'free'}:
            column_labels = ('h', 'k', 'l', 'I', 'si', 'b')
            format_string = '4s 4s 4s 8s 8s 4s'
            file_prefix = False
            file_suffix = False
            zero_line = True
        elif type(hkl_format) in (dict, OrderedDict):
            column_labels = list()
            format_string = str()
            file_prefix = 'COMPOUND_ID          F^2  NDAT 6'
            file_suffix = False
            zero_line = True
            if type(hkl_format) is dict and sys.version_info[0] < 3:
                format_items = hkl_format.iteritems()
            else:
                format_items = hkl_format.items()
            for key, value in format_items:
                column_labels.append(key)
                if int(value) > 0:
                    format_string += str(value) + 's '
                else:
                    format_string += str(abs(int(value))) + 'x '
            column_labels = tuple(column_labels)
            format_string.rstrip(' ')
        else:
            raise TypeError('Format type should be 2, 3, 4, 40, 5, 6,'
                            '"XD", "TONTO", "free" or dict')
        return format_string, column_labels, file_prefix, file_suffix, zero_line

    def make_ball(self, radius=2.0):
        """
        Fill the dataframe with all possible reflections within
        the distance *radius* from the reciprocal space origin.

        :param radius: Maximum distance from the reciprocal space origin
        to placed reflection (in reciprocal Angstrom).
        :type radius: float
        """

        max_index = 25

        # make an initial guess of the hkl ball
        def _make_hkl_ball(i=max_index, _r=radius):
            hkl_grid = np.mgrid[-i:i:2j * i + 1j, -i:i:2j * i + 1j,
                       -i:i:2j * i + 1j]
            h_column = np.concatenate(np.concatenate(hkl_grid[0]))
            k_column = np.concatenate(np.concatenate(hkl_grid[1]))
            l_column = np.concatenate(np.concatenate(hkl_grid[2]))
            hkls = np.matrix((h_column, k_column, l_column)).T
            xyz = np.matrix((self.a_w, self.b_w, self.c_w))
            return hkls[lin.norm(hkls @ xyz, axis=1) <= _r]

        hkl = _make_hkl_ball()

        # increase the ball size until all needed points are in
        previous_length = -1
        while len(hkl) > previous_length and max_index <= self.hkl_limit:
            previous_length = len(hkl)
            max_index = max_index * 2
            hkl = _make_hkl_ball(max_index)

        # create new dataframe using obtained ball of data
        _h, _k, _l = np.vsplit(hkl.T, 3)
        ones = np.ones_like(np.array(_h)[0])
        self.from_dict({'h': np.array(_h)[0], 'k': np.array(_k)[0],
                        'l': np.array(_l)[0], 'I': ones, 'si': ones, 'm': ones})

    def make_stats(self, bins=10, point_group=PG['1'], extinctions=('000',)):
        """
        Analyses dataframe in terms of number of individual, unique and
        theoretically possible reflections, as well as completeness and
        redundancy in given point group.

        Point group is necessary to correctly specify symmetry equivalence
        and it has been described for method :meth:`merge` in detail.
        List of extinctions is necessary to correctly extinct the reference data
        and it has been for method :meth:`extinct` in detail.

        :param bins: Number of individual bins to divide the data into.
        :type bins: int
        :param point_group: Point group used to calculate the statistics.
        :type point_group: kesshou.symmetry.PointGroup
        :param extinctions: Iterable of extinction rules to be applied
        :type extinctions: tuple
        """

        def prepare_base_copy_of_hkl():
            _hkl_base = self.duplicate()
            for _extinction in extinctions:
                _hkl_base.extinct(_extinction)
            return _hkl_base

        hkl_base = prepare_base_copy_of_hkl()

        def prepare_ball_of_hkl(_point_group=PG['1']):
            _hkl_full = hkl_base.duplicate()
            _hkl_full.make_ball(radius=max(self.data['r']))
            _hkl_full.merge(point_group=_point_group)
            for _extinction in extinctions:
                _hkl_full.extinct(_extinction)
            return _hkl_full

        hkl_full = prepare_ball_of_hkl(_point_group=point_group)

        def prepare_merged_hkl(_point_group=PG['1']):
            _hkl_merged_pg1 = hkl_base.duplicate()
            _hkl_merged_pg1.merge(point_group=_point_group)
            return _hkl_merged_pg1

        hkl_merged = prepare_merged_hkl(_point_group=point_group)

        def group_by_resolution(ungrouped_hkl, _bins=bins):
            cube_bins = cubespace(0.0, max(self.data['r']), num=_bins + 1)
            grouped_hkl = ungrouped_hkl.data.groupby(
                pd.cut(ungrouped_hkl.data['r'], cube_bins))
            return grouped_hkl

        grouped_base = group_by_resolution(hkl_base)
        grouped_full = group_by_resolution(hkl_full)
        grouped_merged = group_by_resolution(hkl_merged)

        def make_table_with_stats(_grouped_base, _grouped_full,
                                  _grouped_merged):
            observed = _grouped_base.size()
            independent = _grouped_merged.size()
            theory = _grouped_full.size()
            completeness = independent.div(theory)
            redundancy = observed.div(independent)
            results = pd.concat([observed, independent, theory,
                                 completeness, redundancy], axis=1)
            results.columns = ['Obser', 'Indep', 'Theory', 'Cplt', 'Redund.']
            print(results)

        print('\nStatistics in Selected Point Group')
        make_table_with_stats(grouped_base, grouped_full, grouped_merged)

    def merge(self, point_group=PG['1']):
        """
        Average down each set of redundant reflections present in the dataframe,
        to one reflection.

        The redundancy of reflections is determined using the
        :meth:`find_equivalents` method with appropriate point group.
        Therefore, the merging can be used in different ways depending on given
        point group:

        - For PG['1'], only reflections with exactly the same values
          of h, k and l indices will be merged. Resulting dataframe will not
          contain a duplicate of any reflection.

        - For PG['-1'] the reflections with the same values of h, k and l,
          as well as their Friedel pairs, will be merged. Resulting file will
          be devoid of duplicate reflections as well as any Friedel pairs.

        - For PG['mmm'] all symmetry-equivalent reflections within the "mmm"
          point group will be merged. Please mind that "mmm" point group is
          centrosymmetric, so the Friedel pairs will be merged as well.

        - For PG['mm2'] symmetry-equivalent reflections within the "mmm"
          point group will be merged, but the Friedel pairs will be preserved,
          as the 'mm2' is a inversion-devoid sub-group of the "mmm" point group.

        The procedure will have a different effect on different dataframe keys,
        depending on their "reduce_behaviour" specified in :class:`HklKeys`:

        - Parameters which should be preserved will be kept intact:

            - index *h, k, l*,

            - position in reciprocal space *x, y, z, r* (see :meth:`place`)

            - symmetry equivalence *equiv* (see :meth:`find_equivalents`)

        - Parameters such as intensity *I*, structure factor *F* and
        their uncertainty *si* will be averaged using arithmetic mean.

        - Multiplicity of occurrence *m* will be summed

        - Other parameters which lose their meaning during the merging
        procedure such as batch number *b* will lose their meaning
        and thus will be discarded.

        The merging inevitably removes some information from the dataframe,
        but it can be necessary for some operations. For example, the drawing
        procedures work much better and provide clearer image if multiple points
        occupying the same position in space are reduced to one instance.

        In order to provide an information about equivalence
        for the sake of merging, a *point_group* must be provided (default PG1).
        Point groups and their notation can
        be found within :mod:`kesshou.symmetry` sub-package.
        Please mind that because the symmetry equivalence information is used
        during the merging procedure to determine which points should be merged,
        previously defined values in "equiv" will be overwritten by this method.

        :param point_group: Point Group used to determine symmetry equivalence
        :type point_group: kesshou.symmetry.PointGroup
        """
        self.find_equivalents(point_group=point_group)
        # group the dataframe and obtain all existing keys
        grouped = self.data.groupby('equiv')
        grouped_first = grouped.first().reset_index()
        grouped_mean = grouped.mean().reset_index()
        grouped_sum = grouped.sum().reset_index()
        # for each key apply a necessary reduce operation and add it to data
        data = dict()
        for key in self.data.keys():
            if key in self.keys.reduce_behaviour['keep']:
                data[key] = grouped_first[key]
            elif key in self.keys.reduce_behaviour['add']:
                data[key] = grouped_sum[key]
            elif key in self.keys.reduce_behaviour['average']:
                data[key] = grouped_mean[key]
        self.from_dict(data)

    def place(self):
        """
        Assign reflections their positions in reciprocal space ("x", "y", "z")
        and calculate their distance from origin ("r") in reciprocal Angstrom.
        Save four new keys and their values into the dataframe.
        """
        hkl = self.data.loc[:, ('h', 'k', 'l')].to_numpy()
        abc = np.matrix((self.a_w, self.b_w, self.c_w))
        xyz = hkl @ abc
        self.data['x'] = xyz[:, 0]
        self.data['y'] = xyz[:, 1]
        self.data['z'] = xyz[:, 2]
        self.data['r'] = lin.norm(xyz, axis=1)

    def read(self, hkl_path, hkl_format):
        """
        Read the contents of .hkl file as specified by path and format,
        and store them in the pandas dataframe in `self.data`.
        For a list of all available .hkl formats,
        please refer to the documentation of :func:`interpret_hkl_format`.

        :param hkl_path: Absolute or relative path to the .hkl file.
        :type hkl_path: str
        :param hkl_format: Format of provided .hkl file.
        :type hkl_format: union[int, str, OrderedDict]
        """

        # PREPARE OBJECTS RESPONSIBLE FOR PARSING INPUT
        format_string, column_labels, file_prefix, file_suffix, zero_line = \
            self.interpret_hkl_format(hkl_format)
        self.keys.add(column_labels)

        if format_string is 'XD':
            def parse_line(hkl_line):
                return tuple(hkl_line.strip().split())
        else:
            def parse_line(hkl_line):
                field_struct = struct.Struct(format_string)
                unpack = field_struct.unpack_from
                return tuple(s.decode() for s in unpack(hkl_line.encode()))

        # if free format
        if str(format_string).lower() == 'free':
            self.data = pd.read_csv(filepath_or_buffer=hkl_path, sep=' ',
                                    names=column_labels, header=False)
            return

        # OPEN HKL FILE AND PREPARE CONTAINER
        hkl_file = open(hkl_path, 'r')
        hkl_content = dict()

        # LOAD HKL COLUMN TAGS IF SUPERTYPE IS "XD"
        if format_string is 'XD':
            title_line = hkl_file.readline().strip().split()
            hkl_main_key = 'F' if title_line[1] == 'F' else 'I'
            if title_line[2] != 'NDAT':
                raise KeyError('Loaded hkl file is not of "XD" type.')
            if int(title_line[3]) == -7:
                column_labels = ('h', 'k', 'l', 'b', hkl_main_key, 'si', 'ph')
            elif 6 <= int(title_line[3]) <= 13:
                column_labels = ('h', 'k', 'l', 'b', hkl_main_key, 'si',
                                 't', 'u1', 'u2', 'u3', 'v1', 'v2', 'v3'
                                 )[0:int(title_line[3])]
            else:
                raise KeyError('"NDAT" parameter of loaded hkl file lies'
                               'outside of the expected range')

        # LOAD THE KEYS
        for key in column_labels:
            hkl_content[key] = []
        self.keys.set(column_labels)

        # INTERPRET THE FILE
        hkl_checksum = 0
        for line in hkl_file:
            # IGNORE EMPTY LINES
            if not line.strip():
                continue
            # SAVE EACH LINE TO THE LIST
            hkl_checksum += 1
            for key, value in zip(column_labels, parse_line(line)):
                value = self.keys.get_property(key, 'type')(value)
                hkl_content[key].append(value)
        hkl_file.close()

        # IF IMPERATIVE DATA WAS NOT GIVEN, ADD DEFAULTS
        forgotten_keys = tuple(self.keys.imperatives - set(column_labels))
        for forgotten in forgotten_keys:
            default = self.keys.get_property(forgotten, 'default')
            hkl_content[forgotten] = [default] * hkl_checksum

        # PRODUCE PANDAS DATAFRAME
        self.from_dict(hkl_content)

    def rescale(self, key, factor):
        """
        Multiply all values stored in a column *key* in the dataframe
        by a *factor*. Accepts a single key and a single rescale factor.
        :param key: Column of dataframe to be rescaled by a factor
        :type key: str
        :param factor: A number to multiply all values in column *key* by.
        :type factor: float
        """
        self.data[key] = self.data.apply(lambda row: row[key] * factor, axis=1)

    def transform(self, operations):
        """

        Apply a symmetry operation or list of symmetry operations to transform
        the diffraction pattern.

        If one symmetry operation (3x3 or 4x4 numpy array / matrix) is provided,
        it effectively multiplies the hkl matrix by the operation matrix
        and accordingly alters the self.data dataframe. As a result,
        the length of self.data before and after transformation is the same.

        However, the function behaves slightly counter-intuitively
        if two or more operation matrices are provided. In such case
        the method applies the transformation procedure independently
        for each operation, and then *concatenates* resulting matrices.
        Resulting self.data is len(operations) times longer than the initial.

        The function can use 3x3 or larger (eg. 4x4) matrices, as it selects
        only the upper-left 3x3 segment for the sake of calculations.
        Also, while reconstructing the symmetry of merged reflection file
        it is important to use all symmetry operations, not only generators.

        Single symmetry operations or their lists belonging to certain
        point groups can be imported from :py:mod:`kesshou.symmetry` module.

        :param operations: Iterable of operation matrices to be applied
        :type operations: Union[Tuple[np.ndarray, np.ndarray], np.ndarray]
        """

        def _make_a_list_of_operations():
            _ops = np.array(operations)
            if len(_ops.shape) == 1:
                raise ValueError('Operations should be a list of matrices')
            elif len(_ops.shape) == 2 and _ops.shape[0] == _ops.shape[1]:
                return [_ops]
            elif len(_ops.shape) == 3 and _ops.shape[1] == _ops.shape[2]:
                return operations
            else:
                raise ValueError('Operations should be a list of matrices')
        ops = _make_a_list_of_operations()

        def _build_transformed_dataframe():
            dataframes = list()
            for op in ops:
                df = copy.deepcopy(self.data)
                mat = op[0:3, 0:3]
                hkl = self.data.loc[:, ['h', 'k', 'l']].to_numpy() @ mat
                df['h'] = hkl[:, 0]
                df['k'] = hkl[:, 1]
                df['l'] = hkl[:, 2]
                dataframes.append(df)
            return pd.concat(dataframes, axis=0, ignore_index=True)
        self.data = _build_transformed_dataframe()
        self.place()

    def thin_out(self, target_cplt=1.0):
        """
        Randomly remove reflections from dataframe in order to decrease
        the completeness to *target_cplt* (relatively to initial completeness).
        :param target_cplt: Percentage of data not removed from dataframe
        :type target_cplt: float
        """
        assert 0.0 <= target_cplt <= 1.0, 'target_cplt outside the 0--1 range'
        number_to_delete = int((1 - target_cplt) * len(self))
        indices_to_delete = random.sample(range(0, len(self)), number_to_delete)
        self.data.drop(indices_to_delete, inplace=True)
        self.data.reset_index(drop=True, inplace=True)

    def to_hklres(self, colored='m', master_key='I', path='hkl.res'):
        """
        Export the reflection information from dataframe to .res file,
        so that a software used to visualize .res files can be used
        to visualise a diffraction data in three dimentions.

        :param colored: Which dataframe key, "m" or "b", should be used to
        provide the color to the reflection in final image.
        :type colored: str
        :param master_key: Which dataframe key, "I" or "F", should be used to
        provide the intensity of the reflection.
        :type master_key: str
        :param path: Absolute or relative path where the file should be saved
        :type path: str
        """
        # TODO mercury accepts max of Uiso == 4.999 - fix to this value
        # prepare minima and maxima for scaling purposes
        data_minima, data_maxima = dict(), dict()
        for key in list(self.keys.all):
            data_minima[key] = self.data[key].min()
            data_maxima[key] = self.data[key].max()

        hkl_min = min((data_minima['h'], data_minima['k'], data_minima['l']))
        hkl_max = max((data_maxima['h'], data_maxima['k'], data_maxima['l']))
        scale = 2./max((abs(hkl_max), abs(hkl_min)))

        # print the title line
        a = self.a_r * 1000
        b = self.b_r * 1000
        c = self.c_r * 1000
        al = np.degrees(self.al_r)
        be = np.degrees(self.be_r)
        ga = np.degrees(self.ga_r)
        cell_list = (self.la, a, b, c, al, be, ga)
        hklres_file = open(path, 'w')
        hklres_file.write('TITL hkl visualisation\n')
        hklres_file.write('REM special hkl visualisation file, to be used in'
                          'mercury with hkl.msd style applied\n')
        hklres_file.write('REM reciprocal unit cell has been inflated '
                          'thousandfold for technical reasons\n')
        hklres_file.write('CELL {0:7f} {1:7f} {2:7f} {3:7f} {4:7f} {5:7f} '
                          '{6:7f}\n'.format(*cell_list))
        hklres_file.write('LATT -1\n\n')

        def get_label(integer):
            labels = ('H', 'Li', 'B', 'N', 'F', 'Na', 'Al', 'P', 'Cl')
            return labels[(int(integer)-1) % len(labels)]

        for index, reflection in self.data.iterrows():
            line_pars = dict()
            line_pars['label'] = '{atom}({h_index},{k_index},{l_index})'.format(
                atom=get_label(reflection[colored]),
                h_index=int(reflection['h']),
                k_index=int(reflection['k']),
                l_index=int(reflection['l']))
            line_pars['x'] = float(reflection['h'] * scale)
            line_pars['y'] = float(reflection['k'] * scale)
            line_pars['z'] = float(reflection['l'] * scale)
            line_pars['size'] = np.log(abs(reflection[master_key])+1.0) ** 2 \
                                * 5 * scale ** 2
            hklres_file.write('{label:16}   1 {x: .4f} {y: .4f} {z: .4f} '
                              '11.0000 {size: .12f}\n'.format(**line_pars))
        hklres_file.close()

    def trim(self, limit):
        """
        Remove from dataframes those reflections, which lie further than *limit*
        from the reciprocal space origin point.
        :param limit: Radius of the trimming sphere in reciprocal Angstrom
        :type limit: float
        """
        self.data = self.data.loc[self.data['r'] <= limit]

    def write(self, hkl_path, hkl_format, columns_separator=True):
        """
        Write the contents of dataframe to a .hkl file using specified
        *path* and *format*.
        For a list of all available .hkl formats,
        please refer to the documentation of :func:`interpret_hkl_format`.

        :param hkl_path: Absolute or relative path to the .hkl file.
        :type hkl_path: str
        :param hkl_format: Desired format of .hkl file.
        :type hkl_format: union[int, str, OrderedDict]
        :param columns_separator: should columns be separated using whitespace?
        :type columns_separator: bool
        """

        # PREPARE OBJECTS RESPONSIBLE FOR WRITING OUTPUT
        format_string, column_labels, file_prefix, file_suffix, zero_line = \
            self.interpret_hkl_format(hkl_format)
        if format_string == 'XD':
            format_string, column_labels, file_prefix,\
                file_suffix, zero_line = \
                self.interpret_hkl_format(OrderedDict([('h', 5), ('k', 5),
                                                       ('l', 5), ('b', 5),
                                                       ('I', 10), ('si', 10)]))
        column_sizes, column_formats = [], []
        for column in format_string.split():
            number, letter = int(column[:-1]), column[-1]
            if letter == 'x':
                column_sizes.append(0)
            elif letter == 's':
                column_sizes.append(number)
        for size in column_sizes:
            cs = int(columns_separator)
            column_formats.append('{{:>{}.{}}}'.format(size, size-cs))

        # PREPARE NON-EXISTED, BUT DEMANDED ROWS
        for key in column_labels:
            try:
                self.data[key]
            except KeyError:
                dummy_column = []
                default_value = self.keys.get_property(key, 'default')
                dtype = self.keys.get_property(key, 'dtype')
                for index in range(self.data.shape[0]):
                    dummy_column.append(default_value)
                self.data[key] = pd.Series.from_array(dummy_column, dtype=dtype)

        # WRITE THE FREE FORMAT
        hkl_file = open(hkl_path, 'w')
        if str(hkl_format).lower() == 'free':
            self.data.to_csv(path_or_buf=hkl_path, sep=' ',
                             columns=column_labels, header=False)
            return

        # WRITE PREFIX LINE
        if file_prefix is not False:
            hkl_file.write(file_prefix + '\n')

        # WRITE SELF.DATA CONTENTS
        for index, row in self.data.iterrows():
            # FOR EACH DEMANDED KEY PRINT IT ACCORDING TO FORMATS
            for key, form in zip(column_labels, column_formats):
                hkl_file.write(form.format(
                    str(self.keys.get_property(key, 'type')(row[key]))))
            hkl_file.write('\n')

        # WRITE 0 0 0 LINE
        if zero_line is True:
            for key, form in zip(column_labels, column_formats):
                hkl_file.write(form.format(
                    str(self.keys.get_property(key, 'default'))))

        # WRITE SUFFIX LINE
        if file_suffix is not False:
            hkl_file.write(file_suffix + '\n')
        hkl_file.close()


class HklIo:
    """
    A helper class supporting HklFrame.
    Menages reading and writing hkl files
    into and out of HklFrame's dataframe
    """

    # TODO write hooks for read and write in HklFrame
    # TODO add 'sf' (denoted using letter 'stigma') as alternative sigma for 'F'
    # TODO Fix the documentation using this new object
    # TODO Add hklres generator to this HklIo
    # TODO think about space groups...

    def __init__(self):
        self.__format = 'shelx_4'
        self.keys = HklKeys()
        self.use_separator = True
        self._load_format_dictionaries()

    def _build_format_string(self):
        """
        Prepare an input string for string 'format' method to write hkl data.
        :return: String for str.format() to format hkl data while writing.
        :rtype: str
        """
        built_format_string = str()
        for l, w in zip(self._format_dict['labels'], self._format_dict['widths']):
            built_format_string += ' ' * self.use_separator
            w = abs(w) - int(self.use_separator)
            built_format_string += '{{{0}:>{1}}}'.format(l, w)
        built_format_string += '\n'
        self.__format_string = built_format_string

    @property
    def format(self):
        """
        Return a name of currently used format.
        :return: String with internal representation of hkl format.
        :rtype: str
        """
        return self.__format

    @format.setter
    def format(self, new_format):
        if new_format in self.formats_defined.keys():
            self.__format = new_format
        elif new_format in self.formats_aliases.keys():
            self.__format = self.formats_aliases[new_format]
        else:
            raise KeyError('Unknown hkl format "{}" given'.format(new_format))
        self._build_format_string()

    @property
    def _format_dict(self):
        return self.formats_defined[self.__format]

    @property
    def _format_string(self):
        return self.__format_string

    def _import_custom_format(self, custom_format_string):
        """
        Import format string such as 'h 4 k 4 l 4 I 8 si 8' as format 'custom'.

        :param custom_format_string: string containing alternating data labels
            and column widths (all negative if free format) separated using ' '.
        :type custom_format_string: str
        :return: None
        """
        custom_format_list = custom_format_string.strip(' ').split(' ')
        custom_labels = custom_format_list[0::2]
        custom_widths = [int(width) for width in custom_format_list[1::2]]
        assert len(custom_labels) == len(custom_widths), "Bad custom format"
        self.formats_defined['custom'] = {'labels': custom_labels,
                                          'widths': custom_widths,
                                          'prefix': '',
                                          'suffix': ''}

    @property
    def is_current_format_free(self):
        """
        Return true if currently defined format is free
        :return: True if all format widths are negative; False otherwise.
        :rtype: bool
        """
        return all(width < 0 for width in self._format_dict['widths'])

    def _load_format_dictionaries(self):
        """
        Load dictionaries of defined formats and their aliases from json files.
        :return: None
        """
        current_file_path = Path(__file__).parent.absolute()
        path_of_defined = current_file_path.joinpath('hkl_formats_defined.json')
        path_of_aliases = current_file_path.joinpath('hkl_formats_aliases.json')
        with open(path_of_defined) as file:
            self.formats_defined = json.load(file)
        with open(path_of_aliases) as file:
            self.formats_aliases = json.load(file)

    def _parse_fixed_line(self, line):
        """
        Parse data from a line, where data from each *label* has fixed *width*.
        :param line: string to be parsed based on format dictionary.
        :type line: str
        :return: list of strings extracted from parsed line
        :rtype: list
        """
        slice_end = list(np.cumsum(self._format_dict['widths']))
        slice_beg = [0] + slice_end[:-1]
        parsed = [line[beg:end] for beg, end in zip(slice_beg, slice_end)]
        parsed = np.array(parsed)
        try:
            parsed.astype('float64')
            assert len(parsed) == len(self._format_dict['widths'])
        except (ValueError, AssertionError):
            return None
        return parsed

    def _parse_free_line(self, line):
        """
        Parse data from line, where data from *labels* is separated with space.
        :param line: string to be parsed based on format dictionary.
        :type line: str
        :return: list of strings extracted from parsed line
        :rtype: list
        """
        parsed = np.array(line.strip().split())
        try:
            parsed.astype('float64')
            assert len(parsed) == len(self._format_dict['widths'])
        except (ValueError, AssertionError):
            return None
        return parsed

    def read(self, hkl_path, hkl_format):
        self.format = hkl_format
        self.keys.set(self._format_dict['labels'])
        parse_line = self._parse_free_line if self.is_current_format_free() \
            else self._parse_fixed_line

        def read_file_to_list_of_data():
            list_of_reflections = list()
            with open(hkl_path, 'r') as hkl_file:
                for line in hkl_file.read().splitlines():
                    parsed_line = parse_line(line)
                    if parsed_line is None:
                        continue
                    list_of_reflections.append(parsed_line)
            return np.array(list_of_reflections).transpose()
        array_of_reflections = read_file_to_list_of_data()

        def build_dict_of_reflections(_hkl_array):
            dict_of_data = dict()
            for index, key in enumerate(self._format_dict['labels']):
                key_dtype = self.keys.get_property(key, 'dtype')
                dict_of_data[key] = _hkl_array[index].astype(key_dtype)
                return dict_of_data
        return build_dict_of_reflections(array_of_reflections)
        # TODO returns dict; should call some generator to return HklFrame.data?

    def write(self, hkl_data, hkl_path, hkl_format):
        """hkl_data must be pandas df!"""
        self.format = hkl_format
        with open(hkl_path, 'w') as hkl_file:
            hkl_file.write(self._format_dict['file_prefix'])
            for index, row in hkl_data.iterrows():
                hkl_file.write(self._format_string().format(**row))
            hkl_file.write(self._format_dict['file_suffix'])





if __name__ == '__main__':
    pass
