from kesshou.dataframes import BaseFrame
from kesshou.utility import cubespace, elements_list, is2n, is3n, is4n, is6n
from kesshou.utility import rescale_list_to_range, rescale_list_to_other
from kesshou.symmetry import PG
from pathlib import Path
import copy
import json
import random
import numpy as np
import numpy.linalg as lin
import pandas as pd
import matplotlib.cm


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
        'description': 'Uncertainty of intensity determination',
        'imperative': False,
        'dtype': 'float64',
        'reduce_behaviour': 'average',
        'type': float
    }
    __sf = {
        'default': 0.0,
        'description': 'Uncertainty of structure factor determination',
        'imperative': False,
        'dtype': 'float64',
        'reduce_behaviour': 'average',
        'type': float
    }
    __u = {
        'default': 0.0,
        'description': 'Structure factor to its uncertainty ratio',
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
    defined_keys = {'h', 'k', 'l', 'F', 'I', 'si', 'sf', 'b', 'm', 'la', 'ph',
                    'u', 'r', 't', 'u1', 'u2', 'u3', 'v1', 'v2', 'v3',
                    'x', 'y', 'z', 'equiv'}

    def __init__(self, keys=set()):
        # DEFINE ALL KNOWN KEYS
        self.imperatives = set()
        for key in HklKeys.defined_keys:
            if self.get_property(key, 'imperative'):
                self.imperatives.add(key)
        self.all = set()
        self.add(self.imperatives)
        self.add(keys)

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
    A master object which menages single-crystal diffraction files.
    It utilises other "Hkl" classes to import, store, manipulate and output
    information about single-crystal diffraction patterns.

    HklFrame acts as an container which stores
    the diffraction data (Pandas dataframe, :attr:`table`)
    and elementary crystal cell data (:class:`kesshou.dataframes.Base`).
    Demanding methods belonging to this class are vectorized,
    providing relatively satisfactory performance and high memory capacity.
    HklFrame methods are designed to work in-place, so the work strategy
    is to create a new instance of HklFrame for each reflection dataset,
    manipulate it using methods, eg. :func:`merge` or :func:`trim`, and
    :func:`duplicate` to other object or output using :func:`write` if needed.

    The HklFrame always initiates empty and does not accept any arguments.
    Some of the magic methods, such as :func:`__len__` and :func:`__add__`
    are defined and describe/operate on the :attr:`frame`.
    """

    HKL_LIMIT = 99
    """Highest absolute value of h, k or l index, which can be
    interpreted correctly by current version of the software."""

    def __init__(self):
        """HklFrame constructor"""
        super().__init__()

        self.keys = HklKeys()
        """Object managing keys (column names) of :attr:`table`."""

        self.__la = 0.71069
        """Wavelength of radiation used in experiment."""

        self.table = pd.DataFrame()
        """
        Pandas dataframe containing diffraction data information.
        Each row represents one reflection observation,
        while each column has one piece of information about the reflections.
        For a list of available keys, see :class:`HklKeys`,
        whose instance is used to menage the keys of this table.
        """

    def __add__(self, other):
        """
        Add magic method. Stacks contents of two :attr:`table` while preserving
        meta-information from the first :class:`HklFrame` object.

        :param other: HklFrame to be added to data
        :type other: HklFrame
        :return: HklFrame with concatenated HklFrame.data pandas dataframes
        :rtype: HklFrame
        """
        _copied = self.duplicate()
        _copied.table = pd.concat([self.table, other.table], ignore_index=True)
        return _copied

    def __len__(self):
        """
        Len magic method, number of individual data points in :attr:`table`.
        :return: Number of rows (individual reflections) in `self.data`
        :rtype: int
        """
        return self.table.shape[0]

    def __str__(self):
        """
        Str magic method, provides human-readable representation of data
        :return: Human-readable representation of `self.data`
        :rtype: str
        """
        return self.table.__str__()

    @property
    def la(self):
        """
        Wavelength of radiation used in the diffraction experiment.
        Can be set using popular definitions such as "MoKa" or "CuKb",
        where a and b stand for *alpha* and *beta*.
        Implemented cathode materials include:
        "Ag", "Co", "Cr", "Cu", "Fe", "Mn", "Mo", "Ni", "Pd", "Rh", "Ti", "Zn"
        and have been imported from International Tables
        of Crystallography, Volume C, Table 4.2.4.1, 3rd Edition.
        :return: wavelength of radiation used in experiment
        :rtype: float
        """
        return self.__la

    @la.setter
    def la(self, wavelength):
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

    @property
    def r_lim(self):
        """
        Calculate limiting sphere radius based on :attr:`la`.

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
        df = self.table
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
        the dac-accessible space traced by the tori.

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

        opening_angle_in_radians = np.deg2rad(opening_angle)
        if vector is None:
            l_v = np.array((1.0, 0.0, 0.0))             # vec. parallel to beam
            h = np.dot(lin.inv(self.orientation), l_v)  # calculate hkl vector
            n = h[0] * self.a_w + h[1] * self.b_w + h[2] * self.c_w
        else:                                           # calculate xyz* vector
            n = np.array(vector)
        n = 1 / lin.norm(n) * n                        # normalise the n vector

        self.trim(limit=self.r_lim * np.sin(opening_angle_in_radians))
        # remove reflections perpendicular to disc's normal vector
        xyz = self.table.loc[:, ('x', 'y', 'z')].to_numpy()
        self.table = self.table.loc[~np.isclose(xyz @ n, self.table['r'])]
        # calculate reflection's position in 2D disc reference system "m"
        # m1 / m2 is a coordinate parallel / perpendicular to vector n
        xyz = self.table.loc[:, ('x', 'y', 'z')].to_numpy()
        m1 = np.outer(xyz @ n, n)
        m2 = xyz - m1
        m1 = m1 @ n
        m2 = lin.norm(m2, axis=1)
        # find the middle of two tori, which trace the DAC-limiting shape
        t1 = 1 / 2 * self.r_lim * np.cos(opening_angle_in_radians)
        t2 = 1 / 2 * self.r_lim * np.sin(opening_angle_in_radians)
        # check if points lie in one of two tori making the DAC shape
        in_torus1 = (m1 - t1) ** 2 + (m2 - t2) ** 2 <= (self.r_lim / 2) ** 2
        in_torus2 = (m1 + t1) ** 2 + (m2 - t2) ** 2 <= (self.r_lim / 2) ** 2
        # leave only points which lie in both tori
        self.table = self.table[in_torus1 * in_torus2]
        self.table.reset_index(drop=True, inplace=True)

    def _domain(self, address='hkl'):
        """
        This method limits the reflection data to the ones "living" in address.

        :param address: Address which the reflections must have to be preserved
        :type address: str
        :return: Dataframe containing only reflections with given address
        :rtype: pd.DataFrame
        """

        address = address.lower().replace(' ', '').replace('_', '')
        df = self.table

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

    def duplicate(self):
        """
        Make and return an exact deep copy of this HklFrame.

        :return: A copy of this HklFrame.
        :rtype: HklFrame
        """
        return copy.deepcopy(self)

    def extinct(self, rule='hkl:', point_group=PG['1']):
        """
        Removes from dataframe all reflections which are in a specified
        *domain*, but do not meet the *condition*.
        The *rules* have a format "domain:condition",
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
        should be accessible.

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
            self.table = self.table[~((self.table['h'] == h) &
                                      (self.table['k'] == k) &
                                      (self.table['l'] == l))]
        # reset the indices for other methods to use
        self.table.reset_index(drop=True, inplace=True)

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
        hkl_lim = self.HKL_LIMIT
        self.keys.add(('equiv',))
        self.table['equiv'] = [(-hkl_lim, -hkl_lim, -hkl_lim)] * len(self.table)
        _hkl_matrix = self.table.loc[:, ('h', 'k', 'l')].to_numpy()
        for op in point_group.operations:
            new_hkl = pd.Series(map(tuple, _hkl_matrix @ op[0:3, 0:3]))
            _to_update = self.table['equiv'] < new_hkl
            self.table.loc[_to_update, 'equiv'] = new_hkl.loc[_to_update]

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
        self.table = new_data
        self.extinct('000')
        if not('x' in self.table.columns):
            self.place()
        self._recalculate_structure_factors_and_intensities()

    def fill(self, radius=2.0):
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
        while len(hkl) > previous_length and max_index <= self.HKL_LIMIT:
            previous_length = len(hkl)
            max_index = max_index * 2
            hkl = _make_hkl_ball(max_index)

        # create new dataframe using obtained ball of data
        _h, _k, _l = np.vsplit(hkl.T, 3)
        ones = np.ones_like(np.array(_h)[0])
        self.from_dict({'h': np.array(_h)[0], 'k': np.array(_k)[0],
                        'l': np.array(_l)[0], 'I': ones, 'si': ones, 'm': ones})

    def stats(self, bins=10, point_group=PG['1'], extinctions=('000',)):
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
            _hkl_full.fill(radius=max(self.table['r']))
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
            cube_bins = cubespace(0.0, max(self.table['r']), num=_bins + 1)
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
            return results
        make_table_with_stats(grouped_base, grouped_full, grouped_merged)

    def merge(self, point_group=PG['1']):
        """
        Average down each set of redundant reflections present in the table,
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
        grouped = self.table.groupby('equiv')
        grouped_first = grouped.first().reset_index()
        grouped_mean = grouped.mean().reset_index()
        grouped_sum = grouped.sum().reset_index()
        # for each key apply a necessary reduce operation and add it to data
        data = dict()
        for key in self.table.keys():
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
        hkl = self.table.loc[:, ('h', 'k', 'l')].to_numpy()
        abc = np.matrix((self.a_w, self.b_w, self.c_w))
        xyz = hkl @ abc
        self.table['x'] = xyz[:, 0]
        self.table['y'] = xyz[:, 1]
        self.table['z'] = xyz[:, 2]
        self.table['r'] = lin.norm(xyz, axis=1)
        self.keys.add(('x', 'y', 'z', 'r'))

    def read(self, hkl_path, hkl_format='shelx_4'):
        """
        Read the contents of .hkl file as specified by path and format,
        and store them in the pandas dataframe in `self.data`.
        For a list of all available .hkl formats,
        please refer to :attr:`kesshou.dataframes.HklIo.format`.

        :param hkl_path: Absolute or relative path to the .hkl file.
        :type hkl_path: str
        :param hkl_format: Format of provided .hkl file.
        :type hkl_format: union[int, str, OrderedDict]
        """
        reader = HklReader(hkl_file_path=hkl_path, hkl_file_format=hkl_format)
        dict_of_data = reader.read()
        self.keys.set(dict_of_data.keys())
        forgotten_keys = tuple(self.keys.imperatives - set(dict_of_data.keys()))
        for forgotten in forgotten_keys:
            default = self.keys.get_property(forgotten, 'default')
            length_of_data = max([len(v) for v in dict_of_data.values()])
            dict_of_data[forgotten] = [default] * length_of_data
        self.from_dict(dict_of_data)

    def _recalculate_structure_factors_and_intensities(self):
        """
        Calculate 'I' and 'si' or 'F' and 'sf', depending on which are missing.
        """
        if 'I' in self.keys.all and not('si' in self.keys.all):
            raise KeyError('Intensities "I" are defined, but "si" not.')
        if 'F' in self.keys.all and not('sf' in self.keys.all):
            raise KeyError('Structure factors "F" are defined, but "sf" not.')
        if all(('I' in self.keys.all,
               'si' in self.keys.all,
                not('F' in self.keys.all),
                not('sf' in self.keys.all))):
            self._recalculate_structure_factors_from_intensities()
        if all(('F' in self.keys.all,
               'sf' in self.keys.all,
                not('I' in self.keys.all),
                not('si' in self.keys.all))):
            self._recalculate_intensities_from_structure_factors()

    def _recalculate_structure_factors_from_intensities(self):
        """
        Recalculate the structure factor F and its uncertainty sf.

        Structure factor is calculated as follows:
        *F = signum(I) \* sqrt(abs(I))*.

        Structure factor's uncertainty is calculated as follows:
        *sf = si / (2 \* sqrt(abs(I)))*.
        """
        new_data = copy.deepcopy(self.table)
        signum_of_i = new_data['I'].copy()
        signum_of_i[signum_of_i > 0] = 1
        signum_of_i[signum_of_i < 0] = -1
        absolute_sqrt_of_i = abs(new_data["I"]) ** 0.5
        new_data['F'] = signum_of_i * absolute_sqrt_of_i
        new_data['sf'] = new_data["si"] / (2 * absolute_sqrt_of_i)
        self.table = new_data
        self.keys.add({'F', 'sf'})

    def _recalculate_intensities_from_structure_factors(self):
        """
        Recalculate the intensity I and its uncertainty si.

        Intensity is calculated as follows:
        *I = signum(F) \* F \*\* 2*.

        Intensity's uncertainty is calculated as follows:
        *si = 2 \* sf \* abs(F)*.
        """
        new_data = copy.deepcopy(self.table)
        signum_of_f = new_data['F'].copy()
        signum_of_f[signum_of_f > 0] = 1
        signum_of_f[signum_of_f < 0] = -1
        new_data['I'] = signum_of_f * (abs(new_data['F']) ** 2)
        new_data['si'] = 2 * new_data["sf"] * abs(new_data["F"])
        self.table = new_data
        self.keys.add({'I', 'si'})

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
                df = copy.deepcopy(self.table)
                mat = op[0:3, 0:3]
                hkl = self.table.loc[:, ['h', 'k', 'l']].to_numpy() @ mat
                df['h'] = hkl[:, 0]
                df['k'] = hkl[:, 1]
                df['l'] = hkl[:, 2]
                dataframes.append(df)
            return pd.concat(dataframes, axis=0, ignore_index=True)
        self.table = _build_transformed_dataframe()
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
        self.table.drop(indices_to_delete, inplace=True)
        self.table.reset_index(drop=True, inplace=True)

    def to_res(self, path='hkl.res', colored='m'):
        """
        Export the reflection information from table to .res file,
        so that a software used to visualize .res files can be used
        to visualise a diffraction data in three dimensions.

        :param colored: Which key of dataframe should be visualised using color.
        :type colored: str
        :param path: Absolute or relative path where the file should be saved
        :type path: str
        """
        artist = HklArtist(self)
        artist.write_res(path=path, colored=colored)

    def trim(self, limit):
        """
        Remove from table those reflections, which lie further than *limit*
        from the reciprocal space origin point.
        :param limit: Radius of the trimming sphere in reciprocal Angstrom
        :type limit: float
        """
        self.table = self.table.loc[self.table['r'] <= limit]

    def write(self, hkl_path, hkl_format='shelx_4'):
        """
        Write the contents of dataframe to a .hkl file using specified
        *path* and *format*.
        For a list of all available .hkl formats,
        please refer to :attr:`kesshou.dataframes.HklIo.format`.

        :param hkl_path: Absolute or relative path to the .hkl file.
        :type hkl_path: str
        :param hkl_format: Desired format of .hkl file.
        :type hkl_format: union[int, str, OrderedDict]
        """
        writer = HklWriter(hkl_file_path=hkl_path, hkl_file_format=hkl_format)
        writer.write(hkl_data=self.table)


class HklIo:
    """
    A helper class supporting HklFrame.
    Menages reading and writing hkl files
    into and out of HklFrame's dataframe
    """

    def __init__(self, hkl_file_path, hkl_file_format):
        self.keys = HklKeys()
        self.use_separator = True
        self.file_path = hkl_file_path
        self._load_format_dictionaries()
        self.__format = 'shelx_4'
        self.format = hkl_file_format

    def _build_line_formatter(self):
        """
        Prepare :meth:`line_formatter` to format hkl data while writing.

        :return: String for str.format() to format hkl data while writing.
        :rtype: str
        """
        formatter = str()
        for label, width in zip(self._format_dict['labels'],
                                self._format_dict['widths']):
            width = abs(width)
            text_width = abs(width) - int(self.use_separator)
            formatter += '{{{0}:>{1}.{2}}}'.format(label, width, text_width)
        formatter += '\n'
        self.__line_formatter = formatter

    @property
    def format(self):
        """
        Return a name of currently used hkl file format. Available file formats
        and their aliases are defined internally in .json files and
        have been presented in a table below:

        +----------+----------+------------------------+------+------+--------+
        | Name     | Aliases  | Contents               |Prefix|Suffix| Free   |
        |          |          | (format string)        |      |      | format |
        +==========+==========+========================+======+======+========+
        | free_2   |          | h -4 k -4 l -4 I -8    | NO   | YES  | YES    |
        |          |          | si -8 b -4 la -8       |      | (a)  |        |
        +----------+----------+------------------------+------+------+--------+
        | free_3   |          | h -4 k -4 l -4         | NO   | YES  | YES    |
        |          |          | F -8 sf -8 b -4        |      | (a)  |        |
        +----------+----------+------------------------+------+------+--------+
        | free_4   |          | h -4 k -4 l -4         | NO   | YES  | YES    |
        |          |          | I -8 si -8 b -4        |      | (a)  |        |
        +----------+----------+------------------------+------+------+--------+
        | free_40  | free     | h -4 k -4 l -4         | NO   | YES  | YES    |
        |          |          | I -8 si -8             |      | (a)  |        |
        +----------+----------+------------------------+------+------+--------+
        | free_5   |          | h -4 k -4 l -4         | NO   | YES  | YES    |
        |          |          | I -8 si -8 c -4        |      | (a)  |        |
        +----------+----------+------------------------+------+------+--------+
        | free_6   |          | h -4 k -4 l -4         | NO   | YES  | YES    |
        |          |          | I -8 si -8 m -4        |      | (a)  |        |
        +----------+----------+------------------------+------+------+--------+
        | shelx_2  | 2        | h 4 k 4 l 4 I 8 si 8   | NO   | YES  | NO     |
        |          |          | b 4 la 8               |      | (a)  |        |
        +----------+----------+------------------------+------+------+--------+
        | shelx_3  | 3        | h 4 k 4 l 4            | NO   | YES  | NO     |
        |          |          | F 8 sf 8 b 4           |      | (a)  |        |
        +----------+----------+------------------------+------+------+--------+
        | shelx_4  | 4        | h 4 k 4 l 4            | NO   | YES  | NO     |
        |          |          | I 8 si 8               |      | (a)  |        |
        +----------+----------+------------------------+------+------+--------+
        | shelx_40 | 40       | h 4 k 4 l 4            | NO   | YES  | NO     |
        |          |          | I 8 si 8               |      | (a)  |        |
        +----------+----------+------------------------+------+------+--------+
        | shelx_5  | 5        | h 4 k 4 l 4            | NO   | YES  | NO     |
        |          |          | I 8 si 8 c 4           |      | (a)  |        |
        +----------+----------+------------------------+------+------+--------+
        | shelx_6  | 6        | h 4 k 4 l 4            | NO   | YES  | NO     |
        |          |          | I 8 si 8 m 4           |      | (a)  |        |
        +----------+----------+------------------------+------+------+--------+
        | tonto_F  |          | h -4 k -4 l -4         | YES  | YES  | YES    |
        |          |          | F -8 sf -8             | (b)  | (b)  |        |
        +----------+----------+------------------------+------+------+--------+
        | tonto_I  | tonto    | h -4 k -4 l -4         | YES  | YES  | YES    |
        |          |          | I -8 si -8             | (b)  | (b)  |        |
        +----------+----------+------------------------+------+------+--------+
        | xd_F6    |          | h -4 k -4 l -4 b -3    | YES  | NO   | YES    |
        |          |          | F -13 sf -13           | (c)  |      |        |
        +----------+----------+------------------------+------+------+--------+
        | xd_F7    |          | h -4 k -4 l -4 b -3    | YES  | NO   | YES    |
        |          |          | F -13 sf -13 t -10     | (c)  |      |        |
        +----------+----------+------------------------+------+------+--------+
        | xd_F-7   |          | h -4 k -4 l -4 b -3    | YES  | NO   | YES    |
        |          |          | F -13 sf -13 ph -10    | (c)  |      |        |
        +----------+----------+------------------------+------+------+--------+
        | xd_F13   |          | h -4 k -4 l -4 b -3    | YES  | NO   | YES    |
        |          |          | F -13 sf -13 t -10     | (c)  |      |        |
        |          |          | u1 -10 u2 -10 u3 -10   |      |      |        |
        |          |          | v1 -10 v2 -10 v3 -10   |      |      |        |
        +----------+----------+------------------------+------+------+--------+
        | xd_I6    | xd       | h -4 k -4 l -4 b -3    | YES  | NO   | YES    |
        |          |          | I -13 si -13           | (c)  |      |        |
        +----------+----------+------------------------+------+------+--------+
        | xd_I7    |          | h -4 k -4 l -4 b -3    | YES  | NO   | YES    |
        |          |          | I -13 si -13 t -10     | (c)  |      |        |
        +----------+----------+------------------------+------+------+--------+
        | xd_I-7   |          | h -4 k -4 l -4 b -3    | YES  | NO   | YES    |
        |          |          | I -13 si -13 ph -10    | (c)  |      |        |
        +----------+----------+------------------------+------+------+--------+
        | xd_I13   |          | h -4 k -4 l -4 b -3    | YES  | NO   | YES    |
        |          |          | I -13 si -13 t -10     | (c)  |      |        |
        |          |          | u1 -10 u2 -10 u3 -10   |      |      |        |
        |          |          | v1 -10 v2 -10 v3 -10   |      |      |        |
        +----------+----------+------------------------+------+------+--------+
        | *custom* |          | custom string as above | NO   | NO   | if all |
        |          |          | with keys and widths   |      |      | widths |
        |          |          |                        |      |      | in     |
        |          |          |                        |      |      | format |
        |          |          |                        |      |      | are <0 |
        +----------+----------+------------------------+------+------+--------+

        Three different types of prefix / suffix are supported at the moment:

        - Suffix (a) is a zero-line: a shelx ending line with h = k = l = 0,

        - Prefix and suffix (b) are tonto-characteristic beginning/end of file,

        - Prefix (c) is an xd-characteristic line with info about file content.

        A custom hkl file format can be defined by providing
        a *format string* instead of 'Name'.
        The string should look like the ones in column "contents".
        For the meaning of keys ('I', 'b', 'c' etc.),
        please refer to :class:`HklKeys`.

        :return: Returns a name of currently used format.
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
        self._build_line_formatter()

    @property
    def _format_dict(self):
        """A dictionary with details concerning current format."""
        return self.formats_defined[self.__format]

    @property
    def _line_formatter(self):
        """A string used by string.format() method to print hkl data."""
        return self.__line_formatter

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
        Return true if currently defined format is free,
        i.e. the columns are separated by whitespace.

        :return: True if all format widths are negative; False otherwise.
        :rtype: bool
        """
        return all(width < 0 for width in self._format_dict['widths'])

    def _load_format_dictionaries(self):
        """
        Load dictionaries of defined formats and their aliases from json files.
        """
        current_file_path = Path(__file__).parent.absolute()
        path_of_defined = current_file_path.joinpath('hkl_formats_defined.json')
        path_of_aliases = current_file_path.joinpath('hkl_formats_aliases.json')
        with open(path_of_defined) as file:
            self.formats_defined = json.load(file)
        with open(path_of_aliases) as file:
            self.formats_aliases = json.load(file)


class HklReader(HklIo):
    """
    A helper class for HklFrame,
    Menages reading hkl files and importing data and keys from them
    """

    def __init__(self, hkl_file_path, hkl_file_format):
        super().__init__(hkl_file_path, hkl_file_format)

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

    def read(self):
        """
        Read the contents of file currently pointed by :attr:`hkl_file_path`
        and format :attr:`hkl_file_format` and return them to a dictionary.

        :return: A dictionary containing information read from .hkl file.
        :rtype: dict
        """
        self.keys.set(self._format_dict['labels'])
        parse_line = self._parse_free_line if self.is_current_format_free \
            else self._parse_fixed_line

        def read_file_to_list_of_data():
            list_of_reflections = list()
            with open(self.file_path, 'r') as hkl_file:
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


class HklWriter(HklIo):
    """
    A helper class for HklFrame,
    Menages writing hkl files and exporting data to them
    """
    def __init__(self, hkl_file_path, hkl_file_format):
        super().__init__(hkl_file_path, hkl_file_format)

    def write(self, hkl_data):
        """
        Write data from pandas dataframe `hkl_data` to the file specified 
        at :attr:`hkl_file_path` of format :attr:`hkl_file_format`.
        
        :param hkl_data: Dataframe containing reflection information.
        :type hkl_data: pandas.dataframe
        """"""
        
        hkl_data must be pandas df!"""
        needed_data = hkl_data.loc[:, self._format_dict['labels']].astype(str)
        with open(self.file_path, 'w') as hkl_file:
            hkl_file.write(self._format_dict['prefix'])
            for row_tuple in needed_data.itertuples(index=False):
                row_dict = dict(zip(self._format_dict['labels'], row_tuple))
                hkl_file.write(self._line_formatter.format(**row_dict))
            hkl_file.write(self._format_dict['suffix'])


class HklArtist:
    """
    A class responsible for representing the HklData using either images
    or other files, which can be further visualised.
    """
    COLORMAP = matplotlib.cm.get_cmap('gist_rainbow')
    "Colourmap used to signify the 'colored' property in hkl.res"

    ELEMENT_NAMES = elements_list
    "List of chemical elements used to give reflections in hkl.res a colour"

    MIN_ATOM_DISTANCE = 10.0
    "Minimum available value of reciprocal unit cell length parameter"

    MIN_ATOM_SIZE = 0.00001
    "Maximum available value of U_iso while defining atom sizes"

    MAX_ATOM_SIZE = 4.99999
    "Maximum available value of U_iso while defining atom sizes"

    def __init__(self, hkl_dataframe):
        self.df = hkl_dataframe

    @property
    def maximum_index(self):
        """
        Return largest absolute value of the following: h, k, l, -h, -k, -l.

        :return: Largest absolute hkl index.
        :rtype: int
        """
        maxima = {key: self.df.data[key].max() for key in self.df.data.keys()}
        minima = {key: self.df.data[key].min() for key in self.df.data.keys()}
        return max(maxima['h'], maxima['k'], maxima['l'],
                   -minima['h'], -minima['k'], -minima['l'])

    @property
    def color_dict(self):
        """
        Adictionary containing colours used for current map.
        :return: Dictionary of element_name:rgb_colour pairs.
        :rtype: dict
        """
        linspace = np.linspace([0], [1], len(self.ELEMENT_NAMES))
        elements = self.ELEMENT_NAMES
        return {k: self.COLORMAP(float(v)) for k, v in zip(elements, linspace)}

    @property
    def res_distance_scale(self):
        """
        Define and return a distance scale so that all reciprocal vectors
        a, b and c have length of at least :attr:`RES_MIN_ATOM_DISTANCE`.

        :return: A scale factor to multiply call distances by
        :rtype: float
        """
        min_cell_length = min(self.df.a_r, self.df.b_r, self.df.c_r)
        return self.maximum_index * self.MIN_ATOM_DISTANCE / min_cell_length

    @property
    def res_reciprocal_cell(self):
        """
        Return dictionary of seven unit cell parameters from file's header,
        with distances rescaled to prevent software from showing bonds.

        :return: Dictionary containing la, a, b, c, al, be, ga values.
        :rtype: dict
        """
        return {'la': self.df.la,
                'a': self.df.a_r * self.res_distance_scale,
                'b': self.df.b_r * self.res_distance_scale,
                'c': self.df.c_r * self.res_distance_scale,
                'al': np.rad2deg(self.df.al_r),
                'be': np.rad2deg(self.df.be_r),
                'ga': np.rad2deg(self.df.ga_r)}

    @property
    def res_header(self):
        """
        String containing information which should be find of top of res file -
        title line, file description and unit cell information.

        :return: String containing res file header.
        :rtype: str
        """
        return "TITL Reflection visualisation\n" \
               "REM Special file to be used in mercury with hkl.msd style.\n" \
               "REM Reciprocal unit cell scaled by {sc} to prevent bonding\n" \
               "CELL {la:7f} {a:7f} {b:7f} {c:7f} {al:7f} {be:7f} {ga:7f}\n" \
               "LATT -1\n\n".format(sc=self.res_distance_scale,
                                    **self.res_reciprocal_cell)

    def res_legend(self, colored='m'):
        """
        Return a legend - a least of meanings of colors used in hklres.

        :return: String of res commands containing details of color meaning.
        :rtype: str
        """
        value_range = np.arange(min(self.df.data[colored]),
                                max(self.df.data[colored]))
        elements_range = rescale_list_to_other(list(value_range), elements_list)
        line_string = 'REM {0} = {{v}}: {{e}}, rgba{{c}}.'.format(colored)
        line_list = [line_string.format(v=v, e=e, c=self.color_dict[e])
                     for v, e in zip(value_range, elements_range)]
        return '\n'.join(line_list) + '\n\n'

    @staticmethod
    def res_line(color, h_ind, k_ind, l_ind, x_pos, y_pos, z_pos, u_iso):
        """
        Transform given element symbol, h, k, l indices,
        reflection positions and size into printable line.

        :return: string filled with provided reflection data.
        :rtype: str
        """
        label = '{}({},{},{})'.format(color, h_ind, k_ind, l_ind)
        pos = ' {: 7.5f} {: 7.5f} {: 7.5f}'.format(x_pos, y_pos, z_pos)
        size = ' {: 7.5f}\n'.format(u_iso)
        return '{:16}   1{} 11.0000{}'.format(label, pos, size)

    def write_res(self, path='hkl.res', colored='m'):
        """
        Write the reflection information in .res file according to the data
        passed to the artist object.

        :param colored: Which key of dataframe should be visualised using color.
        :type colored: str
        :param path: Absolute or relative path where the file should be saved
        :type path: str
        """
        color = rescale_list_to_other(self.df.data[colored], elements_list)
        h_ind = self.df.data['h']
        k_ind = self.df.data['k']
        l_ind = self.df.data['l']
        x_pos = (1.0 / self.maximum_index) * self.df.data['h']
        y_pos = (1.0 / self.maximum_index) * self.df.data['k']
        z_pos = (1.0 / self.maximum_index) * self.df.data['l']
        u_iso = rescale_list_to_range(self.df.data['F'],
                                      (self.MIN_ATOM_SIZE,
                                       self.MAX_ATOM_SIZE))
        zipped = zip(color, h_ind, k_ind, l_ind, x_pos, y_pos, z_pos, u_iso)

        file = open(path, 'w')
        file.write(self.res_header)
        file.write(self.res_legend(colored=colored))
        for c, h, k, l, x, y, z, u in zipped:
            file.write(self.res_line(color=c, h_ind=h, k_ind=k, l_ind=l,
                                     x_pos=x, y_pos=y, z_pos=z, u_iso=u))
        file.close()

# TODO method to export the style file using "inspect" module
# TODO get all fixed files to templates


if __name__ == '__main__':
    from kesshou.dataframes import HklFrame
    h1 = HklFrame()
    h1.read('/home/dtchon/_/shelxt2c_mod.hkl', 'tonto_I')
    h1.to_res('/home/dtchon/_/shelxt2c_mod.res', colored='h')

    # TODO Fix the documentation using this new object
    # TODO think about space groups...

    # TODO wrap table/data in getter/setter and make it automatically place,
    # TODO refresh, set keys etc.

    # TODO add point group and extinctions to object, update other methods.
