import copy
import random

import numpy as np
import numpy.linalg as lin
import pandas as pd

import hikari
from hikari.dataframes import BaseFrame
from hikari.symmetry import PG, SG
from hikari.resources import hkl_formats, hkl_aliases, hkl_mercury_style, \
    characteristic_radiation
from hikari.utility import cubespace, chemical_elements, make_abspath, \
    rescale_list_to_range, rescale_list_to_other

pd.options.mode.chained_assignment = 'raise'


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
        'description': 'Observed intensity',
        'imperative': False,
        'dtype': 'float64',
        'reduce_behaviour': 'average',
        'type': float
    }
    __Ic = {
        'default': 1.0,
        'description': 'Calculated intensity',
        'imperative': False,
        'dtype': 'float64',
        'reduce_behaviour': 'average',
        'type': float
    }
    __si = {
        'default': 0.0,
        'description': 'Uncertainty of intensity I determination',
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
    __ze = {
        'default': 0.0,
        'description': 'weighted diff. in observed I - calculated Ic intensity',
        'imperative': False,
        'dtype': 'float64',
        'reduce_behaviour': 'average',
        'type': float
    }
    __ze2 = {
        'default': 0.0,
        'description': 'weighted diff. in obs. I - calculated Ic int., squared',
        'imperative': False,
        'dtype': 'float64',
        'reduce_behaviour': 'average',
        'type': float
    }
    __Icsi = {
        'default': 1.0,
        'description': 'Calculated intensity of reflection div. by exp. sigma',
        'imperative': False,
        'dtype': 'float64',
        'reduce_behaviour': 'average',
        'type': float
    }
    __Iosi = {
        'default': 1.0,
        'description': 'Observed intensity of reflection div. by exp. sigma',
        'imperative': False,
        'dtype': 'float64',
        'reduce_behaviour': 'average',
        'type': float
    }
    __equiv = {
        'default': 0,
        'description': 'integer unique for each set of sym. equiv. reflections',
        'imperative': False,
        'dtype': 'int64',
        'reduce_behaviour': 'keep',
        'type': int
    }
    defined_keys = {'h', 'k', 'l', 'F', 'I', 'si', 'sf', 'b', 'm', 'la', 'ph',
                    'u', 'r', 't', 'u1', 'u2', 'u3', 'v1', 'v2', 'v3',
                    'x', 'y', 'z', 'ze', 'ze2', 'Iosi', 'Icsi', 'equiv'}
    # TODO check if needed and cplt?

    def __init__(self, keys=()):
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
    A master object which manages single-crystal diffraction files.
    It utilises other `Hkl*` classes to import, store, manipulate and output
    information about single-crystal diffraction patterns.

    HklFrame acts as an container which stores
    the diffraction data (Pandas dataframe, :attr:`table`)
    and elementary crystal cell data (:class:`hikari.dataframes.Base`).
    Demanding methods belonging to this class are vectorized,
    providing relatively satisfactory performance and high memory capacity.
    HklFrame methods are designed to work in-place, so the work strategy
    is to create a new instance of HklFrame for each reflection dataset,
    manipulate it using methods, eg. :func:`merge` or :func:`trim`, and
    :func:`copy` to other object or output using :func:`write` if needed.

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
        :param other: HklFrame to be added to data
        :type other: HklFrame
        :return: concatenated :attr:`table` dataframes, with metadata from first
        :rtype: HklFrame
        """
        _copied = self.copy()
        _copied.table = pd.concat([self.table, other.table], ignore_index=True)
        return _copied

    def __len__(self):
        """
        :return: Number of rows (individual reflections) in `self.data`
        :rtype: int
        """
        return self.table.shape[0]

    def __str__(self):
        """
        :return: Human-readable representation of `self.data`
        :rtype: str
        """
        return self.table.__str__()

    @property
    def la(self):
        """
        Wavelength of radiation used in the diffraction experiment.
        Can be set using popular abbreviations such as "MoKa" or "CuKb",
        where *a* and *b* stand for *alpha* and *beta*.
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
        try:
            self.__la = characteristic_radiation[wavelength[:4].lower()]
        except TypeError:
            self.__la = float(wavelength)

    @property
    def r_lim(self):
        """
        :return: Radius of limiting sphere calculated in A-1 based on :attr:`la`
        :rtype: float
        """
        return 2 / self.la

    def _in_dacs(self, opening_angle, vectors):
        oa = np.deg2rad(opening_angle)                # opening angle in radians
        xyz = self.table.loc[:, ('x', 'y', 'z')].to_numpy()
        r = self.table.loc[:, 'r'].to_numpy()
        v = np.array(vectors)
        v = (v.T / lin.norm(v, axis=1)).T              # normalise vectors v
        m1 = np.matmul(v, xyz.T)                       # dist from dac plane "p"
        phi = np.abs(np.arcsin((m1 / r).clip(-1, 1)))  # angle <(plane p, v)
        lim = self.r_lim * np.sin(oa - phi)            # True if in dac
        return r[None, :] < lim

    def dac_trim(self, opening_angle=35.0, vector=None):
        """
        Remove reflections outside the opening_angle DAC-accessible volume.
        Sample/DAC orientation can be supplied either via specifying crystal
        orientation in :class:`hikari.dataframes.BaseFrame`, in
        :attr:`orientation` or providing a xyz\* *vector* perpendicular to the
        dac-accessible disc. For further details refer to `*Tchoń & Makal, IUCrJ
        8, 1006-1017 (2021)* <https://doi.org/10.1107/s2052252521009532>`_.

        :param opening_angle: DAC single opening angle in degrees, default 35.0.
        :type opening_angle: float
        :param vector: Provides information about orientation of crystal
          relative to DAC. If None, current :attr:`orientation` is used instead.
        :type vector: Tuple[float]
        :return: HklFrame containing only reflections in dac-accessible region.
        :rtype: HklFrame
        """
        if vector is None:
            l_v = np.array((1.0, 0.0, 0.0))           # vector parallel to beam
            hkl = lin.inv(self.orientation) @ l_v     # calculate hkl vector
            v = hkl @ self.A_r                        # calculate xyz* vector
        else:
            v = vector
        in_dac = self._in_dacs(opening_angle, vectors=np.array([v, ]))[0]
        self.table = self.table[in_dac]
        self.table.reset_index(drop=True, inplace=True)

    def dacs_count(self, opening_angle=35.0, vectors=np.array((1, 0, 0))):
        """
        Count unique dac-accessible reflections for n crystals placed such that
        vector n is perpendicular to diamond. For details see :meth:`dac_trim`.

        :param opening_angle: DAC single opening angle in degrees, default 35.0.
        :type opening_angle: float
        :param vectors: Array containing rotational axes of available DAC-discs.
        :type vectors: np.array
        :return: Array with numbers of unique reflns in DAC-accessible region.
        :rtype: np.array
        """
        memory_estimate = 26 * len(self) * len(vectors)  # estimate memory use
        cycles_needed = -(memory_estimate // -hikari.MEMORY_SIZE)  # and split
        if cycles_needed > 1:
            vectors_split = np.array_split(vectors, cycles_needed)
            return np.hstack([self.dacs_count(opening_angle, vectors=v)
                              for v in vectors_split])
        else:
            in_dac = self._in_dacs(opening_angle, vectors)
            return np.array([self.table.loc[in_dac[n, :], 'equiv'].nunique()
                             for n in range(vectors.shape[0])])

    def copy(self):
        """
        :return: An exact deep copy of this HklFrame.
        :rtype: HklFrame
        """
        return copy.deepcopy(self)

    def extinct(self, space_group=SG['P1']):
        """
        Removes from dataframe reflections which should be extinct based on
        space :class:`hikari.symmetry.group.Group`. For ref. see ITC-A12.3.5.

        :param space_group: Space group used to extinct the reflections.
        :type space_group: hikari.symmetry.group.Group
        """
        hkls = self.table.loc[:, ['h', 'k', 'l']].to_numpy()
        extinct_flag_list = [o.extincts(hkls) for o in space_group.operations]
        extinct_flag_list_union = np.logical_or.reduce(extinct_flag_list)
        self.table = self.table[~extinct_flag_list_union]
        self.table.reset_index(drop=True, inplace=True)

    def find_equivalents(self, point_group=PG['1']):
        """
        Assign each reflection its symmetry equivalence identifier and store
        it in the `hikari.dataframes.HklFrame.data['equiv']` column.
        The ID is an integer unique for each set of equivalent reflections.

        In order to provide an information about equivalence, a *point_group*
        of reciprocal space must be provided (default PG['1']). Point groups
        and their notation can be found in :mod:`hikari.symmetry` sub-package.

        :param point_group: Point group used to determine symmetry equivalence
        :type point_group: hikari.symmetry.Group
        """
        inc = 10 ** (int(np.log10(self.HKL_LIMIT)) + 2)
        self.keys.add(('equiv',))
        self.table.reset_index(drop=True, inplace=True)
        self.table['equiv'] = -inc**3
        _hkl_matrix = self.table.loc[:, ('h', 'k', 'l')].to_numpy()
        for op in point_group.operations:
            new_equiv = op.transform(_hkl_matrix) @ np.array([inc ** 2, inc, 1])
            _to_update = self.table['equiv'] < new_equiv
            self.table.loc[_to_update, 'equiv'] = new_equiv[_to_update]

    def from_dict(self, dictionary):
        """
        Construct the self.data using information stored in dictionary.
        The dictionary keys must be valid strings, see :class:`HklKeys` for
        a list of valid keys. The dictionary values must be iterable of equal
        size, preferably `numpy.ndarray`.

        :param dictionary: Dictionary with "key - iterable of values" pairs.
        :type dictionary: Dict[str, numpy.ndarray]
        """
        df = pd.DataFrame()
        self.keys.add(dictionary.keys())
        for key, value in dictionary.items():
            typ = self.keys.get_property(key, 'dtype')
            df[key] = pd.Series(value, dtype=typ, name=key)
        self.table = df[(df['h'] != 0) | (df['k'] != 0) | (df['l'] != 0)].copy()
        if not('x' in self.table.columns):
            self.place()
        self._recalculate_structure_factors_and_intensities()

    def fill(self, radius=2.0):
        """
        Fill dataframe with all reflections within *radius* from space origin.

        :param radius: Maximum distance from the reciprocal space origin
            to placed reflection (in reciprocal Angstrom).
        :type radius: float
        """

        max_index = 25

        # make an initial guess of the hkl ball
        def _make_hkl_ball(i=max_index):
            hkl_grid = np.mgrid[-i:i:2j*i+1j, -i:i:2j*i+1j, -i:i:2j*i+1j]
            hkls = np.stack(hkl_grid, -1).reshape(-1, 3)
            xyz = np.matrix((self.a_w, self.b_w, self.c_w))
            return hkls[lin.norm(hkls @ xyz, axis=1) <= radius]
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

    def stats(self, bins=10, space_group=SG['P1']):
        """
        Returns completeness, redundancy, number of all, unique & theoretically
        possible reflections within equal-volume `bins` in given `space group`.

        :param bins: Number of equal-volume bins to divide the data into.
        :type bins: int
        :param space_group: Group used to calculate equivalence and extinctions.
        :type space_group: hikari.symmetry.Group
        :return: String containing table with stats as a function of resolution
        :rtype: str
        """

        point_group = space_group.reciprocate()

        hkl_base = self.copy()
        hkl_base.extinct(space_group)
        hkl_base.find_equivalents(point_group)
        hkl_base.table['_i_to_si'] = hkl_base.table['I'] / hkl_base.table['si']

        hkl_full = self.copy()
        hkl_full.fill(radius=max(self.table['r']))
        hkl_full.merge(point_group)
        hkl_full.extinct(space_group)
        hkl_full.find_equivalents(point_group)

        def group_by_resolution(hkl_):
            bin_limits = cubespace(0.0, max(self.table['r']), num=bins + 1)
            return hkl_.table.groupby(pd.cut(hkl_.table['r'], bin_limits))

        grouped_base = group_by_resolution(hkl_base)
        grouped_full = group_by_resolution(hkl_full)

        observed = grouped_base.size()
        independent = grouped_base['equiv'].nunique()
        theory = grouped_full['equiv'].nunique()
        cpl = independent.div(theory)
        red = observed.div(independent)
        i2si = grouped_base['_i_to_si'].mean()
        out = pd.concat([observed, independent, theory, i2si, cpl, red], axis=1)
        out.columns = ['Obser', 'Indep', 'Theory', 'I/si(I)', 'Cplt', 'Red.']
        return out  # use .reset_index().to_string(index=False) to flatten

    def merge(self, point_group=PG['1']):
        """
        Average down each set of redundant reflections present in the table,
        to one reflection. The redundancy is determined using the
        :meth:`find_equivalents` method with appropriate point group. Therefore,
        the merging can be used in different ways depending on point group:

        - For PG['1'], only reflections with exactly the same h, k and l indices
          will be merged. Resulting dataframe will not contain any duplicates.

        - For PG['-1'] reflections with the same h, k and l as well as their
          Friedel pairs will be merged together to one reflection.

        - For PG['mmm'] all equivalent reflections of "mmm" point group will be
          merged. Since "mmm" is centrosymmetric, Friedel pairs will be merged.

        - For PG['mm2'] symmetry-equivalent reflections within the "mmm"
          point group will be merged, but the Friedel pairs will be preserved.

        The procedure will have a different effect on different dataframe keys,
        depending on their "reduce_behaviour" specified in :class:`HklKeys`.
        Fixed parameters *h, k, l, x, y, z, r* and *equiv* will be preserved;
        Floating points such as intensity *I*, structure factor *F* and their
        uncertainties *si* and *sf* will be averaged using arithmetic mean;
        Multiplicity *m* will be summed; Other parameters which would lose their
        meaning such as batch number *b* will be discarded.

        The merging inevitably removes some information from the dataframe,
        but it can be necessary for some operations. For example, the drawing
        procedures work much better and provide clearer image if multiple points
        occupying the same position in space are reduced to one instance.

        :param point_group: Point Group used to determine symmetry equivalence
        :type point_group: hikari.symmetry.Group
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
        xyz = hkl @ self.A_r
        self.table.loc[:, 'x'] = xyz[:, 0]
        self.table.loc[:, 'y'] = xyz[:, 1]
        self.table.loc[:, 'z'] = xyz[:, 2]
        self.table.loc[:, 'r'] = lin.norm(xyz, axis=1)
        self.keys.add(('x', 'y', 'z', 'r'))

    def calculate_fcf_statistics(self):
        """
        Calculate values of zeta (I - Ic) / si on other stats based on contents
        of fcf files. Save new key and its values into the dataframe.
        """
        ze = (self.table['I'] - self.table['Ic']) / self.table['si']
        self.table['ze'] = ze
        self.table['ze2'] = ze ** 2
        self.table['Iosi'] = self.table['I'] / self.table['si']
        self.table['Icsi'] = self.table['Ic'] / self.table['si']
        self.keys.add(('ze', 'ze2', 'Iosi', 'Icsi'))

    def read(self, hkl_path, hkl_format='shelx_4'):
        """
        Read the contents of .hkl file as specified by path and format,
        and store them in the pandas dataframe in `self.data`.
        For a list of all available .hkl formats,
        please refer to :attr:`hikari.dataframes.HklIo.format`.

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
        point groups can be imported from :py:mod:`hikari.symmetry` module.

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
                hkl = (mat @ self.table.loc[:, ['h', 'k', 'l']].to_numpy().T).T
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
        HklToResConverter(self).convert(path)

    def trim(self, limit):
        """
        Remove reflections further than *limit* from reciprocal space origin.

        :param limit: Radius of the trimming sphere in reciprocal Angstrom
        :type limit: float
        """
        self.table = self.table.loc[self.table['r'] <= limit]

    def write(self, hkl_path, hkl_format='shelx_4'):
        """
        Write the contents of dataframe to a .hkl file using specified
        *path* and *format*.
        For a list of all available .hkl formats,
        please refer to :attr:`hikari.dataframes.HklIo.format`.

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
        self.file_path = make_abspath(hkl_file_path)
        self.formats_defined = hkl_formats
        self.formats_aliases = hkl_aliases
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
        | shelx_fcf| fcf      | h 4 k 4 l 4            | YES  | NO   | NO     |
        |          |          | Ic 14 I 14 si 13       | (*)  |      |        |
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

        Pre/suffixes denoted with (*) are not suported in terms of writing.
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

        if self.is_current_format_free:
            def parse_line(line):
                return self._parse_free_line(line)
        else:
            def parse_line(line):
                return self._parse_fixed_line(line)

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
        """
        needed_data = hkl_data.loc[:, self._format_dict['labels']].astype(str)
        with open(self.file_path, 'w') as hkl_file:
            hkl_file.write(self._format_dict['prefix'])
            for row_tuple in needed_data.itertuples(index=False):
                row_dict = dict(zip(self._format_dict['labels'], row_tuple))
                hkl_file.write(self._line_formatter.format(**row_dict))
            hkl_file.write(self._format_dict['suffix'])


class HklToResConverter:
    """A class responsible for representing hkl data using using .res format"""

    MIN_DISTANCE = 10.0
    ELEMENTS = chemical_elements[:100]
    MIN_U = 0.00001
    MAX_U = 4.99999

    def __init__(self, hkl_dataframe):
        self.df = hkl_dataframe

    @property
    def abc_scale_factor(self):
        return self.largest_absolute_hkl * self.MIN_DISTANCE \
               / min(self.df.a_r, self.df.b_r, self.df.c_r)

    @property
    def largest_absolute_hkl(self):
        return self.df.table[['h', 'k', 'l']].abs().max().max()

    @property
    def x(self):
        return self.df.table['h'] / self.largest_absolute_hkl

    @property
    def y(self):
        return self.df.table['k'] / self.largest_absolute_hkl

    @property
    def z(self):
        return self.df.table['l'] / self.largest_absolute_hkl

    @property
    def u(self):
        if 'F' in self.df.table.columns:
            return rescale_list_to_range(self.df.table['F'],
                                         (self.MIN_U, self.MAX_U))
        else:
            return [1.0] * len(self.df.table)

    @property
    def c(self):
        if 'm' in self.df.table.columns:
            return rescale_list_to_other(self.df.table['m'], self.ELEMENTS)
        else:
            return [self.ELEMENTS[0]] * len(self.df.table)

    @property
    def res_header(self):
        return "TITL Reflection visualisation\n" \
               "REM Special file to be used in mercury with hkl.msd style.\n" \
               "REM Reciprocal unit cell inflated by {sc} to prevent bonds\n" \
               "CELL {la:7f} {a:7f} {b:7f} {c:7f} {al:7f} {be:7f} {ga:7f}\n" \
               "LATT -1\n\n".format(sc=self.abc_scale_factor,
                                    la=self.df.la,
                                    a=self.df.a_r * self.abc_scale_factor,
                                    b=self.df.b_r * self.abc_scale_factor,
                                    c=self.df.c_r * self.abc_scale_factor,
                                    al=np.rad2deg(self.df.al_r),
                                    be=np.rad2deg(self.df.be_r),
                                    ga=np.rad2deg(self.df.ga_r))

    @property
    def atom_list(self):
        cols = [self.c, self.df.table['h'], self.df.table['k'],
                self.df.table['l'], self.x, self.y, self.z, self.u]
        return [list(row) for row in zip(*cols)]

    @staticmethod
    def res_line(_c, _h, _k, _l, _x, _y, _z, _u):
        """
        :return: res line filled based on input element, hkl, xyz and u.
        :rtype: str
        """
        label = f'{_c}({_h},{_k},{_l})'
        pos = f' {_x: 7.5f} {_y: 7.5f} {_z: 7.5f}'
        return f'{label:16}   1{pos} 11.0 {_u: 7.5f}\n'

    def convert(self, path='~'):
        with open(make_abspath(path), 'w') as res:
            res.write(self.res_header)
            for row in self.atom_list:
                res.write(self.res_line(*row))
        with open(make_abspath(path, '../hkl.msd'), 'w') as style:
            style.write(hkl_mercury_style)

    # TODO wrap table/data in getter/setter and make it automatically place,
    # TODO refresh, set keys etc.
