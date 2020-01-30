from collections import OrderedDict
from kesshou.utility import cubespace, is2n, is3n, is4n, is6n
from kesshou.symmetry.pointgroup import *
import copy
import random
import struct
import sys
import numpy as np
import numpy.linalg as lin
import pandas as pd
import matplotlib.cm
import matplotlib.pyplot as plt


class HklCrystal:
    """Simple class for storing crystal information.

    Space          | Direct        | Reciprocal    |
    Scalar         | *_d           | *_r           |
    Vector         | *_v           | *_w           |

    a, b, c        | unit cell lengths in Angstrom / Angstrom^-1
    x, y, z        | normalised unit cell vectors
    al, be, ga     | unit cell angles in radians
    v              | unit cell volume in Angstrom^3 / Angstrom^-3"""

    def __init__(self):
        # DIRECT SPACE SCALARS
        self.__a_d = 1.
        self.__b_d = 1.
        self.__c_d = 1.
        self.__al_d = np.pi / 2
        self.__be_d = np.pi / 2
        self.__ga_d = np.pi / 2
        self.__v_d = 1.
        # DIRECT SPACE VECTORS
        self.__a_v = np.array((1.0, 0.0, 0.0))
        self.__b_v = np.array((0.0, 1.0, 0.0))
        self.__c_v = np.array((0.0, 0.0, 1.0))
        # RECIPROCAL SPACE SCALARS
        self.__a_r = 1.
        self.__b_r = 1.
        self.__c_r = 1.
        self.__al_r = np.pi / 2
        self.__be_r = np.pi / 2
        self.__ga_r = np.pi / 2
        self.__v_r = 1.
        # DIRECT SPACE VECTORS
        self.__a_w = np.array((1.0, 0.0, 0.0))
        self.__b_w = np.array((0.0, 1.0, 0.0))
        self.__c_w = np.array((0.0, 0.0, 1.0))
        # OTHERS
        self.orient_matrix = np.array(((1.0, 0.0, 0.0),
                                       (0.0, 1.0, 0.0),
                                       (0.0, 0.0, 1.0)))

    def edit_cell(self, **parameters):
        """Edit direct space unit cell using dictionary with the following keys:
        a, b, c [in Angstrom] and al, be, ga [in degrees or radians]."""

        # INSERT NEW VALUES OF DIRECT CELL PARAMETERS
        for key, value in parameters.items():
            key = key[:2]
            try:
                value = max((0.01, float(value)))
            except TypeError:
                value = max((0.01, value.nominal_value))
            if key in ('a', 'b', 'c'):
                key = '_HklCrystal__{}_d'.format(key)
                setattr(self, key, value)
            elif key in ('al', 'be', 'ga'):
                key = '_HklCrystal__{}_d'.format(key)
                if value < np.pi:
                    setattr(self, key, value)
                else:
                    setattr(self, key, np.radians(value))

        # RECALCULATE REST OF THE PARAMETERS
        # DIRECT CELL SCALARS
        self.__v_d = self.a_d * self.b_d * self.c_d * \
               (1 - np.cos(self.al_d) ** 2 - np.cos(self.be_d) ** 2
                - np.cos(self.ga_d) ** 2 + 2 * np.cos(self.al_d)
                * np.cos(self.be_d) * np.cos(self.ga_d)) ** 0.5
        # DIRECT CELL VECTORS
        self.__a_v = np.array((self.a_d, 0, 0))
        self.__b_v = np.array((self.b_d * np.cos(self.ga_d),
                        self.b_d * np.sin(self.ga_d), 0))
        self.__c_v = np.array((self.c_d * np.cos(self.be_d),
                        self.c_d * (np.cos(self.al_d) - np.cos(self.be_d)
                        * np.cos(self.ga_d)) / np.sin(self.ga_d),
                        self.v_d / (self.a_d * self.b_d * np.sin(self.ga_d))))
        # RECIPROCAL CELL SCALARS
        self.__a_r = self.b_d * self.c_d * np.sin(self.al_d) / self.v_d
        self.__b_r = self.c_d * self.a_d * np.sin(self.be_d) / self.v_d
        self.__c_r = self.a_d * self.b_d * np.sin(self.ga_d) / self.v_d
        self.__al_r = np.arccos((np.cos(self.be_d) * np.cos(self.ga_d)
                                 - np.cos(self.al_d)) /
                                (np.sin(self.be_d) * np.sin(self.ga_d)))
        self.__be_r = np.arccos((np.cos(self.ga_d) * np.cos(self.al_d)
                                 - np.cos(self.be_d)) /
                                (np.sin(self.ga_d) * np.sin(self.al_d)))
        self.__ga_r = np.arccos((np.cos(self.al_d) * np.cos(self.be_d)
                                 - np.cos(self.ga_d)) /
                                (np.sin(self.al_d) * np.sin(self.be_d)))
        self.__v_r = self.a_r * self.b_r * self.c_r * \
                     (1 - np.cos(self.al_r) ** 2 - np.cos(self.be_r) ** 2
                      - np.cos(self.ga_r) ** 2 + 2 * np.cos(self.al_r)
                      * np.cos(self.be_r) * np.cos(self.ga_r)) ** 0.5
        # RECIPROCAL CELL VECTORS
        self.__a_w = np.cross(self.b_v, self.c_v) / self.v_d
        self.__b_w = np.cross(self.c_v, self.a_v) / self.v_d
        self.__c_w = np.cross(self.a_v, self.b_v) / self.v_d

    def import_from_frame(self, frame):
        """Import necessary crystal parameters from CifFrame"""

        # IMPORT AND CHANGE LATTICE PARAMETERS
        new_parameters = {
            'a': float(frame.data['_cell_length_a']),
            'b': float(frame.data['_cell_length_b']),
            'c': float(frame.data['_cell_length_c']),
            'al': float(frame.data['_cell_angle_alpha']),
            'be': float(frame.data['_cell_angle_beta']),
            'ga': float(frame.data['_cell_angle_gamma'])}
        self.edit_cell(**new_parameters)

        # IMPORT AND CHANGE ORIENTATION MATRIX
        try:
            self.orient_matrix = \
                np.array(((float(frame.data['_diffrn_orient_matrix_UB_11']),
                        float(frame.data['_diffrn_orient_matrix_UB_12']),
                        float(frame.data['_diffrn_orient_matrix_UB_13'])),
                       (float(frame.data['_diffrn_orient_matrix_UB_21']),
                        float(frame.data['_diffrn_orient_matrix_UB_22']),
                        float(frame.data['_diffrn_orient_matrix_UB_23'])),
                       (float(frame.data['_diffrn_orient_matrix_UB_31']),
                        float(frame.data['_diffrn_orient_matrix_UB_32']),
                        float(frame.data['_diffrn_orient_matrix_UB_33']))))
        except KeyError:
            pass

    @property
    def a_d(self):
        return self.__a_d

    @property
    def b_d(self):
        return self.__b_d

    @property
    def c_d(self):
        return self.__c_d

    @property
    def al_d(self):
        return self.__al_d

    @property
    def be_d(self):
        return self.__be_d

    @property
    def ga_d(self):
        return self.__ga_d

    @property
    def v_d(self):
        return self.__v_d

    @property
    def a_v(self):
        return self.__a_v

    @property
    def b_v(self):
        return self.__b_v

    @property
    def c_v(self):
        return self.__c_v

    @property
    def x_v(self):
        return self.__a_v / lin.norm(self.__a_v)

    @property
    def y_v(self):
        return self.__b_v / lin.norm(self.__b_v)

    @property
    def z_v(self):
        return self.__c_v / lin.norm(self.__c_v)

    @property
    def a_r(self):
        return self.__a_r

    @property
    def b_r(self):
        return self.__b_r

    @property
    def c_r(self):
        return self.__c_r

    @property
    def al_r(self):
        return self.__al_r

    @property
    def be_r(self):
        return self.__be_r

    @property
    def ga_r(self):
        return self.__ga_r

    @property
    def v_r(self):
        return self.__v_r

    @property
    def a_w(self):
        return self.__a_w

    @property
    def b_w(self):
        return self.__b_w

    @property
    def c_w(self):
        return self.__c_w

    @property
    def x_w(self):
        return self.__a_w / lin.norm(self.__a_w)

    @property
    def y_w(self):
        return self.__b_w / lin.norm(self.__b_w)

    @property
    def z_w(self):
        return self.__c_w / lin.norm(self.__c_w)


class HklKeys:
    """Object managing .hkl file keys for HklFrame"""
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
        'reduce_behaviour': 'add',
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
    defined_keys = {'h', 'k', 'l', 'F', 'I', 'si', 'b', 'm', 'la', 'ph',
                    'u', 'r', 't', 'u1', 'u2', 'u3', 'v1', 'v2', 'v3',
                    'x', 'y', 'z'}

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


class HklFrame:
    """Single crystal diffraction pattern container"""

    def __init__(self):
        self.crystal = HklCrystal()
        self.data = pd.DataFrame()
        self.keys = HklKeys()
        self.meta = {'wavelength': 0.71069}

    def __len__(self):
        return self.data.shape[0]

    def __add__(self, other):
        c = copy.deepcopy(self)
        c.data = pd.concat([self.data, other.data], ignore_index=True)
        return c

    @property
    def la(self):
        return self.meta['wavelength']

    @property
    def r_lim(self):
        return 2 / self.meta['wavelength']

    def condition(self, equation=''):
        """This method returns part of self.data for which equation is True"""
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

    def domain(self, address='hkl'):
        """This method returns part of self.data which lies within address"""
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

    @staticmethod
    def interpret_hkl_format(hkl_format):
        """Interpret hkl format int / dict / OrderedDict and return
        format_string (for parser) and column labels list"""
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
            raise TypeError(
                'Format type should be 2, 3, 4, 5, 6, "XD", "TONTO" or dict')
        return format_string, column_labels, file_prefix, file_suffix, zero_line

    def calculate_intensity_from_structure_factor(self):
        self.data['I'] = self.data.apply(lambda row: row['F']**2, axis=1)
        self.data['si'] = self.data.apply(lambda row: 2*row['si']*row['F'],
                                          axis=1)

    def from_dict(self, dictionary):
        """Produce pd DataFrame from dictionary of values"""
        new_data = pd.DataFrame()
        self.keys.add(dictionary.keys())
        for key, value in dictionary.items():
            typ = self.keys.get_property(key, 'dtype')
            new_data[key] = pd.Series(value, dtype=typ, name=key)
        self.data = new_data
        self._place()

    def edit_wavelength(self, wavelength):
        """Define wavelength: "CuKa", "Ag", "Polichromatic" or custom in A"""
        sources = {'cr': 2.2896,	'fe': 1.9360,	'co': 1.79,
                   'cu': 1.54056,	'mo': 0.71073,	'zr': 0.69,
                   'ag': 0.56087}
        try:
            self.meta['wavelength'] = sources[wavelength[:2].lower()]
        except TypeError:
            # IF CUSTOM VALUE WAS GIVEN
            self.meta['wavelength'] = float(wavelength)
        except KeyError:
            # NEGATIVE WAVELENGTH SERVES AS SCALE FOR INDIVIDUAL REFLECT. VALUES
            if wavelength[:3] in {'pol', 'Pol', 'POL', 'P'}:
                self.meta['wavelength'] = -1.0

    def _place(self):
        """Assign reflections their positions in reciprocal space (x, y, z)
        and calculate their distance from origin (r) in reciprocal Angstrom"""
        hkl = self.data.loc[:, ('h', 'k', 'l')].to_numpy()
        abc = np.matrix((self.crystal.a_w, self.crystal.b_w, self.crystal.c_w))
        xyz = hkl @ abc
        self.data['x'] = xyz[:, 0]
        self.data['y'] = xyz[:, 1]
        self.data['z'] = xyz[:, 2]
        self.data['r'] = lin.norm(xyz, axis=1)

    def resymmetrify(self, operations=tuple(), merge=True):
        # prepare empty hkl holder
        q = copy.deepcopy(self)
        # For each declared symmetry operation add hkl * operation
        for operation in operations:
            r = copy.deepcopy(self)
            r.transform(operation)
            q = q + r
        if merge is True:
            q.merge()
        # return the result to the dataframe
        self.data = q.data

    def transform(self, matrix):
        """Transform reflection indices using given 3x3 or 4x4 matrix"""
        mat = matrix[0:3, 0:3]
        hkl = self.data.loc[:, ('h', 'k', 'l')].to_numpy()
        hkl = hkl @ mat
        self.data['h'] = hkl[:, 0]
        self.data['k'] = hkl[:, 1]
        self.data['l'] = hkl[:, 2]
        self._place()

    def read(self, hkl_path, hkl_format):
        """Read .hkl file as specified by path and fields
        format: either ordered dictionary with specified fields (minus = ignore)
        or type number"""

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

        # OPEN HKL FILE AND PREPARE CONTAINER
        hkl_file = open(hkl_path, 'r')
        hkl_content = dict()

        # LOAD HKL COLUMN TAGS IF SUPERTYPE IS "XD"
        if format_string is 'XD':
            title_line = hkl_file.readline().strip().split()
            hkl_id = title_line[0]
            hkl_mainkey = 'F' if title_line[1] == 'F' else 'I'
            if title_line[2] != 'NDAT':
                raise KeyError('Loaded hkl file is not of "XD" type.')
            if int(title_line[3]) == -7:
                column_labels = ('h', 'k', 'l', 'b', hkl_mainkey, 'si', 'ph')
            elif 6 <= int(title_line[3]) <= 13:
                column_labels = ('h', 'k', 'l', 'b', hkl_mainkey, 'si',
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
        forgottens = tuple(self.keys.imperatives - set(column_labels))
        for forgotten in forgottens:
            default = self.keys.get_property(forgotten, 'default')
            hkl_content[forgotten] = [default] * hkl_checksum

        # PRODUCE PANDAS DATAFRAME
        self.from_dict(hkl_content)

    def rescale(self, key, factor):
        self.data[key] = self.data.apply(lambda row: row[key] * factor, axis=1)

    def write(self, hkl_path, hkl_format, columns_separator=True):
        """Write .hkl file as specified by path and write_format"""

        # PREPARE OBJECTS RESPONSIBLE FOR WRITING OUTPUT
        format_string, column_labels, file_prefix, file_suffix, zero_line = \
            self.interpret_hkl_format(hkl_format)
        if format_string == 'XD':
            format_string, column_labels, file_prefix, file_suffix, zero_line = \
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

        # WRITE PREFIX LINE
        hkl_file = open(hkl_path, 'w')
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

    def copy(self, empty=False):
        """Make new dataframe which is an exact copy of self.
         Remove all data from self.data if empty is True"""
        new_dataframe = copy.deepcopy(self)
        if empty:
            new_dataframe.extinct()
        return new_dataframe

    def dac(self, opening_angle=40, vector=None):
        """Cut the reflections based on DAC angle (degrees) and
        orientation matrix or perpendicular vector in rec. space (if given).
        Please provide angle from normal vector, not double opening angle."""
        # make the abbreviations for opening angle, self.data and self.crystal
        xl = self.crystal
        oa = np.radians([float(opening_angle)])[0]
        # extinct 000 and make sure
        self.extinct('000')
        assert 'x' in self.data.columns
        # calculate normal vector "n" in reciprocal lattice
        if vector is None:
            l_v = np.array((1.0, 0.0, 0.0))      # vector parallel to beam
            _UB = self.crystal.orient_matrix   # import UB orientation matrix
            h = np.dot(lin.inv(_UB), l_v)        # calculate plane hkl indices
            n = h[0] * xl.a_w + h[1] * xl.b_w + h[2] * xl.c_w # perp. vector
        else:
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

    def thin_out(self, target_completeness=1.0, exponentially=False):
        """Randomly delete reflections to relative desired completeness"""

        if exponentially is False:
            # check whether acceptance is between 0 and 1
            if not 0.0 <= target_completeness <= 1.0:
                raise ValueError('acceptance parameter outside the 0--1 range')

            # create a list of reflections to be deleted
            total_indices = self.data.shape[0]
            number_to_delete = int((1-target_completeness) * total_indices)
            indices_to_delete = \
                random.sample(range(0, total_indices), number_to_delete)

            # delete chosen reflections
            self.data.drop(indices_to_delete, inplace=True)
            self.data.reset_index(drop=True, inplace=True)

        if exponentially is True:
            # the old algorithm, doesn't give exact value of completeness
            falloff = 0.2
            indices_to_delete = []
            for index, reflection in self.data.iterrows():
                v = np.array((reflection['x'], reflection['y'], reflection['z']))
                if random.random() > np.exp(-lin.norm(v) * falloff):
                    indices_to_delete.append(index)
            self.data.drop(indices_to_delete, inplace=True)
            self.data.reset_index(drop=True, inplace=True)
            # TODO vectorize (low priority)

    def trim(self, limit):
        """Cut the reflection further then the limit in A-1"""
        self.data = self.data.loc[self.data['r'] <= limit]

    def merge(self):
        """Average down redundant reflections, e.g. for drawing"""

        # group the dataframe and obtain all existing keys
        grouped = self.data.groupby(['h', 'k', 'l'])
        grouped_first = grouped.first().reset_index()
        grouped_mean = grouped.mean().reset_index()
        grouped_sum = grouped.sum().reset_index()
        keys = self.data.keys()

        # import reduction behaviours for all keys
        _reduction_kept = self.keys.reduce_behaviour['keep']
        _reduction_added = self.keys.reduce_behaviour['add']
        _reduction_averaged = self.keys.reduce_behaviour['average']

        # for each key apply a necessary reflections and add it to data
        data = dict()
        for key in keys:
            if key in _reduction_kept:
                data[key] = grouped_first[key]
            elif key in _reduction_added:
                data[key] = grouped_sum[key]
            elif key in _reduction_averaged:
                data[key] = grouped_mean[key]

        # create dataframe based on obtained data
        self.from_dict(data)

    def overwrite_values(self, other, keys):
        """Take one merged hkl and overwrite its keys to other merged hkl's"""

        # MAKE DEEPCOPIES OF BOTH HKLS AND SORT THEM H-K-L-WISE
        data1 = copy.deepcopy(self.data)
        data2 = copy.deepcopy(other.data)
        data1.sort_values(['h', 'k', 'l'], inplace=True)
        data2.sort_values(['h', 'k', 'l'], inplace=True)
        indices_to_delete = list()

        # GO THROUGH LISTS ELEMENT-WISE UNTIL THE FIRST IS EXHAUSTED
        index1, index2 = 0, 0

        while True:
            try:
                hkl1 = tuple(data1.loc[index1, ['h', 'k', 'l']])
                hkl2 = tuple(data2.loc[index2, ['h', 'k', 'l']])
            except KeyError:
                break
            else:
                if hkl1 == hkl2:
                    for key in keys:
                        data1.at[index1, key] = data2.loc[index2, key]
                    index1 += 1
                    index2 += 1
                    continue
                if hkl1 > hkl2:
                    index2 += 1
                    continue
                if hkl1 < hkl2:
                    indices_to_delete.append(hkl1)
                    index1 += 1
                    continue
                else:
                    raise ValueError('Problem with merging hkl files')

        # DELETE NOT OVERWRITTEN ELEMENTS AND RETURN NEW DATAFRAME
        new_dataframe = copy.deepcopy(self)
        new_dataframe.data.drop(indices_to_delete, inplace=True)
        new_dataframe.data.reset_index(drop=True, inplace=True)
        data1.drop(indices_to_delete, inplace=True)
        data1.reset_index(drop=True, inplace=True)
        self.data = data1
        self._place()

    def draw(self, alpha=False, colored='b', dpi=600, legend=True,
             master_key='I', projection=('h', 'k', 0),
             savepath=False, scale=1.0, showfig=False):
        """Draw a cross-section of reciprocal lattice for given pattern

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
        """

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

    def make_ball(self, radius=2.0, hkl_limit=400):
        """Generate points in a sphere of given radius"""
        # prepare necessary abbreviations of known properties
        a, b, c = self.crystal.a_w, self.crystal.b_w, self.crystal.c_w
        max_index = 25

        # function to make a hkl rectangle
        def _make_hkl_ball(i=max_index, r=radius):
            hkl_grid = np.mgrid[-i:i:2j*i+1j, -i:i:2j*i+1j, -i:i:2j*i+1j]
            _h = np.concatenate(np.concatenate(hkl_grid[0]))
            _k = np.concatenate(np.concatenate(hkl_grid[1]))
            _l = np.concatenate(np.concatenate(hkl_grid[2]))
            _hkl = np.matrix((_h, _k, _l)).T
            return _hkl[lin.norm(_hkl @ np.matrix((a, b, c)), axis=1) <= r]

        # grow the dataframe until all needed points are in
        hkl = _make_hkl_ball()
        previous_length = -1
        while len(hkl) > previous_length and max_index <= hkl_limit:
            previous_length = len(hkl)
            max_index = max_index * 2
            hkl = _make_hkl_ball(max_index)

        # create new dataframe using obtained ball of data
        _h, _k, _l = np.vsplit(hkl.T, 3)
        ones = np.ones_like(np.array(_h)[0])
        data = {'h': np.array(_h)[0],
                'k': np.array(_k)[0],
                'l': np.array(_l)[0],
                'I': ones, 'si': ones, 'm': ones}
        self.from_dict(data)

    def to_hklres(self, colored='m', master_key='I', path='hkl.res'):

        # prepare minima and maxima for scaling purposes
        data_minima, data_maxima = dict(), dict()
        for key in list(self.keys.all):
            data_minima[key] = self.data[key].min()
            data_maxima[key] = self.data[key].max()

        hkl_minimum = min((data_minima['h'], data_minima['k'], data_minima['l']))
        hkl_maximum = max((data_maxima['h'], data_maxima['k'], data_maxima['l']))
        scale = 2./max((abs(hkl_maximum), abs(hkl_minimum)))

        # print the title line
        a = self.crystal.a_r * 1000
        b = self.crystal.b_r * 1000
        c = self.crystal.c_r * 1000
        al = np.degrees(self.crystal.al_r)
        be = np.degrees(self.crystal.be_r)
        ga = np.degrees(self.crystal.ga_r)
        la = self.meta['wavelength']
        cell_list = (la, a, b, c, al, be, ga)
        hklres_file = open(path, 'w')
        hklres_file.write('TITL hkl visualisation\n')
        hklres_file.write('REM special hkl visualisation file, to be used in'
                          'mercury with hkl.msd style applied\n')
        hklres_file.write('REM reciprocal unit cell has been inflated '
                          'thousandfold for technical reasons\n')
        hklres_file.write('CELL {0:7f} {1:7f} {2:7f} {3:7f} {4:7f} {5:7f} '
                          '{6:7f}\n'.format(*cell_list))
        hklres_file.write('LATT -1\n\n')

        # define function for getting good label initial - colour
        # labels = ('H', 'He',
        #          'Li', 'Be', 'B', 'C', 'N', 'O', 'F', 'Ne',
        #          'Na', 'Mg', 'Al', 'Si', 'P', 'S', 'Cl', 'Ar')     # 18c.scale
        labels = ('H', 'Li', 'B', 'N', 'F', 'Na', 'Al', 'P', 'Cl')   # 9c.scale
        # labels = ('H')                                             # all red

        def get_label(integer):
            return labels[(int(integer)-1) % len(labels)]

        # for each reflection write a respective line to .res
        for index, reflection in self.data.iterrows():
            line_pars = dict()
            line_pars['label'] = '{atom}({h},{k},{l})'.format(
                atom=get_label(reflection[colored]),
                h=int(reflection['h']),
                k=int(reflection['k']),
                l=int(reflection['l']))
            line_pars['x'] = float(reflection['h'] * scale)
            line_pars['y'] = float(reflection['k'] * scale)
            line_pars['z'] = float(reflection['l'] * scale)
            line_pars['size'] = np.log(abs(reflection[master_key])+1.0) ** 2 \
                                * 5 * scale ** 2
            hklres_file.write('{label:16}   1 {x: .4f} {y: .4f} {z: .4f} '
                              '11.0000 {size: .12f}\n'.format(**line_pars))

        # close the file
        hklres_file.close()

    def make_stats(self, bins=10, point_group=PG_1, extinctions=(('000', ''),)):
        """This method analyses dataframe in terms of no. of reflections,
        Rint, completeness, redundancy in 'bins' resolution shells."""

        # get max resolution for later calculations and extinct data
        max_resolution = max(self.data['r'])

        # extinct wrong reflections
        for domain, condition in extinctions:
            self.extinct(domain=domain, condition=condition)

        # prepare merged, base, full, resymmetrified merged and
        # resymmetrified unmerged dataframe
        hkl_base = copy.deepcopy(self)
        hkl_full = copy.deepcopy(self)
        hkl_full.make_ball(radius=max_resolution)
        for domain, condition in extinctions:
            hkl_full.extinct(domain=domain, condition=condition)
        hkl_merged = copy.deepcopy(self)
        hkl_merged.merge()
        hkl_resy_u = copy.deepcopy(self)
        hkl_resy_u.resymmetrify(operations=point_group.operations, merge=False)
        hkl_resy_m = copy.deepcopy(hkl_resy_u)
        hkl_resy_m.merge()

        # group the dataframes by resolution
        res_limits = cubespace(0.0, max_resolution, num=bins+1)
        hkl_base_res_bins = pd.cut(hkl_base.data['r'], res_limits)
        hkl_base_by_res = hkl_base.data.groupby(hkl_base_res_bins)
        hkl_full_res_bins = pd.cut(hkl_full.data['r'], res_limits)
        hkl_full_by_res = hkl_full.data.groupby(hkl_full_res_bins)
        hkl_merged_res_bins = pd.cut(hkl_merged.data['r'], res_limits)
        hkl_merged_by_res = hkl_merged.data.groupby(hkl_merged_res_bins)
        hkl_resy_u_res_bins = pd.cut(hkl_resy_u.data['r'], res_limits)
        hkl_resy_u_by_res = hkl_resy_u.data.groupby(hkl_resy_u_res_bins)
        hkl_resy_m_res_bins = pd.cut(hkl_resy_m.data['r'], res_limits)
        hkl_resy_m_by_res = hkl_resy_m.data.groupby(hkl_resy_m_res_bins)

        # count observed, independent and theory reflns for each res. shell
        observed = hkl_base_by_res.size()
        independent = hkl_merged_by_res.size()
        theory = hkl_full_by_res.size()

        # calculate completeness and redundancy in P1
        completeness_nosymm = independent.div(theory)
        redundancy_nosymm = observed.div(independent)

        # count observed and independent after resymmetrization
        observed_after_resy = hkl_resy_u_by_res.size()
        independent_after_resy = hkl_resy_m_by_res.size()

        # calculate completeness and redundancy in point_grpup
        completeness_symm = independent_after_resy.div(theory)
        redundancy_symm = observed_after_resy.div(independent_after_resy)

        # make a final results table
        results = pd.concat([observed, independent, theory,
                             completeness_nosymm, redundancy_nosymm,
                             completeness_symm, redundancy_symm],
                            axis=1)
        results.columns = ['Obser', 'Indep', 'Theory',
                           'Cplt_P1', 'Redu_P1', 'Cplt', 'Redu']
        print(results)
        # TODO last big function to optimise (low priority)

    def _extinct(self, domain='hkl', condition=''):
        """This extinct is faster, but does not work for 'hhl'-type reflns """
        # obtain a copy of dataframe containing only reflections within domain
        dom = self.domain(domain)
        # obtain a part of dataframe which meets the condition
        con = self.condition(condition)
        # obtain a part of domain which does NOT meet the condition
        dom = dom.loc[~dom.isin(con).all(1)]
        # remove from self.data refls in domain which do NOT meet the condition
        self.data = self.data.loc[~self.data.isin(dom).all(1)]
        # reset the indices for other methods to use
        self.data.reset_index(drop=True, inplace=True)

    def extinct(self, domain='hkl', condition='', laue_group=PG_1):
        """Extinct all reflections which meet condition lll:rrr.
        domain describes area which is affected by condition [hkl]
        condition describes reflections which are preserved in domain [none]"""
        # obtain a copy of dataframe containing only reflections within domain
        dom = self.domain(domain)
        # obtain a part of dataframe which meets the condition
        con = self.condition(condition)
        # obtain a part of domain which does NOT meet the condition
        dom = dom.loc[~dom.isin(con).all(1)][['h', 'k', 'l']]
        # make a set with all hkls which should be extinct
        extinct_set = set()
        for op in laue_group.operations:
            extinct_set = extinct_set | {tuple(row) for row in dom.to_numpy() @ op}
        # remove all reflections whose hkls are in extinct set
        for h, k, l in list(extinct_set):
            self.data = self.data[~((self.data['h'] == h) &
                                    (self.data['k'] == k) &
                                    (self.data['l'] == l))]
        # reset the indices for other methods to use
        self.data.reset_index(drop=True, inplace=True)


if __name__ == '__main__':
    pass
    # TODO make resymmetrify and scripts use point group instead of operations
    # TODO add extintion functionality to make_stats module

# TODO 3D call visualise and to pyqtplot (low priority)

# TODO prepare and check installing routines
# TODO manual

# TODO join methods 'transform' and 'resymmetrify'
# TODO resymmetrify uses identity matrix twice
# TODO (redundancy after resymm is one larger than it should be)
