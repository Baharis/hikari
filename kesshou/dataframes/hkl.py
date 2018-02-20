from collections import OrderedDict
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
        self.__a_r = 1 / self.a_d
        self.__b_r = 1 / self.b_d
        self.__c_r = 1 / self.c_d
        self.__al_r = np.pi - self.al_d
        self.__be_r = np.pi - self.be_d
        self.__ga_r = np.pi - self.ga_d
        self.__v_r = 1 / self.v_d
        # RECIPROCAL CELL VECTORS
        self.__a_w = np.array((self.a_r, 0, 0))
        self.__b_w = np.array((self.b_r * np.cos(self.ga_r),
                        self.b_r * np.sin(self.ga_r), 0))
        self.__c_w = np.array((self.c_r * np.cos(self.be_r),
                        self.c_r * (np.cos(self.al_r) - np.cos(self.be_r)
                        * np.cos(self.ga_r)) / np.sin(self.ga_r),
                        self.v_r / (self.a_r * self.b_r * np.sin(self.ga_r))))

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
        'default': 0,
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
        self.meta = {'wavelength': 0.71069}
        self.keys = HklKeys()

    def __len__(self):
        return self.data.shape[0]

    @staticmethod
    def interpret_hkl_format(hkl_format):
        """Interpret hkl format int / dict / OrderedDict and return
        format_string (for parser) and column labels list"""
        if hkl_format == 2:
            column_labels = ('h', 'k', 'l', 'I', 'si', 'b', 'la')
            format_string = '4s 4s 4s 8s 8s 4s 8s'
        elif hkl_format == 3:
            column_labels = ('h', 'k', 'l', 'F', 'si', 'b')
            format_string = '4s 4s 4s 8s 8s 4s'
        elif hkl_format == 4:
            column_labels = ('h', 'k', 'l', 'I', 'si', 'b')
            format_string = '4s 4s 4s 8s 8s 4s'
        elif hkl_format == 5:
            column_labels = ('h', 'k', 'l', 'I', 'si', 'c')
            format_string = '4s 4s 4s 8s 8s 4s'
        elif hkl_format == 6:
            column_labels = ('h', 'k', 'l', 'I', 'si', 'm')
            format_string = '4s 4s 4s 8s 8s 4s'
        elif hkl_format in ('xd', 'Xd', 'xD', 'XD'):
            column_labels, format_string = tuple(), 'XD'
        elif type(hkl_format) in (dict, OrderedDict):
            column_labels = list()
            format_string = str()
            if type(hkl_format) is dict and sys.version_info[0] < 3:
                format_items = hkl_format.iteritems()
            else:
                format_items = hkl_format.items()
            for key, value in format_items:
                column_labels.append(key)
                if int(value) > 0:
                    format_string += value + 's '
                else:
                    format_string += str(abs(int(value))) + 'x '
            column_labels = tuple(column_labels)
            format_string.rstrip(' ')
        else:
            raise TypeError(
                'Format type should be 2, 3, 4, 5, 6, "XD" or dict')
        return format_string, column_labels

    def data_from_dict(self, dictionary):
        """Produce pd DataFrame from dictionary of values"""
        new_data = pd.DataFrame()
        for key, value in dictionary.items():
            typ = self.keys.get_property(key, 'dtype')
            new_data[key] = pd.Series(value, dtype=typ, name=key)
        return new_data

    def edit_wavelength(self, wavelength):
        """Define wavelength [e.g. "CuKa", "Ag", "Individual" or custom in A]"""
        sources = {'Cr': 2.2896,	'Fe': 1.9360,	'Co': 1.79,
                   'Cu': 1.54056,	'Mo': 0.71069,	'Zr': 0.69,
                   'Ag': 0.56087}
        try:
            self.meta['wavelength'] = sources[wavelength[:2]]
        except TypeError:
            # IF CUSTOM VALUE WAS GIVEN
            self.meta['wavelength'] = float(wavelength)
        except KeyError:
            # NEGATIVE WAVELENGTH SERVES AS SCALE FOR INDIVIDUAL REFLECT. VALUES
            if wavelength[:3] in {'ind', 'Ind', 'IND', 'i'}:
                self.meta['wavelength'] = -1.0

    def seek_reflection(self, **param):
        """Find indices of reflections described by dict of keys and values"""
        positions = []
        for index, reflection in self.data.iterrows():
            if all(reflection[key] == param[key] for key in param.keys()):
                positions.append(index)
        return positions

    def drop_zero(self):
        """Delete all reflections with h==h==l==0"""
        indices_to_delete = self.seek_reflection(**{'h': 0, 'k': 0, 'l': 0})
        self.data.drop(indices_to_delete, inplace=True)
        self.data.reset_index(drop=True, inplace=True)

    def place(self):
        """Assign reflections their positions anr calculate dist. from 000"""

        # CREATE EMPTY DATA HOLDERS AND IMPORT CONSTANTS
        _x, _y, _z, _r, _types = list(), list(), list(), list(), dict()
        _a, _b, _c = self.crystal.a_w, self.crystal.b_w, self.crystal.c_w
        _keys = ('x', 'y', 'z', 'r')
        for key in _keys:
            _types[key] = self.keys.get_property(key, 'dtype')

        # POSITION THE REFLECTIONS IN THE RECIPROCAL SPACE
        for index, row in self.data.iterrows():
            _v = row['h'] * _a + row['k'] * _b + row['l'] * _c
            _x.append(_v[0])
            _y.append(_v[1])
            _z.append(_v[2])
            _r.append(lin.norm(_v))

        # ADD KEYS AND APPEND DATA TO DATAFRAME
        self.keys.add(_keys)
        self.data['x'] = pd.Series(_x, index=self.data.index, dtype=_types['x'])
        self.data['y'] = pd.Series(_y, index=self.data.index, dtype=_types['y'])
        self.data['z'] = pd.Series(_z, index=self.data.index, dtype=_types['z'])
        self.data['r'] = pd.Series(_r, index=self.data.index, dtype=_types['r'])

    def transform(self, matrix):
        """Transform reflection pattern using given 3x3 or 4x4 matrix"""
        m_linear = np.array(((matrix[0][0], matrix[0][1], matrix[0][2]),
                             (matrix[1][0], matrix[1][1], matrix[1][2]),
                             (matrix[2][0], matrix[2][1], matrix[2][2])))
        try:
            m_translation = np.array((matrix[0][3], matrix[1][3], matrix[2][3]))
        except IndexError:
            m_translation = np.array((0.0, 0.0, 0.0))
        _x, _y, _z = [], [], []
        for index, row in self.data.iterrows():
            v = np.array((row['x'], row['y'], row['z']))
            v = np.dot(m_linear, v) + m_translation
            _x.append(v[0])
            _y.append(v[1])
            _z.append(v[2])
        self.data['x'] = pd.Series(_x, index=self.data.index)
        self.data['y'] = pd.Series(_y, index=self.data.index)
        self.data['z'] = pd.Series(_z, index=self.data.index)

    def read(self, hkl_path, hkl_format):
        """Read .hkl file as specified by path and fields
        format: either ordered dictionary with specified fields (minus = ignore)
        or type number"""

        # PREPARE OBJECTS RESPONSIBLE FOR PARSING INPUT
        format_string, column_labels = self.interpret_hkl_format(hkl_format)
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
        self.data = self.data_from_dict(hkl_content)

    def write(self, hkl_path, hkl_format):
        """Write .hkl file as specified by path and write_format"""

        # PREPARE OBJECTS RESPONSIBLE FOR WRITING OUTPUT
        first_line = ''
        format_string, column_labels = self.interpret_hkl_format(hkl_format)
        if format_string == 'XD':
            first_line = 'COMPOUND_ID          F^2  NDAT 6'
            format_string, column_labels = self.interpret_hkl_format(
                OrderedDict([('h', 5), ('k', 5), ('l', 5), ('b', 5),
                             ('I', 10), ('si', 10)])
            )
        column_sizes, column_formats = [], []
        for column in format_string.split():
            number, letter = int(column[:-1]), column[-1]
            if letter == 'x':
                column_sizes.append(0)
            elif letter == 's':
                column_sizes.append(number)
        for size in column_sizes:
            column_formats.append('{{:>{}.{}}}'.format(size, size-1))

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

        # WRITE SELF.DATA CONTENTS
        hkl_file = open(hkl_path, 'w')
        for index, row in self.data.iterrows():
            # FOR EACH DEMANDED KEY PRINT IT ACCORDING TO FORMATS
            for key, form in zip(column_labels, column_formats):
                hkl_file.write(form.format(
                    str(self.keys.get_property(key, 'type')(row[key]))))
            hkl_file.write('\n')

        # WRITE LAST DEFAULT 0 0 0 LINE AND CLOSE THE FILE
        for key, form in zip(column_labels, column_formats):
            hkl_file.write(form.format(
                str(self.keys.get_property(key, 'default'))))
        hkl_file.close()

    def dac(self, opening_angle=40):
        """Cut the reflections based on DAC angle (in angles)
        and [h, k, l] indexes of diamond-parallel crystal face"""

        # CALCULATE OPENING ANGLE "oa" AND WAVELENGTH "la"
        oa = np.radians([float(opening_angle)])[0]
        la = self.meta['wavelength']

        # CALCULATE NORMAL VECTOR "_n" IN RECIPROCAL LATTICE
        _r = np.array((1.0, 0.0, 0.0))      # vector parallel to beam
        _UB = self.crystal.orient_matrix    # import UB orientation matrix
        _h = np.dot(lin.inv(_UB), _r)       # calculate plane hkl indices
        _n = _h[0] * self.crystal.a_w + \
             _h[1] * self.crystal.b_w + \
             _h[2] * self.crystal.c_w       # calc. perp. vect. in recipr. space
        _n /= lin.norm(_n)                  # normalise the n vector

        # DELETE REFLECTIONS OUT OF ACCESSIBLE VOLUME
        indices_to_delete = []
        for index, row in self.data.iterrows():
            # CALCULATE RECIPROCAL POSITION VECTOR "_v"
            _v = np.array((row['x'], row['y'], row['z']))
            # CALCULATE PARALLELLED TO NORMAL VECTOR "_p"
            _p = _v - _n * np.dot(_v, _n)
            # EXCLUDE REFLECTIONS PARALLEL TO NORMAL
            if np.isclose([lin.norm(_p)], [0.]):
                indices_to_delete.append(index)
                continue
            # CALCULATE INCLINATION OF "_v" FROM CELL PLANE (PERP. TO NORMAL)
            _ph = np.arccos(
                min((np.dot(_p, _v) / (lin.norm(_p) * lin.norm(_v)), 1)))
            r_max = 2 * np.sin(max(oa-_ph, 0)) / la
            # IF ITS OUT OF ACCESSIBLE VOLUME ADD TO DELETION
            if lin.norm(_v) > r_max:
                indices_to_delete.append(index)
        self.data.drop(indices_to_delete, inplace=True)
        self.data.reset_index(drop=True, inplace=True)

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

    def trim(self, limit):
        """Cut the reflection further then the limit in A-1"""

        indices_to_delete = []
        for index, reflection in self.data.iterrows():
            v = np.array((reflection['x'], reflection['y'], reflection['z']))
            if lin.norm(v) > limit:
                indices_to_delete.append(index)
        self.data.drop(indices_to_delete, inplace=True)
        self.data.reset_index(drop=True, inplace=True)

    def reduce(self):
        """Average down redundant reflections, e.g. for drawing"""

        # SORT THE HKL DATAFRAME AND PREPARE NECESSARY OBJECTS
        self.data = self.data.sort_values(['h', 'k', 'l'])
        superreflections, redundants = dict(), list()
        self.keys.remove(self.keys.reduce_behaviour['discard'])
        for key in self.keys.all:
            superreflections[key] = list()

        # IMPORT REDUCTION BEHAVIOURS
        _reduction_kept = self.keys.reduce_behaviour['keep']
        _reduction_added = self.keys.reduce_behaviour['add']
        _reduction_averaged = self.keys.reduce_behaviour['average']

        # DEFINE HOW TO TREAT REDUNDANT REFLECTIONS
        def average_down_redundant_reflections():
            # DEAL WITH "KEEP"-TYPE KEYS OF REFLECTION
            for key in _reduction_kept:
                superreflections[key].append(redundants[0][key])
            # DEAL WITH "ADD"-TYPE KEYS OF REFLECTION
            for key in _reduction_added:
                _sum = sum([redundant[key] for redundant in redundants])
                superreflections[key].append(_sum)
            # DEAL WITH "AVERAGE"-TYPE KEYS OF REFLECTION
            for key in _reduction_averaged:
                _sum = sum([redundant[key] for redundant in redundants])
                _amount = len(redundants)
                superreflections[key].append(_sum / _amount)

        # ITERATE OVER REFLECTIONS
        for index, reflection in self.data.iterrows():
            # CHECK CURRENT REFLECTION
            try:
                _conditions = [reflection[i] == redundants[0][i] for i in 'hkl']
                reflections_are_redundant = all(_conditions)
            except IndexError:
                # IF ITS THE VERY FIRST REFLECTION, ADD IT TO REDUNDANT ONES
                redundants = [reflection]
                continue
            # IF IT IS, ADD IT TO THE LIST
            if reflections_are_redundant:
                redundants.append(reflection)
            # IF IT ISN'T:
            else:
                average_down_redundant_reflections()
                redundants = [reflection]

        # INCLUDE LAST REDUNDANTS AND CREATE NEW DATAFRAME
        average_down_redundant_reflections()
        self.data = self.data_from_dict(superreflections)

    def calculate_uncertainty(self, master_key):
        """For each reflection calculate u = master_key/sigma(master_key)"""
        uncertainties = []
        for index, row in self.data.iterrows():
            try:
                uncertainty = abs(row[master_key] / row['si'])
            except ZeroDivisionError:
                uncertainty = 0.0
            uncertainties.append(uncertainty)
        self.keys.add(['u'])
        typ = self.keys.get_property('u', 'dtype')
        self.data['u'] = pd.Series(uncertainties, dtype=typ)

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
        self.drop_zero()
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

    def generate_ball(self, radius=2.0):
        """Generate points in a sphere of given radius"""
        a, b, c = self.crystal.a_w, self.crystal.b_w, self.crystal.c_w
        r2 = radius**2

        # function checking whether the point is in sphere
        def _is_in_sphere(h, k, l):
            v = h * a + k * b + l * c
            return np.dot(v, v) <= r2
            # return lin.norm(h * a + k * b + l * c) <= radius # distance method

        # non-negative values generator
        def _positive_integers():
            n = -1
            while True:
                n += 1
                yield n

        # negative values generator
        def _negative_integers():
            n = 0
            while True:
                n -= 1
                yield n

        # prepare container for points
        hkl_content = dict()
        column_labels = ['h', 'k', 'l', 'I', 'si', 'b', 'm']
        for key in column_labels:
            hkl_content[key] = []
        self.keys.set(column_labels)

        def _add_hkl():
            hkl_content['h'].append(h)
            hkl_content['k'].append(k)
            hkl_content['l'].append(l)
            hkl_content['I'].append(1.0)
            hkl_content['si'].append(1.0)
            hkl_content['b'].append(1)
            hkl_content['m'].append(1)

        for h in _positive_integers():
            for k in _positive_integers():
                for l in _positive_integers():
                    if _is_in_sphere(h, k, l):
                        _add_hkl()
                    else:
                        break
                for l in _negative_integers():
                    if _is_in_sphere(h, k, l):
                        _add_hkl()
                    else:
                        break
                if not _is_in_sphere(h, k, 0):
                    break
            for k in _negative_integers():
                for l in _positive_integers():
                    if _is_in_sphere(h, k, l):
                        _add_hkl()
                    else:
                        break
                for l in _negative_integers():
                    if _is_in_sphere(h, k, l):
                        _add_hkl()
                    else:
                        break
                if not _is_in_sphere(h, k, 0):
                    break
            if not _is_in_sphere(h, 0, 0):
                break
        for h in _negative_integers():
            for k in _positive_integers():
                for l in _positive_integers():
                    if _is_in_sphere(h, k, l):
                        _add_hkl()
                    else:
                        break
                for l in _negative_integers():
                    if _is_in_sphere(h, k, l):
                        _add_hkl()
                    else:
                        break
                if not _is_in_sphere(h, k, 0):
                    break
            for k in _negative_integers():
                for l in _positive_integers():
                    if _is_in_sphere(h, k, l):
                        _add_hkl()
                    else:
                        break
                for l in _negative_integers():
                    if _is_in_sphere(h, k, l):
                        _add_hkl()
                    else:
                        break
                if not _is_in_sphere(h, k, 0):
                    break
            if not _is_in_sphere(h, 0, 0):
                break

        # produce new pandas frame with only considered points
        self.data = self.data_from_dict(hkl_content)

    def count_points_in_sphere(self, radius):
        """Calculate amount of points in a sphere of given radius"""

        self.generate_ball(radius)
        return len(self)

    def calculate_statistics(self, radius):
        # TODO forgot to consider radius :)
        pos_reflections = self.count_points_in_sphere(radius)
        copy = self
        copy.drop_zero()
        copy.data = copy.data.sort_values(['r'])
        all_reflections = 0
        obs_reflections = 0
        for index, row in copy.data.iterrows():
            if row['r'] <= radius:
                all_reflections += row['m']
        copy.reduce()
        ind_reflections = 0
        for index, row in copy.data.iterrows():
            if row['r'] <= radius:
                ind_reflections += 1
                try:
                    if row['I']/row['si'] > 2:
                        obs_reflections += row['m']
                except KeyError:
                    if row['F']/row['si'] > 2:
                        obs_reflections += row['m']
        print('RADIUS: ' + str(radius))
        print('all: ' + str(all_reflections))
        print('independent: ' + str(ind_reflections))
        print('observed: ' + str(obs_reflections))
        print('possible: ' + str(pos_reflections))
        print('completeness: ' + str(ind_reflections/pos_reflections))
        print('redundancy: ' + str(all_reflections/ind_reflections) + '\n')


if __name__ == '__main__':
    import copy

    p = HklFrame()
    q = HklFrame()
    p.read('/home/dtchon/git/kesshou/test_data/c1_p1_dt.hkl', 4)
    q.read('/home/dtchon/git/kesshou/test_data/c1_p1_dt_THoa45.hkl', 4)
    from cif import CifFrame, ustrip
    # c = CifFrame()
    # c.read('/home/dtchon/git/kesshou/test_data/exp_353.cif', datablock='exp_353')
    # p.crystal.import_from_frame(c)

    # exp_353 orientation
    #p.crystal.orient_matrix = np.array(((
    #    0.0117244000, 0.0419390000, 0.0270692000),
     #   (0.0571303000, 0.0081278000, 0.0181937000),
     #   (0.0264406000, 0.0361560000, 0.0272830000)))

    #c1_p1_dt (RFpirazCHO) orientation
    p.crystal.orient_matrix = np.array(((
                                            -0.0026344000, 0.0562933000, -0.0106563000),
       (-0.0075935000, -0.0365821000, -0.0163966000),
       (-0.0748129000, 0.0020128000, 0.0016300000)))
    q.crystal.orient_matrix = p.crystal.orient_matrix
    p.crystal.edit_cell(a=7.448, b=8.302, c=28.56, al=90, be=91.12, ga=90)
    q.crystal.edit_cell(a=7.448, b=8.302, c=28.56, al=90, be=91.12, ga=90)
    p.edit_wavelength('Ag')
    q.edit_wavelength('Ag')

    p.place()
    p.drop_zero()
    p.reduce()
    q.place()
    q.drop_zero()

    print('oa', 'experiment', 'theory')
    for angle in range(45, 0, -1):
        p.dac(opening_angle=angle)
        q.dac(opening_angle=angle)
        print(angle, len(p), len(q))


    #p.write('/home/dtchon/git/kesshou/test_data/output.hkl', 'XD')

    #p.crystal.edit_cell(a=11.071, b=12.664, c=16.666)
    #p.place()
    #for radius in np.arange(0.2, 3.0, 0.2):
    #    p.calculate_statistics(radius) #0.625A-1 = 0.8A

# TODO 3D call visualise and to pyqtplot
# TODO some function changes something globally (try to cut some