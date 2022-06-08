from enum import Enum
from collections import OrderedDict, defaultdict
import numpy as np
from hikari.dataframes import BaseFrame
from hikari.resources import Xray_atomic_form_factors
from hikari.utility import split_atom_label, make_abspath

res_instructions = defaultdict(set)


class ResInstructionType(Enum):
    UNIQUE_INT = 1
    UNIQUE_LIST = 2
    UNIQUE_MULTILINE = 3
    REPEATING_SINGLET = -1


class ResFrame(BaseFrame):
    def __init__(self):
        super().__init__()
        self.data = OrderedDict()

    def atomic_form_factor(self, atom, hkl):
        """
        Calculate X-ray atomic form factors for a single atom and a hkl array

        :param atom: Atom/ion name/identifier interpreted by form factor table
        :type atom: str
        :param hkl: A 2D array listing all hkls to consider
        :type hkl: np.array
        :return: A 1D array listing atomic form factors for desired hkls
        :rtype: np.array
        """
        r_star = self.A_r @ hkl.T
        s = Xray_atomic_form_factors.loc[atom]
        sintl2 = (r_star * r_star).sum(axis=0) / 4
        return s['a1'] * np.exp(-s['b1'] * sintl2) + \
               s['a2'] * np.exp(-s['b2'] * sintl2) + \
               s['a3'] * np.exp(-s['b3'] * sintl2) + \
               s['a4'] * np.exp(-s['b4'] * sintl2) + s['c']

    def temperature_factor(self, hkl, u):
        """
        Calculate temperature factor for single u matrix and a hkl array

        :param hkl: A 2D array listing all hkls to consider
        :type hkl: np.array
        :param u: A classical anisotropic displacement parameters matrix
        :type u: np.array
        :return: A 1D array listing temperature factors for desired hkls
        :rtype: np.array
        """
        n_matrix = np.diag([self.a_r, self.b_r, self.c_r])
        hkl_n_u_n_hkl = (hkl @ n_matrix * hkl @ n_matrix @ u.T).sum(axis=-1)
        return np.exp(-2 * np.pi ** 2 * hkl_n_u_n_hkl)

    def form_factor(self, hkl, space_group):
        """
        Calculate form factors based on current structure, hkls, and space group

        :param hkl: A 2D array listing all hkls to consider
        :type hkl: np.array
        :param space_group: Space group describing the internal crystal symmetry
        :type space_group: hikari.symmetry.Group
        :return: A 1D array listing total form factors for desired hkls
        :rtype: np.array
        """
        f = np.zeros_like(hkl.sum(axis=-1), dtype='complex128')
        for k, v in self.data['ATOM'].items():
            atom = split_atom_label(k)[0].title()
            xyz = np.fromiter(map(float, v[1:4]), dtype=np.float)
            occupation = float(v[4]) % 10
            uij_strings = v[5:] if len(v[5:]) == 6 else [v[5]] * 3 + ['0'] * 3
            u11, u22, u33, u23, u13, u12 = map(float, uij_strings)
            u = np.array([[u11, u12, u13], [u12, u22, u23], [u13, u23, u33]])
            for o in space_group.operations:
                xyz_t = (o.tf @ np.array(xyz)).T + o.tl
                u_t = o.tf @ u @ (o**-1).tf
                fa = occupation * self.atomic_form_factor(atom, hkl) * \
                    self.temperature_factor(hkl, u_t)
                f += fa * np.exp(2j * np.pi * xyz_t @ hkl.T)
        return f
    # TODO the factors seem to be evaluated incorrectly, check in the future

    def read(self, path):
        """
        Read data from specified ins/res file and return an OrderedDict

        :param path: Relative or absolute path to the res file to be read
        :type path: str
        :return: None
        :rtype: None
        """

        # SPECIFY DEFAULT VALUES OF NECESSARY KEYS
        self.data['REM'] = ['These comments were found by resins:']
        self.data['ATOM'] = dict()
        self.data['PEAK'] = dict()
        self.data['TITL'] = list()
        self.data['SYMM'] = list()

        # SPECIFY SUPPORTED KEYS AND THEIR TYPES
        key_types = OrderedDict([
            ('TITL', 'multiline'),
            ('CELL', 'listing'),
            ('ZERR', 'listing'),
            ('LATT', 'listing'),
            ('SYMM', 'multiline'),
            ('SFAC', 'listing'),
            ('UNIT', 'listing'),
            ('L.S.', 'listing'),
            ('PLAN', 'listing'),
            ('SIZE', 'listing'),
            ('MORE', 'listing'),
            ('BOND', 'listing'),
            ('CONF', 'listing'),
            ('FMAP', 'listing'),
            ('ACTA', 'listing'),
            ('WGHT', 'listing'),
            ('FVAR', 'listing'),
            ('HKLF', 'listing'),
            ('REM',	'special'),
            ('END',	'special')
        ])

        # READ THE FILE AND JOIN LINES SEPARATED BY '=' SIGN
        with open(make_abspath(path), 'r') as res_file:
            lines = res_file.read().replace('=\n', '').split('\n')
            lines = [line.strip() for line in lines if line.strip()]

        class ResIoStage(Enum):
            TITLE = 0
            PREAMBLE = 1
            ATOMS = 2
            APPENDIX = 3
            END = 4

        # READ THE FILE AND GATHER ALL COMMENTS, ATOMS AND PEAKS
        reading_stage = ResIoStage.PREAMBLE
        lines_to_delete = []
        for index, line in enumerate(lines):
            key = line.split(maxsplit=1)[0]
            values = list(line[len(key):].strip().split())
            if key.upper() == 'FVAR':
                reading_stage = ResIoStage.ATOMS
            elif key.upper() == 'HKLF':
                reading_stage = ResIoStage.APPENDIX
            elif key.upper() == 'END':
                reading_stage = ResIoStage.END
            elif key.upper() == 'REM':
                self.data['REM'].append(line)
                lines_to_delete.append(index)
            elif all((key.upper() not in key_types.keys(),
                      reading_stage is ResIoStage.END,
                      key[0] == 'Q', key[1].isdigit())):
                self.data['PEAK'][key] = values
                lines_to_delete.append(index)
            elif all((key.upper() not in key_types.keys(),
                      reading_stage is ResIoStage.ATOMS,
                      key[0].isalpha())):
                self.data['ATOM'][key] = values
                lines_to_delete.append(index)
            elif key.upper() not in key_types.keys():
                self.data['REM'].append('Line not interpreted: '+line)
                lines_to_delete.append(index)
        for index in reversed(lines_to_delete):
            del lines[index]

        # FOR EACH LINE SPLIT LINE TO KEY AND VALUE AND APPEND SELF
        for line in lines:
            key = line.split(maxsplit=1)[0].upper()
            value = line[len(key):].strip()
            if key_types[key] == 'listing':
                self.data[key] = list(value.split())
            elif key_types[key] == 'string':
                self.data[key] = value
            elif key_types[key] == 'multiline':
                self.data[key].append(value)
            elif key == 'END':
                break

        # BRING REM, ATOM AND PEAK TO THE END AND RETURN DICTIONARY
        self.data.move_to_end('REM')
        self.data.move_to_end('ATOM')
        self.data.move_to_end('PEAK')
