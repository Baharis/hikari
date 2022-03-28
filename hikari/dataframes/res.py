from enum import Enum
from collections import OrderedDict, defaultdict
import numpy as np
from pandas import DataFrame
from hikari.dataframes import BaseFrame
from hikari.resources import Xray_atomic_form_factors

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

    def atomic_form_factor(self, atom, hkl, u):
        s = Xray_atomic_form_factors.loc[atom]
        sintl2 = np.dot(self.A_r @ hkl, self.A_r @ hkl) / 4  # 5.641087
        q = np.exp(-2*np.pi**2 * (hkl.T @ self.G_r @ self.A_d.T @ u @
                                  self.A_d @ self.G_r @ hkl))
        # TODO for some bizarre reason multiplying q by a**2 magically works???
        f = s['a1'] * np.exp(-s['b1'] * sintl2) + \
            s['a2'] * np.exp(-s['b2'] * sintl2) + \
            s['a3'] * np.exp(-s['b3'] * sintl2) + \
            s['a4'] * np.exp(-s['b4'] * sintl2) + s['c']
        return q * f

    def form_factor(self, hkl, space_group):
        f = 0.0
        for k, v in self.data['ATOM'].items():
            atom = k.title()[:2]
            xyz = np.fromiter(map(float, v[1:4]), dtype=np.float)
            percent = float(v[4]) % 10
            u11, u22, u33, u12, u13, u23 = map(float, v[5:11])
            u = np.array([[u11, u12, u13], [u12, u22, u23], [u13, u23, u33]])
            for o in space_group.operations:
                new_xyz = (o.tf @ np.array(xyz)).T + o.tl
                new_u = o.tf @ u @ (o**-1).tf
                f_atom = percent * self.atomic_form_factor(atom, hkl, new_u)
                f += f_atom * np.exp(2 * np.pi * 1j * np.dot(new_xyz, hkl))
        return f
    # TODO imprecise, especially when sintl is large - see far NaCl reflections

    def read(self, path):
        """Read data from specified ins/res file and return an OrderedDict"""

        # SPECIFY DEFAULT VALUES OF NECESSARY KEYS
        self.data['REM'] = ['These comments were found by resins:']
        self.data['ATOM'] = dict()
        self.data['PEAK'] = dict()
        self.data['TITL'] = list()
        self.data['SYMM'] = list()

        # SPECIFY SUPPORTED KEYS AND THEIR TYPES
        key_types = OrderedDict([
            ('TITL',	'multiline'),
            ('CELL',	'listing'),
            ('ZERR',	'listing'),
            ('LATT',	'listing'),
            ('SYMM',	'multiline'),
            ('SFAC',	'listing'),
            ('UNIT',	'listing'),
            ('L.S.',	'listing'),
            ('PLAN',	'listing'),
            ('SIZE',    'listing'),
            ('MORE',	'listing'),
            ('BOND',	'listing'),
            ('CONF',	'listing'),
            ('FMAP',	'listing'),
            ('ACTA',	'listing'),
            ('WGHT',	'listing'),
            ('FVAR',	'listing'),
            ('HKLF',	'listing'),
            ('REM',		'special'),
            ('END',		'special')
        ])

        # READ THE FILE AND JOIN LINES SEPARATED BY '=' SIGN
        with open(path, 'r') as res_file:
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
