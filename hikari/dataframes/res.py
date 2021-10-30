from collections import OrderedDict


class ResFrame:
    def __init__(self, file_path=None):
        self.data = OrderedDict()
        self.meta = dict()
        if file_path is not None:
            self.read(path=file_path)

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
        lines = [line.strip() for line in open(path, 'r') if line.strip()]
        for index, line in enumerate(lines):
            while line[-1:] == '=':
                lines[index] = line[:-1] + lines[index+1]
                del(lines[index+1])

        # READ THE FILE AND GATHER ALL COMMENTS, ATOMS AND PEAKS
        reading_atoms, reading_peaks, lines_to_delete = False, False, []
        for index, line in enumerate(lines):
            key = line.split(maxsplit=1)[0]
            values = list(line[len(key):].strip().split())
            if key.upper() == 'FVAR':
                reading_atoms = True
            elif key.upper() == 'HKLF':
                reading_atoms = False
            elif key.upper() == 'END':
                reading_peaks = True
            elif key.upper() == 'REM':
                self.data['REM'].append(line)
                lines_to_delete.append(index)
            elif all((key.upper() not in key_types.keys(),
                      reading_peaks is True,
                      key[0] == 'Q', key[1].isdigit())):
                self.data['PEAK'][key] = values
                lines_to_delete.append(index)
            elif all((key.upper() not in key_types.keys(),
                      reading_atoms is True, key[0].isalpha())):
                self.data['ATOM'][key] = values
                lines_to_delete.append(index)
            elif key.upper() not in key_types.keys():
                self.data['REM'].append('Line not interpreted: '+line)
                lines_to_delete.append(index)
        for index in reversed(lines_to_delete):
            del (lines[index])

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

        # UPDATE META INFORMATION
        self.meta['comment'] = self.data['REM']
