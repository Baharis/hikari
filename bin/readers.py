from collections import OrderedDict


def read_res(path):
    """Read data from specified ins/res file and return a relevant dictionary"""
    output = OrderedDict()

    # SPECIFY DEFAULT VALUES OF NECESSARY KEYS
    output['REM'] = ['These comments were found by resins:']
    output['ATOM'] = dict()
    output['PEAK'] = dict()
    output['TITL'] = list()
    output['SYMM'] = list()

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
        if line[-1:] == '=':
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
            output['REM'].append(line)
            lines_to_delete.append(index)
        elif all((key.upper() not in key_types.keys(), reading_peaks == True,
        key[0] == 'Q', key[1] in '0123456789')):
            output['PEAK'][key] = values
            lines_to_delete.append(index)
        elif all((key.upper() not in key_types.keys(), reading_atoms == True,
        key[0] in 'abcdefghijklmnopqrstuwxyzABCDEFGHIJKLMNOPQRESUWXYZ')):
            output['ATOM'][key] = values
            lines_to_delete.append(index)
        elif key.upper() not in key_types.keys():
            output['REM'].append('Line not interpreted: '+line)
            lines_to_delete.append(index)
    for index in reversed(lines_to_delete):
        del (lines[index])

    # FOR EACH LINE SPLIT LINE TO KEY AND VALUE AND APPEND SELF
    for line in lines:
        key = line.split(maxsplit=1)[0].upper()
        value = line[len(key):].strip()
        if key_types[key] == 'listing':
            output[key] = list(value.split())
        elif key_types[key] == 'string':
            output[key] = value
        elif key_types[key] == 'multiline':
            output[key].append(value)
        elif key == 'END':
            break

    # BRING REM, ATOM AND PEAK TO THE END AND RETURN DICTIONARY
    output.move_to_end('REM')
    output.move_to_end('ATOM')
    output.move_to_end('PEAK')
    return output


def read_cif(path, data_block='I'):
    """Read data from specified cif file's data block and return a dictionary"""
    output = OrderedDict()

    # READ THE FILE, FIND RELEVANT DATA BLOCK, DELETE REST
    lines = [line.strip() for line in open(path, 'r')]
    lines_to_delete = []
    this_content_is_interesting = False
    for index, line in enumerate(lines):
        if line[:5+len(data_block)] == 'data_' + data_block:
            this_content_is_interesting = True
        elif line[:5] == 'data_':
            this_content_is_interesting = False
        if not this_content_is_interesting:
            lines_to_delete.append(index)
    for index in reversed(lines_to_delete):
        del(lines[index])

    # DELETE BLOCK DEFINITION, ITERATE OVER FILE, TRY TO READ KEY AND VALUE
    del(lines[0])
    index = 0
    while index < len(lines):
        line = lines[index]
        try:
            key = line.split(maxsplit=1)[0]
            value = line[len(key):].strip()
        except IndexError:
            key, value = '', ''

        # IF DEALING WITH EXTERNAL COMMENT
        if '#' in lines[index][:1]:
            pass

        # IF DEALING WITH STANDARD ENTRY
        elif key != 'loop_' and len(value) > 0:
            output[key] = value

        # IF DEALING WITH MULTI-LINE ENTRY
        elif key != 'loop_' and len(key) > 1 and len(value) == 0:
            index += 2
            value = []
            while ';' not in lines[index][:1]:
                value.append(lines[index])
                index += 1
            output[key] = value

        # IF DEALING WITH LOOP
        elif key == 'loop_':
            index += 1
            subkeys, subvalues, subvalues_lines = [], [], []
            while lines[index][:1] == '_':
                subkeys.append(lines[index])
                index += 1
            while len(lines[index].split()) == len(subkeys):
                subvalues_lines.append(list(lines[index].split()))
                index += 1
            subvalues = list(map(list, zip(*subvalues_lines)))
            if not len(subvalues):
                for i in range(len(subkeys)):
                    subvalues.append([])
            output.update(OrderedDict(zip(subkeys, subvalues)))

        index += 1
    return output
