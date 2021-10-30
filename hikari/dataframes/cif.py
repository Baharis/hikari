from collections import OrderedDict


def ustrip(ufloat):
    """Strip uncertainty out of float: 1.23(4) --> 1.23"""
    string = str(ufloat)
    new = ''
    for letter in string:
        if letter is not '(':
            new += letter
        else:
            break
    return new


class CifFrame:
    """

    """
    def __init__(self, file_path=None, file_data_block='I'):
        self.data = OrderedDict()
        self.meta = dict()
        if file_path is not None:
            self.read(path=file_path, datablock=file_data_block)

    def read(self, path, datablock='I'):
        """Read data from specified ins/res file and return an OrderedDict"""

        # SPECIFY SOME META
        self.meta['name'] = datablock
        self.meta['comment'] = str()

        # READ THE FILE, FIND RELEVANT DATA BLOCK, DELETE REST
        lines = [line.strip() for line in open(path, 'r')]
        lines_to_delete = []
        correct_datablock = False
        for index, line in enumerate(lines):
            if line[:5 + len(datablock)] == 'data_' + datablock:
                correct_datablock = True
            elif line[:5] == 'data_':
                correct_datablock = False
            if not correct_datablock:
                lines_to_delete.append(index)
        for index in reversed(lines_to_delete):
            del(lines[index])

        # DELETE BLOCK DEFINITION, ITERATE OVER FILE, TRY TO READ KEY AND VALUE
        try:
            del(lines[0])
        except IndexError:
            pass
        index = 0
        while index in range(len(lines)):
            line = lines[index]
            try:
                key = line.split(maxsplit=1)[0]
                value = line[len(key):].strip()
            except IndexError:
                key, value = '', ''

            # IF IT'S ALREADY SHELX FILE, STOP LOOP
            if key == '_shelx_hkl_file':
                break

            # IF DEALING WITH STANDARD ENTRY
            if key != 'loop_' and len(value) > 0:
                self.data[key] = ustrip(value)

            # IF DEALING WITH MULTI-LINE ENTRY
            elif key != 'loop_' and len(key) > 1 and len(value) == 0:
                index += 1
                value = []
                while ';' not in lines[index][:1]:
                    value.append(lines[index])
                    index += 1
                self.data[key] = value

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
                    try:
                        _ = lines[index]
                    except IndexError:
                        break
                subvalues = list(map(list, zip(*subvalues_lines)))
                if not len(subvalues):
                    for i in range(len(subkeys)):
                        subvalues.append([])
                self.data.update(OrderedDict(zip(subkeys, subvalues)))
            index += 1


if __name__ == '__main__':
    cif = CifFrame(file_path='//test_data/exp_353.cif',
                   file_data_block='exp_353')
    print(len(cif.data.items()))
    for key, value in cif.data.items():
        print(key, '::', value)

# TODO Try using pyCIFrw package to read and write cif information.
