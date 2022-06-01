import re
from collections import OrderedDict
from enum import Enum
from hikari.utility import make_abspath


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


def common_prefix(loop_elements):
    prefix = ''
    for char_index in range(len(loop_elements[0])):
        chars = [loop_element[char_index] for loop_element in loop_elements]
        if len(set(chars)) == 1:
            prefix += chars[0]
        else:
            break
    return prefix.rstrip('_')


class CifFrame:
    """
    A master object which manages cif files. It utilises other `Cif*`
    classes to import and export crystallographic information.

    CifFrame stores all information under `data` attribute inside an ordered
    dictionary. Similarly to other `Frame`s, it is designed to work in-place,
    meaning a `CifFrame` should be first created, and then modified using
    methods such as :func:`read` or :func:`write`, but not chain assignments.

    The `CifFrame` always initiates empty and does not accept any arguments.
    """
    def __init__(self):
        self.data = OrderedDict()
        self.block = ''

    def read(self, path, block='I'):
        """
        Read the contents of .cif file as specified by `path` and data `block`
        and store them in the ordered dictionary in `self.data`.

        :param path: Absolute or relative path to the .cif file.
        :type path: str
        :param block: Name of the cif block to be read (without "data_").
        :type block: str
        """
        reader = CifIO(cif_file_path=path, cif_block_header=block)
        self.data = reader.read()
        self.block = block


class CifIO:
    """
    A helper class supporting CifFrame.
    Menages reading and writing cif files
    into and out of CifFrame's dataframe.
    Based on the IUCr File Syntax version 1.1 Working specification available
    [here](`https://www.iucr.org/resources/cif/spec/version1.1/cifsyntax`)
    """
    WHITESPACE_SUBSTITUTES = {' ': '█'}

    def __init__(self, cif_file_path, cif_block_header):
        self.file_path = make_abspath(cif_file_path)
        self.file_lines = []
        self.data_block_header = cif_block_header
        self.data = OrderedDict()

    class DataBuffer:
        """This class buffers data in temporary dict until flush() is called."""
        def __init__(self, target):
            self.names = []
            self.values = []
            self.target = target

        def parse(self, word):
            """Append the word to names or values based on its first char"""
            if word.startswith('_'):
                self.names.append(word)
            else:
                self.values.append(word)

        def append_to_multiline(self, string):
            """Add the word to values if they're empty, concatenate otherwise"""
            if self.values:
                self.values[-1] = self.values[-1] + '\n' + string
            else:
                self.values.append(string)

        def flush(self):
            """Update the target dict with names and values stored hitherto"""
            d = OrderedDict()
            lv = len(self.values)
            ln = len(self.names)
            if lv == ln:
                d.update({n: v for n, v in zip(self.names, self.values)})
            elif lv % ln == 0 and lv > 0:
                d.update({n: self.values[i::ln]
                          for i, n in enumerate(self.names)})
            else:
                raise IndexError(f'len(values) == {lv} must be a positive '
                                 f'multiple of len(names) == {ln}')
            self.target.update(d)
            self.names = []
            self.values = []

    class ReadingState(Enum):
        """This class stores current cif reading state (eg. inside loop etc.)"""
        default = 0
        loop = 1
        loop_header = 2
        multiline = 3

    @staticmethod
    def remove_outer_quotes(s):
        """
        Return the string without matching outer quotation marks: `'` or `"`
        :param s: the string which should be stripped
        :type s: str
        :return: string `s` without matching outer quotation marks
        :rtype: str
        """
        left = s[0]
        right = s[-1]
        return s.strip(left) if left == right and left in ('"', "'") else s

    def locate_block(self, block_header=''):
        """
        Determine the numbers of first lines of the seeked and next cif block
        :param block_header: name of the cif block which should be localised
        :type block_header: str
        :return: tuple with numbers of 1st-line-to-read and 1st-line-not-to-read
        :rtype: tuple
        """
        data_block_string = 'data_' + block_header
        data_block_start = 1
        data_block_end = None
        for line_number, line in enumerate(self.file_lines):
            if data_block_string in line:
                data_block_start = line_number + 1
                break
        for line_number, line in enumerate(self.file_lines[data_block_start:]):
            if line.startswith('data_'):
                data_block_end = line_number
                break
        return data_block_start, data_block_end

    def parse_lines(self, start, end):
        """
        Read the data from :attr:`~.CifIO.lines` numbered `start` to `end`,
        interpret it, and store it in :attr:`~.CifIO.data` dictionary.
        :param start: number of the first line which data should be read from
        :type start: int
        :param end: number of the first line which should not be read anymore
        :type end: int
        """
        buffer = self.DataBuffer(target=self.data)
        state = self.ReadingState.default
        for line in self.file_lines[start:end]:
            if state is self.ReadingState.loop_header:
                state = self.ReadingState.loop
            if line.startswith('#'):
                continue
            if state is self.ReadingState.multiline and line.startswith(';'):
                buffer.flush()
                state = self.ReadingState.default
                line = line[1:]
            elif state is self.ReadingState.multiline:
                buffer.append_to_multiline(line)
                continue
            elif line.startswith(';'):
                state = self.ReadingState.multiline
                line = line[1:]
            elif line.startswith('loop_'):
                state = self.ReadingState.loop_header
                line = line[5:]
            words = self.split_line(line)
            if not words:
                if state in {self.ReadingState.loop, self.ReadingState.default}:
                    buffer.flush()
                    state = self.ReadingState.default
                elif state is self.ReadingState.multiline:
                    buffer.append_to_multiline('')
                continue
            if words[0].startswith('_') and state is self.ReadingState.default:
                buffer.flush()
            for word in words:
                buffer.parse(word)
            if not words and state is self.ReadingState.loop:
                pass
        buffer.flush()

    def read(self):
        """
        Read the contents of cif currently pointed by :attr:`~.CifIO.file_path`
        and block :attr:`~.CifIO.data_block_header` and return them to a dict.
        :return: A dictionary containing information read from .cif file.
        :rtype: dict
        """
        with open(self.file_path, 'r') as cif_file:
            self.file_lines = cif_file.read().splitlines()
        block_start, block_end = self.locate_block(self.data_block_header)
        self.parse_lines(start=block_start, end=block_end)
        return self.data

    def split_line(self, line):
        """
        Split line into words, keeping words inside quotation marks together.
        :param line: line to be split based on whitespace into words
        :type line: str
        :return: list of words obtained from splitting
        :rtype: list
        """
        substituted_line = self.substitute_whitespace_in_quotes(line)
        words = []
        for word in substituted_line.strip().split():
            word = self.substitute_whitespace_in_quotes(word, reverse=True)
            words.append(self.remove_outer_quotes(word))
        return words

    def substitute_whitespace_in_quotes(self, string, reverse=False):
        """
        Substitute whitespace between matching quotation marks with substitutes.
        :param string: text in which whitespace will be substituted
        :type string: str
        :param reverse: if True, change the substitutes back to whitespace
        :type reverse: bool
        :return: string where whitespace/substitutes inside quotes were changed
        :rtype: str
        """
        # see: https://stackoverflow.com/q/46967465/, https://regex101.com/
        split_string = []
        matching_quotes_regex = r"""(["'])((?:\\\1|(?:(?!\1)).)*)(\1)"""
        for i, m in enumerate(re.split(matching_quotes_regex, string)):
            if i % 4 == 2:
                for ws, sub in self.WHITESPACE_SUBSTITUTES.items():
                    m = m.replace(sub, ws) if reverse else m.replace(ws, sub)
            split_string.append(m)
        return ''.join(split_string)


if __name__ == '__main__':
    cifio = CifIO(cif_file_path='~/x/HiPHAR/anders_script/rfpirazB_100K_SXD.cif',
                  cif_block_header='rfpirazB_100K_SXD')
    cifio.read()
    for k, v in cifio.data.items():
        print(f'{k} :: {repr(v)}')

# TODO Try using pyCIFrw package to read and write cif information.
