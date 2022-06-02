import re
from collections import OrderedDict
from enum import Enum
from functools import lru_cache
from itertools import chain, zip_longest

from hikari.utility import make_abspath


class CifBlock(OrderedDict):
    """
    CifBlock object handles all data inside an individual block of Cif file.
    It is a subclass of an `OrderedDict` and, as such, features a lot
    of similarities with python dictionary while preserving item order.
    Individual Cif items can be accessed or assigned using a dict-like syntax.
    """
    def __init__(self, *args):
        super().__init__(*args)

    def get_as_type(self, key, typ):
        """
        Get `self[key]` and convert it (element-wise for lists) to `typ`
        :param key: key associated with accessed element
        :type key: str
        :param typ: type/method applied to value or every element of value
        :return: converted value of `self[key]` or `default`
        :rtype: Union[list, str]
        """
        value = self[key]
        if value and typ:
            if isinstance(value, str):
                value = typ(value)
            elif isinstance(value, list):
                value = map(typ, value)
            else:
                raise TypeError(f'Unknown value type of {value}: {type(value)}')
        return value


class CifFrame(OrderedDict):
    """
    A master object which manages cif files. It utilises other `Cif*` classes
    to manage multiple :class:`CifBlock`s with crystallographic information.
    It is a subclass of an `OrderedDict` and, as such, features a lot
    of similarities with python dictionary while preserving item order.
    Individual Cif blocks and items within them can be accessed or assigned
    using a single- or nested- dict-like syntax.

    Similarly to other `Frame`s, `CifFrame` is designed to work in-place,
    meaning it should be first created, and only then accessed using
    methods such as :func:`read` or :func:`write`, but not chain assignments.

    Unlike OrderedDict, CifBlock always initiates empty and does not accept
    any parameters at creation.
    """

    def read(self, path):
        """
        Read the contents of .cif file specified by the `path` parameter.
        Store each found block as a {block_name: CifBlock} pair.

        :param path: Absolute or relative path to the .cif file.
        :type path: str
        """
        reader = CifReader(cif_file_path=path)
        self.update(reader.read())


class CifIO:
    """
    A base class for `CifRead` and `CifWrite`. This class and its inheritors
    base on the IUCr File Syntax version 1.1 Working specification available
    [here](`https://www.iucr.org/resources/cif/spec/version1.1/cifsyntax`)
    """
    MATCHING_QUOTES_REGEX = re.compile(r"""(["'])((?:\\\1|(?!\1).)*)(\1)""")
    WHITESPACE_SUBSTITUTES = {' ': '█', '\t': '▄'}

    def __init__(self, cif_file_path):
        self.file_path = make_abspath(cif_file_path)
        self.file_lines = []
        self.data = OrderedDict()


class CifReader(CifIO):
    """A helper class managing reading cif files into a `CifFrame`."""

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
            elif ln == 0:
                raise IndexError(f'Orphan values found while flushing: '
                                 f'{self.values}')
            elif lv % ln == 0 and lv > 0:
                d.update({n: self.values[i::ln]
                          for i, n in enumerate(self.names)})
            else:
                raise IndexError(f'len(values) == {lv} must be a positive '
                                 f'multiple of len(names) == {ln}')
            self.target.update(d)
            self.__init__(target=self.target)

    @property
    def blocks(self):
        """A dictionary of all blocks names and their positions in cif file."""
        return self._blocks(lines=tuple(self.file_lines))

    @lru_cache(maxsize=1)
    def _blocks(self, lines):
        return OrderedDict({l[5:]: i for i, l in enumerate(lines)
                            if l.startswith('data_')})

    class ReadingState(Enum):
        """This class stores current cif reading state (eg. inside loop etc.)"""
        default = 0
        loop = 1
        loop_header = 2
        multiline = 3

    def parse_lines(self, start, end):
        """
        Read the data from :attr:`~.CifIO.lines` numbered `start` to `end`,
        interpret it, and return it as an instance of an `OrderedDict`.
        :param start: number of the first line which data should be read from
        :type start: int
        :param end: number of the first line which should not be read anymore
        :type end: int
        :return: ordered dictionary with name: value pairs for all parsed lines
        :rtype: OrderedDict
        """
        parsed_data = OrderedDict()
        buffer = self.DataBuffer(target=parsed_data)
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
                buffer.flush()
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
        return parsed_data

    def split_line(self, line):
        """
        Split line into words, keeping words inside quotation marks together.
        :param line: line to be split based on whitespace into words
        :type line: str
        :return: list of words obtained from splitting
        :rtype: list
        """
        substituted_line = self.substitute_quoted_whitespace(line)
        return [self.revert_whitespace(word) for word
                in substituted_line.strip().split()]

    def read(self):
        """
        Read the contents of cif currently pointed by :attr:`~.CifIO.file_path`
        and block :attr:`~.CifIO.data_block_header` and return them to a dict.
        :return: A dictionary containing information read from .cif file.
        :rtype: dict
        """
        with open(self.file_path, 'r') as cif_file:
            self.file_lines = cif_file.read().splitlines()
        block_names = self.blocks.keys()
        block_starts = [v + 1 for v in self.blocks.values()]
        block_ends = list(block_starts)[1:] + [None]
        for n, s, e in zip(block_names, block_starts, block_ends):
            self.data[n] = CifBlock(self.parse_lines(s, e))
        return self.data

    def revert_whitespace(self, string):
        """
        Change the substitute characters in supplied `string` back to whitespace
        :param string: text in which whitespace will be reverted
        :type string: str
        :return: string where substitutes were changed back to whitespace
        :rtype: str
        """
        new_string = string
        for ws, sub in self.WHITESPACE_SUBSTITUTES.items():
            new_string = new_string.replace(sub, ws)
        return new_string

    def substitute_quoted_whitespace(self, string):
        """
        Substitute whitespace between matching quotation marks with substitutes
        and remove the outer quotation marks
        :param string: text in which whitespace will be substituted
        :type string: str
        :return: string where whitespace inside quotes were substituted
        :rtype: str
        """
        # see: https://stackoverflow.com/q/46967465/, https://regex101.com/
        split_by_quotes = self.MATCHING_QUOTES_REGEX.split(string)[::2]
        quoted = split_by_quotes[::2]
        unquoted = split_by_quotes[1::2]
        for ws, sub in self.WHITESPACE_SUBSTITUTES.items():
            unquoted = [w.replace(ws, sub) for w in unquoted]
        split_by_quotes = chain.from_iterable(zip_longest(unquoted, quoted))
        return ''.join([w for w in split_by_quotes if w is not None])


if __name__ == '__main__':
    c = CifFrame()
    c.read(path='~/x/HiPHAR/anders_script/rfpirazB_100K_SXD.cif')
    for k, v in c['rfpirazB_100K_SXD'].items():
        print(f'{k} :: {repr(v)}')
