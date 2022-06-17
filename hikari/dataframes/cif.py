import abc
import re
import pathlib
import tempfile

from collections import OrderedDict
from enum import Enum
from functools import lru_cache
from typing import List, Union, Dict, TextIO

from hikari.resources import cif_core_dict
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

    def get_as_type(self, key, typ, default=None):
        """
        Get value of `self[key]` converted to `typ`. If value is a list,
        convert its contents element-wise.

        :param key: key associated with accessed element
        :type key: str
        :param typ: type/function applied to a value or its every element
        :type typ: Callable
        :param default: if given, return it on KeyError
        :type default: Any
        :return: converted value of `self[key]` or `default`
        :rtype: Union[List, str]
        """
        try:
            value = self[key]
        except KeyError as e:
            if default is not None:
                value = default
            else:
                raise e
        else:
            if value and typ:
                if isinstance(value, str):
                    value = typ(value)
                elif isinstance(value, list):
                    value = list(map(typ, value))
                else:
                    raise TypeError(f'Unknown value type'
                                    f'of {value}: {type(value)}')
        return value

    def read(self, path: str, block: str) -> None:
        """
        Read the contents of .cif file specified by the `path` parameter, but
        access and store only the `block` data block in self.

        :param path: Absolute or relative path to the .cif file.
        :type path: str
        :param block: Name of the cif data block to be accessed
        :type block: str
        """
        reader = CifReader(cif_file_path=path)
        self.update(reader.read()[block])

    def write(self, path: str) -> None:
        """
        Write the contents of `CifBlock` to the .cif file specified
        by the `path` parameter, using 'hikari' as block name.

        :param path: Absolute or relative path to the .cif file.

        """
        writer = CifWriter(cif_file_path=path)
        writer.write(cif_frame=CifFrame({'hikari': self}))


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

    def read(self, path: str) -> None:
        """
        Read the contents of .cif file specified by the `path` parameter.
        Store each found block as a {block_name: CifBlock} pair.

        :param path: Absolute or relative path to the .cif file.
        """
        reader = CifReader(cif_file_path=path)
        self.update(reader.read())

    def write(self, path: str) -> None:
        """
        Write the contents of `CifFrame` to the .cif file specified
        by the `path` parameter.

        :param path: Absolute or relative path to the .cif file.
        """
        writer = CifWriter(cif_file_path=path)
        writer.write(cif_frame=self)


class CifValidator(OrderedDict):
    """
    This object reads an appropriate cif core dictionary and uses it in order to
    format or validate all entries passing through it.

    The `CifValidator` contains all keys from core cif dictionary. In order
    to access individual values, use `.get()` instead of bracket notation.
    """

    def __init__(self):
        super().__init__()
        with tempfile.TemporaryDirectory() as temp_dir:
            temp_dic_path = str(pathlib.Path(temp_dir) / 'cif_core.dic')
            with open(temp_dic_path, 'w+') as f:
                f.write(cif_core_dict)
            reader = CifReader(cif_file_path=temp_dic_path, validate=False)
            self.update(reader.read())

    def __contains__(self, item):
        try:
            _ = self.get(item)
        except KeyError:
            return False
        else:
            return True

    def get(self, key, default=None):
        key, _key = (key[1:], key) if key.startswith('_') else (key, '_' + key)
        value = OrderedDict()
        try:
            value = self[key]
        except KeyError as e:
            for self_key, self_value in self.items():
                if key.startswith(self_key):
                    if _key in self_value.get('_name', []):
                        value = self_value
            if not value:
                value = default
        return value

    def get__category(self, key, default=None):
        value = self.get(key)
        if value is not None:
            _category = value.get('_category', default)
        else:
            _category = default
        return _category

    def get__list(self, key, default=None):
        value = self.get(key)
        if value is not None:
            got = value.get('_list')
            _list = True if got == 'yes' else False if got == 'no' else default
        else:
            _list = default
        return _list


class CifIOBuffer(abc.ABC):
    def __init__(self, target):
        self.names = []
        self.values = []

    @abc.abstractmethod
    def add(self, data):
        pass

    @abc.abstractmethod
    def flush(self):
        pass


class CifIO(abc.ABC):
    """
    A base class for `CifRead` and `CifWrite`. This class and its inheritors
    base on the IUCr File Syntax version 1.1 Working specification available
    [here](`https://www.iucr.org/resources/cif/spec/version1.1/cifsyntax`)
    """
    COMMENT_REGEX = \
        re.compile(r"(?<=\s)(#.*)(?=$)|(?<=^)(#.*)(?=$)")
    MATCHING_QUOTES_REGEX = \
        re.compile(r"(?<!\b)([\"'])((?:\\\1|(?!\1\s).)*.)(\1)(?!\b)")
    MATCHING_OUTER_QUOTES_REGEX = \
        re.compile(r"(?<=^)([\"'])((?:\\\1|(?!\1\s).)*.)(\1)(?=$)")
    WHITESPACE_SUBSTITUTES = {' ': '█', '\t': '▄'}

    def __init__(self, cif_file_path, validate=True):
        self.file_path = make_abspath(cif_file_path)
        self.file_lines = []
        self.validate = validate


class CifReaderBuffer(CifIOBuffer):
    """Buffer for reading data from cif file into `CifReader`"""

    def __init__(self, target):
        super().__init__(target=target)
        self.target: OrderedDict = target
        self.multilines = []

    def add(self, word):
        """Append the word to names or values based on its first char"""
        if word.startswith('_'):
            if self.values:
                self.flush()
            self.names.append(word)
        else:
            self.values.append(CifReader.release_quote(word))

    def flush(self):
        """Update the target dict with names and values stored hitherto"""
        d = OrderedDict()
        lv = len(self.values)
        ln = len(self.names)
        if lv == ln == 0:
            pass
        elif ln == 0:
            raise IndexError(f'Orphan values found while '
                             f'flushing buffer: {self.values}')
        elif lv % ln == 0:
            d.update({n: self.values[i::ln] for i, n in
                      enumerate(self.names)})
        else:
            raise IndexError(
                f'len(values) == {lv} % len(names) == {ln} mus'
                f't be zero: {self.values} % {self.names}')
        self.target.update(d)
        self.__init__(target=self.target)

    def initiate_multiline(self):
        self.multilines = []

    def append_to_multiline(self, string):
        """Add the word to values if they're empty, concatenate otherwise"""
        self.multilines.append(string)

    def terminate_multiline(self):
        self.values.append('\n'.join(self.multilines))


class CifReader(CifIO):
    """A helper class managing reading cif files into `CifFrame`s."""

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

    def format_dictionary(self, parsed_dict_: Dict[str, List[str]]) \
            -> Dict[str, Union[str, List[str]]]:
        """
        Reformat a dictionary of parsed data so that the format of every name
        and value agrees with the cif core dictionary stored in `CifValidator`.

        :param parsed_dict_: Dictionary with data pairs
        :return: Data dictionary with correctly formatted data names and values
        """

        def item_value_should_be_a_list(k_, v_):
            is_listable = cif_core_validator.get__list(k_) \
                if self.validate else False
            is_long = len(v_) > 1
            is_a_validator_name_field = not self.validate and k_ == '_name'
            return is_listable or is_long or is_a_validator_name_field

        new_dict = OrderedDict()
        for k, v in parsed_dict_.items():
            if item_value_should_be_a_list(k, v):
                new_dict[k] = v
            else:
                new_dict[k] = v[0]
        return new_dict

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
        buffer = CifReaderBuffer(target=parsed_data)
        state = self.ReadingState.default
        for line in self.file_lines[start:end]:
            line = self.protect_quotes(line)
            line = self.strip_comments(line)
            if state is self.ReadingState.loop_header:
                state = self.ReadingState.loop
            if line.startswith(';') and state != self.ReadingState.multiline:
                buffer.initiate_multiline()
                state = self.ReadingState.multiline
                line = line[1:]
            elif line.startswith(';') and state is self.ReadingState.multiline:
                buffer.terminate_multiline()
                state = self.ReadingState.default
                continue
            if state is self.ReadingState.multiline:
                buffer.append_to_multiline(line)
                continue
            elif line.lstrip().startswith('loop_'):
                buffer.flush()
                state = self.ReadingState.loop_header
                line = line.lstrip()[5:]
            words = line.strip().split()
            if not words and self.ReadingState.multiline:
                buffer.append_to_multiline(line)
                continue
            if words[0].startswith('_') and state is self.ReadingState.default:
                buffer.flush()
            for word in words:
                buffer.add(word)
            if not words and state is self.ReadingState.loop:
                pass
        buffer.flush()
        formatted_data = self.format_dictionary(parsed_data)
        return formatted_data

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
        block_starts = self.blocks.values()
        block_ends = list(block_starts)[1:] + [None]
        read_data = OrderedDict()
        for n, s, e in zip(block_names, block_starts, block_ends):
            read_data[n] = CifBlock(self.parse_lines(s + 1, e))
        return read_data

    @classmethod
    def protect_quotes(cls, string):
        """
        Substitute whitespace between matching quotation marks with substitutes
        and remove the outer quotation marks

        :param string: text in which whitespace will be substituted
        :type string: str
        :return: string where whitespace inside quotes were substituted
        :rtype: str
        """
        # see: https://stackoverflow.com/q/46967465/, https://regex101.com/
        split_by_quotes = cls.MATCHING_QUOTES_REGEX.split(string)
        quoted = split_by_quotes[2::4]
        for ws, sub in cls.WHITESPACE_SUBSTITUTES.items():
            quoted = [w.replace(ws, sub) for w in quoted]
        split_by_quotes[2::4] = quoted
        return ''.join(split_by_quotes)

    @classmethod
    def release_quote(cls, string):
        """
        Change the substitute characters in supplied `string` back
        to whitespace, remove matching outer quotation marks, and return string

        :param string: text where whitespace will be reverted and quotes removed
        :type string: str
        :return: modified output string
        :rtype: str
        """
        new_str = ''.join(cls.MATCHING_OUTER_QUOTES_REGEX.split(string)[::2])
        for ws, sub in cls.WHITESPACE_SUBSTITUTES.items():
            new_str = new_str.replace(sub, ws)
        return new_str

    @classmethod
    def strip_comments(cls, string):
        """
        Remove everything following "#" at the start of line or after whitespace

        :param string: string where comments should be removed
        :type string: str
        :return: string with comments removed
        :rtype: str
        """
        return cls.COMMENT_REGEX.split(string)[0]

# TODO: Looks like reader has some problems with reading olex2 files. check


class CifWriterBuffer(CifIOBuffer):
    """Buffer for writing data from `CifReader` into cif file """

    MAX_NAME_LENGTH = 33
    MAX_LINE_LENGTH = 80
    MIN_STEP_LENGTH = 2
    WHITESPACE = {' ', '\t', '\n'}

    def __init__(self, target):
        super().__init__(target=target)
        self.target: TextIO = target
        self.data = OrderedDict()
        self.current__category = ''
        self.current__list = False
        self.current_len = 0

    def add(self, data: tuple):
        k_, v_ = data
        k__category = cif_core_validator.get__category(k_)
        k__list = cif_core_validator.get__list(k_) or isinstance(v_, list)
        v_len = len(v_) if isinstance(v_, list) else 0

        if all([self.current__category == k__category,
                self.current__list == k__list,
                self.current_len == v_len]):
            self.names.append(k_)
            self.values.append(v_)
        else:
            self.flush()
            self.names = [k_]
            self.values = [v_]
            self.current__category = k__category
            self.current__list = k__list
            self.current_len = v_len

    def flush(self):
        s = '\n'
        if self.current__list is True:
            s += self.format_table()
        else:
            for n_, v_ in zip(self.names, self.values):
                s += self.format_line(n_, v_) + '\n'
        self.target.write(s)

    def format_line(self, k, v):
        name_string = f'{k:<{self.MAX_NAME_LENGTH}}'
        step_string = ' ' * self.MIN_STEP_LENGTH
        value_string = self.enquote(v)
        if len(name_string + step_string + value_string) > self.MAX_LINE_LENGTH:
            step_string = '\n '
            value_string = self.enquote(v, force=True)
        if value_string.startswith(';'):
            step_string = '\n'
        return name_string + step_string + value_string

    def format_table(self):
        column_widths = [max(map(len, v)) for v in self.values]
        if sum(column_widths) + len(column_widths) >= self.MAX_LINE_LENGTH:
            pass  # TODO: break long loop tables rows into multiple
        formatted_string = '_loop\n'
        for name in self.names:
            formatted_string += f' {name}\n'
        for value_row in list(map(list, zip(*self.values))):
            enquoted_value_row = map(self.enquote, value_row)
            formatted_string += f' {" ".join(enquoted_value_row)}\n'
        return formatted_string

    def enquote(self, text, force=False):
        if any(whitespace in text for whitespace in self.WHITESPACE) or force:
            if '\n' in text:
                quoted = f';{text}\n;'
            elif "'" not in text:
                quoted = f"'{text}'"
            elif '"' not in text:
                quoted = f'"{text}"'
            elif '\n;' not in text:
                quoted = f';{text}\n;'
            else:
                raise ValueError(f'Unable to quote text: {text}')
        else:
            quoted = text
        return quoted


class CifWriter(CifIO):
    """A helper class managing writing `CifFrame`s into cif files"""

    def write(self, cif_frame):
        with open(self.file_path, 'w') as cif_file:
            buffer = CifWriterBuffer(target=cif_file)
            first_block = True
            for block_name, block in cif_frame.items():
                if not first_block:
                    cif_file.write('\n\n')
                cif_file.write(f'data_{block_name}')
                for data in block.items():
                    buffer.add(data)
                buffer.flush()
                first_block = False




cif_core_validator = CifValidator()


if __name__ == '__main__':
    c = CifFrame()
    c.read('/home/dtchon/x/HiPHAR/anders_script/RFpirazB_cplt100.fractional.cif1')
    c.write('/home/dtchon/x/HiPHAR/anders_script/out.cif')


