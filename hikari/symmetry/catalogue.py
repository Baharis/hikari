from copy import deepcopy
import pickle
import re
from typing import List

import numpy as np
import pandas as pd

from hikari.utility.typing import PathLike
from hikari.symmetry.group import Group


# ~~~~~~~~~~~~~~~~~~~~~~~~~~ CATALOGUE KEY REGISTRY ~~~~~~~~~~~~~~~~~~~~~~~~~~ #


class GroupCatalogueKeyRegistrar(type):
    """Metaclass for `GroupCatalogueKey`s, registers these that define `name`"""
    REGISTRY = {}

    def __new__(mcs, name, bases, attrs):
        new_cls = type.__new__(mcs, name, bases, attrs)
        if hasattr(new_cls, 'name') and new_cls.name:
            mcs.REGISTRY[new_cls.name] = new_cls
        return new_cls

    @property
    def accessors(self):
        """Lists `GroupCatalogueKey`s whose accessor priority is not 0"""
        return [k for k, v in self.REGISTRY.items() if v.accessor_priority]


class GroupCatalogueKey(metaclass=GroupCatalogueKeyRegistrar):
    """Base Class for every GroupCatalogueKey."""
    name: str = ''                    # if not empty, how key will be registered
    accessor_priority: float = 0.     # if >0 used to access groups (inc. order)
    dependencies: list = []           # other keys needed to construct this key

    @classmethod
    def construct(cls, table: pd.DataFrame) -> pd.Series:
        """Abstract method to implement if key might have to be constructed"""


# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ CATALOGUE KEYS ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ #


class GroupCatalogueKeyNC(GroupCatalogueKey):
    """Unique group identification string composed of group number:setting"""
    name = 'n:c'


class GroupCatalogueKeyHM(GroupCatalogueKey):
    """Full international Hermann-Mauguin name split with `_` with :setting"""
    name = 'HM'


class GroupCatalogueKeyHall(GroupCatalogueKey):
    name = 'Hall'


class GroupCatalogueKeyGroup(GroupCatalogueKey):
    name = 'group'
    dependencies = [GroupCatalogueKeyHall]

    @classmethod
    def construct(cls, table: pd.DataFrame) -> pd.Series:
        return pd.Series(Group.from_hall_symbol(h) for h in table['Hall'])


class GroupCatalogueKeyNumber(GroupCatalogueKey):
    name = 'number'
    dependencies = [GroupCatalogueKeyNC]

    @classmethod
    def construct(cls, table: pd.DataFrame) -> pd.Series:
        print(table)
        return (table['n:c'] + ':').str.split(':', expand=True).iloc[:, 0].astype(int)


class GroupCatalogueKeySetting(GroupCatalogueKey):
    name = 'setting'
    dependencies = [GroupCatalogueKeyNC]

    @classmethod
    def construct(cls, table: pd.DataFrame) -> pd.Series:
        return table['n:c'].str.split(':', expand=True).iloc[:, -1].fillna('').astype(str)


class GroupCatalogueKeyHMShort(GroupCatalogueKey):
    """Shortened `HM` symbol where all underscores were simply removed"""
    name = 'HM-short'
    dependencies = [GroupCatalogueKeyHM]

    @classmethod
    def construct(cls, table: pd.DataFrame) -> pd.Series:
        return table['HM'].str.replace('_', '')


class GroupCatalogueKeyHMSimple(GroupCatalogueKey):
    """`HM-short` without setting and with `1` removed for monoclinic system"""
    name = 'HM-simple'
    dependencies = [GroupCatalogueKeyHM, GroupCatalogueKeyHMShort, GroupCatalogueKeyGroup]

    @classmethod
    def construct(cls, table: pd.DataFrame) -> pd.Series:
        has_colon = table['HM'].str.contains(':')
        gsm = Group.System.monoclinic
        monoclinic = np.array([g.system == gsm for g in table['group']])
        simple = table['HM-short'].copy()
        simple.loc[has_colon] = simple.loc[has_colon].str.replace(r'\:.', '', regex=True)
        simple.loc[monoclinic] = table.loc[monoclinic, 'HM'].str.replace('_1', '')
        return simple.str.replace('_', '')


class GroupCatalogueKeyHMNumbered(GroupCatalogueKey):
    """Nicely-formatted name with `number: HM-short` to be used in GUIs"""
    name = 'HM-numbered'
    dependencies = [GroupCatalogueKeyNumber, GroupCatalogueKeyHMShort]

    @classmethod
    def construct(cls, table: pd.DataFrame) -> pd.Series:
        return table['number'].astype(str).str.cat(table['HM-short'], sep=': ')


class GroupCatalogueKeyStandard(GroupCatalogueKey):
    name = 'standard'
    dependencies = [GroupCatalogueKeyNumber]

    @classmethod
    def construct(cls, table: pd.DataFrame) -> pd.Series:
        return table['number'].rolling(2).var().ne(0)


def _resolve_construct_order(keys: List[GroupCatalogueKey]) -> List[GroupCatalogueKey]:
    """
    Return `GroupCatalogueKey`s in an order that warrants that
    key's dependencies are constructed before it
    """
    unordered = deepcopy(keys)
    ordered = []
    def find_constructable(unordered_: List[GroupCatalogueKey]) -> GroupCatalogueKey:
        for key_i, key in enumerate(unordered_):
            if all(d in ordered for d in key.dependencies):
                return unordered_.pop(key_i)
            raise RuntimeError('Circular dependency when creating `GroupCatalogue`')
    while len(ordered) < len(keys):
        ordered.append(find_constructable(unordered))
    return ordered


# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ CATALOGUE CLASS ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ #


class GroupCatalogue:
    """
    Manage generating and mappings of point and space groups.
    Relies on a built-in pandas DataFrame `table` to store all the information.
    Individual columns are named & generated based on `GroupCatalogueKey` data.

    Some notes on the uniqueness of columns pairwise for accessing:

    - Column `Hall` has 3 groups appear twice due to inconsistency of HM names:
      `c_2_2_-1ac`, `a_2_2_-1ab`, and `a_2_2_-1ab`.
    - There are no overlaps between `HM` and `Hall` column names
    - There are no overlaps between `HM-short` and `Hall` column names
    - There are no overlaps between `HM-simple` and `Hall` column names
    - There are no overlaps between `HM-simple` and `HM` column names
    - There are 345 overlaps between `HM-simple` and `HM-short` column names
    """

    REST_COL_FORMAT = {'n:c': '7.7s', 'HM-short': '9.9s', 'Hall': '14.14s'}

    def __init__(self, table: pd.DataFrame) -> None:
        for key in _resolve_construct_order(list(GroupCatalogueKey.REGISTRY.values())):
            if key.name not in table:
                table[key.name] = key.construct(table)
        self.table: pd.DataFrame = table

    @classmethod
    def _from_ssv(cls, ssv_path: PathLike) -> 'GroupCatalogue':
        """Generate a `GroupCatalogue` from a .ssv file"""
        with open(ssv_path, 'r') as ssv_file:
            table = pd.read_csv(ssv_file, comment='#', sep=r'\s+')
        return cls(table)

    @classmethod
    def _from_pickle(cls, pickle_path: PathLike) -> 'GroupCatalogue':
        """Load a `GroupCatalogue` from a .pickle file"""
        with open(pickle_path, 'br') as pickle_file:
            new = pickle.load(pickle_file)
        if isinstance(new, cls):
            return new
        else:
            raise TypeError(f'Loaded pickle is not instance of {cls}')

    def _as_rest_table(self, txt_path: PathLike) -> None:
        """Write an instant of `GroupCatalogue` to a .txt file as ReST table"""
        r = self.REST_COL_FORMAT
        field_len_re = re.compile(r'^(\d*)(?:\.\d+)?[a-z]*$')
        field_lens = [int(field_len_re.match(v).group(1)) for v in r.values()]
        sep = '+-' + '-+-'.join(['-' * f for f in field_lens]) + '-+'
        fmt = '| ' + ' | '.join([f'{{:{f}.{f}s}}' for f in field_lens]) + ' |'
        lines = [fmt.format(*r.keys())]
        for t in self.table[list(r.keys())].itertuples(index=False):
            lines.append(fmt.format(*t))
        rest = sep + '\n' + ('\n' + sep + '\n').join(lines) + '\n' + sep
        with open(txt_path, 'w') as txt_file:
            txt_file.write(rest)

    def _to_pickle(self, pickle_path: PathLike) -> None:
        """Development only; Dump `self` to a .pickle file"""
        with open(pickle_path, 'bw') as pickle_file:
            pickle.dump(self, file=pickle_file)

    @property
    def standard(self) -> 'GroupCatalogue':
        """A subset of current catalogue with standard-setting groups only"""
        standard = deepcopy(self.table[self.table['standard']]).reset_index()
        return self.__class__(standard)
