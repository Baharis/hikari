import warnings
from copy import deepcopy
from functools import reduce
from operator import or_
import pickle
import re
from typing import List, Union

import numpy as np
import pandas as pd
from pandas.core.interchange.dataframe_protocol import DataFrame

from hikari.utility.typing import PathLike
from hikari.symmetry.group import Group


# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ CATALOG KEYS ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ #


class GroupCatalogKey:
    """Base Class for every following `GroupCatalogKey`. Each named
    `GroupCatalogKey` represents a single column in the `GroupCatalog` table."""
    name: str = ''                    # if not empty, how key will be registered
    accessor_priority: float = 0.     # if >0 used to access groups (inc. order)
    dependencies: list = []           # other keys needed to construct this key

    @classmethod
    def construct(cls, table: pd.DataFrame) -> pd.Series:
        """Abstract method to implement only if key might have to be constructed"""


class GroupCatalogKeyNC(GroupCatalogKey):
    """Unique group identification string composed of group number:setting"""
    name = 'n:c'


class GroupCatalogKeyNumber(GroupCatalogKey):
    """Integer assigned to each groups in ITC A, shared by all settings of a group"""
    name = 'number'
    accessor_priority = 350.
    dependencies = [GroupCatalogKeyNC]

    @classmethod
    def construct(cls, table: pd.DataFrame) -> pd.Series:
        return (table['n:c'] + ':').str.split(':', expand=True).iloc[:, 0].astype(int)


class GroupCatalogKeySetting(GroupCatalogKey):
    """Unique setting symbol (typically axis dir.) distinguishes same-number groups"""
    name = 'setting'
    dependencies = [GroupCatalogKeyNC]

    @classmethod
    def construct(cls, table: pd.DataFrame) -> pd.Series:
        return table['n:c'].str.split(':', expand=True).iloc[:, -1].fillna('').astype(str)


class GroupCatalogKeyHM(GroupCatalogKey):
    """Full international Hermann-Mauguin name split with `_` with :setting"""
    name = 'HM'
    accessor_priority = 150.


class GroupCatalogKeyHall(GroupCatalogKey):
    """Full Hall symbol of given group used to recreate it"""
    name = 'Hall'
    accessor_priority = 250.


class GroupCatalogKeyGroup(GroupCatalogKey):
    """Abstract base class for point and space group columns / constructors"""
    name = 'group'
    dependencies = [GroupCatalogKeyNumber, GroupCatalogKeyHM, GroupCatalogKeyHall]


class GroupCatalogKeyPointGroup(GroupCatalogKeyGroup):
    """Point and space groups are identical, so must recreate full Hall symbol"""

    @classmethod
    def construct(cls, table: pd.DataFrame) -> pd.Series:
        groups = []
        for n, name, symbol in zip(table['number'], table['HM'], table['Hall']):
            symbol = 'p_'.join(re.fullmatch(r'(-?)(.+)', symbol).groups())
            group = Group.from_hall_symbol(symbol)
            group.name = name
            group.number = n
            groups.append(group)
        return pd.Series(groups)


class GroupCatalogKeySpaceGroup(GroupCatalogKeyGroup):
    """Creates and names `hikari.symmetry.Group` object based on Hall symbol"""

    @classmethod
    def construct(cls, table: pd.DataFrame) -> pd.Series:
        groups = []
        for n, name, symbol in zip(table['number'], table['HM'], table['Hall']):
            group = Group.from_hall_symbol(symbol)
            group.name = name
            group.number = n
            groups.append(group)
        return pd.Series(groups)


class GroupCatalogKeyHMShort(GroupCatalogKey):
    """Shortened `HM` symbol where all underscores were simply removed"""
    name = 'HM_short'
    accessor_priority = 140.
    dependencies = [GroupCatalogKeyHM]

    @classmethod
    def construct(cls, table: pd.DataFrame) -> pd.Series:
        return table['HM'].str.replace('_', '')


class GroupCatalogKeyHMSimple(GroupCatalogKey):
    """`HM_short` without setting and with `1` removed for monoclinic system"""
    name = 'HM_simple'
    accessor_priority = 130.
    dependencies = [GroupCatalogKeyHM, GroupCatalogKeyHMShort, GroupCatalogKeyGroup]

    @classmethod
    def construct(cls, table: pd.DataFrame) -> pd.Series:
        has_colon = table['HM'].str.contains(':')
        gsm = Group.System.monoclinic
        monoclinic = np.array([g.system == gsm for g in table['group']])
        simple = table['HM_short'].copy()
        simple.loc[has_colon] = simple.loc[has_colon].str.replace(r'\:.', '', regex=True)
        simple.loc[monoclinic] = '_' + table.loc[monoclinic, 'HM']  # needed for PG
        simple.loc[monoclinic] = simple.loc[monoclinic].str.replace('_1', '')
        return simple.str.replace('_', '')


class GroupCatalogKeyHMNumbered(GroupCatalogKey):
    """Nicely-formatted name with `number: HM_short` to be used in GUIs"""
    name = 'HM_numbered'
    dependencies = [GroupCatalogKeyNumber, GroupCatalogKeyHMShort]

    @classmethod
    def construct(cls, table: pd.DataFrame) -> pd.Series:
        return table['number'].astype(str).str.cat(table['HM_short'], sep=': ')


class GroupCatalogKeyStandard(GroupCatalogKey):
    """Boolean column with True if setting is standard (first with given number)"""
    name = 'standard'
    dependencies = [GroupCatalogKeyNumber]

    @classmethod
    def construct(cls, table: pd.DataFrame) -> pd.Series:
        return table['number'].rolling(2).var().ne(0)


# ~~~~~~~~~~~~~~~~~~~~~~~~~ CATALOG HELPER FUNCTIONS ~~~~~~~~~~~~~~~~~~~~~~~~~ #


def _resolve_construct_order(keys: List[GroupCatalogKey]) -> List[GroupCatalogKey]:
    """
    Return `GroupCatalogKey`s in an order that warrants that
    key's dependencies are constructed before it
    """
    unordered = deepcopy(keys)
    ordered = []
    def find_constructable(unordered_: List[GroupCatalogKey]) -> GroupCatalogKey:
        for key_i, key in enumerate(unordered_):
            if all(any(issubclass(o, d) for o in ordered) for d in key.dependencies):
                return unordered_.pop(key_i)
        raise RuntimeError('Circular dependency when creating `GroupCatalog`')
    while len(ordered) < len(keys):
        ordered.append(find_constructable(unordered))
    return ordered


# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ CATALOG CLASS ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ #


class GroupCatalog:
    """
    Manage generating and mappings of point and space groups.
    Relies on a built-in pandas DataFrame `table` to store all the information.
    Individual columns are named & generated based on `GroupCatalogKey` data.

    Some notes on the uniqueness of columns pairwise for accessing:

    - Column `Hall` has 3 groups appear twice due to inconsistency of HM names:
      `c_2_2_-1ac`, `a_2_2_-1ab`, and `b_2_2_-1ab`.
    - There are no overlaps between `HM` and `Hall` column names
    - There are no overlaps between `HM_short` and `Hall` column names
    - There are no overlaps between `HM_simple` and `Hall` column names
    - There are no overlaps between `HM_simple` and `HM` column names
    - There are 345 overlaps between `HM_simple` and `HM_short` column names
    """
    KEYS: List[GroupCatalogKey] = [
        GroupCatalogKeyNC,
        GroupCatalogKeyNumber,
        GroupCatalogKeySetting,
        GroupCatalogKeyHM,
        GroupCatalogKeyHall,
        GroupCatalogKeyHMShort,
        GroupCatalogKeyHMSimple,
        GroupCatalogKeyHMNumbered,
        GroupCatalogKeyStandard,
    ]
    REST_COL_FORMAT = {'n:c': '7.7s', 'HM_short': '9.9s', 'Hall': '14.14s'}

    def __init__(self, table: pd.DataFrame) -> None:
        for key in _resolve_construct_order(self.KEYS):
            if key.name not in table:
                table[key.name] = key.construct(table)
        self.table: pd.DataFrame = table

    @classmethod
    def from_ssv(cls, ssv_path: PathLike) -> 'GroupCatalog':
        """Generate a `GroupCatalog` from a .ssv file"""
        with open(ssv_path, 'r') as ssv_file:
            table = pd.read_csv(ssv_file, comment='#', sep=r'\s+')
        return cls(table)

    @classmethod
    def from_pickle(cls, pickle_path: PathLike) -> 'GroupCatalog':
        """Load a `GroupCatalog` from a .pickle file"""
        with open(pickle_path, 'br') as pickle_file:
            new = pickle.load(pickle_file)
        if isinstance(new, cls):
            return new
        else:
            raise TypeError(f'Loaded pickle is not instance of {cls}')

    def as_rest_table(self, txt_path: PathLike) -> None:
        """Write an instant of `GroupCatalog` to a .txt file as ReST table"""
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

    def to_pickle(self, pickle_path: PathLike) -> None:
        """Development only; Dump `self` to a .pickle file"""
        with open(pickle_path, 'bw') as pickle_file:
            pickle.dump(self, file=pickle_file)

    @property
    def accessors(self) -> List['GroupCatalogKey']:
        """Lists `cls.KEYS whose accessor priority is not 0 in decreasing order"""
        accessors = [k for k in self.KEYS if k.accessor_priority]
        return sorted(accessors, key=lambda a: a.accessor_priority, reverse=True)

    @property
    def standard(self) -> 'GroupCatalog':
        """A subset of current catalog with standard-setting groups only"""
        standard = deepcopy(self.table[self.table['standard']]).reset_index()
        return self.__class__(standard)

    # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~ SMART GETTERS ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ #

    def _get_by_key(self, key: Union[str, int]) -> DataFrame:
        """Iterate over accessors; whenever key is in accessor, return matching group"""
        matching = []
        for accessor in self.accessors:
            matching.append(self.table[self.table.loc[:, accessor.name] == key])
        return pd.concat(matching, axis=0)

    def _get_by_kwargs(self, **kwargs) -> DataFrame:
        """Return the first group that matches all queries specified in kwargs"""
        masks = []
        for key, value in kwargs.items():
            masks.append(self.table[key].eq(value))
        return self.table[reduce(or_, masks)]

    def get(self, key: Union[str, int] = None, **kwargs) -> Union[Group, None]:
        got1 = self._get_by_key(key) if key else pd.DataFrame()
        got2 = self._get_by_kwargs(**kwargs) if kwargs else pd.DataFrame()
        if len(got1) == 0 and len(got2) == 0:
            return
        elif len(got1) > 0 and len(got2) == 0:
            got = got1
        elif len(got1) == 0 and len(got2) > 0:
            got = got2
        else:
            got = got1.merge(got2, how='inner', on=['n:c'])
        got.drop_duplicates(subset='n:c', inplace=True, ignore_index=True)
        if len(got) > 1:
            matches = list(got['n:c'])
            warnings.warn(f'get({key=}, {kwargs=}) yielded multiple {matches=}. '
                          f'Returning first result in standard setting, if possible')
            std = self.standard.table
            for nc in matches:
                if nc in std['n:c']:
                    return deepcopy(std.loc[std.loc['n:c'] == nc]['group'][0])
        return deepcopy(got['group'].iloc[0] if len(got) else None)

    def __getitem__(self, item: Union[str, int]) -> Group:
        got = self.get(key=item)
        if not got:
            raise KeyError(f'Unknown key: {item}')
        return got


class PointGroupCatalog(GroupCatalog):
    """Subclass of `GroupCatalog` specialized to hold space groups"""
    KEYS = GroupCatalog.KEYS + [GroupCatalogKeyPointGroup]


class SpaceGroupCatalog(GroupCatalog):
    """Subclass of `GroupCatalog` specialized to hold space groups"""
    KEYS = GroupCatalog.KEYS + [GroupCatalogKeySpaceGroup]

