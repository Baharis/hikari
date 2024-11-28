import warnings
from collections.abc import Buffer
from copy import deepcopy
from functools import reduce
from operator import and_
from pathlib import Path
import pickle
import re
from typing import List, Union

import numpy as np
import pandas as pd

from hikari.resources import point_groups_pickle, space_groups_pickle
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
    name = 'n_c'


class GroupCatalogKeyNumber(GroupCatalogKey):
    """Integer assigned to each groups in ITC A, shared by all settings of a group"""
    name = 'number'
    accessor_priority = 350.
    dependencies = [GroupCatalogKeyNC]

    @classmethod
    def construct(cls, table: pd.DataFrame) -> pd.Series:
        return (table['n_c'] + ':').str.split(':', expand=True).iloc[:, 0].astype(int)


class GroupCatalogKeySetting(GroupCatalogKey):
    """Unique setting symbol (typically axis dir.) distinguishes same-number groups"""
    name = 'setting'
    dependencies = [GroupCatalogKeyNC]

    @classmethod
    def construct(cls, table: pd.DataFrame) -> pd.Series:
        return table['n_c'].str.split(':', expand=True).iloc[:, -1].fillna('').astype(str)


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
            # symbol = 'p_'.join(re.fullmatch(r'(-?)(.+)', symbol).groups())
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


# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ CATALOG WARNING ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ #


class AmbiguousGroupAccessorWarning(UserWarning):
    """Raised if the accessors provided to get the `Group` match >1 `Group`"""


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
    REST_COL_FORMAT = {
        'n_c': '7.7s',
        'Hall': '16.16s',  # 14 is the longest
        'HM': '16.16s',  # 12 is the longest
        'HM_short': '11.11s',  # 9 is the longest
        'HM_simple': '11.11s',  # 7 is the longest

    }

    def __init__(self, table: pd.DataFrame) -> None:
        for key in _resolve_construct_order(self.KEYS):
            if key.name not in table:
                table[key.name] = key.construct(table)
        self.table: pd.DataFrame = table

    def __eq__(self, other):
        if isinstance(other, self.__class__):
            return self.table.equals(other.table)
        return NotImplemented

    @classmethod
    def from_bytes(cls, bytes_: Buffer):
        """Load directly from bytes object. Drastically speeds up library load"""
        new = pickle.loads(bytes_)
        if isinstance(new, cls):
            return new
        else:
            raise TypeError(f'Loaded pickle is not instance of {cls}')

    @classmethod
    def from_pickle(cls, pickle_path: PathLike) -> 'GroupCatalog':
        """Load from pickled bytes file. Drastically speeds up library load"""
        with open(pickle_path, 'br') as pickled_bytes:
            return cls.from_bytes(pickled_bytes.read())

    @classmethod
    def from_ssv(cls, ssv_path: PathLike) -> 'GroupCatalog':
        """New from a .ssv file. Used for development or custom `GroupCatalog`s."""
        with open(ssv_path, 'r') as ssv_file:
            table = pd.read_csv(ssv_file, comment='#', sep=r'\s+')
        return cls(table)

    def to_rest_table(self, txt_path: PathLike) -> None:
        """Generate a .txt file with ReST table containing `GroupCatalog` elements."""
        r = self.REST_COL_FORMAT
        field_len_re = re.compile(r'^(\d*)(?:\.\d+)?[a-z]*$')
        field_lens = [int(field_len_re.match(v).group(1)) for v in r.values()]
        sep = '+-' + '-+-'.join(['-' * f for f in field_lens]) + '-+'
        fmt = '| ' + ' | '.join([f'{{:{f}.{f}s}}' for f in field_lens]) + ' |'
        lines = [fmt.format(*r.keys())]
        underscore_escape = re.compile(r'_(?=[-(])')
        for t in self.table[list(r.keys())].itertuples(index=False):
            t = [underscore_escape.sub(repl=r'\_', string=tt) for tt in t]
            lines.append(fmt.format(*t))
        rest = sep + '\n' + ('\n' + sep + '\n').join(lines) + '\n' + sep
        with open(txt_path, 'w') as txt_file:
            txt_file.write(rest)

    def to_pickle(self, pickle_path: PathLike) -> None:
        r"""
        Dump a `GroupCatalog` to a .pkl file.
        Used for development purposes or saving custom `GroupCatalog`s.
        Loading pickles is drastically faster than generating Catalog from ssv.

        .. code-block:: python

            PG = PointGroupCatalog.from_ssv(r'hikari\resources\Hall_symbols_PG.dat')
            SG = SpaceGroupCatalog.from_ssv(r'hikari\resources\Hall_symbols_SG.dat')

        """
        with open(pickle_path, 'bw') as pickle_file:
            pickle.dump(self, file=pickle_file)  # noqa - this code is OK

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

    def _get_by_key(self, key: Union[str, int]) -> pd.DataFrame:
        """Iterate over accessors; whenever key is in accessor, return matching group"""
        matching = []
        for accessor in self.accessors:
            matching.append(self.table[self.table.loc[:, accessor.name] == key])
        return pd.concat(matching, axis=0)

    def _get_by_kwargs(self, **kwargs) -> pd.DataFrame:
        """Return the first group that matches all queries specified in kwargs"""
        masks = []
        for key, value in kwargs.items():
            masks.append(self.table[key].eq(value))
        return self.table[reduce(and_, masks)]

    def get(self, key: Union[str, int] = None, **kwargs) -> Union[Group, None]:
        """Get first `Group` matching provided anonymous&known accessors or None"""
        got1 = deepcopy(self._get_by_key(key)) if key else pd.DataFrame()
        got2 = deepcopy(self._get_by_kwargs(**kwargs)) if kwargs else pd.DataFrame()
        if len(got1) == 0 and len(got2) == 0:
            return
        elif len(got1) > 0 and len(got2) == 0:
            got = got1
        elif len(got1) == 0 and len(got2) > 0:
            got = got2
        else:
            got = got1.merge(got2, how='inner', on=['n_c'])
        got.drop_duplicates(subset='n_c', inplace=True, ignore_index=True)
        got.sort_values(by='standard', kind='stable', ascending=False, inplace=True)
        first_got = got['group'].iloc[0] if len(got) else None
        if len(got) > 1:
            matches = list(got['HM_numbered'])
            msg = f'get({key=}, {kwargs=}) yielded multiple {matches=}. '\
                  f'Returning first in standard setting, if possible: {first_got.name}'
            warnings.warn(msg, AmbiguousGroupAccessorWarning)
        return deepcopy(first_got)

    def __getitem__(self, item: Union[str, int]) -> Group:
        """Get first `Group` matching provided anonymous accessor or raise"""
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


def regenerate_group_catalog_pickles():
    r"""
    This function replaces current `resources/*_group.pickle`s.
    It should be run from hikari's source directory with hikari imported
    as module whenever any changes to `GroupCatalog` class are made.
    """
    pg = PointGroupCatalog.from_ssv(Path('resources/Hall_symbols_PG.dat'))
    sg = SpaceGroupCatalog.from_ssv(Path('resources/Hall_symbols_SG.dat'))
    pg.to_pickle(Path('resources/point_groups.pickle'))
    sg.to_pickle(Path('resources/space_groups.pickle'))


PG = PointGroupCatalog.from_bytes(point_groups_pickle)
SG = SpaceGroupCatalog.from_bytes(space_groups_pickle)
