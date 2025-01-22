import warnings
from copy import deepcopy
from functools import reduce
import json
from operator import and_
from pathlib import Path
import re
from typing import Any, Union

import numpy as np
import pandas as pd

from hikari.resources import (point_groups_json, space_groups_json,
                              point_groups_dataframe, space_groups_dataframe)
from hikari.symmetry.operations import BoundedOperation
from hikari.utility.typing import PathLike
from hikari.symmetry.group import Group


# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ CATALOG KEYS ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ #


class GroupCatalogKey:
    """Base Class for every following `GroupCatalogKey`. Each named
    `GroupCatalogKey` represents a single column in the `GroupCatalog` table."""
    name: str = ''                  # if not empty, how key will be registered
    accessor_priority: float = 0.   # if > 0 used to access groups (inc. order)
    dtype: Union[str, Any] = 'str'  # used for type-casting at column creation
    dependencies: list = []         # other keys needed to construct this key

    @classmethod
    def construct(cls, table: pd.DataFrame) -> pd.Series:
        """Abstract method to implement only if key might have to be constructed"""


class GroupCatalogKeyCompulsory(GroupCatalogKey):
    """Subclass of `GroupCatalogKey`, specifies key must exist in input table"""

    @classmethod
    def construct(cls, table: pd.DataFrame):
        if cls.name not in table:
            raise KeyError('Compulsory `GroupCatalogKey` missing: ' + cls.name)


class GroupCatalogKeyNC(GroupCatalogKeyCompulsory):
    """Unique group identification string composed of group number:setting"""
    name = 'n_c'


class GroupCatalogKeyNumber(GroupCatalogKey):
    """Integer assigned to each groups in ITC A, shared by all settings of a group"""
    name = 'number'
    accessor_priority = 350.
    dependencies = [GroupCatalogKeyNC]
    dtype = 'int16'

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


class GroupCatalogKeyHM(GroupCatalogKeyCompulsory):
    """Full international Hermann-Mauguin name split with `_` with :setting"""
    name = 'HM'
    accessor_priority = 150.


class GroupCatalogKeyHall(GroupCatalogKeyCompulsory):
    """Full Hall symbol of given group used to recreate it"""
    name = 'Hall'
    accessor_priority = 250.


class GroupCatalogKeyGroup(GroupCatalogKey):
    """Creates and names `hikari.symmetry.Group` object based on Hall symbol"""
    name = 'group'
    dtype = Group
    dependencies = [GroupCatalogKeyNumber, GroupCatalogKeyHM, GroupCatalogKeyHall]

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
    dtype = bool

    @classmethod
    def construct(cls, table: pd.DataFrame) -> pd.Series:
        return table['number'].rolling(2).var().ne(0)


# ~~~~~~~~~~~~~~~~~~~~~~~~~ CATALOG HELPER FUNCTIONS ~~~~~~~~~~~~~~~~~~~~~~~~~ #


def _resolve_construct_order(keys: list[GroupCatalogKey]) -> list[GroupCatalogKey]:
    """
    Return `GroupCatalogKey`s in an order that warrants that
    key's dependencies are constructed before it
    """
    unordered = deepcopy(keys)
    ordered = []
    def find_constructable(unordered_: list[GroupCatalogKey]) -> GroupCatalogKey:
        for key_i, key in enumerate(unordered_):
            if all(any(issubclass(o, d) for o in ordered) for d in key.dependencies):
                return unordered_.pop(key_i)
        raise RuntimeError('Circular dependency when creating `GroupCatalog`')
    while len(ordered) < len(keys):
        ordered.append(find_constructable(unordered))
    return ordered


# ~~~~~~~~~~~~~~~~~~~~~~~ CATALOG JSON ENCODER/DECODER ~~~~~~~~~~~~~~~~~~~~~~~ #


class GroupCatalogJSONEncoder(json.JSONEncoder):
    """Handles serialization of `GroupCatalog` to a `.json` file"""
    def default(self, gc: 'GroupCatalog') -> dict:
        if not isinstance(gc, GroupCatalog):
            return super().default(gc)  # raises TypeError for other input
        records = []
        for row in gc.table.itertuples(index=False):
            record = dict(row._asdict())  # noqa - It should be protected
            record['group'] = {
                'generators': [g.code for g in record['group'].generators],
                'operations': [o.code for o in record['group'].operations],
            }
            records.append(record)
        return {'_type': gc.__class__.__name__, 'table': records}


class GroupCatalogJSONDecoder(json.JSONDecoder):
    """Handles deserialization of `GroupCatalog` from a `.json` file"""
    def __init__(self, *args, **kwargs) -> None:
        json.JSONDecoder.__init__(self, object_hook=self.object_hook, *args, **kwargs)

    def object_hook(self, obj) -> Any:  # noqa - signature must match
        if '_type' not in obj:
            return obj
        _type = obj['_type']
        if _type == 'GroupCatalog':
            records = []
            for record in obj['table']:
                group = record['group']
                group = Group.from_generators_operations(
                    generators=[BoundedOperation.from_code(c) for c in group['generators']],
                    operations=[BoundedOperation.from_code(c) for c in group['operations']])
                group.name = record['HM']
                group.number = record['number']
                record.update({'group': group})
                records.append(record)
            return GroupCatalog(table=pd.DataFrame.from_records(records))
        return obj


# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ CATALOG WARNING ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ #


class AmbiguousGroupAccessorWarning(UserWarning):
    """Raised if the accessors provided to get the `Group` match >1 `Group`"""


# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ CATALOG CLASS ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ #


class GroupCatalog:
    """
    Manage generating and mappings of point and space groups.
    Relies on a built-in pandas DataFrame `table` to store all the information.
    Individual columns are named & generated based on `GroupCatalogKey` data.

    For

    Some notes on the uniqueness of columns pairwise for accessing:

    - Column `Hall` has 3 groups appear twice due to inconsistency of HM names:
      `c_2_2_-1ac`, `a_2_2_-1ab`, and `b_2_2_-1ab`.
    - There are no overlaps between `HM` and `Hall` column names
    - There are no overlaps between `HM_short` and `Hall` column names
    - There are no overlaps between `HM_simple` and `Hall` column names
    - There are no overlaps between `HM_simple` and `HM` column names
    - There are 345 overlaps between `HM_simple` and `HM_short` column names
    """
    KEYS: list[GroupCatalogKey] = [
        GroupCatalogKeyNC,
        GroupCatalogKeyNumber,
        GroupCatalogKeySetting,
        GroupCatalogKeyHM,
        GroupCatalogKeyHall,
        GroupCatalogKeyGroup,
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
        if 'n_c' not in table and table.index.name == 'n_c':
            table.reset_index(inplace=True)
        for key in _resolve_construct_order(self.KEYS):
            if key.name not in table:
                table[key.name] = key.construct(table)
            gck = {k.name: k for k in self.KEYS}[key.name]
            if not (table[key.name].dtype == gck.dtype or (
                    isinstance(gck.dtype, type) and
                    all(isinstance(v, gck.dtype) for v in table[key.name]))):
                table[key.name] = table[key.name].astype(gck.dtype)
        self.table: pd.DataFrame = table

    def __eq__(self, other):
        if isinstance(other, self.__class__):
            return self.table.equals(other.table)
        return NotImplemented

    def __len__(self) -> int:
        return len(self.table)

    @classmethod
    def from_json(cls, text: str) -> 'GroupCatalog':
        """Load from a json-formatted string"""
        return json.loads(text, cls=GroupCatalogJSONDecoder)

    def to_json(self, json_path: PathLike) -> None:
        with open(json_path, 'w') as json_file:
            # noinspection PyTypeChecker
            json.dump(obj=self, fp=json_file, cls=GroupCatalogJSONEncoder,
                      ensure_ascii=False, indent=4)

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

    @property
    def accessors(self) -> list['GroupCatalogKey']:
        """Lists `cls.KEYS whose accessor priority is not 0 in decreasing order"""
        accessors = [k for k in self.KEYS if k.accessor_priority]
        return sorted(accessors, key=lambda a: a.accessor_priority, reverse=True)

    @property
    def standard(self) -> 'GroupCatalog':
        """A subset of current catalog with standard-setting groups only"""
        standard = deepcopy(self.table[self.table['standard']]).reset_index(drop=True)
        return self.__class__(standard)

    # ~~~~~~~~~~~~~~~~~~~~ DUCK-TYPING DICT-LIKE INTERFACE ~~~~~~~~~~~~~~~~~~~ #

    def keys(self) -> list[str]:
        return list(self.table['n_c'])

    def values(self) -> list[Group]:
        return list(self.table['group'])

    def items(self) -> list[tuple[Union[int, str], Group]]:
        return [(k, v) for k, v in zip(self.table['n_c'], self.table['group'])]

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
            msg = (f'get(key={key}, kwargs={kwargs}) yielded multiple matches: '
                   f'{matches}. Returning first in standard setting, '
                   f'if possible: {first_got.name}')
            warnings.warn(msg, AmbiguousGroupAccessorWarning)
        return deepcopy(first_got)

    def __getitem__(self, item: Union[str, int]) -> Group:
        """Get first `Group` matching provided anonymous accessor or raise"""
        got = self.get(key=item)
        if not got:
            raise KeyError(f'Unknown key: {item}')
        return got


def regenerate_group_catalog_jsons():
    r"""
    This function regenerates current `resources/*_group.json`s from `.wsv`s.
    It should be run from hikari's parent directory with hikari imported
    as module whenever any changes to `GroupCatalog` class are made.
    """
    pg = GroupCatalog(point_groups_dataframe)
    sg = GroupCatalog(space_groups_dataframe)
    pg.to_json(Path('hikari/resources/point_groups.json'))
    sg.to_json(Path('hikari/resources/space_groups.json'))


# ~~~~~~~~~~~~~~~~~~~~~~~~ PRE-DEFINED GROUPS CATALOGS ~~~~~~~~~~~~~~~~~~~~~~~~ #


PG = GroupCatalog.from_json(point_groups_json)
r"""
Since hikari's groups do not carry information about lattice translations,
hikari does not differentiate between point groups and space groups.
As a result, all groups are instances of the same class `Group`,
and are stored in almost identical `GroupCatalog`s.

This pre-defined `GroupCatalog` `PG` holds and provides access
to all pre-defined point groups available in hikari.
The point groups are named and generated based on a wsv file in resources.
Individual groups can be accessed in two different ways:

  - using a dict-like anonymous accessor: `PG[32]`, `PG['m-3m']`;
  - using `get()` method: `PG.get(number=32)`, `PG.get(HM_short='m-3m')`;

The `get()` method can use a combination of various keywords to unambiguously
localize the correct space group. The following keyword arguments can be used:

  - `n_c` - unique "number:setting" string identifying each group
  - `number` - index of the point group as given in ICT A (integer)\*
  - `setting` - arbitrary string declaring group setting
  - `HM` - Full underscore-delimited Hermann-Mauguin symbol\*
  - `HM_short` - Hermann-Mauguin symbol with underscores removed\*
  - `HM_simple` - Short Hermann-Mauguin symbol without setting information\*
  - `HM_numbered` - "number: Short Hermann-Mauguin symbol" string
  - `Hall` - Full Hall symbol\*
  - `standard` - True for groups in standard setting only

The keywords marked with "*" are called "accessors" and can be used to `get`
the group using bracket notation `get[accessor]`. If `get` finds ambiguity,
If `get` finds ambiguity, for example `PG['2/m']` matches all settings
of point group #5: `1_2/m_1` (std), `1_1_2/m`, and `2/m_1_1`,
the user is warned and the first group in std setting, if possible, is returned.

The following table lists all possible values of selected point `Group`
keywords. Please mind that in the raw docstring all `\` should be ignored.

+---------+------------------+------------------+-------------+-------------+
| n_c     | Hall             | HM               | HM_short    | HM_simple   |
+---------+------------------+------------------+-------------+-------------+
| 1       | 1                | 1                | 1           | 1           |
+---------+------------------+------------------+-------------+-------------+
| 2       | -1               | -1               | -1          | -1          |
+---------+------------------+------------------+-------------+-------------+
| 3:b     | 2y               | 1_2_1            | 121         | 2           |
+---------+------------------+------------------+-------------+-------------+
| 3:c     | 2                | 1_1_2            | 112         | 2           |
+---------+------------------+------------------+-------------+-------------+
| 3:a     | 2x               | 2_1_1            | 211         | 2           |
+---------+------------------+------------------+-------------+-------------+
| 4:b     | -2y              | 1_m_1            | 1m1         | m           |
+---------+------------------+------------------+-------------+-------------+
| 4:c     | -2               | 1_1_m            | 11m         | m           |
+---------+------------------+------------------+-------------+-------------+
| 4:a     | -2x              | m_1_1            | m11         | m           |
+---------+------------------+------------------+-------------+-------------+
| 5:b     | -2y              | 1_2/m_1          | 12/m1       | 2/m         |
+---------+------------------+------------------+-------------+-------------+
| 5:c     | -2               | 1_1_2/m          | 112/m       | 2/m         |
+---------+------------------+------------------+-------------+-------------+
| 5:a     | -2x              | 2/m_1_1          | 2/m11       | 2/m         |
+---------+------------------+------------------+-------------+-------------+
| 6       | 2_2              | 2_2_2            | 222         | 222         |
+---------+------------------+------------------+-------------+-------------+
| 7       | 2\_-2            | m_m_2            | mm2         | mm2         |
+---------+------------------+------------------+-------------+-------------+
| 8       | -2_2             | m_m_m            | mmm         | mmm         |
+---------+------------------+------------------+-------------+-------------+
| 9       | 4                | 4                | 4           | 4           |
+---------+------------------+------------------+-------------+-------------+
| 10      | -4               | -4               | -4          | -4          |
+---------+------------------+------------------+-------------+-------------+
| 11      | -4               | 4/m              | 4/m         | 4/m         |
+---------+------------------+------------------+-------------+-------------+
| 12      | 4_2              | 4_2_2            | 422         | 422         |
+---------+------------------+------------------+-------------+-------------+
| 13      | 4\_-2            | 4_m_m            | 4mm         | 4mm         |
+---------+------------------+------------------+-------------+-------------+
| 14      | -4_2             | -4_2_m           | -42m        | -42m        |
+---------+------------------+------------------+-------------+-------------+
| 15      | -4_2             | 4/m_m_m          | 4/mmm       | 4/mmm       |
+---------+------------------+------------------+-------------+-------------+
| 16      | 3                | 3                | 3           | 3           |
+---------+------------------+------------------+-------------+-------------+
| 17      | -3               | -3               | -3          | -3          |
+---------+------------------+------------------+-------------+-------------+
| 18      | 3_2              | 3_1_2            | 312         | 312         |
+---------+------------------+------------------+-------------+-------------+
| 19      | 3\_-2"           | 3_m_1            | 3m1         | 3m1         |
+---------+------------------+------------------+-------------+-------------+
| 20      | -3_2             | -3_1_m           | -31m        | -31m        |
+---------+------------------+------------------+-------------+-------------+
| 21      | 6                | 6                | 6           | 6           |
+---------+------------------+------------------+-------------+-------------+
| 22      | -6               | -6               | -6          | -6          |
+---------+------------------+------------------+-------------+-------------+
| 23      | -6               | 6/m              | 6/m         | 6/m         |
+---------+------------------+------------------+-------------+-------------+
| 24      | 6_2              | 6_2_2            | 622         | 622         |
+---------+------------------+------------------+-------------+-------------+
| 25      | 6\_-2            | 6_m_m            | 6mm         | 6mm         |
+---------+------------------+------------------+-------------+-------------+
| 26      | -6_2             | -6_m_2           | -6m2        | -6m2        |
+---------+------------------+------------------+-------------+-------------+
| 27      | -6_2             | 6/m_m_m          | 6/mmm       | 6/mmm       |
+---------+------------------+------------------+-------------+-------------+
| 28      | 2_2_3            | 2_3              | 23          | 23          |
+---------+------------------+------------------+-------------+-------------+
| 29      | -2_2_3           | m\_-3            | m-3         | m-3         |
+---------+------------------+------------------+-------------+-------------+
| 30      | 4_2_3            | 4_3_2            | 432         | 432         |
+---------+------------------+------------------+-------------+-------------+
| 31      | -4_2_3           | -4_3_m           | -43m        | -43m        |
+---------+------------------+------------------+-------------+-------------+
| 32      | -4_2_3           | m\_-3_m          | m-3m        | m-3m        |
+---------+------------------+------------------+-------------+-------------+
"""


SG = GroupCatalog.from_json(space_groups_json)
r"""
Since hikari's groups do not carry information about lattice translations,
hikari does not differentiate between point groups and space groups.
As a result, all groups are instances of the same class `Group`,
and are stored in almost identical `GroupCatalog`s.

This pre-defined `SpaceGroupCatalog` `SG` holds and provides access
to all pre-defined space groups available in hikari.
The space groups are named and generated based on a wsv file in resources.
Individual groups can be accessed in two different ways:

  - using a dict-like anonymous accessor: `PG[62]`, `PG['Pnma']`;
  - using `get()` method: `PG.get(number=62)`, `PG.get(HM_short='Pnma')`;

The `get()` method can use a combination of various keywords to unambiguously
localize the correct space group. The following keyword arguments can be used:

  - `n_c` - unique "number:setting" string identifying each group
  - `number` - index of the point group as given in ICT A (integer)\*
  - `setting` - arbitrary string declaring group setting
  - `HM` - Full underscore-delimited Hermann-Mauguin symbol\*
  - `HM_short` - Hermann-Mauguin symbol with underscores removed\*
  - `HM_simple` - Short Hermann-Mauguin symbol without setting information\*
  - `HM_numbered` - "number: Short Hermann-Mauguin symbol" string
  - `Hall` - Full Hall symbol\*
  - `standard` - True for groups in standard setting only

The keywords marked with "*" are called "accessors" and can be used to `get`
the group using bracket notation `get[accessor]`.
If `get` finds ambiguity, for example `SG['P2/m']` matches all settings
of space group #10: `P_1_2/m_1` (std), `P_1_1_2/m`, and `P_2/m_1_1`,
the user is warned and the first group in std setting, if possible, is returned.

The following table lists all possible values of selected space `Group`
keywords. Please mind that in the raw docstring all `\` should be ignored.

+---------+-----------------+-----------------+--------------+--------------+
| n_c     | Hall            | HM              | HM_short     | HM_simple    |
+---------+-----------------+-----------------+--------------+--------------+
| 1       | p_1             | P_1             | P1           | P1           |
+---------+-----------------+-----------------+--------------+--------------+
| 2       | -p_1            | P\_-1           | P-1          | P-1          |
+---------+-----------------+-----------------+--------------+--------------+
| 3:b     | p_2y            | P_1_2_1         | P121         | P2           |
+---------+-----------------+-----------------+--------------+--------------+
| 3:c     | p_2             | P_1_1_2         | P112         | P2           |
+---------+-----------------+-----------------+--------------+--------------+
| 3:a     | p_2x            | P_2_1_1         | P211         | P2           |
+---------+-----------------+-----------------+--------------+--------------+
| 4:b     | p_2yb           | P_1_21_1        | P1211        | P21          |
+---------+-----------------+-----------------+--------------+--------------+
| 4:c     | p_2c            | P_1_1_21        | P1121        | P21          |
+---------+-----------------+-----------------+--------------+--------------+
| 4:a     | p_2xa           | P_21_1_1        | P2111        | P21          |
+---------+-----------------+-----------------+--------------+--------------+
| 5:b1    | c_2y            | C_1_2_1         | C121         | C2           |
+---------+-----------------+-----------------+--------------+--------------+
| 5:b2    | a_2y            | A_1_2_1         | A121         | A2           |
+---------+-----------------+-----------------+--------------+--------------+
| 5:b3    | i_2y            | I_1_2_1         | I121         | I2           |
+---------+-----------------+-----------------+--------------+--------------+
| 5:c1    | a_2             | A_1_1_2         | A112         | A2           |
+---------+-----------------+-----------------+--------------+--------------+
| 5:c2    | b_2             | B_1_1_2         | B112         | B2           |
+---------+-----------------+-----------------+--------------+--------------+
| 5:c3    | i_2             | I_1_1_2         | I112         | I2           |
+---------+-----------------+-----------------+--------------+--------------+
| 5:a1    | b_2x            | B_2_1_1         | B211         | B2           |
+---------+-----------------+-----------------+--------------+--------------+
| 5:a2    | c_2x            | C_2_1_1         | C211         | C2           |
+---------+-----------------+-----------------+--------------+--------------+
| 5:a3    | i_2x            | I_2_1_1         | I211         | I2           |
+---------+-----------------+-----------------+--------------+--------------+
| 6:b     | p\_-2y          | P_1_m_1         | P1m1         | Pm           |
+---------+-----------------+-----------------+--------------+--------------+
| 6:c     | p\_-2           | P_1_1_m         | P11m         | Pm           |
+---------+-----------------+-----------------+--------------+--------------+
| 6:a     | p\_-2x          | P_m_1_1         | Pm11         | Pm           |
+---------+-----------------+-----------------+--------------+--------------+
| 7:b1    | p\_-2yc         | P_1_c_1         | P1c1         | Pc           |
+---------+-----------------+-----------------+--------------+--------------+
| 7:b2    | p\_-2yac        | P_1_n_1         | P1n1         | Pn           |
+---------+-----------------+-----------------+--------------+--------------+
| 7:b3    | p\_-2ya         | P_1_a_1         | P1a1         | Pa           |
+---------+-----------------+-----------------+--------------+--------------+
| 7:c1    | p\_-2a          | P_1_1_a         | P11a         | Pa           |
+---------+-----------------+-----------------+--------------+--------------+
| 7:c2    | p\_-2ab         | P_1_1_n         | P11n         | Pn           |
+---------+-----------------+-----------------+--------------+--------------+
| 7:c3    | p\_-2b          | P_1_1_b         | P11b         | Pb           |
+---------+-----------------+-----------------+--------------+--------------+
| 7:a1    | p\_-2xb         | P_b_1_1         | Pb11         | Pb           |
+---------+-----------------+-----------------+--------------+--------------+
| 7:a2    | p\_-2xbc        | P_n_1_1         | Pn11         | Pn           |
+---------+-----------------+-----------------+--------------+--------------+
| 7:a3    | p\_-2xc         | P_c_1_1         | Pc11         | Pc           |
+---------+-----------------+-----------------+--------------+--------------+
| 8:b1    | c\_-2y          | C_1_m_1         | C1m1         | Cm           |
+---------+-----------------+-----------------+--------------+--------------+
| 8:b2    | a\_-2y          | A_1_m_1         | A1m1         | Am           |
+---------+-----------------+-----------------+--------------+--------------+
| 8:b3    | i\_-2y          | I_1_m_1         | I1m1         | Im           |
+---------+-----------------+-----------------+--------------+--------------+
| 8:c1    | a\_-2           | A_1_1_m         | A11m         | Am           |
+---------+-----------------+-----------------+--------------+--------------+
| 8:c2    | b\_-2           | B_1_1_m         | B11m         | Bm           |
+---------+-----------------+-----------------+--------------+--------------+
| 8:c3    | i\_-2           | I_1_1_m         | I11m         | Im           |
+---------+-----------------+-----------------+--------------+--------------+
| 8:a1    | b\_-2x          | B_m_1_1         | Bm11         | Bm           |
+---------+-----------------+-----------------+--------------+--------------+
| 8:a2    | c\_-2x          | C_m_1_1         | Cm11         | Cm           |
+---------+-----------------+-----------------+--------------+--------------+
| 8:a3    | i\_-2x          | I_m_1_1         | Im11         | Im           |
+---------+-----------------+-----------------+--------------+--------------+
| 9:b1    | c\_-2yc         | C_1_c_1         | C1c1         | Cc           |
+---------+-----------------+-----------------+--------------+--------------+
| 9:b2    | a\_-2yab        | A_1_n_1         | A1n1         | An           |
+---------+-----------------+-----------------+--------------+--------------+
| 9:b3    | i\_-2ya         | I_1_a_1         | I1a1         | Ia           |
+---------+-----------------+-----------------+--------------+--------------+
| 9:-b1   | a\_-2ya         | A_1_a_1         | A1a1         | Aa           |
+---------+-----------------+-----------------+--------------+--------------+
| 9:-b2   | c\_-2yac        | C_1_n_1         | C1n1         | Cn           |
+---------+-----------------+-----------------+--------------+--------------+
| 9:-b3   | i\_-2yc         | I_1_c_1         | I1c1         | Ic           |
+---------+-----------------+-----------------+--------------+--------------+
| 9:c1    | a\_-2a          | A_1_1_a         | A11a         | Aa           |
+---------+-----------------+-----------------+--------------+--------------+
| 9:c2    | b\_-2ab         | B_1_1_n         | B11n         | Bn           |
+---------+-----------------+-----------------+--------------+--------------+
| 9:c3    | i\_-2b          | I_1_1_b         | I11b         | Ib           |
+---------+-----------------+-----------------+--------------+--------------+
| 9:-c1   | b\_-2b          | B_1_1_b         | B11b         | Bb           |
+---------+-----------------+-----------------+--------------+--------------+
| 9:-c2   | a\_-2ab         | A_1_1_n         | A11n         | An           |
+---------+-----------------+-----------------+--------------+--------------+
| 9:-c3   | i\_-2a          | I_1_1_a         | I11a         | Ia           |
+---------+-----------------+-----------------+--------------+--------------+
| 9:a1    | b\_-2xb         | B_b_1_1         | Bb11         | Bb           |
+---------+-----------------+-----------------+--------------+--------------+
| 9:a2    | c\_-2xac        | C_n_1_1         | Cn11         | Cn           |
+---------+-----------------+-----------------+--------------+--------------+
| 9:a3    | i\_-2xc         | I_c_1_1         | Ic11         | Ic           |
+---------+-----------------+-----------------+--------------+--------------+
| 9:-a1   | c\_-2xc         | C_c_1_1         | Cc11         | Cc           |
+---------+-----------------+-----------------+--------------+--------------+
| 9:-a2   | b\_-2xab        | B_n_1_1         | Bn11         | Bn           |
+---------+-----------------+-----------------+--------------+--------------+
| 9:-a3   | i\_-2xb         | I_b_1_1         | Ib11         | Ib           |
+---------+-----------------+-----------------+--------------+--------------+
| 10:b    | -p_2y           | P_1_2/m_1       | P12/m1       | P2/m         |
+---------+-----------------+-----------------+--------------+--------------+
| 10:c    | -p_2            | P_1_1_2/m       | P112/m       | P2/m         |
+---------+-----------------+-----------------+--------------+--------------+
| 10:a    | -p_2x           | P_2/m_1_1       | P2/m11       | P2/m         |
+---------+-----------------+-----------------+--------------+--------------+
| 11:b    | -p_2yb          | P_1_21/m_1      | P121/m1      | P21/m        |
+---------+-----------------+-----------------+--------------+--------------+
| 11:c    | -p_2c           | P_1_1_21/m      | P1121/m      | P21/m        |
+---------+-----------------+-----------------+--------------+--------------+
| 11:a    | -p_2xa          | P_21/m_1_1      | P21/m11      | P21/m        |
+---------+-----------------+-----------------+--------------+--------------+
| 12:b1   | -c_2y           | C_1_2/m_1       | C12/m1       | C2/m         |
+---------+-----------------+-----------------+--------------+--------------+
| 12:b2   | -a_2y           | A_1_2/m_1       | A12/m1       | A2/m         |
+---------+-----------------+-----------------+--------------+--------------+
| 12:b3   | -i_2y           | I_1_2/m_1       | I12/m1       | I2/m         |
+---------+-----------------+-----------------+--------------+--------------+
| 12:c1   | -a_2            | A_1_1_2/m       | A112/m       | A2/m         |
+---------+-----------------+-----------------+--------------+--------------+
| 12:c2   | -b_2            | B_1_1_2/m       | B112/m       | B2/m         |
+---------+-----------------+-----------------+--------------+--------------+
| 12:c3   | -i_2            | I_1_1_2/m       | I112/m       | I2/m         |
+---------+-----------------+-----------------+--------------+--------------+
| 12:a1   | -b_2x           | B_2/m_1_1       | B2/m11       | B2/m         |
+---------+-----------------+-----------------+--------------+--------------+
| 12:a2   | -c_2x           | C_2/m_1_1       | C2/m11       | C2/m         |
+---------+-----------------+-----------------+--------------+--------------+
| 12:a3   | -i_2x           | I_2/m_1_1       | I2/m11       | I2/m         |
+---------+-----------------+-----------------+--------------+--------------+
| 13:b1   | -p_2yc          | P_1_2/c_1       | P12/c1       | P2/c         |
+---------+-----------------+-----------------+--------------+--------------+
| 13:b2   | -p_2yac         | P_1_2/n_1       | P12/n1       | P2/n         |
+---------+-----------------+-----------------+--------------+--------------+
| 13:b3   | -p_2ya          | P_1_2/a_1       | P12/a1       | P2/a         |
+---------+-----------------+-----------------+--------------+--------------+
| 13:c1   | -p_2a           | P_1_1_2/a       | P112/a       | P2/a         |
+---------+-----------------+-----------------+--------------+--------------+
| 13:c2   | -p_2ab          | P_1_1_2/n       | P112/n       | P2/n         |
+---------+-----------------+-----------------+--------------+--------------+
| 13:c3   | -p_2b           | P_1_1_2/b       | P112/b       | P2/b         |
+---------+-----------------+-----------------+--------------+--------------+
| 13:a1   | -p_2xb          | P_2/b_1_1       | P2/b11       | P2/b         |
+---------+-----------------+-----------------+--------------+--------------+
| 13:a2   | -p_2xbc         | P_2/n_1_1       | P2/n11       | P2/n         |
+---------+-----------------+-----------------+--------------+--------------+
| 13:a3   | -p_2xc          | P_2/c_1_1       | P2/c11       | P2/c         |
+---------+-----------------+-----------------+--------------+--------------+
| 14:b1   | -p_2ybc         | P_1_21/c_1      | P121/c1      | P21/c        |
+---------+-----------------+-----------------+--------------+--------------+
| 14:b2   | -p_2yn          | P_1_21/n_1      | P121/n1      | P21/n        |
+---------+-----------------+-----------------+--------------+--------------+
| 14:b3   | -p_2yab         | P_1_21/a_1      | P121/a1      | P21/a        |
+---------+-----------------+-----------------+--------------+--------------+
| 14:c1   | -p_2ac          | P_1_1_21/a      | P1121/a      | P21/a        |
+---------+-----------------+-----------------+--------------+--------------+
| 14:c2   | -p_2n           | P_1_1_21/n      | P1121/n      | P21/n        |
+---------+-----------------+-----------------+--------------+--------------+
| 14:c3   | -p_2bc          | P_1_1_21/b      | P1121/b      | P21/b        |
+---------+-----------------+-----------------+--------------+--------------+
| 14:a1   | -p_2xab         | P_21/b_1_1      | P21/b11      | P21/b        |
+---------+-----------------+-----------------+--------------+--------------+
| 14:a2   | -p_2xn          | P_21/n_1_1      | P21/n11      | P21/n        |
+---------+-----------------+-----------------+--------------+--------------+
| 14:a3   | -p_2xac         | P_21/c_1_1      | P21/c11      | P21/c        |
+---------+-----------------+-----------------+--------------+--------------+
| 15:b1   | -c_2yc          | C_1_2/c_1       | C12/c1       | C2/c         |
+---------+-----------------+-----------------+--------------+--------------+
| 15:b2   | -a_2yab         | A_1_2/n_1       | A12/n1       | A2/n         |
+---------+-----------------+-----------------+--------------+--------------+
| 15:b3   | -i_2ya          | I_1_2/a_1       | I12/a1       | I2/a         |
+---------+-----------------+-----------------+--------------+--------------+
| 15:-b1  | -a_2ya          | A_1_2/a_1       | A12/a1       | A2/a         |
+---------+-----------------+-----------------+--------------+--------------+
| 15:-b2  | -c_2yac         | C_1_2/n_1       | C12/n1       | C2/n         |
+---------+-----------------+-----------------+--------------+--------------+
| 15:-b3  | -i_2yc          | I_1_2/c_1       | I12/c1       | I2/c         |
+---------+-----------------+-----------------+--------------+--------------+
| 15:c1   | -a_2a           | A_1_1_2/a       | A112/a       | A2/a         |
+---------+-----------------+-----------------+--------------+--------------+
| 15:c2   | -b_2ab          | B_1_1_2/n       | B112/n       | B2/n         |
+---------+-----------------+-----------------+--------------+--------------+
| 15:c3   | -i_2b           | I_1_1_2/b       | I112/b       | I2/b         |
+---------+-----------------+-----------------+--------------+--------------+
| 15:-c1  | -b_2b           | B_1_1_2/b       | B112/b       | B2/b         |
+---------+-----------------+-----------------+--------------+--------------+
| 15:-c2  | -a_2ab          | A_1_1_2/n       | A112/n       | A2/n         |
+---------+-----------------+-----------------+--------------+--------------+
| 15:-c3  | -i_2a           | I_1_1_2/a       | I112/a       | I2/a         |
+---------+-----------------+-----------------+--------------+--------------+
| 15:a1   | -b_2xb          | B_2/b_1_1       | B2/b11       | B2/b         |
+---------+-----------------+-----------------+--------------+--------------+
| 15:a2   | -c_2xac         | C_2/n_1_1       | C2/n11       | C2/n         |
+---------+-----------------+-----------------+--------------+--------------+
| 15:a3   | -i_2xc          | I_2/c_1_1       | I2/c11       | I2/c         |
+---------+-----------------+-----------------+--------------+--------------+
| 15:-a1  | -c_2xc          | C_2/c_1_1       | C2/c11       | C2/c         |
+---------+-----------------+-----------------+--------------+--------------+
| 15:-a2  | -b_2xab         | B_2/n_1_1       | B2/n11       | B2/n         |
+---------+-----------------+-----------------+--------------+--------------+
| 15:-a3  | -i_2xb          | I_2/b_1_1       | I2/b11       | I2/b         |
+---------+-----------------+-----------------+--------------+--------------+
| 16      | p_2_2           | P_2_2_2         | P222         | P222         |
+---------+-----------------+-----------------+--------------+--------------+
| 17      | p_2c_2          | P_2_2_21        | P2221        | P2221        |
+---------+-----------------+-----------------+--------------+--------------+
| 17:cab  | p_2a_2a         | P_21_2_2        | P2122        | P2122        |
+---------+-----------------+-----------------+--------------+--------------+
| 17:bca  | p_2_2b          | P_2_21_2        | P2212        | P2212        |
+---------+-----------------+-----------------+--------------+--------------+
| 18      | p_2_2ab         | P_21_21_2       | P21212       | P21212       |
+---------+-----------------+-----------------+--------------+--------------+
| 18:cab  | p_2bc_2         | P_2_21_21       | P22121       | P22121       |
+---------+-----------------+-----------------+--------------+--------------+
| 18:bca  | p_2ac_2ac       | P_21_2_21       | P21221       | P21221       |
+---------+-----------------+-----------------+--------------+--------------+
| 19      | p_2ac_2ab       | P_21_21_21      | P212121      | P212121      |
+---------+-----------------+-----------------+--------------+--------------+
| 20      | c_2c_2          | C_2_2_21        | C2221        | C2221        |
+---------+-----------------+-----------------+--------------+--------------+
| 20:cab  | a_2a_2a         | A_21_2_2        | A2122        | A2122        |
+---------+-----------------+-----------------+--------------+--------------+
| 20:bca  | b_2_2b          | B_2_21_2        | B2212        | B2212        |
+---------+-----------------+-----------------+--------------+--------------+
| 21      | c_2_2           | C_2_2_2         | C222         | C222         |
+---------+-----------------+-----------------+--------------+--------------+
| 21:cab  | a_2_2           | A_2_2_2         | A222         | A222         |
+---------+-----------------+-----------------+--------------+--------------+
| 21:bca  | b_2_2           | B_2_2_2         | B222         | B222         |
+---------+-----------------+-----------------+--------------+--------------+
| 22      | f_2_2           | F_2_2_2         | F222         | F222         |
+---------+-----------------+-----------------+--------------+--------------+
| 23      | i_2_2           | I_2_2_2         | I222         | I222         |
+---------+-----------------+-----------------+--------------+--------------+
| 24      | i_2b_2c         | I_21_21_21      | I212121      | I212121      |
+---------+-----------------+-----------------+--------------+--------------+
| 25      | p_2\_-2         | P_m_m_2         | Pmm2         | Pmm2         |
+---------+-----------------+-----------------+--------------+--------------+
| 25:cab  | p\_-2_2         | P_2_m_m         | P2mm         | P2mm         |
+---------+-----------------+-----------------+--------------+--------------+
| 25:bca  | p\_-2\_-2       | P_m_2_m         | Pm2m         | Pm2m         |
+---------+-----------------+-----------------+--------------+--------------+
| 26      | p_2c\_-2        | P_m_c_21        | Pmc21        | Pmc21        |
+---------+-----------------+-----------------+--------------+--------------+
| 26:ba-c | p_2c\_-2c       | P_c_m_21        | Pcm21        | Pcm21        |
+---------+-----------------+-----------------+--------------+--------------+
| 26:cab  | p\_-2a_2a       | P_21_m_a        | P21ma        | P21ma        |
+---------+-----------------+-----------------+--------------+--------------+
| 26:-cba | p\_-2_2a        | P_21_a_m        | P21am        | P21am        |
+---------+-----------------+-----------------+--------------+--------------+
| 26:bca  | p\_-2\_-2b      | P_b_21_m        | Pb21m        | Pb21m        |
+---------+-----------------+-----------------+--------------+--------------+
| 26:a-cb | p\_-2b\_-2      | P_m_21_b        | Pm21b        | Pm21b        |
+---------+-----------------+-----------------+--------------+--------------+
| 27      | p_2\_-2c        | P_c_c_2         | Pcc2         | Pcc2         |
+---------+-----------------+-----------------+--------------+--------------+
| 27:cab  | p\_-2a_2        | P_2_a_a         | P2aa         | P2aa         |
+---------+-----------------+-----------------+--------------+--------------+
| 27:bca  | p\_-2b\_-2b     | P_b_2_b         | Pb2b         | Pb2b         |
+---------+-----------------+-----------------+--------------+--------------+
| 28      | p_2\_-2a        | P_m_a_2         | Pma2         | Pma2         |
+---------+-----------------+-----------------+--------------+--------------+
| 28:ba-c | p_2\_-2b        | P_b_m_2         | Pbm2         | Pbm2         |
+---------+-----------------+-----------------+--------------+--------------+
| 28:cab  | p\_-2b_2        | P_2_m_b         | P2mb         | P2mb         |
+---------+-----------------+-----------------+--------------+--------------+
| 28:-cba | p\_-2c_2        | P_2_c_m         | P2cm         | P2cm         |
+---------+-----------------+-----------------+--------------+--------------+
| 28:bca  | p\_-2c\_-2c     | P_c_2_m         | Pc2m         | Pc2m         |
+---------+-----------------+-----------------+--------------+--------------+
| 28:a-cb | p\_-2a\_-2a     | P_m_2_a         | Pm2a         | Pm2a         |
+---------+-----------------+-----------------+--------------+--------------+
| 29      | p_2c\_-2ac      | P_c_a_21        | Pca21        | Pca21        |
+---------+-----------------+-----------------+--------------+--------------+
| 29:ba-c | p_2c\_-2b       | P_b_c_21        | Pbc21        | Pbc21        |
+---------+-----------------+-----------------+--------------+--------------+
| 29:cab  | p\_-2b_2a       | P_21_a_b        | P21ab        | P21ab        |
+---------+-----------------+-----------------+--------------+--------------+
| 29:-cba | p\_-2ac_2a      | P_21_c_a        | P21ca        | P21ca        |
+---------+-----------------+-----------------+--------------+--------------+
| 29:bca  | p\_-2bc\_-2c    | P_c_21_b        | Pc21b        | Pc21b        |
+---------+-----------------+-----------------+--------------+--------------+
| 29:a-cb | p\_-2a\_-2ab    | P_b_21_a        | Pb21a        | Pb21a        |
+---------+-----------------+-----------------+--------------+--------------+
| 30      | p_2\_-2bc       | P_n_c_2         | Pnc2         | Pnc2         |
+---------+-----------------+-----------------+--------------+--------------+
| 30:ba-c | p_2\_-2ac       | P_c_n_2         | Pcn2         | Pcn2         |
+---------+-----------------+-----------------+--------------+--------------+
| 30:cab  | p\_-2ac_2       | P_2_n_a         | P2na         | P2na         |
+---------+-----------------+-----------------+--------------+--------------+
| 30:-cba | p\_-2ab_2       | P_2_a_n         | P2an         | P2an         |
+---------+-----------------+-----------------+--------------+--------------+
| 30:bca  | p\_-2ab\_-2ab   | P_b_2_n         | Pb2n         | Pb2n         |
+---------+-----------------+-----------------+--------------+--------------+
| 30:a-cb | p\_-2bc\_-2bc   | P_n_2_b         | Pn2b         | Pn2b         |
+---------+-----------------+-----------------+--------------+--------------+
| 31      | p_2ac\_-2       | P_m_n_21        | Pmn21        | Pmn21        |
+---------+-----------------+-----------------+--------------+--------------+
| 31:ba-c | p_2bc\_-2bc     | P_n_m_21        | Pnm21        | Pnm21        |
+---------+-----------------+-----------------+--------------+--------------+
| 31:cab  | p\_-2ab_2ab     | P_21_m_n        | P21mn        | P21mn        |
+---------+-----------------+-----------------+--------------+--------------+
| 31:-cba | p\_-2_2ac       | P_21_n_m        | P21nm        | P21nm        |
+---------+-----------------+-----------------+--------------+--------------+
| 31:bca  | p\_-2\_-2bc     | P_n_21_m        | Pn21m        | Pn21m        |
+---------+-----------------+-----------------+--------------+--------------+
| 31:a-cb | p\_-2ab\_-2     | P_m_21_n        | Pm21n        | Pm21n        |
+---------+-----------------+-----------------+--------------+--------------+
| 32      | p_2\_-2ab       | P_b_a_2         | Pba2         | Pba2         |
+---------+-----------------+-----------------+--------------+--------------+
| 32:cab  | p\_-2bc_2       | P_2_c_b         | P2cb         | P2cb         |
+---------+-----------------+-----------------+--------------+--------------+
| 32:bca  | p\_-2ac\_-2ac   | P_c_2_a         | Pc2a         | Pc2a         |
+---------+-----------------+-----------------+--------------+--------------+
| 33      | p_2c\_-2n       | P_n_a_21        | Pna21        | Pna21        |
+---------+-----------------+-----------------+--------------+--------------+
| 33:ba-c | p_2c\_-2ab      | P_b_n_21        | Pbn21        | Pbn21        |
+---------+-----------------+-----------------+--------------+--------------+
| 33:cab  | p\_-2bc_2a      | P_21_n_b        | P21nb        | P21nb        |
+---------+-----------------+-----------------+--------------+--------------+
| 33:-cba | p\_-2n_2a       | P_21_c_n        | P21cn        | P21cn        |
+---------+-----------------+-----------------+--------------+--------------+
| 33:bca  | p\_-2n\_-2ac    | P_c_21_n        | Pc21n        | Pc21n        |
+---------+-----------------+-----------------+--------------+--------------+
| 33:a-cb | p\_-2ac\_-2n    | P_n_21_a        | Pn21a        | Pn21a        |
+---------+-----------------+-----------------+--------------+--------------+
| 34      | p_2\_-2n        | P_n_n_2         | Pnn2         | Pnn2         |
+---------+-----------------+-----------------+--------------+--------------+
| 34:cab  | p\_-2n_2        | P_2_n_n         | P2nn         | P2nn         |
+---------+-----------------+-----------------+--------------+--------------+
| 34:bca  | p\_-2n\_-2n     | P_n_2_n         | Pn2n         | Pn2n         |
+---------+-----------------+-----------------+--------------+--------------+
| 35      | c_2\_-2         | C_m_m_2         | Cmm2         | Cmm2         |
+---------+-----------------+-----------------+--------------+--------------+
| 35:cab  | a\_-2_2         | A_2_m_m         | A2mm         | A2mm         |
+---------+-----------------+-----------------+--------------+--------------+
| 35:bca  | b\_-2\_-2       | B_m_2_m         | Bm2m         | Bm2m         |
+---------+-----------------+-----------------+--------------+--------------+
| 36      | c_2c\_-2        | C_m_c_21        | Cmc21        | Cmc21        |
+---------+-----------------+-----------------+--------------+--------------+
| 36:ba-c | c_2c\_-2c       | C_c_m_21        | Ccm21        | Ccm21        |
+---------+-----------------+-----------------+--------------+--------------+
| 36:cab  | a\_-2a_2a       | A_21_m_a        | A21ma        | A21ma        |
+---------+-----------------+-----------------+--------------+--------------+
| 36:-cba | a\_-2_2a        | A_21_a_m        | A21am        | A21am        |
+---------+-----------------+-----------------+--------------+--------------+
| 36:bca  | b\_-2\_-2b      | B_b_21_m        | Bb21m        | Bb21m        |
+---------+-----------------+-----------------+--------------+--------------+
| 36:a-cb | b\_-2b\_-2      | B_m_21_b        | Bm21b        | Bm21b        |
+---------+-----------------+-----------------+--------------+--------------+
| 37      | c_2\_-2c        | C_c_c_2         | Ccc2         | Ccc2         |
+---------+-----------------+-----------------+--------------+--------------+
| 37:cab  | a\_-2a_2        | A_2_a_a         | A2aa         | A2aa         |
+---------+-----------------+-----------------+--------------+--------------+
| 37:bca  | b\_-2b\_-2b     | B_b_2_b         | Bb2b         | Bb2b         |
+---------+-----------------+-----------------+--------------+--------------+
| 38      | a_2\_-2         | A_m_m_2         | Amm2         | Amm2         |
+---------+-----------------+-----------------+--------------+--------------+
| 38:ba-c | b_2\_-2         | B_m_m_2         | Bmm2         | Bmm2         |
+---------+-----------------+-----------------+--------------+--------------+
| 38:cab  | b\_-2_2         | B_2_m_m         | B2mm         | B2mm         |
+---------+-----------------+-----------------+--------------+--------------+
| 38:-cba | c\_-2_2         | C_2_m_m         | C2mm         | C2mm         |
+---------+-----------------+-----------------+--------------+--------------+
| 38:bca  | c\_-2\_-2       | C_m_2_m         | Cm2m         | Cm2m         |
+---------+-----------------+-----------------+--------------+--------------+
| 38:a-cb | a\_-2\_-2       | A_m_2_m         | Am2m         | Am2m         |
+---------+-----------------+-----------------+--------------+--------------+
| 39      | a_2\_-2b        | A_e_m_2         | Aem2         | Aem2         |
+---------+-----------------+-----------------+--------------+--------------+
| 39:ba-c | b_2\_-2a        | B_m_a_2         | Bma2         | Bma2         |
+---------+-----------------+-----------------+--------------+--------------+
| 39:cab  | b\_-2a_2        | B_2_c_m         | B2cm         | B2cm         |
+---------+-----------------+-----------------+--------------+--------------+
| 39:-cba | c\_-2a_2        | C_2_m_b         | C2mb         | C2mb         |
+---------+-----------------+-----------------+--------------+--------------+
| 39:bca  | c\_-2a\_-2a     | C_m_2_a         | Cm2a         | Cm2a         |
+---------+-----------------+-----------------+--------------+--------------+
| 39:a-cb | a\_-2b\_-2b     | A_c_2_m         | Ac2m         | Ac2m         |
+---------+-----------------+-----------------+--------------+--------------+
| 40      | a_2\_-2a        | A_m_a_2         | Ama2         | Ama2         |
+---------+-----------------+-----------------+--------------+--------------+
| 40:ba-c | b_2\_-2b        | B_b_m_2         | Bbm2         | Bbm2         |
+---------+-----------------+-----------------+--------------+--------------+
| 40:cab  | b\_-2b_2        | B_2_m_b         | B2mb         | B2mb         |
+---------+-----------------+-----------------+--------------+--------------+
| 40:-cba | c\_-2c_2        | C_2_c_m         | C2cm         | C2cm         |
+---------+-----------------+-----------------+--------------+--------------+
| 40:bca  | c\_-2c\_-2c     | C_c_2_m         | Cc2m         | Cc2m         |
+---------+-----------------+-----------------+--------------+--------------+
| 40:a-cb | a\_-2a\_-2a     | A_m_2_a         | Am2a         | Am2a         |
+---------+-----------------+-----------------+--------------+--------------+
| 41      | a_2\_-2ab       | A_e_a_2         | Aea2         | Aea2         |
+---------+-----------------+-----------------+--------------+--------------+
| 41:ba-c | b_2\_-2ab       | B_b_a_2         | Bba2         | Bba2         |
+---------+-----------------+-----------------+--------------+--------------+
| 41:cab  | b\_-2ab_2       | B_2_c_b         | B2cb         | B2cb         |
+---------+-----------------+-----------------+--------------+--------------+
| 41:-cba | c\_-2ac_2       | C_2_c_b         | C2cb         | C2cb         |
+---------+-----------------+-----------------+--------------+--------------+
| 41:bca  | c\_-2ac\_-2ac   | C_c_2_a         | Cc2a         | Cc2a         |
+---------+-----------------+-----------------+--------------+--------------+
| 41:a-cb | a\_-2ab\_-2ab   | A_c_2_a         | Ac2a         | Ac2a         |
+---------+-----------------+-----------------+--------------+--------------+
| 42      | f_2\_-2         | F_m_m_2         | Fmm2         | Fmm2         |
+---------+-----------------+-----------------+--------------+--------------+
| 42:cab  | f\_-2_2         | F_2_m_m         | F2mm         | F2mm         |
+---------+-----------------+-----------------+--------------+--------------+
| 42:bca  | f\_-2\_-2       | F_m_2_m         | Fm2m         | Fm2m         |
+---------+-----------------+-----------------+--------------+--------------+
| 43      | f_2\_-2d        | F_d_d_2         | Fdd2         | Fdd2         |
+---------+-----------------+-----------------+--------------+--------------+
| 43:cab  | f\_-2d_2        | F_2_d_d         | F2dd         | F2dd         |
+---------+-----------------+-----------------+--------------+--------------+
| 43:bca  | f\_-2d\_-2d     | F_d_2_d         | Fd2d         | Fd2d         |
+---------+-----------------+-----------------+--------------+--------------+
| 44      | i_2\_-2         | I_m_m_2         | Imm2         | Imm2         |
+---------+-----------------+-----------------+--------------+--------------+
| 44:cab  | i\_-2_2         | I_2_m_m         | I2mm         | I2mm         |
+---------+-----------------+-----------------+--------------+--------------+
| 44:bca  | i\_-2\_-2       | I_m_2_m         | Im2m         | Im2m         |
+---------+-----------------+-----------------+--------------+--------------+
| 45      | i_2\_-2c        | I_b_a_2         | Iba2         | Iba2         |
+---------+-----------------+-----------------+--------------+--------------+
| 45:cab  | i\_-2a_2        | I_2_c_b         | I2cb         | I2cb         |
+---------+-----------------+-----------------+--------------+--------------+
| 45:bca  | i\_-2b\_-2b     | I_c_2_a         | Ic2a         | Ic2a         |
+---------+-----------------+-----------------+--------------+--------------+
| 46      | i_2\_-2a        | I_m_a_2         | Ima2         | Ima2         |
+---------+-----------------+-----------------+--------------+--------------+
| 46:ba-c | i_2\_-2b        | I_b_m_2         | Ibm2         | Ibm2         |
+---------+-----------------+-----------------+--------------+--------------+
| 46:cab  | i\_-2b_2        | I_2_m_b         | I2mb         | I2mb         |
+---------+-----------------+-----------------+--------------+--------------+
| 46:-cba | i\_-2c_2        | I_2_c_m         | I2cm         | I2cm         |
+---------+-----------------+-----------------+--------------+--------------+
| 46:bca  | i\_-2c\_-2c     | I_c_2_m         | Ic2m         | Ic2m         |
+---------+-----------------+-----------------+--------------+--------------+
| 46:a-cb | i\_-2a\_-2a     | I_m_2_a         | Im2a         | Im2a         |
+---------+-----------------+-----------------+--------------+--------------+
| 47      | -p_2_2          | P_m_m_m         | Pmmm         | Pmmm         |
+---------+-----------------+-----------------+--------------+--------------+
| 48:1    | p_2_2\_-1n      | P_n_n_n:1       | Pnnn:1       | Pnnn         |
+---------+-----------------+-----------------+--------------+--------------+
| 48:2    | -p_2ab_2bc      | P_n_n_n:2       | Pnnn:2       | Pnnn         |
+---------+-----------------+-----------------+--------------+--------------+
| 49      | -p_2_2c         | P_c_c_m         | Pccm         | Pccm         |
+---------+-----------------+-----------------+--------------+--------------+
| 49:cab  | -p_2a_2         | P_m_a_a         | Pmaa         | Pmaa         |
+---------+-----------------+-----------------+--------------+--------------+
| 49:bca  | -p_2b_2b        | P_b_m_b         | Pbmb         | Pbmb         |
+---------+-----------------+-----------------+--------------+--------------+
| 50:1    | p_2_2\_-1ab     | P_b_a_n:1       | Pban:1       | Pban         |
+---------+-----------------+-----------------+--------------+--------------+
| 50:2    | -p_2ab_2b       | P_b_a_n:2       | Pban:2       | Pban         |
+---------+-----------------+-----------------+--------------+--------------+
| 50:1cab | p_2_2\_-1bc     | P_n_c_b:1       | Pncb:1       | Pncb         |
+---------+-----------------+-----------------+--------------+--------------+
| 50:2cab | -p_2b_2bc       | P_n_c_b:2       | Pncb:2       | Pncb         |
+---------+-----------------+-----------------+--------------+--------------+
| 50:1bca | p_2_2\_-1ac     | P_c_n_a:1       | Pcna:1       | Pcna         |
+---------+-----------------+-----------------+--------------+--------------+
| 50:2bca | -p_2a_2c        | P_c_n_a:2       | Pcna:2       | Pcna         |
+---------+-----------------+-----------------+--------------+--------------+
| 51      | -p_2a_2a        | P_m_m_a         | Pmma         | Pmma         |
+---------+-----------------+-----------------+--------------+--------------+
| 51:ba-c | -p_2b_2         | P_m_m_b         | Pmmb         | Pmmb         |
+---------+-----------------+-----------------+--------------+--------------+
| 51:cab  | -p_2_2b         | P_b_m_m         | Pbmm         | Pbmm         |
+---------+-----------------+-----------------+--------------+--------------+
| 51:-cba | -p_2c_2c        | P_c_m_m         | Pcmm         | Pcmm         |
+---------+-----------------+-----------------+--------------+--------------+
| 51:bca  | -p_2c_2         | P_m_c_m         | Pmcm         | Pmcm         |
+---------+-----------------+-----------------+--------------+--------------+
| 51:a-cb | -p_2_2a         | P_m_a_m         | Pmam         | Pmam         |
+---------+-----------------+-----------------+--------------+--------------+
| 52      | -p_2a_2bc       | P_n_n_a         | Pnna         | Pnna         |
+---------+-----------------+-----------------+--------------+--------------+
| 52:ba-c | -p_2b_2n        | P_n_n_b         | Pnnb         | Pnnb         |
+---------+-----------------+-----------------+--------------+--------------+
| 52:cab  | -p_2n_2b        | P_b_n_n         | Pbnn         | Pbnn         |
+---------+-----------------+-----------------+--------------+--------------+
| 52:-cba | -p_2ab_2c       | P_c_n_n         | Pcnn         | Pcnn         |
+---------+-----------------+-----------------+--------------+--------------+
| 52:bca  | -p_2ab_2n       | P_n_c_n         | Pncn         | Pncn         |
+---------+-----------------+-----------------+--------------+--------------+
| 52:a-cb | -p_2n_2bc       | P_n_a_n         | Pnan         | Pnan         |
+---------+-----------------+-----------------+--------------+--------------+
| 53      | -p_2ac_2        | P_m_n_a         | Pmna         | Pmna         |
+---------+-----------------+-----------------+--------------+--------------+
| 53:ba-c | -p_2bc_2bc      | P_n_m_b         | Pnmb         | Pnmb         |
+---------+-----------------+-----------------+--------------+--------------+
| 53:cab  | -p_2ab_2ab      | P_b_m_n         | Pbmn         | Pbmn         |
+---------+-----------------+-----------------+--------------+--------------+
| 53:-cba | -p_2_2ac        | P_c_n_m         | Pcnm         | Pcnm         |
+---------+-----------------+-----------------+--------------+--------------+
| 53:bca  | -p_2_2bc        | P_n_c_m         | Pncm         | Pncm         |
+---------+-----------------+-----------------+--------------+--------------+
| 53:a-cb | -p_2ab_2        | P_m_a_n         | Pman         | Pman         |
+---------+-----------------+-----------------+--------------+--------------+
| 54      | -p_2a_2ac       | P_c_c_a         | Pcca         | Pcca         |
+---------+-----------------+-----------------+--------------+--------------+
| 54:ba-c | -p_2b_2c        | P_c_c_b         | Pccb         | Pccb         |
+---------+-----------------+-----------------+--------------+--------------+
| 54:cab  | -p_2a_2b        | P_b_a_a         | Pbaa         | Pbaa         |
+---------+-----------------+-----------------+--------------+--------------+
| 54:-cba | -p_2ac_2c       | P_c_a_a         | Pcaa         | Pcaa         |
+---------+-----------------+-----------------+--------------+--------------+
| 54:bca  | -p_2bc_2b       | P_b_c_b         | Pbcb         | Pbcb         |
+---------+-----------------+-----------------+--------------+--------------+
| 54:a-cb | -p_2b_2ab       | P_b_a_b         | Pbab         | Pbab         |
+---------+-----------------+-----------------+--------------+--------------+
| 55      | -p_2_2ab        | P_b_a_m         | Pbam         | Pbam         |
+---------+-----------------+-----------------+--------------+--------------+
| 55:cab  | -p_2bc_2        | P_m_c_b         | Pmcb         | Pmcb         |
+---------+-----------------+-----------------+--------------+--------------+
| 55:bca  | -p_2ac_2ac      | P_c_m_a         | Pcma         | Pcma         |
+---------+-----------------+-----------------+--------------+--------------+
| 56      | -p_2ab_2ac      | P_c_c_n         | Pccn         | Pccn         |
+---------+-----------------+-----------------+--------------+--------------+
| 56:cab  | -p_2ac_2bc      | P_n_a_a         | Pnaa         | Pnaa         |
+---------+-----------------+-----------------+--------------+--------------+
| 56:bca  | -p_2bc_2ab      | P_b_n_b         | Pbnb         | Pbnb         |
+---------+-----------------+-----------------+--------------+--------------+
| 57      | -p_2c_2b        | P_b_c_m         | Pbcm         | Pbcm         |
+---------+-----------------+-----------------+--------------+--------------+
| 57:ba-c | -p_2c_2ac       | P_c_a_m         | Pcam         | Pcam         |
+---------+-----------------+-----------------+--------------+--------------+
| 57:cab  | -p_2ac_2a       | P_m_c_a         | Pmca         | Pmca         |
+---------+-----------------+-----------------+--------------+--------------+
| 57:-cba | -p_2b_2a        | P_m_a_b         | Pmab         | Pmab         |
+---------+-----------------+-----------------+--------------+--------------+
| 57:bca  | -p_2a_2ab       | P_b_m_a         | Pbma         | Pbma         |
+---------+-----------------+-----------------+--------------+--------------+
| 57:a-cb | -p_2bc_2c       | P_c_m_b         | Pcmb         | Pcmb         |
+---------+-----------------+-----------------+--------------+--------------+
| 58      | -p_2_2n         | P_n_n_m         | Pnnm         | Pnnm         |
+---------+-----------------+-----------------+--------------+--------------+
| 58:cab  | -p_2n_2         | P_m_n_n         | Pmnn         | Pmnn         |
+---------+-----------------+-----------------+--------------+--------------+
| 58:bca  | -p_2n_2n        | P_n_m_n         | Pnmn         | Pnmn         |
+---------+-----------------+-----------------+--------------+--------------+
| 59:1    | p_2_2ab\_-1ab   | P_m_m_n:1       | Pmmn:1       | Pmmn         |
+---------+-----------------+-----------------+--------------+--------------+
| 59:2    | -p_2ab_2a       | P_m_m_n:2       | Pmmn:2       | Pmmn         |
+---------+-----------------+-----------------+--------------+--------------+
| 59:1cab | p_2bc_2\_-1bc   | P_n_m_m:1       | Pnmm:1       | Pnmm         |
+---------+-----------------+-----------------+--------------+--------------+
| 59:2cab | -p_2c_2bc       | P_n_m_m:2       | Pnmm:2       | Pnmm         |
+---------+-----------------+-----------------+--------------+--------------+
| 59:1bca | p_2ac_2ac\_-1ac | P_m_n_m:1       | Pmnm:1       | Pmnm         |
+---------+-----------------+-----------------+--------------+--------------+
| 59:2bca | -p_2c_2a        | P_m_n_m:2       | Pmnm:2       | Pmnm         |
+---------+-----------------+-----------------+--------------+--------------+
| 60      | -p_2n_2ab       | P_b_c_n         | Pbcn         | Pbcn         |
+---------+-----------------+-----------------+--------------+--------------+
| 60:ba-c | -p_2n_2c        | P_c_a_n         | Pcan         | Pcan         |
+---------+-----------------+-----------------+--------------+--------------+
| 60:cab  | -p_2a_2n        | P_n_c_a         | Pnca         | Pnca         |
+---------+-----------------+-----------------+--------------+--------------+
| 60:-cba | -p_2bc_2n       | P_n_a_b         | Pnab         | Pnab         |
+---------+-----------------+-----------------+--------------+--------------+
| 60:bca  | -p_2ac_2b       | P_b_n_a         | Pbna         | Pbna         |
+---------+-----------------+-----------------+--------------+--------------+
| 60:a-cb | -p_2b_2ac       | P_c_n_b         | Pcnb         | Pcnb         |
+---------+-----------------+-----------------+--------------+--------------+
| 61      | -p_2ac_2ab      | P_b_c_a         | Pbca         | Pbca         |
+---------+-----------------+-----------------+--------------+--------------+
| 61:ba-c | -p_2bc_2ac      | P_c_a_b         | Pcab         | Pcab         |
+---------+-----------------+-----------------+--------------+--------------+
| 62      | -p_2ac_2n       | P_n_m_a         | Pnma         | Pnma         |
+---------+-----------------+-----------------+--------------+--------------+
| 62:ba-c | -p_2bc_2a       | P_m_n_b         | Pmnb         | Pmnb         |
+---------+-----------------+-----------------+--------------+--------------+
| 62:cab  | -p_2c_2ab       | P_b_n_m         | Pbnm         | Pbnm         |
+---------+-----------------+-----------------+--------------+--------------+
| 62:-cba | -p_2n_2ac       | P_c_m_n         | Pcmn         | Pcmn         |
+---------+-----------------+-----------------+--------------+--------------+
| 62:bca  | -p_2n_2a        | P_m_c_n         | Pmcn         | Pmcn         |
+---------+-----------------+-----------------+--------------+--------------+
| 62:a-cb | -p_2c_2n        | P_n_a_m         | Pnam         | Pnam         |
+---------+-----------------+-----------------+--------------+--------------+
| 63      | -c_2c_2         | C_m_c_m         | Cmcm         | Cmcm         |
+---------+-----------------+-----------------+--------------+--------------+
| 63:ba-c | -c_2c_2c        | C_c_m_m         | Ccmm         | Ccmm         |
+---------+-----------------+-----------------+--------------+--------------+
| 63:cab  | -a_2a_2a        | A_m_m_a         | Amma         | Amma         |
+---------+-----------------+-----------------+--------------+--------------+
| 63:-cba | -a_2_2a         | A_m_a_m         | Amam         | Amam         |
+---------+-----------------+-----------------+--------------+--------------+
| 63:bca  | -b_2_2b         | B_b_m_m         | Bbmm         | Bbmm         |
+---------+-----------------+-----------------+--------------+--------------+
| 63:a-cb | -b_2b_2         | B_m_m_b         | Bmmb         | Bmmb         |
+---------+-----------------+-----------------+--------------+--------------+
| 64      | -c_2ac_2        | C_m_c_e         | Cmce         | Cmce         |
+---------+-----------------+-----------------+--------------+--------------+
| 64:ba-c | -c_2ac_2ac      | C_c_m_b         | Ccmb         | Ccmb         |
+---------+-----------------+-----------------+--------------+--------------+
| 64:cab  | -a_2ab_2ab      | A_b_m_a         | Abma         | Abma         |
+---------+-----------------+-----------------+--------------+--------------+
| 64:-cba | -a_2_2ab        | A_c_a_m         | Acam         | Acam         |
+---------+-----------------+-----------------+--------------+--------------+
| 64:bca  | -b_2_2ab        | B_b_c_m         | Bbcm         | Bbcm         |
+---------+-----------------+-----------------+--------------+--------------+
| 64:a-cb | -b_2ab_2        | B_m_a_b         | Bmab         | Bmab         |
+---------+-----------------+-----------------+--------------+--------------+
| 65      | -c_2_2          | C_m_m_m         | Cmmm         | Cmmm         |
+---------+-----------------+-----------------+--------------+--------------+
| 65:cab  | -a_2_2          | A_m_m_m         | Ammm         | Ammm         |
+---------+-----------------+-----------------+--------------+--------------+
| 65:bca  | -b_2_2          | B_m_m_m         | Bmmm         | Bmmm         |
+---------+-----------------+-----------------+--------------+--------------+
| 66      | -c_2_2c         | C_c_c_m         | Cccm         | Cccm         |
+---------+-----------------+-----------------+--------------+--------------+
| 66:cab  | -a_2a_2         | A_m_a_a         | Amaa         | Amaa         |
+---------+-----------------+-----------------+--------------+--------------+
| 66:bca  | -b_2b_2b        | B_b_m_b         | Bbmb         | Bbmb         |
+---------+-----------------+-----------------+--------------+--------------+
| 67      | -c_2a_2         | C_m_m_e         | Cmme         | Cmme         |
+---------+-----------------+-----------------+--------------+--------------+
| 67:ba-c | -c_2a_2a        | C_m_m_b         | Cmmb         | Cmmb         |
+---------+-----------------+-----------------+--------------+--------------+
| 67:cab  | -a_2b_2b        | A_b_m_m         | Abmm         | Abmm         |
+---------+-----------------+-----------------+--------------+--------------+
| 67:-cba | -a_2_2b         | A_c_m_m         | Acmm         | Acmm         |
+---------+-----------------+-----------------+--------------+--------------+
| 67:bca  | -b_2_2a         | B_m_c_m         | Bmcm         | Bmcm         |
+---------+-----------------+-----------------+--------------+--------------+
| 67:a-cb | -b_2a_2         | B_m_a_m         | Bmam         | Bmam         |
+---------+-----------------+-----------------+--------------+--------------+
| 68:1    | c_2_2\_-1ac     | C_c_c_e:1       | Ccce:1       | Ccce         |
+---------+-----------------+-----------------+--------------+--------------+
| 68:2    | -c_2a_2ac       | C_c_c_e:2       | Ccce:2       | Ccce         |
+---------+-----------------+-----------------+--------------+--------------+
| 68:1ba- | c_2_2\_-1ac     | C_c_c_b:1       | Cccb:1       | Cccb         |
+---------+-----------------+-----------------+--------------+--------------+
| 68:2ba- | -c_2a_2c        | C_c_c_b:2       | Cccb:2       | Cccb         |
+---------+-----------------+-----------------+--------------+--------------+
| 68:1cab | a_2_2\_-1ab     | A_b_a_a:1       | Abaa:1       | Abaa         |
+---------+-----------------+-----------------+--------------+--------------+
| 68:2cab | -a_2a_2b        | A_b_a_a:2       | Abaa:2       | Abaa         |
+---------+-----------------+-----------------+--------------+--------------+
| 68:1-cb | a_2_2\_-1ab     | A_c_a_a:1       | Acaa:1       | Acaa         |
+---------+-----------------+-----------------+--------------+--------------+
| 68:2-cb | -a_2ab_2b       | A_c_a_a:2       | Acaa:2       | Acaa         |
+---------+-----------------+-----------------+--------------+--------------+
| 68:1bca | b_2_2\_-1ab     | B_b_c_b:1       | Bbcb:1       | Bbcb         |
+---------+-----------------+-----------------+--------------+--------------+
| 68:2bca | -b_2ab_2b       | B_b_c_b:2       | Bbcb:2       | Bbcb         |
+---------+-----------------+-----------------+--------------+--------------+
| 68:1a-c | b_2_2\_-1ab     | B_b_a_b:1       | Bbab:1       | Bbab         |
+---------+-----------------+-----------------+--------------+--------------+
| 68:2a-c | -b_2b_2ab       | B_b_a_b:2       | Bbab:2       | Bbab         |
+---------+-----------------+-----------------+--------------+--------------+
| 69      | -f_2_2          | F_m_m_m         | Fmmm         | Fmmm         |
+---------+-----------------+-----------------+--------------+--------------+
| 70:1    | f_2_2\_-1d      | F_d_d_d:1       | Fddd:1       | Fddd         |
+---------+-----------------+-----------------+--------------+--------------+
| 70:2    | -f_2uv_2vw      | F_d_d_d:2       | Fddd:2       | Fddd         |
+---------+-----------------+-----------------+--------------+--------------+
| 71      | -i_2_2          | I_m_m_m         | Immm         | Immm         |
+---------+-----------------+-----------------+--------------+--------------+
| 72      | -i_2_2c         | I_b_a_m         | Ibam         | Ibam         |
+---------+-----------------+-----------------+--------------+--------------+
| 72:cab  | -i_2a_2         | I_m_c_b         | Imcb         | Imcb         |
+---------+-----------------+-----------------+--------------+--------------+
| 72:bca  | -i_2b_2b        | I_c_m_a         | Icma         | Icma         |
+---------+-----------------+-----------------+--------------+--------------+
| 73      | -i_2b_2c        | I_b_c_a         | Ibca         | Ibca         |
+---------+-----------------+-----------------+--------------+--------------+
| 73:ba-c | -i_2a_2b        | I_c_a_b         | Icab         | Icab         |
+---------+-----------------+-----------------+--------------+--------------+
| 74      | -i_2b_2         | I_m_m_a         | Imma         | Imma         |
+---------+-----------------+-----------------+--------------+--------------+
| 74:ba-c | -i_2a_2a        | I_m_m_b         | Immb         | Immb         |
+---------+-----------------+-----------------+--------------+--------------+
| 74:cab  | -i_2c_2c        | I_b_m_m         | Ibmm         | Ibmm         |
+---------+-----------------+-----------------+--------------+--------------+
| 74:-cba | -i_2_2b         | I_c_m_m         | Icmm         | Icmm         |
+---------+-----------------+-----------------+--------------+--------------+
| 74:bca  | -i_2_2a         | I_m_c_m         | Imcm         | Imcm         |
+---------+-----------------+-----------------+--------------+--------------+
| 74:a-cb | -i_2c_2         | I_m_a_m         | Imam         | Imam         |
+---------+-----------------+-----------------+--------------+--------------+
| 75      | p_4             | P_4             | P4           | P4           |
+---------+-----------------+-----------------+--------------+--------------+
| 76      | p_4w            | P_41            | P41          | P41          |
+---------+-----------------+-----------------+--------------+--------------+
| 77      | p_4c            | P_42            | P42          | P42          |
+---------+-----------------+-----------------+--------------+--------------+
| 78      | p_4cw           | P_43            | P43          | P43          |
+---------+-----------------+-----------------+--------------+--------------+
| 79      | i_4             | I_4             | I4           | I4           |
+---------+-----------------+-----------------+--------------+--------------+
| 80      | i_4bw           | I_41            | I41          | I41          |
+---------+-----------------+-----------------+--------------+--------------+
| 81      | p\_-4           | P\_-4           | P-4          | P-4          |
+---------+-----------------+-----------------+--------------+--------------+
| 82      | i\_-4           | I\_-4           | I-4          | I-4          |
+---------+-----------------+-----------------+--------------+--------------+
| 83      | -p_4            | P_4/m           | P4/m         | P4/m         |
+---------+-----------------+-----------------+--------------+--------------+
| 84      | -p_4c           | P_42/m          | P42/m        | P42/m        |
+---------+-----------------+-----------------+--------------+--------------+
| 85:1    | p_4ab\_-1ab     | P_4/n:1         | P4/n:1       | P4/n         |
+---------+-----------------+-----------------+--------------+--------------+
| 85:2    | -p_4a           | P_4/n:2         | P4/n:2       | P4/n         |
+---------+-----------------+-----------------+--------------+--------------+
| 86:1    | p_4n\_-1n       | P_42/n:1        | P42/n:1      | P42/n        |
+---------+-----------------+-----------------+--------------+--------------+
| 86:2    | -p_4bc          | P_42/n:2        | P42/n:2      | P42/n        |
+---------+-----------------+-----------------+--------------+--------------+
| 87      | -i_4            | I_4/m           | I4/m         | I4/m         |
+---------+-----------------+-----------------+--------------+--------------+
| 88:1    | i_4bw\_-1bw     | I_41/a:1        | I41/a:1      | I41/a        |
+---------+-----------------+-----------------+--------------+--------------+
| 88:2    | -i_4ad          | I_41/a:2        | I41/a:2      | I41/a        |
+---------+-----------------+-----------------+--------------+--------------+
| 89      | p_4_2           | P_4_2_2         | P422         | P422         |
+---------+-----------------+-----------------+--------------+--------------+
| 90      | p_4ab_2ab       | P_4_21_2        | P4212        | P4212        |
+---------+-----------------+-----------------+--------------+--------------+
| 91      | p_4w_2c         | P_41_2_2        | P4122        | P4122        |
+---------+-----------------+-----------------+--------------+--------------+
| 92      | p_4abw_2nw      | P_41_21_2       | P41212       | P41212       |
+---------+-----------------+-----------------+--------------+--------------+
| 93      | p_4c_2          | P_42_2_2        | P4222        | P4222        |
+---------+-----------------+-----------------+--------------+--------------+
| 94      | p_4n_2n         | P_42_21_2       | P42212       | P42212       |
+---------+-----------------+-----------------+--------------+--------------+
| 95      | p_4cw_2c        | P_43_2_2        | P4322        | P4322        |
+---------+-----------------+-----------------+--------------+--------------+
| 96      | p_4nw_2abw      | P_43_21_2       | P43212       | P43212       |
+---------+-----------------+-----------------+--------------+--------------+
| 97      | i_4_2           | I_4_2_2         | I422         | I422         |
+---------+-----------------+-----------------+--------------+--------------+
| 98      | i_4bw_2bw       | I_41_2_2        | I4122        | I4122        |
+---------+-----------------+-----------------+--------------+--------------+
| 99      | p_4\_-2         | P_4_m_m         | P4mm         | P4mm         |
+---------+-----------------+-----------------+--------------+--------------+
| 100     | p_4\_-2ab       | P_4_b_m         | P4bm         | P4bm         |
+---------+-----------------+-----------------+--------------+--------------+
| 101     | p_4c\_-2c       | P_42_c_m        | P42cm        | P42cm        |
+---------+-----------------+-----------------+--------------+--------------+
| 102     | p_4n\_-2n       | P_42_n_m        | P42nm        | P42nm        |
+---------+-----------------+-----------------+--------------+--------------+
| 103     | p_4\_-2c        | P_4_c_c         | P4cc         | P4cc         |
+---------+-----------------+-----------------+--------------+--------------+
| 104     | p_4\_-2n        | P_4_n_c         | P4nc         | P4nc         |
+---------+-----------------+-----------------+--------------+--------------+
| 105     | p_4c\_-2        | P_42_m_c        | P42mc        | P42mc        |
+---------+-----------------+-----------------+--------------+--------------+
| 106     | p_4c\_-2ab      | P_42_b_c        | P42bc        | P42bc        |
+---------+-----------------+-----------------+--------------+--------------+
| 107     | i_4\_-2         | I_4_m_m         | I4mm         | I4mm         |
+---------+-----------------+-----------------+--------------+--------------+
| 108     | i_4\_-2c        | I_4_c_m         | I4cm         | I4cm         |
+---------+-----------------+-----------------+--------------+--------------+
| 109     | i_4bw\_-2       | I_41_m_d        | I41md        | I41md        |
+---------+-----------------+-----------------+--------------+--------------+
| 110     | i_4bw\_-2c      | I_41_c_d        | I41cd        | I41cd        |
+---------+-----------------+-----------------+--------------+--------------+
| 111     | p\_-4_2         | P\_-4_2_m       | P-42m        | P-42m        |
+---------+-----------------+-----------------+--------------+--------------+
| 112     | p\_-4_2c        | P\_-4_2_c       | P-42c        | P-42c        |
+---------+-----------------+-----------------+--------------+--------------+
| 113     | p\_-4_2ab       | P\_-4_21_m      | P-421m       | P-421m       |
+---------+-----------------+-----------------+--------------+--------------+
| 114     | p\_-4_2n        | P\_-4_21_c      | P-421c       | P-421c       |
+---------+-----------------+-----------------+--------------+--------------+
| 115     | p\_-4\_-2       | P\_-4_m_2       | P-4m2        | P-4m2        |
+---------+-----------------+-----------------+--------------+--------------+
| 116     | p\_-4\_-2c      | P\_-4_c_2       | P-4c2        | P-4c2        |
+---------+-----------------+-----------------+--------------+--------------+
| 117     | p\_-4\_-2ab     | P\_-4_b_2       | P-4b2        | P-4b2        |
+---------+-----------------+-----------------+--------------+--------------+
| 118     | p\_-4\_-2n      | P\_-4_n_2       | P-4n2        | P-4n2        |
+---------+-----------------+-----------------+--------------+--------------+
| 119     | i\_-4\_-2       | I\_-4_m_2       | I-4m2        | I-4m2        |
+---------+-----------------+-----------------+--------------+--------------+
| 120     | i\_-4\_-2c      | I\_-4_c_2       | I-4c2        | I-4c2        |
+---------+-----------------+-----------------+--------------+--------------+
| 121     | i\_-4_2         | I\_-4_2_m       | I-42m        | I-42m        |
+---------+-----------------+-----------------+--------------+--------------+
| 122     | i\_-4_2bw       | I\_-4_2_d       | I-42d        | I-42d        |
+---------+-----------------+-----------------+--------------+--------------+
| 123     | -p_4_2          | P_4/m_m_m       | P4/mmm       | P4/mmm       |
+---------+-----------------+-----------------+--------------+--------------+
| 124     | -p_4_2c         | P_4/m_c_c       | P4/mcc       | P4/mcc       |
+---------+-----------------+-----------------+--------------+--------------+
| 125:1   | p_4_2\_-1ab     | P_4/n_b_m:1     | P4/nbm:1     | P4/nbm       |
+---------+-----------------+-----------------+--------------+--------------+
| 125:2   | -p_4a_2b        | P_4/n_b_m:2     | P4/nbm:2     | P4/nbm       |
+---------+-----------------+-----------------+--------------+--------------+
| 126:1   | p_4_2\_-1n      | P_4/n_n_c:1     | P4/nnc:1     | P4/nnc       |
+---------+-----------------+-----------------+--------------+--------------+
| 126:2   | -p_4a_2bc       | P_4/n_n_c:2     | P4/nnc:2     | P4/nnc       |
+---------+-----------------+-----------------+--------------+--------------+
| 127     | -p_4_2ab        | P_4/m_b_m       | P4/mbm       | P4/mbm       |
+---------+-----------------+-----------------+--------------+--------------+
| 128     | -p_4_2n         | P_4/m_n_c       | P4/mnc       | P4/mnc       |
+---------+-----------------+-----------------+--------------+--------------+
| 129:1   | p_4ab_2ab\_-1ab | P_4/n_m_m:1     | P4/nmm:1     | P4/nmm       |
+---------+-----------------+-----------------+--------------+--------------+
| 129:2   | -p_4a_2a        | P_4/n_m_m:2     | P4/nmm:2     | P4/nmm       |
+---------+-----------------+-----------------+--------------+--------------+
| 130:1   | p_4ab_2n\_-1ab  | P_4/n_c_c:1     | P4/ncc:1     | P4/ncc       |
+---------+-----------------+-----------------+--------------+--------------+
| 130:2   | -p_4a_2ac       | P_4/n_c_c:2     | P4/ncc:2     | P4/ncc       |
+---------+-----------------+-----------------+--------------+--------------+
| 131     | -p_4c_2         | P_42/m_m_c      | P42/mmc      | P42/mmc      |
+---------+-----------------+-----------------+--------------+--------------+
| 132     | -p_4c_2c        | P_42/m_c_m      | P42/mcm      | P42/mcm      |
+---------+-----------------+-----------------+--------------+--------------+
| 133:1   | p_4n_2c\_-1n    | P_42/n_b_c:1    | P42/nbc:1    | P42/nbc      |
+---------+-----------------+-----------------+--------------+--------------+
| 133:2   | -p_4ac_2b       | P_42/n_b_c:2    | P42/nbc:2    | P42/nbc      |
+---------+-----------------+-----------------+--------------+--------------+
| 134:1   | p_4n_2\_-1n     | P_42/n_n_m:1    | P42/nnm:1    | P42/nnm      |
+---------+-----------------+-----------------+--------------+--------------+
| 134:2   | -p_4ac_2bc      | P_42/n_n_m:2    | P42/nnm:2    | P42/nnm      |
+---------+-----------------+-----------------+--------------+--------------+
| 135     | -p_4c_2ab       | P_42/m_b_c      | P42/mbc      | P42/mbc      |
+---------+-----------------+-----------------+--------------+--------------+
| 136     | -p_4n_2n        | P_42/m_n_m      | P42/mnm      | P42/mnm      |
+---------+-----------------+-----------------+--------------+--------------+
| 137:1   | p_4n_2n\_-1n    | P_42/n_m_c:1    | P42/nmc:1    | P42/nmc      |
+---------+-----------------+-----------------+--------------+--------------+
| 137:2   | -p_4ac_2a       | P_42/n_m_c:2    | P42/nmc:2    | P42/nmc      |
+---------+-----------------+-----------------+--------------+--------------+
| 138:1   | p_4n_2ab\_-1n   | P_42/n_c_m:1    | P42/ncm:1    | P42/ncm      |
+---------+-----------------+-----------------+--------------+--------------+
| 138:2   | -p_4ac_2ac      | P_42/n_c_m:2    | P42/ncm:2    | P42/ncm      |
+---------+-----------------+-----------------+--------------+--------------+
| 139     | -i_4_2          | I_4/m_m_m       | I4/mmm       | I4/mmm       |
+---------+-----------------+-----------------+--------------+--------------+
| 140     | -i_4_2c         | I_4/m_c_m       | I4/mcm       | I4/mcm       |
+---------+-----------------+-----------------+--------------+--------------+
| 141:1   | i_4bw_2bw\_-1bw | I_41/a_m_d:1    | I41/amd:1    | I41/amd      |
+---------+-----------------+-----------------+--------------+--------------+
| 141:2   | -i_4bd_2        | I_41/a_m_d:2    | I41/amd:2    | I41/amd      |
+---------+-----------------+-----------------+--------------+--------------+
| 142:1   | i_4bw_2aw\_-1bw | I_41/a_c_d:1    | I41/acd:1    | I41/acd      |
+---------+-----------------+-----------------+--------------+--------------+
| 142:2   | -i_4bd_2c       | I_41/a_c_d:2    | I41/acd:2    | I41/acd      |
+---------+-----------------+-----------------+--------------+--------------+
| 143     | p_3             | P_3             | P3           | P3           |
+---------+-----------------+-----------------+--------------+--------------+
| 144     | p_31            | P_31            | P31          | P31          |
+---------+-----------------+-----------------+--------------+--------------+
| 145     | p_32            | P_32            | P32          | P32          |
+---------+-----------------+-----------------+--------------+--------------+
| 146:h   | r_3             | R_3:h           | R3:h         | R3           |
+---------+-----------------+-----------------+--------------+--------------+
| 146:r   | p_3*            | R_3:r           | R3:r         | R3           |
+---------+-----------------+-----------------+--------------+--------------+
| 147     | -p_3            | P\_-3           | P-3          | P-3          |
+---------+-----------------+-----------------+--------------+--------------+
| 148:h   | -r_3            | R\_-3:h         | R-3:h        | R-3          |
+---------+-----------------+-----------------+--------------+--------------+
| 148:r   | -p_3*           | R\_-3:r         | R-3:r        | R-3          |
+---------+-----------------+-----------------+--------------+--------------+
| 149     | p_3_2           | P_3_1_2         | P312         | P312         |
+---------+-----------------+-----------------+--------------+--------------+
| 150     | p_3_2"          | P_3_2_1         | P321         | P321         |
+---------+-----------------+-----------------+--------------+--------------+
| 151     | p_31_2\_(0_0_4) | P_31_1_2        | P3112        | P3112        |
+---------+-----------------+-----------------+--------------+--------------+
| 152     | p_31_2"         | P_31_2_1        | P3121        | P3121        |
+---------+-----------------+-----------------+--------------+--------------+
| 153     | p_32_2\_(0_0_2) | P_32_1_2        | P3212        | P3212        |
+---------+-----------------+-----------------+--------------+--------------+
| 154     | p_32_2"         | P_32_2_1        | P3221        | P3221        |
+---------+-----------------+-----------------+--------------+--------------+
| 155:h   | r_3_2"          | R_3_2:h         | R32:h        | R32          |
+---------+-----------------+-----------------+--------------+--------------+
| 155:r   | p_3*_2          | R_3_2:r         | R32:r        | R32          |
+---------+-----------------+-----------------+--------------+--------------+
| 156     | p_3\_-2"        | P_3_m_1         | P3m1         | P3m1         |
+---------+-----------------+-----------------+--------------+--------------+
| 157     | p_3\_-2         | P_3_1_m         | P31m         | P31m         |
+---------+-----------------+-----------------+--------------+--------------+
| 158     | p_3\_-2"c       | P_3_c_1         | P3c1         | P3c1         |
+---------+-----------------+-----------------+--------------+--------------+
| 159     | p_3\_-2c        | P_3_1_c         | P31c         | P31c         |
+---------+-----------------+-----------------+--------------+--------------+
| 160:h   | r_3\_-2"        | R_3_m:h         | R3m:h        | R3m          |
+---------+-----------------+-----------------+--------------+--------------+
| 160:r   | p_3*\_-2        | R_3_m:r         | R3m:r        | R3m          |
+---------+-----------------+-----------------+--------------+--------------+
| 161:h   | r_3\_-2"c       | R_3_c:h         | R3c:h        | R3c          |
+---------+-----------------+-----------------+--------------+--------------+
| 161:r   | p_3*\_-2n       | R_3_c:r         | R3c:r        | R3c          |
+---------+-----------------+-----------------+--------------+--------------+
| 162     | -p_3_2          | P\_-3_1_m       | P-31m        | P-31m        |
+---------+-----------------+-----------------+--------------+--------------+
| 163     | -p_3_2c         | P\_-3_1_c       | P-31c        | P-31c        |
+---------+-----------------+-----------------+--------------+--------------+
| 164     | -p_3_2"         | P\_-3_m_1       | P-3m1        | P-3m1        |
+---------+-----------------+-----------------+--------------+--------------+
| 165     | -p_3_2"c        | P\_-3_c_1       | P-3c1        | P-3c1        |
+---------+-----------------+-----------------+--------------+--------------+
| 166:h   | -r_3_2"         | R\_-3_m:h       | R-3m:h       | R-3m         |
+---------+-----------------+-----------------+--------------+--------------+
| 166:r   | -p_3*_2         | R\_-3_m:r       | R-3m:r       | R-3m         |
+---------+-----------------+-----------------+--------------+--------------+
| 167:h   | -r_3_2"c        | R\_-3_c:h       | R-3c:h       | R-3c         |
+---------+-----------------+-----------------+--------------+--------------+
| 167:r   | -p_3*_2n        | R\_-3_c:r       | R-3c:r       | R-3c         |
+---------+-----------------+-----------------+--------------+--------------+
| 168     | p_6             | P_6             | P6           | P6           |
+---------+-----------------+-----------------+--------------+--------------+
| 169     | p_61            | P_61            | P61          | P61          |
+---------+-----------------+-----------------+--------------+--------------+
| 170     | p_65            | P_65            | P65          | P65          |
+---------+-----------------+-----------------+--------------+--------------+
| 171     | p_62            | P_62            | P62          | P62          |
+---------+-----------------+-----------------+--------------+--------------+
| 172     | p_64            | P_64            | P64          | P64          |
+---------+-----------------+-----------------+--------------+--------------+
| 173     | p_6c            | P_63            | P63          | P63          |
+---------+-----------------+-----------------+--------------+--------------+
| 174     | p\_-6           | P\_-6           | P-6          | P-6          |
+---------+-----------------+-----------------+--------------+--------------+
| 175     | -p_6            | P_6/m           | P6/m         | P6/m         |
+---------+-----------------+-----------------+--------------+--------------+
| 176     | -p_6c           | P_63/m          | P63/m        | P63/m        |
+---------+-----------------+-----------------+--------------+--------------+
| 177     | p_6_2           | P_6_2_2         | P622         | P622         |
+---------+-----------------+-----------------+--------------+--------------+
| 178     | p_61_2\_(0_0_5) | P_61_2_2        | P6122        | P6122        |
+---------+-----------------+-----------------+--------------+--------------+
| 179     | p_65_2\_(0_0_1) | P_65_2_2        | P6522        | P6522        |
+---------+-----------------+-----------------+--------------+--------------+
| 180     | p_62_2\_(0_0_4) | P_62_2_2        | P6222        | P6222        |
+---------+-----------------+-----------------+--------------+--------------+
| 181     | p_64_2\_(0_0_2) | P_64_2_2        | P6422        | P6422        |
+---------+-----------------+-----------------+--------------+--------------+
| 182     | p_6c_2c         | P_63_2_2        | P6322        | P6322        |
+---------+-----------------+-----------------+--------------+--------------+
| 183     | p_6\_-2         | P_6_m_m         | P6mm         | P6mm         |
+---------+-----------------+-----------------+--------------+--------------+
| 184     | p_6\_-2c        | P_6_c_c         | P6cc         | P6cc         |
+---------+-----------------+-----------------+--------------+--------------+
| 185     | p_6c\_-2        | P_63_c_m        | P63cm        | P63cm        |
+---------+-----------------+-----------------+--------------+--------------+
| 186     | p_6c\_-2c       | P_63_m_c        | P63mc        | P63mc        |
+---------+-----------------+-----------------+--------------+--------------+
| 187     | p\_-6_2         | P\_-6_m_2       | P-6m2        | P-6m2        |
+---------+-----------------+-----------------+--------------+--------------+
| 188     | p\_-6c_2        | P\_-6_c_2       | P-6c2        | P-6c2        |
+---------+-----------------+-----------------+--------------+--------------+
| 189     | p\_-6\_-2       | P\_-6_2_m       | P-62m        | P-62m        |
+---------+-----------------+-----------------+--------------+--------------+
| 190     | p\_-6c\_-2c     | P\_-6_2_c       | P-62c        | P-62c        |
+---------+-----------------+-----------------+--------------+--------------+
| 191     | -p_6_2          | P_6/m_m_m       | P6/mmm       | P6/mmm       |
+---------+-----------------+-----------------+--------------+--------------+
| 192     | -p_6_2c         | P_6/m_c_c       | P6/mcc       | P6/mcc       |
+---------+-----------------+-----------------+--------------+--------------+
| 193     | -p_6c_2         | P_63/m_c_m      | P63/mcm      | P63/mcm      |
+---------+-----------------+-----------------+--------------+--------------+
| 194     | -p_6c_2c        | P_63/m_m_c      | P63/mmc      | P63/mmc      |
+---------+-----------------+-----------------+--------------+--------------+
| 195     | p_2_2_3         | P_2_3           | P23          | P23          |
+---------+-----------------+-----------------+--------------+--------------+
| 196     | f_2_2_3         | F_2_3           | F23          | F23          |
+---------+-----------------+-----------------+--------------+--------------+
| 197     | i_2_2_3         | I_2_3           | I23          | I23          |
+---------+-----------------+-----------------+--------------+--------------+
| 198     | p_2ac_2ab_3     | P_21_3          | P213         | P213         |
+---------+-----------------+-----------------+--------------+--------------+
| 199     | i_2b_2c_3       | I_21_3          | I213         | I213         |
+---------+-----------------+-----------------+--------------+--------------+
| 200     | -p_2_2_3        | P_m\_-3         | Pm-3         | Pm-3         |
+---------+-----------------+-----------------+--------------+--------------+
| 201:1   | p_2_2_3\_-1n    | P_n\_-3:1       | Pn-3:1       | Pn-3         |
+---------+-----------------+-----------------+--------------+--------------+
| 201:2   | -p_2ab_2bc_3    | P_n\_-3:2       | Pn-3:2       | Pn-3         |
+---------+-----------------+-----------------+--------------+--------------+
| 202     | -f_2_2_3        | F_m\_-3         | Fm-3         | Fm-3         |
+---------+-----------------+-----------------+--------------+--------------+
| 203:1   | f_2_2_3\_-1d    | F_d\_-3:1       | Fd-3:1       | Fd-3         |
+---------+-----------------+-----------------+--------------+--------------+
| 203:2   | -f_2uv_2vw_3    | F_d\_-3:2       | Fd-3:2       | Fd-3         |
+---------+-----------------+-----------------+--------------+--------------+
| 204     | -i_2_2_3        | I_m\_-3         | Im-3         | Im-3         |
+---------+-----------------+-----------------+--------------+--------------+
| 205     | -p_2ac_2ab_3    | P_a\_-3         | Pa-3         | Pa-3         |
+---------+-----------------+-----------------+--------------+--------------+
| 206     | -i_2b_2c_3      | I_a\_-3         | Ia-3         | Ia-3         |
+---------+-----------------+-----------------+--------------+--------------+
| 207     | p_4_2_3         | P_4_3_2         | P432         | P432         |
+---------+-----------------+-----------------+--------------+--------------+
| 208     | p_4n_2_3        | P_42_3_2        | P4232        | P4232        |
+---------+-----------------+-----------------+--------------+--------------+
| 209     | f_4_2_3         | F_4_3_2         | F432         | F432         |
+---------+-----------------+-----------------+--------------+--------------+
| 210     | f_4d_2_3        | F_41_3_2        | F4132        | F4132        |
+---------+-----------------+-----------------+--------------+--------------+
| 211     | i_4_2_3         | I_4_3_2         | I432         | I432         |
+---------+-----------------+-----------------+--------------+--------------+
| 212     | p_4acd_2ab_3    | P_43_3_2        | P4332        | P4332        |
+---------+-----------------+-----------------+--------------+--------------+
| 213     | p_4bd_2ab_3     | P_41_3_2        | P4132        | P4132        |
+---------+-----------------+-----------------+--------------+--------------+
| 214     | i_4bd_2c_3      | I_41_3_2        | I4132        | I4132        |
+---------+-----------------+-----------------+--------------+--------------+
| 215     | p\_-4_2_3       | P\_-4_3_m       | P-43m        | P-43m        |
+---------+-----------------+-----------------+--------------+--------------+
| 216     | f\_-4_2_3       | F\_-4_3_m       | F-43m        | F-43m        |
+---------+-----------------+-----------------+--------------+--------------+
| 217     | i\_-4_2_3       | I\_-4_3_m       | I-43m        | I-43m        |
+---------+-----------------+-----------------+--------------+--------------+
| 218     | p\_-4n_2_3      | P\_-4_3_n       | P-43n        | P-43n        |
+---------+-----------------+-----------------+--------------+--------------+
| 219     | f\_-4a_2_3      | F\_-4_3_c       | F-43c        | F-43c        |
+---------+-----------------+-----------------+--------------+--------------+
| 220     | i\_-4bd_2c_3    | I\_-4_3_d       | I-43d        | I-43d        |
+---------+-----------------+-----------------+--------------+--------------+
| 221     | -p_4_2_3        | P_m\_-3_m       | Pm-3m        | Pm-3m        |
+---------+-----------------+-----------------+--------------+--------------+
| 222:1   | p_4_2_3\_-1n    | P_n\_-3_n:1     | Pn-3n:1      | Pn-3n        |
+---------+-----------------+-----------------+--------------+--------------+
| 222:2   | -p_4a_2bc_3     | P_n\_-3_n:2     | Pn-3n:2      | Pn-3n        |
+---------+-----------------+-----------------+--------------+--------------+
| 223     | -p_4n_2_3       | P_m\_-3_n       | Pm-3n        | Pm-3n        |
+---------+-----------------+-----------------+--------------+--------------+
| 224:1   | p_4n_2_3\_-1n   | P_n\_-3_m:1     | Pn-3m:1      | Pn-3m        |
+---------+-----------------+-----------------+--------------+--------------+
| 224:2   | -p_4bc_2bc_3    | P_n\_-3_m:2     | Pn-3m:2      | Pn-3m        |
+---------+-----------------+-----------------+--------------+--------------+
| 225     | -f_4_2_3        | F_m\_-3_m       | Fm-3m        | Fm-3m        |
+---------+-----------------+-----------------+--------------+--------------+
| 226     | -f_4a_2_3       | F_m\_-3_c       | Fm-3c        | Fm-3c        |
+---------+-----------------+-----------------+--------------+--------------+
| 227:1   | f_4d_2_3\_-1d   | F_d\_-3_m:1     | Fd-3m:1      | Fd-3m        |
+---------+-----------------+-----------------+--------------+--------------+
| 227:2   | -f_4vw_2vw_3    | F_d\_-3_m:2     | Fd-3m:2      | Fd-3m        |
+---------+-----------------+-----------------+--------------+--------------+
| 228:1   | f_4d_2_3\_-1ad  | F_d\_-3_c:1     | Fd-3c:1      | Fd-3c        |
+---------+-----------------+-----------------+--------------+--------------+
| 228:2   | -f_4ud_2vw_3    | F_d\_-3_c:2     | Fd-3c:2      | Fd-3c        |
+---------+-----------------+-----------------+--------------+--------------+
| 229     | -i_4_2_3        | I_m\_-3_m       | Im-3m        | Im-3m        |
+---------+-----------------+-----------------+--------------+--------------+
| 230     | -i_4bd_2c_3     | I_a\_-3_d       | Ia-3d        | Ia-3d        |
+---------+-----------------+-----------------+--------------+--------------+
"""
