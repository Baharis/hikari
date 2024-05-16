from copy import deepcopy
import pickle
import re

import numpy as np
import pandas as pd

from hikari.utility.typing import PathLike
from hikari.symmetry.group import Group


class GroupCatalogue:
    """Manage generating and mappings of point and space groups"""

    REST_COL_FORMAT = {'n:c': '7.7s', 'HM-short': '9.9s',
                       #'Schoenflies': '7.7s',
                       'Hall': '14.14s'}

    def __init__(self, table: pd.DataFrame) -> 'GroupCatalogue':
        has_colon = table['n:c'].str.contains(':')
        table['n:c0'] = table['n:c']
        table.loc[~has_colon, 'n:c0'] = table.loc[~has_colon, 'n:c'] + ':0'
        table[['number', 'setting']] = table['n:c'].str.split(':', expand=True)
        table['number'] = table['number'].astype(int)
        if 'group' not in table.columns:
            table['group'] = [Group.from_hall_symbol(h) for h in table['Hall']]
        if 'HM-short' not in table.columns:
            table['HM-short'] = table['HM'].str.replace('_', '')
        if 'HM-simple' not in table.columns:
            colon = table['HM'].str.contains(':')
            gsm = Group.System.monoclinic
            m = np.array([g.system == gsm for g in table['group']])
            table['HM-simple'] = table['HM-short']
            table.loc[colon, 'HM-simple'] = \
                table.loc[colon, 'HM-simple'].str.replace(r'\:.', '', regex=True)
            table.loc[m, 'HM-simple'] = table.loc[m, 'HM'].str.replace('_1', '')
            table['HM-simple'] = table['HM-simple'].str.replace('_', '')
        if 'standard' not in table.columns:
            table['standard'] = table['number'].rolling(2).var() != 0
        self.table: pd.DataFrame = table

    @classmethod
    def _from_ssv(cls, ssv_path: PathLike) -> 'GroupCatalogue':
        """Generate a `GroupCatalogue` from a .ssv file"""
        with open(ssv_path, 'r') as ssv_file:
            table = pd.read_csv(ssv_file, comment='#', delim_whitespace=True)
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
        field_len_re = re.compile('^(\d*)(?:\.\d+)?[a-z]*$')
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
        """Dump an instance of self to a .pickle file"""
        with open(pickle_path, 'bw') as pickle_file:
            pickle.dump(self, file=pickle_file)

    @property
    def standard(self) -> 'GroupCatalogue':
        standard = deepcopy(self.table[self.table['standard']]).reset_index()
        return self.__class__(standard)

    def export_dict(self, col: str) -> dict:
        exported = self.table[[col, 'group']]
        return {r[0]: r[1] for r in exported.itertuples(index=False)}
