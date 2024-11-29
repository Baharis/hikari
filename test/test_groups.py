import unittest

from hikari.resources import (point_groups_dataframe, space_groups_dataframe,
                              point_groups_json, space_groups_json)
from hikari.symmetry import Group, PG, SymmOp, SG
from hikari.symmetry.catalog import GroupCatalog

import pandas as pd


class TestGroup(unittest.TestCase):
    def test_group_initializes(self):
        sg230_generator_codes = [
            'x,y,z', '-x+1/2,-y,z+1/2', '-x,y+1/2,-z+1/2', 'z,x,y',
            'y+3/4,x+1/4,-z+1/4', '-x,-y,-z', 'x+1/2,y+1/2,z+1/2'
        ]
        sg230_generators = [SymmOp.from_code(c) for c in sg230_generator_codes]
        _ = Group(*sg230_generators)

    def test_point_groups_pickle_consistent_with_wsv(self):
        pgc1 = GroupCatalog(point_groups_dataframe)
        pgc2 = GroupCatalog.from_json(point_groups_json)
        pd.testing.assert_frame_equal(pgc1.table, pgc2.table)
        pd.testing.assert_frame_equal(pgc1.table, PG.table)

    def test_space_groups_pickle_consistent_with_wsv(self):
        sgc1 = GroupCatalog(space_groups_dataframe)
        sgc2 = GroupCatalog.from_json(space_groups_json)
        pd.testing.assert_frame_equal(sgc1.table, sgc2.table)
        pd.testing.assert_frame_equal(sgc1.table, SG.table)
