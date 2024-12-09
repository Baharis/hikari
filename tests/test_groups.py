import os
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


class TestPredefinedGroupCatalogs(unittest.TestCase):
    def test_group_catalogs_length(self):
        self.assertGreaterEqual(len(PG), 43)
        self.assertGreaterEqual(len(SG), 530)

    def test_group_catalogs_standard(self):
        self.assertEqual(len(PG.standard), 32)
        self.assertEqual(len(SG.standard), 230)

    def test_group_catalogs_dict_interface(self):
        for gc in [PG, SG]:
            keys = list(gc.keys())
            values = list(gc.values())
            items = list(gc.items())
            for k1, v1, (k2, v2) in zip(keys, values, items):
                self.assertEqual(k1, k2)
                self.assertIs(v1, v2)

    @unittest.skipIf(os.environ.get('HIKARI_TESTS_RAPID'), reason='Rapid test mode')
    def test_point_group_catalogs_consistent_wsv_json_pg(self):
        pgc1 = GroupCatalog(point_groups_dataframe)
        pgc2 = GroupCatalog.from_json(point_groups_json)
        pd.testing.assert_frame_equal(pgc1.table, pgc2.table)
        pd.testing.assert_frame_equal(pgc1.table, PG.table)

    @unittest.skipIf(os.environ.get('HIKARI_TESTS_RAPID'), reason='Rapid test mode')
    def test_space_group_catalogs_consistent_wsv_json_sg(self):
        sgc1 = GroupCatalog(space_groups_dataframe)
        sgc2 = GroupCatalog.from_json(space_groups_json)
        pd.testing.assert_frame_equal(sgc1.table, sgc2.table)
        pd.testing.assert_frame_equal(sgc1.table, SG.table)
