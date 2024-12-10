import os
import tempfile
import unittest
import warnings

from hikari.resources import (point_groups_dataframe, space_groups_dataframe,
                              point_groups_json, space_groups_json)
from hikari.symmetry import Group, PG, SymmOp, SG
from hikari.symmetry.catalog import GroupCatalog, AmbiguousGroupAccessorWarning

import pandas as pd


class TestGroup(unittest.TestCase):
    def test_group_initializes(self):
        sg230_generator_codes = [
            'x,y,z', '-x+1/2,-y,z+1/2', '-x,y+1/2,-z+1/2', 'z,x,y',
            'y+3/4,x+1/4,-z+1/4', '-x,-y,-z', 'x+1/2,y+1/2,z+1/2'
        ]
        sg230_generators = [SymmOp.from_code(c) for c in sg230_generator_codes]
        _ = Group(*sg230_generators)


class TestPointGroupCatalog(unittest.TestCase):
    catalogue_object: GroupCatalog = PG
    catalogue_dataframe: pd.DataFrame = point_groups_dataframe
    catalogue_json_text: str = point_groups_json
    catalogue_length: int = 43
    catalogue_standards: int = 32
    catalogue_sample_HM_simple: str = '2/m'

    def test_group_catalogs_length(self):
        self.assertGreaterEqual(len(self.catalogue_object), b=self.catalogue_length)

    def test_group_catalogs_standard(self):
        self.assertEqual(len(self.catalogue_object.standard),
                         second=self.catalogue_standards)

    def test_group_catalogs_dict_interface(self):
        keys = list(self.catalogue_object.keys())
        values = list(self.catalogue_object.values())
        items = list(self.catalogue_object.items())
        for k1, v1, (k2, v2) in zip(keys, values, items):
            self.assertEqual(k1, k2)
            self.assertIs(v1, v2)

    def test_group_catalogs_json(self):
        json_file = tempfile.NamedTemporaryFile(mode='w+', delete=False)
        try:
            self.catalogue_object.to_json(json_file.name)
            json_file.flush()
            json_file.seek(0)
            json_text = json_file.read()
            self.assertEqual(json_text, self.catalogue_json_text)
        finally:
            json_file.close()
            os.unlink(json_file.name)

    def test_group_catalogs_rest(self):
        """No simple mechanism to access instance docstring - test write only"""
        with tempfile.NamedTemporaryFile(mode='w+', delete=False) as rest_file:
            self.catalogue_object.to_rest_table(rest_file.name)

    def test_group_catalogs_getitem(self):
        with self.assertWarns(AmbiguousGroupAccessorWarning):  # smart
            g1 = self.catalogue_object[self.catalogue_sample_HM_simple]
        with warnings.catch_warnings(record=True) as caught_warnings:  # explicit
            g2 = self.catalogue_object.get(
                HM_simple=self.catalogue_sample_HM_simple,
                standard=True)
            self.assertEqual(len(caught_warnings), 0)
        self.assertEqual(g1, g2)

    @unittest.skipIf(os.environ.get('HIKARI_TESTS_RAPID'), reason='Rapid test mode')
    def test_group_catalogs_consistency(self):
        gc1 = GroupCatalog(self.catalogue_dataframe)
        gc2 = GroupCatalog.from_json(self.catalogue_json_text)
        gc3 = self.catalogue_object
        self.assertEqual(gc1, gc2)
        self.assertEqual(gc1, gc3)


class TestSpaceGroupCatalog(TestPointGroupCatalog):
    catalogue_object: GroupCatalog = SG
    catalogue_dataframe: pd.DataFrame = space_groups_dataframe
    catalogue_json_text: str = space_groups_json
    catalogue_length: int = 530
    catalogue_standards: int = 230
    catalogue_sample_HM_simple: str = 'P2/m'
