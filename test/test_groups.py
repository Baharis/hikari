import unittest

from hikari.resources import (hall_symbols_pg_table, hall_symbols_sg_table,
                              point_groups_pickle, space_groups_pickle)
from hikari.symmetry import Group, PG, SymmOp, SG
from hikari.symmetry.catalog import PointGroupCatalog, SpaceGroupCatalog


class TestGroup(unittest.TestCase):
    def test_group_initializes(self):
        sg230_generator_codes = [
            'x,y,z', '-x+1/2,-y,z+1/2', '-x,y+1/2,-z+1/2', 'z,x,y',
            'y+3/4,x+1/4,-z+1/4', '-x,-y,-z', 'x+1/2,y+1/2,z+1/2'
        ]
        sg230_generators = [SymmOp.from_code(c) for c in sg230_generator_codes]
        _ = Group(*sg230_generators)

    def test_point_groups_pickle_consistent_with_wsv(self):
        pgc1 = PointGroupCatalog(hall_symbols_pg_table)
        pgc2 = PointGroupCatalog.from_bytes(point_groups_pickle)
        all_keys = list({*pgc1.keys(), *pgc2.keys()})
        for key in all_keys:
            self.assertEqual(pgc1.get(n_c=key), pgc2.get(n_c=key))

    def test_space_groups_pickle_consistent_with_wsv(self):
        sgc1 = SpaceGroupCatalog(hall_symbols_sg_table)
        sgc2 = SpaceGroupCatalog.from_bytes(space_groups_pickle)
        all_keys = list({*sgc1.keys(), *sgc2.keys()})
        for key in all_keys:
            self.assertEqual(sgc1.get(n_c=key), sgc2.get(n_c=key))
