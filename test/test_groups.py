import unittest

from hikari.resources import hall_symbols_table, point_groups_json,\
    space_groups_json
from hikari.symmetry import Group, PG, SymmOp, SG
from hikari.symmetry.group import _unpack_group_dictionary_from_json  # noqa


class TestGroup(unittest.TestCase):
    def test_group_initializes(self):
        sg230_generator_codes = [
            'x,y,z', '-x+1/2,-y,z+1/2', '-x,y+1/2,-z+1/2', 'z,x,y',
            'y+3/4,x+1/4,-z+1/4', '-x,-y,-z', 'x+1/2,y+1/2,z+1/2'
        ]
        sg230_generators = [SymmOp.from_code(c) for c in sg230_generator_codes]
        _ = Group(*sg230_generators)

    def test_point_groups_consistent_with_pickles(self):
        json_pg_dict = _unpack_group_dictionary_from_json(point_groups_json)
        pickled_pg_dict = PG
        all_pg_keys = list({*json_pg_dict.keys(), *pickled_pg_dict.keys()})
        for pg_key in all_pg_keys:
            json_pg = json_pg_dict[pg_key]
            pickled_pg = pickled_pg_dict[pg_key]
            self.assertEqual(json_pg, pickled_pg)

    def test_space_groups_consistent_with_pickles(self):
        json_sg_dict = _unpack_group_dictionary_from_json(space_groups_json)
        pickled_sg_dict = SG
        all_sg_keys = list({*json_sg_dict.keys(), *pickled_sg_dict.keys()})
        for sg_key in all_sg_keys:
            json_sg = json_sg_dict[sg_key]
            pickled_sg = pickled_sg_dict[sg_key]
            self.assertEqual(json_sg, pickled_sg)

    def test_space_groups_picked_match_hall(self):
        for i in range(1, 231):
            sg_index_regex = rf'^{i}(?::[\da-z-]+)?$'
            mask = hall_symbols_table.index.str.contains(sg_index_regex)
            sg_index_first = list(mask).index(True)
            hall_symbol = hall_symbols_table['Hall entry'].iloc[sg_index_first]
            print(str(i) + ' ' + hall_symbol)
            self.assertEqual(SG[i], Group.from_hall_symbol(hall_symbol))
