from pkgutil import get_data as get
import json
import pickle


def _load_json(filename):
    return json.loads(get(__name__, filename).decode('utf-8'))


def _load_pickle(filename):
    return pickle.loads(get(__name__, filename))


potency_map_template = get(__name__, 'cplt_map_template.gnu').decode('utf-8')
point_groups_json = _load_json('point_groups.json')
space_groups_json = _load_json('space_groups.json')
point_groups_pickle = get(__name__, 'point_groups.pickle')
space_groups_pickle = get(__name__, 'space_groups.pickle')
point_groups_dictionary = _load_pickle('point_groups.pickle')
space_groups_dictionary = _load_pickle('space_groups.pickle')
hkl_formats = _load_json('hkl_formats_defined.json')
hkl_aliases = _load_json('hkl_formats_aliases.json')
