from pkgutil import get_data as get
import json


def _load_json(filename):
    return json.loads(get(__name__, filename).decode('utf-8'))


potency_map_template = get(__name__, 'cplt_map_template.gnu').decode('utf-8')
point_groups_json = get(__name__, 'point_groups.json').decode('utf-8')
space_groups_json = get(__name__, 'space_groups.json').decode('utf-8')
point_groups_pickle = get(__name__, 'point_groups.pickle')
space_groups_pickle = get(__name__, 'space_groups.pickle')
hkl_formats = _load_json('hkl_formats_defined.json')
hkl_aliases = _load_json('hkl_formats_aliases.json')
