from pkgutil import get_data as get

potency_map_template = get(__name__, 'cplt_map_template.gnu').decode('utf-8')
point_groups_json = get(__name__, 'point_groups.json').decode('utf-8')
space_groups_json = get(__name__, 'space_groups.json').decode('utf-8')
