from pkgutil import get_data as get

potency_map_template = get(__name__, 'cplt_map_template.gnu').decode('utf-8')
