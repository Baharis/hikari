from pkgutil import get_data as get
import io
import json
import pickle
import pandas as pd


def _load_indexed_csv(filename):
    return pd.read_csv(io.StringIO(get(__name__, filename).decode('utf-8')), 
                       comment='#', index_col=0)


def _load_json(filename):
    return json.loads(get(__name__, filename).decode('utf-8'))


def _load_pickle(filename):
    return pickle.loads(get(__name__, filename))


def _save_pickle(data, filename):
    pickle.dump(data, open(filename, 'wb'), protocol=4)


gnuplot_angular_heatmap_template = \
    get(__name__, 'gnuplot_angular_heatmap_template.gnu').decode('utf-8')
point_groups_json = _load_json('point_groups.json')
space_groups_json = _load_json('space_groups.json')
point_groups_pickle = get(__name__, 'point_groups.pickle')
space_groups_pickle = get(__name__, 'space_groups.pickle')
hkl_formats = _load_json('hkl_formats_defined.json')
hkl_aliases = _load_json('hkl_formats_aliases.json')
hkl_mercury_style = get(__name__, 'hkl.msd').decode('utf-8')
characteristic_radiation = _load_json('characteristic_radiation.json')
cif_core_dict = get(__name__, 'cif_core_2.4.5.dic').decode('utf-8')
Xray_atomic_form_factors = _load_indexed_csv('Xray_atomic_form_factors.csv')
