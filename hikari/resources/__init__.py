from __future__ import annotations

import io
import json

import pandas as pd
from importlib.resources import open_text, open_binary


def _load_bytes(resource_name: str) -> bytes:
    with open_binary(__name__, resource_name) as f:
        return f.read()


def _load_indexed_csv(resource_name: str) -> pd.DataFrame:
    with open_text(__name__, resource_name, encoding='utf-8') as f:
        s = io.StringIO(f.read())
    return pd.read_csv(s, comment='#', index_col=0)


def _load_indexed_wsv(resource_name: str) -> pd.DataFrame:
    with open_text(__name__, resource_name, encoding='utf-8') as f:
        s = io.StringIO(f.read())
    return pd.read_csv(s, comment='#', index_col=0, sep=r'\s+')


def _load_json(resource_name: str) -> dict:
    with open_text(__name__, resource_name, encoding='utf-8') as f:
        return json.load(f)


def _load_text(resource_name: str) -> str:
    with open_text(__name__, resource_name, encoding='utf-8') as f:
        return f.read()


gnuplot_angular_heatmap_template = _load_text('gnuplot_angular_heatmap_template.gnu')
point_groups_pickle = _load_bytes('point_groups.pickle')
space_groups_pickle = _load_bytes('space_groups.pickle')
hkl_formats = _load_json('hkl_formats_defined.json')
hkl_aliases = _load_json('hkl_formats_aliases.json')
hkl_mercury_style = _load_text('hkl.msd')
characteristic_radiation = _load_json('characteristic_radiation.json')
cif_core_dict = _load_text('cif_core_2.4.5.dic')
Xray_atomic_form_factors = _load_indexed_csv('Xray_atomic_form_factors.csv')
hall_symbols_pg_table = _load_indexed_wsv('Hall_symbols_PG.wsv')
hall_symbols_sg_table = _load_indexed_wsv('Hall_symbols_SG.wsv')
