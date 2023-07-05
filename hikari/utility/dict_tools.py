"""
This file contains tools to work with dictionaries used in the package.
"""


def dict_union(*dicts: dict) -> dict:
    """Return a union of dicts; if conflicts, later dicts take precedence."""
    return dict(i for d in dicts for i in d.items())
