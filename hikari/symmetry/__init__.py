"""
Module containing symmetry operations and point groups
to be used by other objects.
"""
import pickle
from .operations import SymmOp
from .group import Group
from hikari.resources import point_groups_pickle
from hikari.resources import space_groups_pickle
PG = pickle.loads(point_groups_pickle)
SG = pickle.loads(space_groups_pickle)
